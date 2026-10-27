#include "core/system/ins_loc_system.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <rclcpp/rclcpp.hpp>
#include <glog/logging.h>
#include <std_msgs/msg/string.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <cgi430_interfaces/msg/latitude.hpp>
#include <cgi430_interfaces/msg/longitude.hpp>
#include <cgi430_interfaces/msg/altitude.hpp>
#include <cgi430_interfaces/msg/attitude.hpp>
#include <cgi430_interfaces/msg/earth_velocity.hpp>
#include <cgi430_interfaces/msg/position_sigma.hpp>
#include <cgi430_interfaces/msg/attitude_sigma.hpp>
#include <cgi430_interfaces/msg/velocity_sigma.hpp>
#include <cgi430_interfaces/msg/ins_status.hpp>
#include "core/navigation/ins_navigation.h"
#include "core/navigation/ins_cloud.h"
#include "core/lio/multi_lidar_fusion.h"
#include "core/lio/pointcloud_preprocess.h"
#include "core/system/sany_localization_output.h"
#include "utils/functional_safety_heartbeat.h"
#include "wrapper/ros_utils.h"

namespace lightning {
namespace {
double Steady() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
double Positive(const YAML::Node& n, const char* key, double fallback) {
    const double v = n[key] ? n[key].as<double>() : fallback;
    if (!std::isfinite(v) || v <= 0) throw std::invalid_argument(std::string("invalid ") + key);
    return v;
}
std::filesystem::path RelativeTo(const std::string& yaml, const std::string& path) {
    std::filesystem::path p(path);
    return p.is_absolute() ? p : std::filesystem::path(yaml).parent_path()/p;
}
SE3 ReadExtrinsic(const YAML::Node& n) {
    if (!n || !n["confirmed"].as<bool>()) throw std::invalid_argument("cloud.rear_from_primary must be calibrated and confirmed");
    const auto t = n["translation_m"].as<std::vector<double>>();
    const auto q = n["quaternion_xyzw"].as<std::vector<double>>();
    if (t.size()!=3 || q.size()!=4) throw std::invalid_argument("invalid rear_from_primary dimensions");
    const Vec3d trans(t[0],t[1],t[2]);
    Quatd quat(q[3],q[0],q[1],q[2]);
    if (!trans.allFinite() || !quat.coeffs().allFinite() || std::abs(quat.norm()-1)>1e-3)
        throw std::invalid_argument("invalid rear_from_primary transform");
    quat.normalize(); return SE3(quat,trans);
}
}

struct InsLocSystem::Impl {
    rclcpp::Node::SharedPtr node;  ///< INS 在线定位的 ROS 节点。
    std::vector<rclcpp::SubscriptionBase::SharedPtr> subscriptions;  ///< 各导航字段的 ROS 订阅器，持有以维持订阅生命周期。
    std::unique_ptr<functional_safety::HeartbeatPublisher> heartbeat;  ///< 独立发布 INS 定位节点的功能安全心跳。
    rclcpp::TimerBase::SharedPtr timer;  ///< 周期检查健康状态与待投影点云的定时器。
    std::unique_ptr<tf2_ros::TransformBroadcaster> tf;  ///< 可选的地图到后轴 TF 发布器。
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub;  ///< 发布后轴位姿 PoseStamped。
    rclcpp::Publisher<lightning::msg::VehiclePose>::SharedPtr vehicle_pub;  ///< 发布车辆位姿与速度。
    rclcpp::Publisher<lightning::msg::LocalizationStatus>::SharedPtr status_pub;  ///< 发布定位运行状态。
    rclcpp::Publisher<lightning::msg::FaultStatus>::SharedPtr fault_pub;  ///< 发布导航质量失效的故障状态。
    rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr path_pub;  ///< 发布后轴可视化轨迹。
    rclcpp::Publisher<std_msgs::msg::String>::SharedPtr diagnostic_pub;  ///< 发布 INS 导航质量与点云处理诊断。
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr inv_pub /**< 发布后轴坐标系下的投影点云。 */, map_pub /**< 发布地图坐标系下的投影点云。 */;
    ins::Sample sample;  ///< 最近收到的各类 INS 导航字段与质量信息。
    ins::QualityPolicy policy;  ///< 导航字段时效、标准差和恢复门控参数。
    std::unique_ptr<ins::GeoReference> reference;  ///< 从大地坐标到配置局部地图的坐标参考。
    std::unique_ptr<ins::RecoveryGate> gate;  ///< 连续有效导航样本的恢复门控。
    std::unique_ptr<ins::PoseBuffer> poses;  ///< 为点云扫描时刻插值提供的后轴位姿历史。
    std::array<double,5> consumed{};  ///< 各主导航字段最近已消费的时间戳（s）。
    std::string reason = "awaiting_navigation" /**< 当前导航不可用、恢复或有效状态的原因。 */, frame = "map" /**< 业务地图输出坐标系名称。 */, rear_frame = "rear_axle" /**< 后轴坐标系名称。 */;
    double offset = 0 /**< CAN 接收时间戳到 ROS 时间域的校正偏移（s）。 */, geoid_separation = 0 /**< 将设备高度换算为椭球高使用的大地水准面差（m）。 */, latest_published = 0 /**< 最近记录发布后轴位姿的时间戳（s）。 */, cloud_wait = 0.3 /**< 等待导航位姿覆盖点云扫描区间的超时（s）。 */;
    double last_path_stamp = 0 /**< 最近一次可视化路径采样的时间戳（s）。 */, max_radius = 20000 /**< 允许发布的位姿到局部地图原点最大距离（m）。 */;
    bool clouds_enabled = false /**< 是否启用融合点云的后轴与地图投影发布。 */, ever_good = false /**< 导航是否曾通过质量与恢复门控。 */, recording_failed = false /**< 审计记录发生写入失败后的故障锁存标志。 */;
    std::size_t cloud_drops = 0 /**< 点云无法获得有效导航覆盖而丢弃的次数。 */, cloud_count = 0 /**< 控制地图点云降频发布的输入帧计数。 */;
    int map_every = 1;  ///< 地图点云发布的帧间隔。
    std::ofstream raw_log /**< 原始 CAN 导航字段的审计文件流。 */, solutions /**< 导航转换与质量判定结果的审计文件流。 */, tum /**< 实际发布的后轴位姿 TUM 文件流。 */, events /**< 导航状态变化事件记录文件流。 */;
    std::filesystem::path audit_directory;  ///< 本次运行的审计产物目录。
    nav_msgs::msg::Path path;  ///< 用于发布的后轴可视化路径。
    std::atomic<std::uint32_t> pose_seq{0} /**< 车辆位姿消息的通信序列号。 */, fault_seq{0} /**< 故障状态消息的通信序列号。 */;
    MultiLidarConfig lidar_config;  ///< 多雷达组帧及固定外参配置。
    MultiLidarFrameAssembler assembler;  ///< 将各雷达点云组成主雷达坐标系融合帧。
    SelfPointFilterConfig self_filter;  ///< 投影前的车体自反射过滤配置。
    std::map<int,std::shared_ptr<PointCloudPreprocess>> preprocessors;  ///< 按雷达 ID 保存的点云预处理器。
    SE3 rear_from_primary;  ///< 将主雷达点坐标变换到后轴坐标的固定外参。
    double cloud_leaf = 0;  ///< 投影点云的体素降采样分辨率（m）。
    struct Pending { FusedLidarFrame frame /**< 等待导航位姿覆盖扫描区间的融合点云帧。 */; double arrival /**< 该待处理点云到达的单调时钟时刻（s）。 */; };
    std::deque<Pending> pending;  ///< 等待导航历史覆盖扫描区间的融合点云队列。

    void Reject(const std::string& why) {
        // 质量中断也切断插值历史，恢复后不能跨越这次失效拼接点云。
        gate->Reject(); poses->Clear();
        if (reason != why) {
            events << std::setprecision(16) << node->now().seconds() << ',' << why << '\n';
            LOG(WARNING) << "INS_ONLY unavailable: " << why;
        }
        reason = why;
        heartbeat->SetState(ever_good ? functional_safety::NodeState::kFault : functional_safety::NodeState::kNotReady);
    }

    void Update(ins::Field field, const cgi430_interfaces::msg::CanFrame& f) {
        const double stamp = ToSec(f.header.stamp) + offset;
        sample.valid[field] = f.valid && std::isfinite(stamp) && stamp > 0 && stamp > sample.time[field];
        if (!sample.valid[field]) Reject("invalid_or_nonmonotonic_field_"+std::to_string(field));
        if (stamp > sample.time[field]) sample.time[field] = stamp;
        sample.arrival[field] = Steady();
        raw_log << std::setprecision(16) << ToSec(f.header.stamp) << ',' << stamp << ','
                << static_cast<int>(field) << ',' << f.can_id << ',' << f.valid << ',' << static_cast<int>(f.dlc);
        for (auto v : f.data) raw_log << ',' << static_cast<int>(v);
        raw_log << '\n';
    }

    template<class M, class F>
    void Subscribe(const std::string& topic, ins::Field field, F extract) {
        subscriptions.push_back(node->create_subscription<M>(topic, rclcpp::SensorDataQoS().keep_last(100),
            [this,field,extract](const std::shared_ptr<M> msg) {
                Update(field,msg->frame); extract(*msg); OnNavigation();
            }));
    }

    void LogSolution(const ins::Pose* pose, bool published, const std::string& why) {
        solutions << std::setprecision(16) << sample.Stamp() << ',' << sample.latitude << ','
                  << sample.longitude << ',' << sample.altitude << ',' << sample.heading << ','
                  << sample.pitch << ',' << sample.roll << ',' << sample.system_state << ','
                  << sample.satellite_status << ',' << sample.differential_age;
        for (const auto* v : {&sample.velocity, &sample.position_sigma, &sample.attitude_sigma, &sample.velocity_sigma})
            for (int i=0;i<3;++i) solutions << ',' << (*v)[i];
        for (double t : sample.time) solutions << ',' << t;
        if (pose) {
            const auto& p=pose->map_from_rear.translation();
            const auto q=pose->map_from_rear.unit_quaternion();
            solutions << ',' << p.x() << ',' << p.y() << ',' << p.z() << ','
                      << q.x() << ',' << q.y() << ',' << q.z() << ',' << q.w();
            for (int i=0;i<3;++i) solutions << ',' << pose->map_velocity[i];
            solutions << ',' << pose->SignedSpeed();
        } else for (int i=0;i<11;++i) solutions << ",nan";
        solutions << ',' << published << ',' << why << '\n';
    }

    /**
     * @brief 在字段更新后评估导航解，驱动恢复门控与业务输出。
     * @par 主要流程
     * 1. <b>检查质量与字段推进</b>：拒绝审计故障，检查质量；主导航字段未全部更新时返回。
     * 2. <b>转换坐标并记录拒绝结果</b>：修正高度后转换到地图坐标，拒绝质量失效或超出场地范围的解。
     * 3. <b>推进恢复门控与位姿历史</b>：保存带有效标志的位姿，恢复样本不足时停止本次发布。
     * 4. <b>发布导航并处理待用点云</b>：记录通过结果，调用 Publish()，再尝试消耗点云队列。
     * @note 阶段名称对应函数体；质量失效解可进入审计转换，但不会作为有效位姿发布。
     */
    void OnNavigation() {
        // 阶段 1：检查质量与字段推进。
        if (recording_failed) { Reject("recording_failed"); return; }
        // [ins-new-navigation-fields]
        const std::string bad = ins::CheckSample(sample,policy,node->now().seconds(),Steady());
        if (!bad.empty()) Reject(bad);
        // Position, attitude and velocity must each advance. Quality messages
        // can be lower rate only within the same explicit skew/age bounds.
        for (std::size_t i=0;i<consumed.size();++i) if (sample.time[i] <= consumed[i]) return;
        for (std::size_t i=0;i<consumed.size();++i) consumed[i] = sample.time[i];
        // [ins-new-navigation-fields]
        // 阶段 2：转换坐标并记录拒绝结果。
        ins::Sample converted = sample;
        converted.altitude += geoid_separation;  // h = H + N; explicit datum only.
        ins::Pose pose;
        try {
            if (!std::isfinite(sample.heading) || !std::isfinite(sample.pitch) || !std::isfinite(sample.roll) ||
                !sample.velocity.allFinite()) throw std::invalid_argument("nonfinite navigation");
            pose = ins::Convert(converted,*reference);
        } catch (const std::exception&) {
            Reject("coordinate_conversion"); LogSolution(nullptr,false,reason); return;
        }
        // Rejected solutions remain available for evaluation, never for business output.
        if (!bad.empty()) { LogSolution(&pose,false,bad); return; }
        if (pose.map_from_rear.translation().norm() > max_radius) {
            Reject("outside_site_radius"); LogSolution(&pose,false,reason); return;
        }
        // 阶段 3：推进恢复门控与位姿历史。
        // 恢复中的位姿也入历史，但 global_valid=false 使点云不能用它们跨区间插值。
        // [ins-navigation-recovery]
        pose.global_valid = gate->Observe(pose.stamp,true);
        if (!poses->Add(pose)) { Reject("nonmonotonic_pose"); LogSolution(&pose,false,reason); return; }
        if (!gate->Good()) { reason="recovering"; LogSolution(&pose,false,reason); return; }
        // 阶段 4：发布导航并处理待用点云。
        reason="good";
        LogSolution(&pose,true,reason);
        Publish(pose);
        // [ins-navigation-recovery]
        DrainClouds();
    }

    /**
     * @brief 将通过门控的后轴导航解发布为业务位姿和速度。
     * @par 主要流程
     * 1. <b>检查门控并构造消息</b>：拒绝时间不递增或恢复门控未通过的位姿，填充消息与通信头。
     * 2. <b>发布并更新工作状态</b>：发送业务消息、尝试写 TUM，并推进心跳状态和工作计数。
     * 3. <b>更新 TF 与轨迹</b>：按开关发送 TF，按时间间隔更新有界可视化路径。
     * @note 阶段名称对应函数体；日志写入故障由 Health() 周期检查，不是发布前的落盘事务。
     */
    void Publish(const ins::Pose& p) {
        // 阶段 1：检查门控并构造消息。
        // [ins-business-publish]
        if (p.stamp <= latest_published || !gate->Good()) return;
        auto vehicle=sany_output::MakeVehiclePoseMessage(p.map_from_rear,p.SignedSpeed(),p.stamp,frame);
        auto pose=sany_output::MakePoseMessage(vehicle);
        vehicle.comm_header.source_id=heartbeat->NodeId();
        vehicle.comm_header.seq=functional_safety::HeartbeatPublisher::Advance(pose_seq);
        vehicle.comm_header.stamp_us=functional_safety::HeartbeatPublisher::UnixMicrosecondsNow();
        // 阶段 2：发布并更新工作状态。
        pose_pub->publish(pose); vehicle_pub->publish(vehicle);
        sany_output::WriteTumPoseLine(tum,pose,latest_published);
        heartbeat->SetState(functional_safety::NodeState::kRunning); heartbeat->RecordWork();
        ever_good=true;
        // [ins-business-publish]
        // 阶段 3：更新 TF 与轨迹。
        if (tf) {
            geometry_msgs::msg::TransformStamped m;
            m.header=pose.header; m.child_frame_id=rear_frame;
            m.transform.translation.x=pose.pose.position.x; m.transform.translation.y=pose.pose.position.y;
            m.transform.translation.z=pose.pose.position.z; m.transform.rotation=pose.pose.orientation;
            tf->sendTransform(m);
        }
        if (p.stamp-last_path_stamp >= 0.1) {
            path.header=pose.header; path.poses.push_back(pose); last_path_stamp=p.stamp;
            if (path.poses.size()>500) path.poses.erase(path.poses.begin());
        }
    }

    /**
     * @brief 校验并预处理单雷达点云，组帧后等待导航位姿覆盖扫描区间。
     * @par 主要流程
     * 1. <b>校验点布局</b>：检查紧密单行、本机字节序及 Livox 字段类型和边界。
     * 2. <b>预处理并检查点值</b>：转换点云，拒绝空云、非有限坐标或无效点时间。
     * 3. <b>组帧入队并尝试发布</b>：添加雷达帧，对融合帧做自反射过滤，限制等待队列后调用 DrainClouds()。
     * @note 阶段名称对应函数体；入队不代表导航覆盖充分，实际发布还需通过去畸变检查。
     */
    void Cloud(int id, const sensor_msgs::msg::PointCloud2::SharedPtr& msg) {
        // 阶段 1：校验点布局。
        // The reused Livox iterator expects one packed, native-endian scan.
        if (msg->height!=1 || !msg->point_step || msg->is_bigendian ||
            static_cast<std::uint64_t>(msg->width)*msg->point_step!=msg->row_step ||
            msg->data.size()!=msg->row_step) { ++cloud_drops; return; }
        for (const auto& field : std::vector<std::pair<std::string,std::uint8_t>>{
                {"x",7},{"y",7},{"z",7},{"intensity",7},{"tag",2},{"line",2},{"timestamp",8}}) {
            auto it=std::find_if(msg->fields.begin(),msg->fields.end(),[&](const auto& f){return f.name==field.first;});
            const std::uint32_t size=field.second==8 ? 8 : (field.second==7 ? 4 : 1);
            if (it==msg->fields.end() || it->datatype!=field.second || it->count!=1 ||
                static_cast<std::uint64_t>(it->offset)+size>msg->point_step) { ++cloud_drops; return; }
        }
        // 阶段 2：预处理并检查点值。
        CloudPtr cloud(new PointCloudType);
        try { preprocessors.at(id)->Process(msg,cloud); }
        catch (const std::exception& e) { ++cloud_drops; LOG(WARNING) << "INS cloud input: " << e.what(); return; }
        if (cloud->empty() || std::any_of(cloud->begin(),cloud->end(),[](const auto& p){
            return !p.getVector3fMap().allFinite() || !std::isfinite(p.time) || p.time<0;
        })) { ++cloud_drops; return; }
        // 阶段 3：组帧入队并尝试发布。
        // Individual secondary clouds have not yet been rotated into primary axes.
        if (!assembler.AddCloud(id,ToSec(msg->header.stamp),cloud)) { ++cloud_drops; return; }
        FusedLidarFrame ready;
        while (assembler.PopReady(ready)) {
            FilterSelfPoints(*ready.cloud,self_filter,rear_from_primary.rotationMatrix());
            pending.push_back({std::move(ready),Steady()});
        }
        while (pending.size()>8) { pending.pop_front(); ++cloud_drops; }
        DrainClouds();
    }

    /**
     * @brief 按队首顺序等待导航覆盖，丢弃不可用帧并发布可去畸变的点云。
     * @par 主要流程
     * 1. <b>复核导航质量</b>：质量中断时重置恢复门控和插值历史。
     * 2. <b>检查队首可处理性</b>：丢弃等待超时帧；门控或扫描末端覆盖不足时保留队列并退出循环。
     * 3. <b>去畸变与降采样</b>：拒绝缺雷达或无法插值的帧，成功后按配置降采样。
     * 4. <b>发布并移出队首</b>：发布后轴点云，按帧间隔发布地图点云，再移除本帧。
     * @note 阶段名称对应函数体；阶段 2–4 循环执行，此函数不会阻塞等待导航数据到达。
     */
    void DrainClouds() {
        // 阶段 1：复核导航质量。
        const auto bad=ins::CheckSample(sample,policy,node->now().seconds(),Steady());
        if (!bad.empty()) Reject(bad);
        while (!pending.empty()) {
            // 阶段 2：检查队首可处理性。
            auto& p=pending.front();
            const auto& stats=p.frame.stats;
            if (Steady()-p.arrival>cloud_wait) { ++cloud_drops; pending.pop_front(); continue; }
            if (!gate->Good() || poses->Latest()<stats.end_time) break;
            // 阶段 3：去畸变与降采样。
            CloudPtr deskewed; SE3 map_from_rear;
            if (stats.partial || !ins::Deskew(p.frame.cloud,stats.begin_time,stats.end_time,
                                            rear_from_primary,*poses,deskewed,map_from_rear)) {
                ++cloud_drops; pending.pop_front(); continue;
            }
            if (cloud_leaf>0) deskewed=DownsamplePreservingSource(deskewed,cloud_leaf);
            // 阶段 4：发布并移出队首。
            inv_pub->publish(sany_output::MakeCloudMessage(deskewed,stats.begin_time,stats.end_time,SE3(),rear_frame));
            if (cloud_count++ % map_every==0)
                map_pub->publish(sany_output::MakeCloudMessage(deskewed,stats.begin_time,stats.end_time,map_from_rear,frame));
            pending.pop_front();
        }
    }

    /**
     * @brief 周期检查导航与审计状态，发布健康信息并推进等待中的点云。
     * @par 主要流程
     * 1. <b>检查有效性并发布定位状态</b>：复核导航时效，刷新审计流并锁存写入故障。
     * 2. <b>发布故障与诊断快照</b>：根据恢复门控和等待队列统计构造故障及诊断消息。
     * 3. <b>发布轨迹并推进点云</b>：门控通过且路径非空时发布路径，再尝试处理点云队列。
     * @note 阶段名称对应函数体；检查发现故障时会切断插值历史。
     */
    void Health() {
        // 阶段 1：检查有效性并发布定位状态。
        // 日志故障在 flush 后检查并锁存；这不是每次 publish 前的写盘事务。
        // [ins-health-status]
        const auto bad=ins::CheckSample(sample,policy,node->now().seconds(),Steady());
        if (!bad.empty()) Reject(bad);
        raw_log.flush(); solutions.flush(); tum.flush(); events.flush();
        if (!raw_log || !solutions || !tum || !events) { recording_failed=true; Reject("recording_failed"); }
        lightning::msg::LocalizationStatus s;
        s.header.stamp=node->now();
        s.status=gate->Good() ? s.STATUS_GOOD : (ever_good ? s.STATUS_FAIL : s.STATUS_INITIALIZING);
        status_pub->publish(s);
        // [ins-health-status]
        // 阶段 2：发布故障与诊断快照。
        lightning::msg::FaultStatus f;
        f.header=s.header; f.comm_header.source_id=heartbeat->NodeId();
        f.comm_header.seq=functional_safety::HeartbeatPublisher::Advance(fault_seq);
        f.comm_header.stamp_us=functional_safety::HeartbeatPublisher::UnixMicrosecondsNow();
        f.level=gate->Good() ? f.LEVEL_NO_FAULT : f.LEVEL_P0;
        f.fault_type=gate->Good() ? 0 : 2; f.description=gate->Good() ? "" : "ins_only: "+reason;
        fault_pub->publish(f);
        std_msgs::msg::String diagnostic;
        diagnostic.data="mode=ins_only reason="+reason+" cloud_drops="+std::to_string(cloud_drops)+
                        " incomplete_scan_drops="+std::to_string(assembler.InsufficientLidarDropCount())+
                        " pending_clouds="+std::to_string(pending.size());
        diagnostic_pub->publish(diagnostic);
        // 阶段 3：发布轨迹并推进点云。
        if (gate->Good() && !path.poses.empty()) path_pub->publish(path);
        DrainClouds();
    }

    /**
     * @brief 按开关配置纯 INS 模式的雷达组帧与投影输入。
     * @par 主要流程
     * 1. <b>检查开关与投影参数</b>：未启用时返回，否则验证已确认的外参和发布策略。
     * 2. <b>载入雷达配置并初始化组帧</b>：保存配置快照，检查多雷达、自反射及 Livox 预处理参数。
     * 3. <b>创建预处理器与点云订阅</b>：为每台雷达配置独立预处理器，注册带雷达 ID 的回调。
     * @note 阶段名称对应函数体；配置异常交由外层 Init() 捕获。
     */
    void InitClouds(const YAML::Node& config,const std::string& yaml) {
        // 阶段 1：检查开关与投影参数。
        const auto cloud=config["cloud"];
        clouds_enabled=cloud && cloud["enabled"].as<bool>();
        if (!clouds_enabled) return;
        rear_from_primary=ReadExtrinsic(cloud["rear_from_primary"]);
        cloud_wait=Positive(cloud,"max_wait_sec",0.3);
        map_every=cloud["map_publish_every"] ? cloud["map_publish_every"].as<int>() : 1;
        cloud_leaf=cloud["voxel_leaf_m"] ? cloud["voxel_leaf_m"].as<double>() : 0;
        if (map_every<1 || !std::isfinite(cloud_leaf) || cloud_leaf<0) throw std::invalid_argument("invalid cloud output options");
        // 阶段 2：载入雷达配置并初始化组帧。
        const YAML::Node base=cloud["lidar_config"] ? YAML::LoadFile(RelativeTo(yaml,cloud["lidar_config"].as<std::string>()).string()) : config;
        std::ofstream snapshot(audit_directory/"lidar_config.yaml");
        snapshot << base;
        if (!snapshot) throw std::runtime_error("cannot save lidar configuration snapshot");
        std::string error;
        if (!LoadMultiLidarConfig(base,lidar_config,&error) || !lidar_config.enabled ||
            !LoadSelfPointFilterConfig(base,self_filter,&error)) throw std::invalid_argument("multi-lidar config: "+error);
        assembler.Reset(lidar_config);
        const auto c=base["fasterlio"];
        if (c["lidar_type"].as<int>()!=1 || c["point_filter_num"].as<int>()<1 ||
            c["scan_line"].as<int>()<1 || c["scan_line"].as<int>()>256 ||
            !std::isfinite(c["blind"].as<double>()) || c["blind"].as<double>()<0)
            throw std::invalid_argument("ins_only cloud input requires valid Livox PointCloud2 settings");
        const double point_time_scale=Positive(c,"livox_point_time_scale",1);
        // 阶段 3：创建预处理器与点云订阅。
        for (const auto& lidar:lidar_config.lidars) {
            auto pre=std::make_shared<PointCloudPreprocess>();
            pre->Set(static_cast<LidarType>(c["lidar_type"].as<int>()),c["blind"].as<double>(),c["point_filter_num"].as<int>());
            pre->NumScans()=c["scan_line"].as<int>(); pre->TimeScale()=c["time_scale"].as<double>();
            pre->LivoxPointTimeScale()=point_time_scale;
            if (base["roi"]) pre->SetHeightROI(base["roi"]["height_max"].as<float>(),base["roi"]["height_min"].as<float>());
            preprocessors[lidar.id]=pre;
            subscriptions.push_back(node->create_subscription<sensor_msgs::msg::PointCloud2>(lidar.lidar_topic,
                rclcpp::SensorDataQoS().keep_last(4),[this,id=lidar.id](const sensor_msgs::msg::PointCloud2::SharedPtr msg){Cloud(id,msg);}));
        }
    }
};

InsLocSystem::InsLocSystem() : impl_(std::make_unique<Impl>()) {}
InsLocSystem::~InsLocSystem() = default;
bool InsLocSystem::Init(const std::string& yaml, const std::string& published_tum) {
    try {
        // 阶段 1：校验坐标参考。
        auto& x=*impl_;
        const auto config=YAML::LoadFile(yaml), c=config["ins_only"], geo=c["georeference"];
        // [ins-fixed-reference]
        if (config["output"] && config["output"]["fixed_map_transform"] &&
            config["output"]["fixed_map_transform"]["enabled"].as<bool>())
            throw std::invalid_argument("ins_only already outputs fixed ENU; disable legacy output.fixed_map_transform");
        if (!c || !geo["confirmed"].as<bool>()) throw std::invalid_argument("ins_only.georeference must be confirmed");
        if (geo["datum"].as<std::string>()!="WGS84") throw std::invalid_argument("this version requires verified WGS84 input; no implicit datum conversion");
        if (c["output_reference"].as<std::string>()!="rear_axle" || c["attitude_reference"].as<std::string>()!="rear_body")
            throw std::invalid_argument("CGI output must be configured for rear axle and rear body");
        const auto origin=geo["origin_llh"].as<std::vector<double>>();
        if (origin.size()!=3) throw std::invalid_argument("origin_llh requires latitude, longitude, ellipsoid height");
        x.reference=std::make_unique<ins::GeoReference>(origin[0],origin[1],origin[2]);
        // [ins-fixed-reference]
        const auto height=geo["input_height"].as<std::string>();
        if (height=="orthometric") x.geoid_separation=geo["geoid_separation_m"].as<double>();
        else if (height!="ellipsoidal") throw std::invalid_argument("input_height must be verified ellipsoidal or orthometric");
        if (!std::isfinite(x.geoid_separation)) throw std::invalid_argument("invalid geoid separation");
        // 阶段 2：配置质量门控。
        x.max_radius=Positive(geo,"max_radius_m",20000);
        x.offset=c["receive_time_offset_sec"] ? c["receive_time_offset_sec"].as<double>() : 0;
        if (!std::isfinite(x.offset)) throw std::invalid_argument("invalid time offset");
        const auto quality=c["quality"];
        x.policy.max_age=Positive(quality,"max_age_sec",0.15);
        x.policy.max_skew=Positive(quality,"max_skew_sec",0.02);
        x.policy.max_future=Positive(quality,"max_future_sec",0.02);
        x.policy.max_differential_age=Positive(quality,"max_differential_age_sec",2);
        x.policy.max_horizontal_sigma=Positive(quality,"max_horizontal_sigma_m",0.1);
        x.policy.max_vertical_sigma=Positive(quality,"max_vertical_sigma_m",0.2);
        x.policy.max_attitude_sigma=Positive(quality,"max_attitude_sigma_deg",0.5);
        x.policy.max_velocity_sigma=Positive(quality,"max_velocity_sigma_mps",0.2);
        x.policy.recovery_samples=quality["recovery_samples"] ? quality["recovery_samples"].as<int>() : 10;
        if (x.policy.recovery_samples<1 || x.policy.max_skew>x.policy.max_age) throw std::invalid_argument("invalid recovery/skew policy");
        x.gate=std::make_unique<ins::RecoveryGate>(x.policy.recovery_samples);
        x.poses=std::make_unique<ins::PoseBuffer>(2,Positive(c,"pose_interpolation_max_gap_sec",0.05));
        // 阶段 3：创建节点与心跳。
        x.node=std::make_shared<rclcpp::Node>("lightning_slam");
        if (c["use_sim_time"]) x.node->set_parameter(rclcpp::Parameter("use_sim_time",c["use_sim_time"].as<bool>()));
        x.frame=config["output"] && config["output"]["map_frame"] ? config["output"]["map_frame"].as<std::string>() : "map";
        if (x.frame.empty() || x.frame==x.rear_frame) throw std::invalid_argument("invalid map frame");
        functional_safety::HeartbeatConfig hc;
        hc.node_id=6; hc.topic="/diagnostics/heartbeat/lightning_slam";
        x.heartbeat=std::make_unique<functional_safety::HeartbeatPublisher>(*x.node,hc);
        // 阶段 4：打开审计文件。
        const std::filesystem::path dir=std::filesystem::path(c["audit_root"] ? c["audit_root"].as<std::string>() : "runs")/
            ("ins_only_"+std::to_string(functional_safety::HeartbeatPublisher::UnixMicrosecondsNow()));
        std::filesystem::create_directories(dir);
        x.audit_directory=dir;
        std::ofstream snapshot(dir/"config.yaml"); snapshot << config;
        if (!snapshot) throw std::runtime_error("cannot save configuration snapshot");
        x.raw_log.open(dir/"can_frames.csv"); x.solutions.open(dir/"solutions.csv"); x.events.open(dir/"events.csv");
        const auto tp=published_tum.empty() ? dir/"published_rear_axle.tum" : std::filesystem::path(published_tum);
        if (std::filesystem::exists(tp)) throw std::invalid_argument("refusing to overwrite published TUM");
        if (!tp.parent_path().empty()) std::filesystem::create_directories(tp.parent_path());
        x.tum.open(tp);
        if (!x.raw_log || !x.solutions || !x.events || !x.tum) throw std::runtime_error("cannot open audit outputs");
        x.raw_log << "raw_stamp,corrected_stamp,field,can_id,valid,dlc,b0,b1,b2,b3,b4,b5,b6,b7\n";
        x.solutions << "stamp,latitude,longitude,altitude,heading,pitch,roll,system_state,satellite_status,differential_age,"
                    "ve,vn,vu,sigma_e,sigma_n,sigma_u,sigma_heading,sigma_pitch,sigma_roll,sigma_ve,sigma_vn,sigma_vu,"
                    "t_lat,t_lon,t_alt,t_att,t_vel,t_pos_sigma,t_att_sigma,t_status,t_vel_sigma,"
                    "x,y,z,qx,qy,qz,qw,vx,vy,vz,signed_speed,published,reason\n";
        x.events << "stamp,reason\n";
        // 阶段 5：创建输出通道。
        const auto qos=rclcpp::QoS(1000).best_effort();
        x.pose_pub=x.node->create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic",qos);
        x.vehicle_pub=x.node->create_publisher<lightning::msg::VehiclePose>("/localization/pose_vel",qos);
        x.status_pub=x.node->create_publisher<lightning::msg::LocalizationStatus>("/localization/loc_status",rclcpp::QoS(10).best_effort());
        x.fault_pub=x.node->create_publisher<lightning::msg::FaultStatus>("/localization/fault_status",rclcpp::QoS(10).best_effort());
        x.path_pub=x.node->create_publisher<nav_msgs::msg::Path>("/localization/path",rclcpp::QoS(1).best_effort());
        x.diagnostic_pub=x.node->create_publisher<std_msgs::msg::String>("/localization/ins_diagnostics",rclcpp::QoS(10).best_effort());
        x.inv_pub=x.node->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInv",rclcpp::QoS(1).best_effort());
        x.map_pub=x.node->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInL",rclcpp::QoS(1).best_effort());
        if (config["system"]["pub_tf"] && config["system"]["pub_tf"].as<bool>()) x.tf=std::make_unique<tf2_ros::TransformBroadcaster>(x.node);
        // 阶段 6：注册导航与点云输入。
        const auto ns=c["topic_prefix"] ? c["topic_prefix"].as<std::string>() : "/cgi430";
        if (ns.empty() || ns.front()!='/' || ns.back()=='/') throw std::invalid_argument("invalid CGI topic_prefix");
        namespace ci=cgi430_interfaces::msg;
        x.Subscribe<ci::Latitude>(ns+"/position/latitude",ins::Latitude,[&x](const auto& m){x.sample.latitude=m.latitude_deg;});
        x.Subscribe<ci::Longitude>(ns+"/position/longitude",ins::Longitude,[&x](const auto& m){x.sample.longitude=m.longitude_deg;});
        x.Subscribe<ci::Altitude>(ns+"/position/altitude",ins::Altitude,[&x](const auto& m){x.sample.altitude=m.altitude_m;});
        x.Subscribe<ci::Attitude>(ns+"/attitude",ins::Attitude,[&x](const auto& m){x.sample.heading=m.heading_deg;x.sample.pitch=m.pitch_deg;x.sample.roll=m.roll_deg;});
        x.Subscribe<ci::EarthVelocity>(ns+"/velocity",ins::Velocity,[&x](const auto& m){x.sample.velocity=Vec3d(m.east_mps,m.north_mps,m.up_mps);});
        x.Subscribe<ci::PositionSigma>(ns+"/position/sigma",ins::PositionSigma,[&x](const auto& m){x.sample.position_sigma=Vec3d(m.east_m,m.north_m,m.up_m);});
        x.Subscribe<ci::AttitudeSigma>(ns+"/attitude/sigma",ins::AttitudeSigma,[&x](const auto& m){x.sample.attitude_sigma=Vec3d(m.heading_deg,m.pitch_deg,m.roll_deg);});
        x.Subscribe<ci::VelocitySigma>(ns+"/velocity/sigma",ins::VelocitySigma,[&x](const auto& m){x.sample.velocity_sigma=Vec3d(m.east_mps,m.north_mps,m.up_mps);});
        x.Subscribe<ci::InsStatus>(ns+"/ins/status",ins::Status,[&x](const auto& m){x.sample.system_state=m.system_state;x.sample.satellite_status=m.satellite_status;x.sample.differential_age=m.differential_age_s;});
        x.InitClouds(config,yaml);
        // 阶段 7：启动健康检查。
        x.timer=x.node->create_wall_timer(std::chrono::milliseconds(50),[&x]{x.Health();});
        LOG(INFO) << "INS_ONLY initialized, fixed ENU origin=" << origin[0] << ',' << origin[1]
                  << " audit=" << dir << " source=CGI_rear_axle no_LIO_no_map_matching";
        return true;
    } catch (const std::exception& e) { LOG(ERROR) << "ins_only initialization: " << e.what(); return false; }
}
void InsLocSystem::Spin() { rclcpp::spin(impl_->node); }
}  // namespace lightning

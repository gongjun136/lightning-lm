#pragma once

#include <pcl/registration/icp.h>
#include <chrono>
#include <deque>
#include <future>
#include <fstream>
#include <iostream>
#include <pcl/kdtree/kdtree_flann.h>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <thread>

#include "common/nav_state.h"
#include "common/timed_pose.h"
#include "core/localization/global_relocalizer.h"
#include "core/localization/localization_result.h"
#include "core/localization/point_to_plane_registration.h"
#include "core/maps/tiled_map.h"
#include "utils/compute_profiling.h"

#include "pclomp/ndt_omp_impl.hpp"

namespace lightning::ui {
class PangolinWindow;
}
namespace lightning::diagnostics {
class MatchFailureRecorder;
}

namespace lightning::loc {

/// 激光定位对外接口类
class LidarLoc {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using UL = std::unique_lock<std::mutex>;

    /// 定位方法
    enum class LocMethod {
        NDT_OMP,  // OMP版本NDT
    };

    struct Options {
        Options() {}

        bool display_realtime_cloud_ = false;          ///< 是否显示实时点云
        bool debug_ = false;                           ///< 是否使用单测模式
        LocMethod match_method_ = LocMethod::NDT_OMP;  ///< 匹配方式
        bool try_self_extrap_ = false;                 ///< 是否尝试自己的外推pose
        bool with_height_ = true;                      ///< 建图期间是否带有高度约束？
        bool force_2d_ = true;                         ///< 强制在2D空间
        float min_init_confidence_ = 0.1;              ///< 初始化时要求的最小分值
        float min_tracking_confidence_ = 0.1;  ///< 正常跟踪接受匹配结果的最低置信度。
        int relocalization_lost_frame_threshold_ = 5;  ///< 触发失锁状态的连续匹配失败帧数。
        bool init_with_fp_ = true;                     ///< 是否使用功能点进行初始化
        bool enable_parking_static_ = false;           ///< 是否在静止时输出固定位置
        bool enable_icp_adjust_ = false;               ///< 是否使用icp调整ndt匹配结果提高定位精度

        /// 点云过滤
        // float filter_z_min_ = -1.0;
        float filter_z_max_ = 30.0;
        // float filter_intensity_min_ = 10.0;
        // float filter_intensity_max_ = 100.0;

        /// 地图配置
        TiledMap::Options map_option_ = TiledMap::Options();

        // 动态图层更新相关
        bool update_dynamic_cloud_ = true;     ///< 是否使用定位后的点云来更新动态图层
        double update_kf_dis_ = 5.0;           ///< 每隔多少米更新一次
        double update_kf_time_ = 10.0;         ///< 每隔多少时间更新一次
        double update_lidar_loc_score_ = 2.2;  ///< 更新时激光定位匹配分值阈值
        double lidar_loc_odom_th_ = 0.3;       ///< 激光两帧匹配结果与对应lidarodom结果差的阈值，超过则认为lidarodom异常

        double max_update_cache_dis_ = 30.0;  ///< 更新动态图层的缓冲距离
        std::string recover_pose_path_ = "./data/recover_pose.txt";  ///< 定位恢复位姿记录文件的路径。
        bool enable_relocalization_map_consistency_ = true;  ///< 是否对重定位候选执行地图一致性验证。
        double relocalization_bounds_margin_ = 0.5;  ///< 地图空间包围范围的外扩余量（m）。
        double relocalization_nearest_neighbor_distance_ = 0.5;  ///< 计算地图重合率使用的最大近邻距离（m）。
        double relocalization_min_inside_xy_ratio_ = 0.90;  ///< 接受候选所需的最小水平范围内点比例。
        double relocalization_min_overlap_ratio_ = 0.20;  ///< 接受候选所需的最小静态地图重合率。
        double relocalization_precheck_min_overlap_ratio_ = 0.20;  ///< 候选预检查的最小地图重合率。
        double relocalization_min_gravity_alignment_cos_ = 0.95;  ///< 接受候选所需的最小重力方向一致性余弦。
        int relocalization_map_consistency_max_points_ = 50000;  ///< 精配后地图一致性检查的采样点数上限。
        int relocalization_precheck_max_points_ = 50000;  ///< 候选预检查的采样点数上限。
        int relocalization_validation_workers_ = 2;  ///< 候选验证使用的并行工作线程数。
        int relocalization_max_refinement_candidates_ = 0;  ///< 一次重定位允许精配的候选数配置。
        int relocalization_confirmation_count_ = 2;  ///< 接受候选所需的连续确认次数。
        double relocalization_confirmation_max_translation_ = 1.0;  ///< 连续确认之间允许的最大平移差（m）。
        double relocalization_confirmation_max_rotation_deg_ = 5.0;  ///< 连续确认之间允许的最大旋转差（deg）。
        double relocalization_confirmation_max_interval_ = 1.0;  ///< 保留的确认间隔配置（s），当前仅读取并校验，不作为候选接受门控。
        std::string relocalization_refinement_backend_ = "ndt";  ///< 候选精配后端名称。
        PointToPlaneRegistration::Options relocalization_plane_icp_options_;  ///< 点面 ICP 候选精配参数。
        std::string relocalization_debug_dir_;  ///< 候选验证的调试产物输出目录。
    };

    struct MatchStats {
        double confidence = 0.0;  ///< 本次匹配得到的置信度。
        int iterations = 0;  ///< 本次匹配使用的迭代次数。
        bool success = false;  ///< 本次地图定位是否成功。
        int active_map_chunks = 0;  ///< 本次匹配加载的地图块数量。
        bool relocalization_attempted = false;  ///< 本帧是否尝试全局重定位。
        bool relocalization_candidate_found = false;  ///< 本次重定位是否检索到候选。
        bool relocalization_accepted = false;  ///< 本次重定位是否接受候选。
        int relocalization_candidate_id = -1;  ///< 本次重定位候选的数据库 ID。
        double relocalization_score = 0.0;  ///< 本次重定位候选检索分值。
        int relocalization_candidate_count = 0;  ///< 本次检索返回的候选数量。
        int relocalization_candidates_prechecked = 0;  ///< 本次完成预检查的候选数量。
        int relocalization_confirmation_count = 0;  ///< 当前候选连续确认的次数。
        int relocalization_query_submap_size = 0;  ///< 本次查询使用的子地图帧数。
        double relocalization_search_time_ms = 0.0;  ///< 本次候选检索耗时（ms）。
        std::string relocalization_reason;  ///< 本次重定位的接受或拒绝原因。
        bool map_consistency_evaluated = false;  ///< 是否执行地图一致性验证。
        bool map_consistency_passed = false;  ///< 是否通过地图一致性验证。
        std::size_t map_consistency_points = 0;  ///< 地图一致性验证使用的点数。
        double map_inside_xy_ratio = 0.0;  ///< 候选点落在地图水平边界内的比例。
        double map_inside_xyz_ratio = 0.0;  ///< 候选点落在地图三维边界内的比例。
        double map_overlap_ratio = 0.0;  ///< 候选点与静态地图近邻重合的比例。
        double map_gravity_alignment_cos = 0.0;  ///< 候选重力方向与地图竖直方向的一致性余弦。
    };

    explicit LidarLoc(Options options = Options());
    virtual ~LidarLoc();

    /**
     * @brief 读取地图定位配置，准备地图与重定位后端，最后启动地图更新线程。
     * @return 必需配置或地图准备失败返回 false；后端不可用是否致命取决于地图归一化要求。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>读取预算并检查地图元数据</b>：读取 YAML，验证计算预算；启用地图高度归一化时要求地图元数据可读。
     * 2. <b>读取跟踪与重定位门控</b>：设置地图加载范围、匹配阈值与确认条件；后端专属配置可覆盖共用重定位参数。
     * 3. <b>读取候选细化参数</b>：装配平面 ICP 参数及后端专属覆盖；此处读取配置，不执行配准。
     * 4. <b>校验参数并补齐定位选项</b>：验证细化和地图一致性参数，再读取功能点、静止处理、ICP 微调等定位开关。
     * 5. <b>加载地图并准备匹配目标</b>：设置动态图层策略，加载地图索引与可选验证缓存；有功能点时加载首个功能点附近的地图。
     * 6. <b>初始化全局重定位后端</b>：按配置创建 BTC 或 SOLiD；后端初始化失败时，归一化地图要求返回 false，否则保留 NDT 定位能力。
     * 7. <b>加入恢复位姿并启动地图线程</b>：存在恢复文件时将其读作功能点，再启动地图更新线程；恢复位姿尚未通过本帧匹配验证。
     */
    bool Init(const std::string& config_path);

    /// 处理 lidar odom 消息
    bool ProcessLO(const NavState& state);

    /**
     * @brief 对定位点云做基础检查，运行初始化或跟踪并记录本帧耗时。
     * @return 点数不足返回 false；返回 true 仅表示已调用 Align()，定位有效性应另查结果状态。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查定位点云</b>：要求指针非空，并拒绝少于 50 个点的输入。
     * 2. <b>准备帧统计并运行匹配</b>：重置帧统计、保存输入点云，交给 Align() 处理初始化、静止或跟踪分支。
     * 3. <b>记录耗时并结束处理</b>：启用统计时记录匹配结果及各步骤耗时；函数随后返回 true，不以此替代定位有效标志。
     */
    bool ProcessCloud(CloudPtr cloud_input);

    /// 处理DR状态
    bool ProcessDR(const NavState& state);

    /// 获取位姿
    NavState GetState();

    /// 设置当前定位状态标志位
    void SetInitRltState();

    /**
     * 尝试在另一个位置进行定位
     * @param input
     * @param pose
     * @return
     */
    bool TryOtherSolution(CloudPtr input, SE3& pose);

    /// 使用功能点初始化
    bool InitWithFP(CloudPtr input, const SE3& fp_pose);

    /// 更新全局地图
    bool UpdateGlobalMap();

    /// @brief 定位是否已经成功初始化
    bool LocInited();

    /// @brief 激光定位重置接口
    void ResetLastPose(const SE3& last_pose);

    /**
     * @brief 以地图系预测位姿启动局部配准，并检查收敛、数值和分数。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>准备目标与初值</b>：保存输入初值并清空匹配统计；目标地图不可用时直接返回 false。
     * 2. <b>执行 NDT 并检查基础结果</b>：持有匹配锁，选择分辨率和点数预算；用收敛、有限数值与分数阈值判定基础成功。
     * 3. <b>执行可选 ICP 微调</b>：只在已初始化且 NDT 成功时执行；ICP 收敛且平移、转角修正分别不超过 0.05 m、0.05 rad 才采用。
     * 4. <b>写回候选位姿与统计</b>：数值有效时写回候选位姿并记录结果；返回 false 不保证 pose 保持原值，也不表示已通过 Align() 的全部门控。
     * @param pose 输入地图系初值；输出配准后的位姿，返回 false 时也可能被修改。
     * @param confidence 输出 NDT 分数。
     * @param input 输入点云
     * @param output 输出点云
     * @param use_rough_res 是否使用粗分辨率
     * @return 局部配准通过基础检查返回 true；Align() 仍须检查里程计增量和可选地图一致性。
     * @see @ref loc_read_acceptance "局部配准通过后还要验证什么"
     */
    bool Localize(SE3& pose, double& confidence, CloudPtr input, CloudPtr output, bool use_rough_res = false);

    /// 设置UI
    void SetUI(std::shared_ptr<ui::PangolinWindow> ui) { ui_ = ui; }

    /// 设置init pose
    void SetInitialPose(SE3 init_pose);

    /// Clear the current map alignment and restart global BTC initialization.
    /// This is also useful when an external health monitor detects localization
    /// failure before the internal consecutive-match threshold is reached.
    void RequestGlobalRelocalization();

    /// 获取定位结果
    LocalizationResult GetLocalizationResult() {
        UL lock(result_mutex_);
        return localization_result_;
    }

    MatchStats GetLastMatchStats() const;

    void Finish();

    /// 激光定位是否认为LO有效
    bool LidarLocThinkLOReliable() { return lo_reliable_; }

   private:
    friend struct LidarLocRegressionAccess;
    // 内部函数  ==========================================================================
    /**
     * @brief 按初始化或跟踪状态处理点云，决定是否接受地图位姿或退回 LO 外推。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>对齐帧时间与里程计</b>：使用去畸变结束时刻，分别为本帧关联 LO 和 DR 位姿。
     * 2. <b>处理静止保持</b>：已初始化且 parking_ 为真时沿用上一地图位姿，更新结果后直接返回。
     * 3. <b>执行首次定位或重定位</b>：未初始化时尝试外部初值、可选功能点和全局候选；该分支结束本帧，不继续跟踪流程。
     * 4. <b>构造跟踪初值</b>：将 LO 相对运动接到上次地图位姿上；自身外推只用于后续平滑性比较。
     * 5. <b>配准并限制单帧修正</b>：加载附近地图并调用 Localize()，保留原始解，再将修正缩放为 0.1 并应用可选二维限制。
     * 6. <b>验证匹配并决定恢复</b>：先用原始配准解检查 LO 增量，再做可选地图一致性验证；失败退回 LO 初值，连续失败可触发重定位。
     * 7. <b>更新位姿与定位状态</b>：成对更新地图位姿与里程计时刻，写入 GOOD、FOLLOWING_DR 或 INITIALIZING 及相关状态。
     * 8. <b>更新动态图层与恢复记录</b>：按成功、距离或时间、分数条件更新动态图层，维护历史，并在文件可打开时保存恢复位姿。
     * @param input 已通过 ProcessCloud() 基础检查的定位点云。
     * @see @ref loc_read_prior "地图匹配初值", @ref loc_read_recovery "拒绝匹配后的状态"
     */
    void Align(const CloudPtr& input);

    /**
     * 寻找当前帧对应的LO相对位姿
     * @param timestamp
     * @return
     */
    bool AssignLOPose(double timestamp);

    /**
     * 寻找当前帧对应的DR相对位姿
     * @param timestamp
     * @return
     */
    bool AssignDRPose(double timestamp);

    /**
     * 检查车辆是否静止
     * @param timestamp
     * @return
     */
    bool CheckStatic(double timestamp);

    /**
     * 更新自身状态
     * @param input
     */
    void UpdateState(const CloudPtr& input);

    /**
     * 更新地图
     */
    void UpdateMapThread();

    /**
     * 使用网格搜索best yaw
     */
    bool YawSearch(SE3& pose, double& confidence, CloudPtr input, CloudPtr output);

    bool CheckLidarOdomValid(const SE3& current_pose_esti, double& delta_posi);
    /**
     * @brief 细化全局候选并做跨帧确认，通过后提交地图位姿和 GOOD 状态。
     * @details 候选循环内执行阶段 4–6；未通过细化或地图检查时继续下一候选。
     * @return 仅最终接受返回 true；无候选、验证失败或仍待确认均返回 false。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>取得待确认或新检索候选</b>：要求后端就绪且有当前 LO；后续帧优先沿用待确认的 map→odom 关系，否则向后端累积并检索。
     * 2. <b>检查检索结果并建立统计</b>：无返回结果时已结束；有结果后记录检索状态，未被后端接受或候选为空时返回 false。
     * 3. <b>预检地图并限制细化数量</b>：启用检查时并行预检并等待完成；筛选后在细化上限为正时按重叠率和检索分数排序裁剪。
     * 4. <b>逐个细化并验证地图</b>：为候选加载目标，执行 NDT 或平面 ICP，再验证地图一致性；失败继续下一候选。
     * 5. <b>累计跨帧一致确认</b>：比较 map→odom 关系，更新连续确认计数；确认不足时保存候选并返回 false。
     * 6. <b>提交接受的地图位姿</b>：确认达标后写入位姿、GOOD 状态与配对里程计，清空待确认状态并返回 true；PGO 重置由上层处理。
     * 7. <b>记录候选全部失败</b>：仅候选列表耗尽时到达；记录失败，必要时丢弃复用的待确认候选。
     * @see @ref loc_read_confirmation "全局候选的连续确认"
     */
    bool TryGlobalRelocalization(const CloudPtr& input);
    struct MapConsistencyResult {
        bool evaluated = false;  ///< 是否执行了地图一致性检查。
        bool passed = false;  ///< 候选是否通过地图一致性检查。
        std::size_t points = 0;  ///< 一致性检查实际使用的点数。
        double inside_xy_ratio = 0.0;  ///< 落在地图水平包围范围内的点比例。
        double inside_xyz_ratio = 0.0;  ///< 落在地图三维包围范围内的点比例。
        double overlap_ratio = 0.0;  ///< 在配置距离内找到静态地图近邻的点比例。
        double gravity_alignment_cos = 0.0;  ///< 候选重力方向与地图竖直方向的一致性余弦。
    };
    bool BuildRelocalizationMapCache();
    MapConsistencyResult EvaluateRelocalizationMapConsistency(
        const CloudPtr& input, const SE3& pose, std::size_t worker_index,
        double min_overlap_ratio, std::size_t maximum_points) const;
    void ApplyMapConsistencyResult(const MapConsistencyResult& result);
    bool ValidateRelocalizationMapConsistency(const CloudPtr& input, const SE3& pose);
    void SaveRelocalizationBirdseye(const CloudPtr& static_map, const CloudPtr& scan_world,
                                    const Vec3d& map_min, const Vec3d& map_max,
                                    const MatchStats& stats);

    struct FrameProfiling {
        profiling::TimingSample assign_pose;  ///< 本帧分派 LIO/DR 位姿的耗时。
        profiling::TimingSample map_load;  ///< 本帧加载活动地图的耗时。
        profiling::TimingSample ndt_align;  ///< 本帧 NDT 匹配的耗时。
        profiling::TimingSample icp_adjust;  ///< 本帧 ICP 调整的耗时。
        profiling::TimingSample state_update;  ///< 本帧更新定位状态的耗时。
        profiling::TimingSample dynamic_map;  ///< 本帧更新动态图层的耗时。
        profiling::TimingSample recover_pose_io;  ///< 本帧记录恢复位姿的 I/O 耗时。
        profiling::TimingSample relocalization_poll;  ///< 本帧检查重定位任务结果的耗时。
        int ndt_calls = 0;  ///< 本帧调用 NDT 的次数。
        std::size_t ndt_input_points = 0;  ///< 本帧送入 NDT 的点数。
    } frame_profiling_;

    // 成员变量  ==========================================================================
    Options options_;  ///< 跟踪、地图更新与重定位配置。
    int ndt_threads_ = 4;  ///< NDT 匹配使用的并行线程上限。
    std::size_t ndt_max_points_ = 0;  ///< NDT 输入点数上限，0 表示不限制。

    std::mutex match_mutex_;  ///< 锁定pcl_ndt指针

    using NDTType = pclomp::NormalDistributionsTransform<PointType, PointType>;
    NDTType::Ptr pcl_ndt_ = nullptr;  ///< 正常跟踪使用的 NDT 匹配器。
    NDTType::Ptr pcl_ndt_rough_ = nullptr;  ///< 粗分辨率

    using ICPType = pcl::IterativeClosestPoint<PointType, PointType>;
    ICPType::Ptr pcl_icp_ = nullptr;  ///< 可选的匹配结果 ICP 调整器。

    CloudPtr current_scan_ = nullptr;                   ///< 当前扫描
    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  ///< 可选的定位点云与匹配状态显示窗口。
    SE3 last_loc_pose_;  ///< 最近一次保存的地图定位位姿。

    std::mutex initial_pose_mutex_;  ///< 初始定位锁
    bool initial_pose_set_ = false;  ///< 定位是否被手动设置
    SE3 initial_pose_;               ///< 手动设置的初始位姿
    bool loc_inited_ = false;        ///< 定位是否初始化成功

    double current_timestamp_ = 0;  ///< 本次输入的时间戳
    double last_timestamp_ = 0;     ///< 上次输入的时间戳

    bool last_lo_pose_set_ = false;  ///< 是否已有上一帧 LIO 相对位姿。
    bool current_lo_pose_set_ = false;  ///< 当前帧是否找到对应 LIO 位姿。
    SE3 last_lo_pose_;     ///< 上一次相对位置，相对位置来自LO
    SE3 current_lo_pose_;  ///< 本次的LO相对位置

    bool last_dr_pose_set_ = false;  ///< 是否已有上一帧 DR 相对位姿。
    bool current_dr_pose_set_ = false;  ///< 当前帧是否找到对应 DR 位姿。
    SE3 last_dr_pose_;     ///< 上一次相对位置，相对位置来自DR
    SE3 current_dr_pose_;  ///< 本次的DR相对位置
    bool parking_ = false;  ///< 当前地图匹配是否处于静止保持状态。

    double try_other_guess_trans_th_ = 0.3;               ///< 在初始估计相差多少时，尝试其他的解
    double try_other_guess_rot_th_ = 0.5 * M_PI / 180.0;  ///< 在初始估计相差多少时，尝试其他的解
    // double low_vel_th_ = 1.0;                             // 低速阈值m/s
    // double update_cache_dis_ = 0;                         // 动态图层的更新缓冲距离

    std::deque<double> ave_scores_;  ///< 近期匹配的分值情况

    Vec3d current_vel_b_ = Vec3d::Zero();  ///< 本次的车体系下的速度
    Vec3d current_vel_ = Vec3d::Zero();    ///< 本次dr的速度

    std::deque<TimedPose> lidar_loc_pose_queue_;  ///< lidar odom pose 轨迹

    bool lo_reliable_ = true;  ///< 当前 LIO 相对运动是否可信。
    int lo_reliable_cnt_ = 0;  ///< LIO 可靠性判定的累计计数。

    SE3 last_abs_pose_;            ///< 上一次绝对定位，绝对定位来自于地图匹配得到的位姿
    TimedPose last_dyn_upd_pose_;  ///< 上次更新动态图层时使用的位姿
    SE3 current_abs_pose_;         ///< 本次的绝对定位
    bool last_abs_pose_set_ = false;  ///< 是否已有上一帧地图绝对位姿。
    double current_score_ = 1e5;  ///< 最近地图匹配的原始分值。
    MatchStats last_match_stats_;  ///< 最近一帧地图匹配与重定位的诊断汇总。
    int match_fail_count_ = 0;  ///< 连续地图匹配失败的次数。
    int static_count_ = 0;  ///< 静止判定的累计次数。

    int rtk_reset_cnt_ = 0;  ///< RTK重置计数

    std::mutex result_mutex_;  ///< 保护最近输出的定位结果。
    LocalizationResult localization_result_;  ///< 输出结果

    // 相对运动观测队列
    std::mutex lo_pose_mutex_;  ///< 保护 LIO 状态历史队列。
    std::deque<NavState> lo_pose_queue_;  ///< 按时间保存的 LIO 状态，用于匹配对应扫描。

    std::mutex dr_pose_mutex_;  ///< 保护 DR 状态历史队列。
    std::deque<NavState> dr_pose_queue_;  ///< 按时间保存的 DR 状态，用于位姿猜测。

    /// 功能点初始化的记录
    std::vector<SE3> fp_init_fail_pose_vec_;
    double fp_last_tried_time_ = 0;  ///< 最近尝试功能点初始化的时间戳（s）。

    bool update_map_quit_ = false;  ///< 请求地图更新线程退出的标志。
    std::thread update_map_thread_;            ///< 地图更新
    std::shared_ptr<TiledMap> map_ = nullptr;  ///< 按位姿加载活动块并管理动态层的分块地图。
    std::unique_ptr<GlobalRelocalizer> global_relocalizer_;  ///< 当前选择的全局候选检索后端。
    std::string relocalization_backend_name_ = "btc";  ///< 全局检索后端名称，如 BTC 或 SOLiD。
    std::shared_ptr<PointToPlaneRegistration> relocalization_plane_registration_;  ///< 重定位候选的点面精配验证器。
    CloudPtr relocalization_static_map_;  ///< 重定位一致性验证使用的静态地图点云。
    std::unique_ptr<diagnostics::MatchFailureRecorder> match_failure_recorder_;  ///< 有界保存被拒绝的原始匹配证据。
    Vec3d relocalization_map_min_ = Vec3d::Zero();  ///< 静态地图空间包围盒下界（m）。
    Vec3d relocalization_map_max_ = Vec3d::Zero();  ///< 静态地图空间包围盒上界（m）。
    std::vector<std::unique_ptr<pcl::KdTreeFLANN<PointType>>> relocalization_kdtrees_;  ///< 按验证工作线程保存的静态地图最近邻索引。
    struct PendingRelocalization {
        bool valid = false;  ///< 当前待确认的重定位候选是否有效。
        SE3 T_map_odom;  ///< 候选对应的里程计坐标系到地图坐标系变换。
        SE3 pose;  ///< 候选精配后的地图定位位姿。
        int candidate_id = -1;  ///< 待确认候选的数据库 ID。
        int query_submap_size = 0;  ///< 待确认候选使用的查询子地图帧数。
        double score = 0.0;  ///< 待确认候选的检索分值。
        double ndt_confidence = 0.0;  ///< 待确认候选精配后的 NDT 置信度。
        double timestamp = 0.0;  ///< 当前候选对应的传感器时间戳（s）。
        int confirmation_count = 0;  ///< 当前候选已通过的连续确认次数。
    } pending_relocalization_;
    bool relocalization_event_active_ = false;  ///< 是否已发出本轮重定位开始事件。
    double map_height_ = 0;  ///< 地图高度约束使用的参考高度（m）。
    int relocalization_debug_index_ = 0;  ///< 重定位调试产物的编号。

    bool has_set_pose_ = false;  ///< 退出保存的抑制标志；当前仅默认初始化，SetInitialPose() 修改的是 initial_pose_set_。

    std::ofstream recover_pose_out_;  ///< 记录定位恢复位姿的输出文件流。
};

}  // namespace lightning::loc

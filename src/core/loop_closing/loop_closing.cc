//
// Created by xiang on 25-4-21.
//

#include "core/loop_closing/loop_closing.h"
#include "common/keyframe.h"
#include "common/loop_candidate.h"
#include "utils/pointcloud_utils.h"

#include <pcl/common/transforms.h>
#include <pcl/registration/ndt.h>

#include "core/opti_algo/algo_select.h"
#include "core/lightning_math.hpp"
#include "core/robust_kernel/cauchy.h"
#include "core/types/edge_se3.h"
#include "core/types/edge_se3_height_prior.h"
#include "core/types/vertex_se3.h"
#include "io/yaml_io.h"

namespace lightning {

LoopClosing::~LoopClosing() {
    if (options_.online_mode_) {
        kf_thread_.Quit();
    }
}

void LoopClosing::Init(const std::string yaml_path) {
    // 1. 初始化Miao图优化器
    miao::OptimizerConfig config(miao::AlgorithmType::LEVENBERG_MARQUARDT,
                                 miao::LinearSolverType::LINEAR_SOLVER_SPARSE_EIGEN, false);
    config.incremental_mode_ = true;                  // 启用增量优化模式，支持动态添加约束边
    optimizer_ = miao::SetupOptimizer<6, 3>(config);  // 6自由度SE3位姿，3自由度点云

    // 2. 初始化运动约束信息矩阵（相邻关键帧间的约束权重）
    info_motion_.setIdentity();
    // 平移权重：1/(trans_noise^2)，控制位置约束强度
    info_motion_.block<3, 3>(0, 0) =
        Mat3d::Identity() * 1.0 / (options_.motion_trans_noise_ * options_.motion_trans_noise_);
    // 旋转权重：1/(rot_noise^2)，控制姿态约束强度
    info_motion_.block<3, 3>(3, 3) =
        Mat3d::Identity() * 1.0 / (options_.motion_rot_noise_ * options_.motion_rot_noise_);

    // 3. 初始化回环约束信息矩阵（回环检测约束的权重）
    info_loops_.setIdentity();
    // 回环平移权重：通常比运动约束权重小，允许更大的误差
    info_loops_.block<3, 3>(0, 0) = Mat3d::Identity() * 1.0 / (options_.loop_trans_noise_ * options_.loop_trans_noise_);
    // 回环旋转权重：与运动约束相当或略小
    info_loops_.block<3, 3>(3, 3) = Mat3d::Identity() * 1.0 / (options_.loop_rot_noise_ * options_.loop_rot_noise_);

    // 4. 从YAML配置文件加载回环检测参数
    if (!yaml_path.empty()) {
        YAML_IO yaml(yaml_path);

        options_.loop_kf_gap_ = yaml.GetValue<int>("loop_closing", "loop_kf_gap");          // 回环检测间隔
        options_.min_id_interval_ = yaml.GetValue<int>("loop_closing", "min_id_interval");  // 最小ID间隔
        options_.closest_id_th_ = yaml.GetValue<int>("loop_closing", "closest_id_th");      // 最近ID阈值
        options_.max_range_ = yaml.GetValue<double>("loop_closing", "max_range");           // 最大搜索距离
        options_.ndt_score_th_ = yaml.GetValue<double>("loop_closing", "ndt_score_th");     // NDT配准分数阈值
        options_.with_height_ = yaml.GetValue<bool>("loop_closing", "with_height");         // 是否启用高度约束

        auto extrinT = yaml.GetValue<std::vector<double>>("fasterlio", "extrinsic_T");
        auto extrinR = yaml.GetValue<std::vector<double>>("fasterlio", "extrinsic_R");
        T_imu_lidar_ = SE3(Eigen::Quaterniond(math::MatFromArray<double>(extrinR)).normalized(),
                           math::VecFromArray<double>(extrinT));
    }

    // 5. 在线模式下启动异步处理线程（离线模式跳过此步骤）
    if (options_.online_mode_) {
        LOG(INFO) << "loop closing module is running in online mode";
        kf_thread_.SetProcFunc([this](Keyframe::Ptr kf) { HandleKF(kf); });  // 设置关键帧处理回调
        kf_thread_.SetName("handle loop closure");                           // 设置线程名称
        kf_thread_.Start();                                                  // 启动处理线程
    }
}

void LoopClosing::AddKF(Keyframe::Ptr kf) {
    if (options_.online_mode_) {
        // 在线模式下不阻塞LIO前端，将关键帧交给后台线程异步检测回环。
        kf_thread_.AddMessage(kf);
    } else {
        // 离线模式通常按数据顺序回放，直接同步处理便于保证流程确定性。
        HandleKF(kf);
    }
}

void LoopClosing::HandleKF(Keyframe::Ptr kf) {
    // 同一个关键帧可能被外部重复投递，直接跳过避免重复加边和重复优化。
    if (kf == last_kf_) {
        return;
    }

    // 记录当前关键帧，并加入历史序列；后续候选搜索和位姿图优化都依赖这个序列。
    cur_kf_ = kf;
    all_keyframes_.emplace_back(kf);

    // 先按关键帧ID间隔和空间距离，从历史关键帧中筛选可能形成回环的候选。
    DetectLoopCandidates();

    if (options_.verbose_) {
        LOG(INFO) << "lc: get kf " << cur_kf_->GetID() << " candi: " << candidates_.size();
    }

    // 对候选帧做点云匹配，估计当前帧与历史帧之间的回环相对位姿。
    ComputeLoopCandidates();

    // 将里程计约束和有效回环约束加入位姿图，并更新关键帧优化位姿。
    PoseOptimization();

    // 标记本帧已处理，用于下一次调用时去重。
    last_kf_ = kf;
}

void LoopClosing::DetectLoopCandidates() {
    // 每次只针对当前关键帧重新生成候选集合，避免沿用上一帧的候选结果。
    candidates_.clear();

    // all_keyframes_按关键帧创建顺序保存；下面依赖这个顺序做近邻ID剪枝。
    auto& kfs_mapping = all_keyframes_;
    Keyframe::Ptr check_first = nullptr;  // 记录第一个候选回环帧，用于跳过连续相近的候选帧

    // 第一帧没有历史轨迹可回环，只记录为上一次尝试回环的关键帧。
    if (last_loop_kf_ == nullptr) {
        last_loop_kf_ = cur_kf_;
        return;
    }

    // 控制回环检测频率：距离上一次触发候选搜索太近时跳过，减少NDT匹配和图优化开销。
    if (last_loop_kf_ && (cur_kf_->GetID() - last_loop_kf_->GetID()) <= options_.loop_kf_gap_) {
        LOG(INFO) << "skip because last loop kf: " << last_loop_kf_->GetID();
        return;
    }

    // 遍历所有历史关键帧，寻找回环候选
    for (auto kf : kfs_mapping) {
        // 1. 跳过连续相近的候选帧，避免冗余
        //    一段历史轨迹里相邻关键帧往往描述同一片区域，只保留间隔足够大的候选。
        if (check_first != nullptr && abs(int(kf->GetID() - check_first->GetID())) <= options_.min_id_interval_) {
            continue;
        }

        // 2. 跳过与当前帧过于接近的关键帧，避免短期回环
        //    当前帧附近的历史帧通常只是连续里程计约束，不应作为回环约束重复加入图优化。
        if (abs(int(kf->GetID() - cur_kf_->GetID())) < options_.closest_id_th_) {
            break;  // 轨迹按ID排序，后续帧会更近，直接退出
        }

        // 3. 计算空间距离，只考虑XY平面距离（忽略高度差）
        //    候选搜索用优化位姿，因为它已经融合了此前回环/位姿图优化的全局一致性。
        Vec3d dt = kf->GetOptPose().translation() - cur_kf_->GetOptPose().translation();
        double t2d = dt.head<2>().norm();  // x-y平面距离
        double range_th = options_.max_range_;

        // 4. 如果空间距离在阈值内，添加为回环候选
        if (t2d < range_th) {
            LoopCandidate c(kf->GetID(), cur_kf_->GetID());
            // 计算初始相对位姿变换（LIO坐标系下的变换）
            // 后续NDT会继续精配准，这里提供前端里程计给出的相对位姿初值。
            c.Tij_ = kf->GetLIOPose().inverse() * cur_kf_->GetLIOPose();

            candidates_.emplace_back(c);
            check_first = kf;  // 记录第一个候选帧，用于后续间隔控制
        }
    }

    if (!candidates_.empty()) {
        // 只有本轮确实找到候选时才更新时间戳，避免无候选帧抑制后续回环搜索。
        last_loop_kf_ = cur_kf_;
    }

    if (options_.verbose_ && !candidates_.empty()) {
        LOG(INFO) << "lc candi: " << candidates_.size();
    }
}

/**
 * @brief 计算并筛选当前关键帧的回环候选。
 *
 * @details 对 DetectLoopCandidates() 生成的候选逐个执行点云配准，
 *          并根据 NDT 分数阈值保留可信回环约束。筛选完成后，
 *          candidates_ 中只保留后续位姿图优化可使用的成功候选。
 */
void LoopClosing::ComputeLoopCandidates() {
    if (candidates_.empty()) {
        return;
    }

    // 对每个候选执行NDT配准，更新候选的相对位姿和匹配分数。
    std::for_each(candidates_.begin(), candidates_.end(), [this](LoopCandidate& c) { ComputeForCandidate(c); });

    // 仅保留分数超过阈值的候选，避免低置信度回环进入位姿图。
    std::vector<LoopCandidate> succ_candidates;
    for (const auto& lc : candidates_) {
        // LOG(INFO) << "candi " << lc.idx1_ << ", " << lc.idx2_ << " s: " << lc.ndt_score_;
        if (lc.ndt_score_ > options_.ndt_score_th_) {
            succ_candidates.emplace_back(lc);
        }
    }

    if (options_.verbose_) {
        LOG(INFO) << "success: " << succ_candidates.size() << "/" << candidates_.size();
    }

    candidates_.swap(succ_candidates);
}

/**
 * @brief 对单个回环候选执行NDT配准。
 *
 * @param c 待计算的回环候选。函数会更新其 NDT 匹配分数 ndt_score_
 *          以及从历史关键帧到当前关键帧的相对位姿 Tij_。
 *
 * @details 以候选历史关键帧附近的优化位姿构建目标子图，
 *          将当前关键帧点云作为源点云进行多分辨率NDT配准，
 *          最后把配准得到的当前帧世界位姿转换为候选相对位姿。
 */
void LoopClosing::ComputeForCandidate(lightning::LoopCandidate& c) {
    // LOG(INFO) << "aligning " << c.idx1_ << " with " << c.idx2_;
    const int submap_idx_range = 40;
    auto kf1 = all_keyframes_.at(c.idx1_), kf2 = all_keyframes_.at(c.idx2_);

    // 构建给定关键帧附近的稀疏子图；每4帧取一帧以控制NDT目标点云规模。
    auto build_submap = [this](int given_id, bool build_in_world) -> CloudPtr {
        CloudPtr submap(new PointCloudType);
        for (int idx = -submap_idx_range; idx < submap_idx_range; idx += 4) {
            int id = idx + given_id;
            if (id < 0 || id >= all_keyframes_.size()) {
                continue;
            }

            auto kf = all_keyframes_[id];
            CloudPtr cloud = kf->GetCloud();

            // RemoveGround(cloud, 0.1);

            if (cloud->empty()) {
                continue;
            }

            // 关键帧点云保存在Lidar系，优化位姿是IMU到世界系，因此投影时需要补上T_imu_lidar_。
            SE3 Twl = kf->GetOptPose() * T_imu_lidar_;

            if (!build_in_world) {
                // 需要局部子图时，以given_id关键帧IMU系为参考系重新表达邻近Lidar点云。
                Twl = all_keyframes_.at(given_id)->GetOptPose().inverse() * Twl;
            }

            CloudPtr cloud_trans(new PointCloudType);
            pcl::transformPointCloud(*cloud, *cloud_trans, Twl.matrix());

            *submap += *cloud_trans;
        }
        return submap;
    };

    auto submap_kf1 = build_submap(kf1->GetID(), true);

    // 当前关键帧点云作为源点云，后续会配准到历史子图所在的世界系。
    CloudPtr submap_kf2 = kf2->GetCloud();

    if (submap_kf1->empty() || submap_kf2->empty()) {
        // 点云为空时无法形成可信回环约束，置零分数供上层筛除。
        c.ndt_score_ = 0;
        return;
    }

    // 使用当前关键帧Lidar到世界的位姿作为NDT初值，加快收敛并降低局部极值风险。
    Mat4f Tw2 = (kf2->GetOptPose() * T_imu_lidar_).matrix().cast<float>();

    /// 不同分辨率下的匹配
    CloudPtr output(new PointCloudType);
    std::vector<double> res{10.0, 5.0, 2.0, 1.0};

    CloudPtr rough_map1, rough_map2;
    // 多分辨率NDT由粗到细迭代，上一层结果作为下一层初值。
    // TODO: 这里最好不要用PCL的NDT，后续可替换为自研实现。
    for (auto& r : res) {
        pcl::NormalDistributionsTransform<PointType, PointType> ndt;
        ndt.setTransformationEpsilon(0.05);
        ndt.setStepSize(0.7);
        ndt.setMaximumIterations(40);

        ndt.setResolution(r);
        rough_map1 = VoxelGrid(submap_kf1, r * 0.1);
        rough_map2 = VoxelGrid(submap_kf2, r * 0.1);
        ndt.setInputTarget(rough_map1);
        ndt.setInputSource(rough_map2);

        ndt.align(*output, Tw2);
        Tw2 = ndt.getFinalTransformation();

        c.ndt_score_ = ndt.getTransformationProbability();
    }

    Mat4d T = Tw2.cast<double>();
    Quatd q(T.block<3, 3>(0, 0));
    q.normalize();  // NDT矩阵数值误差可能破坏正交性，转四元数后归一化。
    Vec3d t = T.block<3, 1>(0, 3);

    // NDT估计的是当前Lidar帧到世界系的位姿，图优化边使用IMU位姿，需要转回T_wi。
    SE3 T_w_lidar2(q, t);
    SE3 T_w_imu2 = T_w_lidar2 * T_imu_lidar_.inverse();
    c.Tij_ = kf1->GetOptPose().inverse() * T_w_imu2;

    // pcl::io::savePCDFileBinaryCompressed(
    //     "./data/lc_" + std::to_string(c.idx1_) + "_" + std::to_string(c.idx2_) + "_out.pcd", *output);
    // pcl::io::savePCDFileBinaryCompressed(
    //     "./data/lc_" + std::to_string(c.idx1_) + "_" + std::to_string(c.idx2_) + "_tgt.pcd", *rough_map1);
}
// TODO，这里只是位姿图优化，为考虑点云特征
void LoopClosing::PoseOptimization() {
    auto v = std::make_shared<miao::VertexSE3>();
    v->SetId(cur_kf_->GetID());
    v->SetEstimate(cur_kf_->GetOptPose());

    optimizer_->AddVertex(v);
    kf_vert_.emplace_back(v);

    /// 上一个关键帧的运动约束
    /// TODO 3D激光最好是跟前面多个帧都有关联
    // 连续运动约束：与前面2个关键帧建立运动约束
    for (int i = 1; i < 3; i++) {
        int id = cur_kf_->GetID() - i;
        if (id >= 0) {
            auto last_kf = all_keyframes_[id];
            auto e = std::make_shared<miao::EdgeSE3>();
            e->SetVertex(0, optimizer_->GetVertex(last_kf->GetID()));
            e->SetVertex(1, v);

            // 计算两帧之间的相对运动变换
            SE3 motion = last_kf->GetLIOPose().inverse() * cur_kf_->GetLIOPose();
            e->SetMeasurement(motion);
            e->SetInformation(info_motion_);
            optimizer_->AddEdge(e);
        }
    }

    // 高度约束：防止Z轴漂移，保持高度稳定
    if (options_.with_height_) {
        auto e = std::make_shared<miao::EdgeHeightPrior>();
        e->SetVertex(0, v);
        e->SetMeasurement(0);  // 期望高度为0
        e->SetInformation(Mat1d::Identity() * 1.0 / (options_.height_noise_ * options_.height_noise_));
        optimizer_->AddEdge(e);
    }

    // 回环约束：基于检测到的回环候选帧添加约束
    for (auto& c : candidates_) {
        auto e = std::make_shared<miao::EdgeSE3>();
        e->SetVertex(0, optimizer_->GetVertex(c.idx1_));
        e->SetVertex(1, optimizer_->GetVertex(c.idx2_));
        e->SetMeasurement(c.Tij_);  // 回环相对变换
        e->SetInformation(info_loops_);

        // 添加Cauchy鲁棒核函数，降低异常值影响
        auto rk = std::make_shared<miao::RobustKernelCauchy>();
        rk->SetDelta(options_.rk_loop_th_);
        e->SetRobustKernel(rk);

        optimizer_->AddEdge(e);
        edge_loops_.emplace_back(e);  // 保存回环边引用
    }

    if (optimizer_->GetEdges().empty()) {
        return;
    }

    if (candidates_.empty()) {
        return;
    }

    optimizer_->InitializeOptimization();
    optimizer_->SetVerbose(false);

    optimizer_->Optimize(20);

    /// remove outliers
    int cnt_outliers = 0;
    for (auto& e : edge_loops_) {
        if (e->GetRobustKernel() == nullptr) {
            continue;
        }

        if (e->Chi2() > e->GetRobustKernel()->Delta()) {
            e->SetLevel(1);
            cnt_outliers++;
        } else {
            e->SetRobustKernel(nullptr);
        }
    }

    if (options_.verbose_) {
        LOG(INFO) << "loop outliers: " << cnt_outliers << "/" << edge_loops_.size();
    }

    /// get results
    for (auto& vert : kf_vert_) {
        SE3 pose = vert->Estimate();
        all_keyframes_[vert->GetId()]->SetOptPose(pose);
    }

    if (loop_cb_) {
        loop_cb_();
    }

    LOG(INFO) << "optimize finished, loops: " << edge_loops_.size();

    // LOG(INFO) << "lc: cur kf " << cur_kf_->GetID() << ", opt: " << cur_kf_->GetOptPose().translation().transpose()
    //           << ", lio: " << cur_kf_->GetLIOPose().translation().transpose();
}

}  // namespace lightning

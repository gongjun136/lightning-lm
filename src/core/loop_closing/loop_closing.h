//
// Created by xiang on 25-4-21.
//

#ifndef LIGHTNING_LOOP_CLOSING_H
#define LIGHTNING_LOOP_CLOSING_H

#include "common/keyframe.h"
#include "common/loop_candidate.h"
#include "utils/async_message_process.h"

#include "core/graph/optimizer.h"
#include "core/types/edge_se3.h"

namespace lightning {

/**
 * @brief 基于NDT匹配和位姿图优化的回环检测模块。
 *
 * LoopClosing接收LIO前端产生的关键帧，按关键帧间隔搜索历史候选帧，
 * 使用NDT估计回环约束，再将里程计约束和回环约束加入位姿图做增量优化。
 * 在线模式下关键帧通过异步线程处理，离线模式下AddKF()会同步执行完整流程。
 */
class LoopClosing {
   public:
    /**
     * @brief 回环检测与位姿图优化参数。
     */
    struct Options {
        Options() {}

        bool verbose_ = true;       ///< 是否输出调试信息。
        bool online_mode_ = false;  ///< 是否使用在线异步处理模式。

        int loop_kf_gap_ = 20;       ///< 两次回环检测之间至少间隔的关键帧数量。
        int min_id_interval_ = 20;   ///< 候选关键帧之间的最小ID间隔。
        int closest_id_th_ = 50;     ///< 历史关键帧与当前关键帧的最小ID间隔。
        double max_range_ = 30.0;    ///< 候选关键帧允许的最大平面距离。
        double ndt_score_th_ = 1.0;  ///< NDT匹配成功的分数阈值。

        double motion_trans_noise_ = 0.1;               ///< 相邻关键帧运动约束的平移噪声。
        double motion_rot_noise_ = 3.0 * M_PI / 180.0;  ///< 相邻关键帧运动约束的旋转噪声。

        double loop_trans_noise_ = 0.2;               ///< 回环约束的平移噪声。
        double loop_rot_noise_ = 3.0 * M_PI / 180.0;  ///< 回环约束的旋转噪声。

        double rk_loop_th_ = 5.2 / 5;  ///< 回环边的鲁棒核阈值。

        bool with_height_ = true;     ///< 是否在位姿图中加入高度先验约束。
        double height_noise_ = 0.1;   ///< 高度先验约束噪声。
    };

    /**
     * @brief 构造回环检测模块。
     * @param options 回环检测和优化参数。
     */
    LoopClosing(Options options = Options()) { options_ = options; }

    /**
     * @brief 析构模块；在线模式下会停止异步关键帧处理线程。
     */
    ~LoopClosing();

    /**
     * @brief 初始化优化器、信息矩阵和在线处理线程。
     * @param yaml_path YAML配置文件路径；为空时使用Options默认值。
     */
    void Init(const std::string yaml_path);

    /**
     * @brief 添加一个待处理关键帧。
     * @param kf LIO前端生成的关键帧。
     */
    void AddKF(Keyframe::Ptr kf);

    /**
     * @brief 回环优化完成后的通知回调类型。
     */
    using LoopClosedCallback = std::function<void()>;

    /**
     * @brief 设置回环闭合回调。
     * @param cb 当检测到有效回环并完成优化后调用的函数。
     */
    void SetLoopClosedCB(LoopClosedCallback cb) { loop_cb_ = cb; }

   protected:
    /**
     * @brief 处理单个关键帧，包含候选检测、候选匹配和位姿图优化。
     * @param kf 待处理关键帧。
     */
    void HandleKF(Keyframe::Ptr kf);

    /**
     * @brief 从历史关键帧中筛选当前关键帧的回环候选。
     */
    void DetectLoopCandidates();

    /**
     * @brief 对回环候选执行匹配并保留成功候选。
     */
    void ComputeLoopCandidates();

    /**
     * @brief 计算单个回环候选的相对位姿约束。
     * @param c 待计算的回环候选。
     */
    void ComputeForCandidate(LoopCandidate& c);

    /**
     * @brief 根据运动边和回环边执行位姿图优化。
     */
    void PoseOptimization();

    Options options_;  ///< 当前回环检测参数。

    Keyframe::Ptr last_kf_ = nullptr;           ///< 上一次处理的关键帧，用于避免重复处理。
    Keyframe::Ptr last_loop_kf_ = nullptr;      ///< 上一次成功触发回环检测的关键帧。
    Keyframe::Ptr cur_kf_ = nullptr;            ///< 当前正在处理的关键帧。
    std::vector<Keyframe::Ptr> all_keyframes_;  ///< 所有已处理的关键帧。
    std::vector<LoopCandidate> candidates_;     ///< 当前关键帧的回环候选。

    AsyncMessageProcess<Keyframe::Ptr> kf_thread_;  ///< 在线模式下的异步关键帧处理线程。

    std::shared_ptr<miao::Optimizer> optimizer_ = nullptr;  ///< Miao图优化器实例。

    Mat6d info_motion_ = Mat6d::Identity();  ///< 相邻关键帧运动约束的信息矩阵。
    Mat6d info_loops_ = Mat6d::Identity();   ///< 回环约束的信息矩阵。

    std::vector<std::shared_ptr<miao::VertexSE3>> kf_vert_;  ///< 位姿图中的关键帧SE3顶点。
    std::vector<std::shared_ptr<miao::EdgeSE3>> edge_loops_;  ///< 回环检测生成的SE3约束边。

    LoopClosedCallback loop_cb_;  ///< 回环优化完成后的通知回调。
};

}  // namespace lightning

#endif  // LIGHTNING_LOOP_CLOSING_H

#ifndef LIGHTNING_BACKEND_BACKEND_PIPELINE_H
#define LIGHTNING_BACKEND_BACKEND_PIPELINE_H

#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "core/backend/btc_loop_detector.h"
#include "core/backend/hierarchical_bundle_adjustment.h"

namespace lightning::backend {

enum class BackendMode { kLegacy, kBaBtcHba, kDisabled };

BackendMode ReadBackendMode(const std::string& yaml_path);
const char* BackendModeName(BackendMode mode);

struct BackendPipelineOptions {
    bool online_mode = false;
    bool verbose = true;
    int pose_graph_iterations = 20;
    int pose_graph_outlier_iterations = 10;
    double motion_translation_noise = 0.10;
    double motion_rotation_noise_deg = 3.0;
    double loop_translation_noise = 0.20;
    double loop_rotation_noise_deg = 3.0;
    double robust_kernel_delta = 1.0;
    double loop_outlier_chi2 = 25.0;
    BundleAdjustmentOptions local_ba;
    BtcLoopDetectorOptions btc;
    HierarchicalBundleAdjustmentOptions hba;
};

struct BackendRuntimeSummary {
    std::size_t keyframes = 0;
    std::size_t local_ba_attempts = 0;
    std::size_t local_ba_accepted = 0;
    std::size_t btc_descriptors = 0;
    std::size_t btc_candidates = 0;
    std::size_t loops_accepted = 0;
    std::size_t loops_applied = 0;
    std::size_t loops_graph_inliers = 0;
    std::size_t hba_runs = 0;
    std::size_t hba_accepted = 0;
    double local_ba_time_ms = 0.0;
    double btc_time_ms = 0.0;
    double pose_graph_time_ms = 0.0;
    double hba_time_ms = 0.0;
};

class BackendPipeline {
   public:
    BackendPipeline() = default;
    ~BackendPipeline();

    bool Init(const std::string& yaml_path, bool online_mode);
    void AddKeyframe(const Keyframe::Ptr& keyframe);
    void RequestGlobalOptimization(const std::string& reason = "manual");
    void WaitUntilIdle(bool force_global_optimization);
    void Shutdown();

    using OptimizedCallback = std::function<void()>;
    void SetOptimizedCallback(OptimizedCallback callback);

    SE3 GetMapToOdom() const;
    std::vector<BtcLoopResult> GetLoopResults() const;
    BackendRuntimeSummary GetRuntimeSummary() const;
    bool SaveDiagnostics(const std::string& directory,
                         const map_frame::Metadata* map_metadata = nullptr) const;
    bool SaveRelocalizationDatabase(
        const std::string& directory,
        const map_frame::Metadata* map_metadata = nullptr) const;
    const BackendPipelineOptions& GetOptions() const { return options_; }

   private:
    struct LoopConstraint {
        BtcLoopResult detection;
        SE3 T_history_imu_current_imu;
        bool graph_inlier = true;
    };

    void HandleKeyframe(const Keyframe::Ptr& keyframe);
    void WorkerLoop();
    void HbaLoop();
    void RequestHba(bool rerun_pose_graph, const std::string& reason);
    bool OptimizePoseGraph();
    void UpdateMapToOdomAndNotify();

    BackendPipelineOptions options_;
    SE3 T_imu_lidar_;
    std::unique_ptr<VoxelBundleAdjuster> local_ba_;
    std::unique_ptr<BtcLoopDetector> btc_;
    std::unique_ptr<HierarchicalBundleAdjuster> hba_;

    mutable std::mutex data_mutex_;
    std::vector<Keyframe::Ptr> keyframes_;
    std::vector<BtcLoopResult> loop_results_;
    std::vector<LoopConstraint> loop_constraints_;
    BackendRuntimeSummary runtime_;
    SE3 T_map_odom_;
    OptimizedCallback optimized_callback_;

    std::mutex optimization_mutex_;

    std::mutex queue_mutex_;
    std::condition_variable queue_cv_;
    std::condition_variable queue_idle_cv_;
    std::deque<Keyframe::Ptr> queue_;
    std::thread worker_thread_;
    bool worker_processing_ = false;
    bool worker_stop_ = false;

    std::mutex hba_mutex_;
    std::condition_variable hba_cv_;
    std::condition_variable hba_idle_cv_;
    std::thread hba_thread_;
    bool hba_requested_ = false;
    bool hba_rerun_pose_graph_ = false;
    bool hba_running_ = false;
    bool hba_stop_ = false;
    std::string hba_reason_;

    bool initialized_ = false;
    std::size_t last_local_ba_keyframe_count_ = 0;
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_BACKEND_PIPELINE_H

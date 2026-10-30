#pragma once

/**
 * @file
 * @brief 在线定位的统一计算预算及 YAML、环境变量加载接口。
 * @see lightning::compute::ComputeBudget, lightning::compute::LoadComputeBudget()
 */

#include <string>

namespace YAML {
class Node;
}

namespace lightning::compute {

constexpr const char* kLioThreadsEnv = "LIGHTNING_LM_LIO_THREADS";  ///< 覆盖 LIO 与预处理线程数。
constexpr const char* kNdtThreadsEnv = "LIGHTNING_LM_NDT_THREADS";  ///< 覆盖 NDT 地图匹配线程数。
constexpr const char* kNdtMaxPointsEnv = "LIGHTNING_LM_NDT_MAX_POINTS";  ///< 覆盖跟踪阶段的 NDT 输入点数上限。
constexpr const char* kSolidIcpWorkersEnv = "LIGHTNING_LM_SOLID_ICP_WORKERS";  ///< 覆盖 SOLiD 候选粗配准线程数。
constexpr const char* kSolidWorkerNiceEnv = "LIGHTNING_LM_SOLID_WORKER_NICE";  ///< 覆盖 SOLiD 工作线程的 nice 下限。
constexpr const char* kSolidCpuAffinityEnv = "LIGHTNING_LM_SOLID_CPU_AFFINITY";  ///< 覆盖 SOLiD 工作线程的 CPU 列表。

/**
 * @brief 统一描述 LIO、NDT 和 SOLiD 在线定位的线程、调度与点数预算。
 * @details 调用方可先填入模块原有默认值，再通过 LoadComputeBudget() 覆盖配置。
 * 本结构仅保存预算；线程创建、调度应用和点云采样由使用预算的模块完成。
 * @see LoadComputeBudget(), @ref sensor_pipeline "传感器处理与算力预算"
 */
struct ComputeBudget {
    int lio_threads = 1;  ///< LIO 与预处理并行线程数，有效范围 [1, 128]。
    int ndt_threads = 4;  ///< NDT 地图匹配并行线程数，有效范围 [1, 128]。
    int solid_icp_workers = 8;  ///< SOLiD 候选粗配准工作线程数，有效范围 [1, 128]。
    int solid_worker_nice = 0;  ///< SOLiD 工作线程的 nice 下限 [0, 19]；0 表示不主动调整 nice。
    std::string solid_cpu_affinity;  ///< CPU 编号/闭区间列表，如 "1-3,5"；空值表示不施加额外亲和性限制。
    int ndt_max_points = 0;  ///< 跟踪阶段 NDT 输入点数上限，0 表示不限制，否则为 [100, 1000000]；不用于全局初始化。
};

/**
 * @brief 以现有预算为默认值，读取 YAML 和环境覆盖，校验成功后整体写回。
 * @par 主要流程
 * 1. <b>加载候选配置</b>：复制现有预算，读取可选的 compute_budget 配置及整数环境覆盖。
 * 2. <b>校验最终预算</b>：应用 CPU 列表环境覆盖，检查点数、线程、nice 和 CPU 列表。
 * 3. <b>提交预算</b>：全部通过后才更新 budget，失败或异常时保留原值。
 * @param root YAML 根节点，compute_budget 为可选的配置映射。
 * @param budget 输入为调用方默认预算，成功时被最终结果整体替换。
 * @param error 可选的失败原因输出；成功时不清空其已有内容。
 * @return 配置读取和范围校验成功返回 true；检测到转换或校验错误返回 false。
 * @throws YAML::Exception YAML 节点索引等未被内部转换处理捕获的结构错误可能向外传播。
 * @note 优先级为环境变量、YAML、传入默认值。整数环境变量未设置或为空时不覆盖；
 * CPU 列表环境变量设置为空字符串时会清除已有列表。每次调用重新读取环境变量。
 * CPU 列表仅校验语法和编号范围，实际允许的 CPU 集合与调度应用由线程配置阶段检查。
 * @see ComputeBudget, lightning::threading::ParseCpuList()
 */
bool LoadComputeBudget(const YAML::Node& root, ComputeBudget& budget,
                       std::string* error = nullptr);

}  // namespace lightning::compute

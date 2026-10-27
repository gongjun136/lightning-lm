#pragma once
#include <memory>
#include <string>

namespace lightning {
/** CGI rear-axle navigation, strict publication gating and INS-driven LiDAR deskew.
 * Owns a single-threaded ROS executor's node; heartbeat has its own worker.
 * Init validates the immutable site configuration and opens audit files.
 * Invalid configuration returns false; Spin runs until ROS shutdown.
 * See @ref ins_only_operation for input clocks, references and failure semantics.
 */
class InsLocSystem {
 public:
    InsLocSystem();
    ~InsLocSystem();
    /**
     * @brief 校验纯 INS 定位配置，连接导航输入、业务输出和审计文件。
     * @par 主要流程
     * 1. <b>校验坐标参考</b>：确认 WGS84、后轴参考点、ENU 原点及高度基准。
     * 2. <b>配置质量门控</b>：读取时间偏移、质量阈值、恢复样本数和位姿插值间隔。
     * 3. <b>创建节点与心跳</b>：设置 ROS 时间和地图坐标系，启动独立心跳线程。
     * 4. <b>打开审计文件</b>：保存配置快照，打开原始字段、导航结果、事件和 TUM 文件。
     * 5. <b>创建输出通道</b>：创建位姿、状态、轨迹、点云发布器及可选 TF。
     * 6. <b>注册导航与点云输入</b>：订阅 CGI 导航字段，按开关初始化雷达输入。
     * 7. <b>启动健康检查</b>：注册周期检查定时器，记录初始化完成并返回。
     * @note 阶段名称对应函数体；ROS 回调由 Spin() 调度，初始化成功不代表导航质量已通过。
     * @return 成功返回 true；捕获到标准异常时记录错误并返回 false，已创建的文件不会回滚删除。
     */
    bool Init(const std::string& yaml_path, const std::string& published_tum = "");
    void Spin();
 private:
    struct Impl;
    std::unique_ptr<Impl> impl_;  ///< INS 输入、质量门控、点云投影与业务发布的实现对象。
};
}  // namespace lightning

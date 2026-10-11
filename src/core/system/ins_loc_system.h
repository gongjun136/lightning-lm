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
    bool Init(const std::string& yaml_path, const std::string& published_tum = "");
    void Spin();
 private:
    struct Impl;
    std::unique_ptr<Impl> impl_;  ///< INS 输入、质量门控、点云投影与业务发布的实现对象。
};
}  // namespace lightning

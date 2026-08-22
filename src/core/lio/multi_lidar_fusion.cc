#include "core/lio/multi_lidar_fusion.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <tuple>
#include <utility>

namespace lightning {
namespace {

void SetError(std::string* error, const std::string& message) {
    if (error) {
        *error = message;
    }
}

std::string InferImuTopic(std::string lidar_topic) {
    const std::string token = "/lidar_";
    const auto pos = lidar_topic.find(token);
    if (pos != std::string::npos) {
        lidar_topic.replace(pos, token.size(), "/imu_");
    }
    return lidar_topic;
}

bool ParseId(const std::string& key, int& id) {
    std::string digits;
    for (const char c : key) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits.push_back(c);
        }
    }
    if (digits.empty()) {
        return false;
    }
    id = std::stoi(digits);
    return true;
}

bool ReadTransform(const YAML::Node& node, Mat3d& R, Vec3d& t, std::string* error) {
    if (!node) {
        SetError(error, "missing lidar transform");
        return false;
    }
    std::vector<double> values;
    try {
        values = node.as<std::vector<double>>();
    } catch (const YAML::Exception& e) {
        SetError(error, e.what());
        return false;
    }
    if (values.size() != 12 && values.size() != 16) {
        SetError(error, "lidar transform must contain 12 or 16 numbers");
        return false;
    }
    if (!std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); })) {
        SetError(error, "lidar transform contains a non-finite value");
        return false;
    }
    if (values.size() == 16 &&
        (std::abs(values[12]) > 1e-9 || std::abs(values[13]) > 1e-9 || std::abs(values[14]) > 1e-9 ||
         std::abs(values[15] - 1.0) > 1e-9)) {
        SetError(error, "lidar transform bottom row must be [0, 0, 0, 1]");
        return false;
    }
    R << values[0], values[1], values[2], values[4], values[5], values[6], values[8], values[9], values[10];
    t << values[3], values[7], values[11];
    const double orthogonality_error = (R.transpose() * R - Mat3d::Identity()).cwiseAbs().maxCoeff();
    if (orthogonality_error > 1e-3 || std::abs(R.determinant() - 1.0) > 1e-3) {
        SetError(error, "lidar rotation is not a valid SO(3) matrix");
        return false;
    }
    return true;
}

bool ValidateConfig(const MultiLidarConfig& config, std::string* error) {
    if (!config.enabled) {
        return true;
    }
    if (config.frame_period <= 0.0 || config.match_tolerance <= 0.0 ||
        config.match_tolerance >= 0.5 * config.frame_period || config.reorder_window < config.frame_period) {
        SetError(error, "invalid multi_lidar timing parameters");
        return false;
    }
    if (config.online_extrinsic_estimation) {
        SetError(error, "online multi-lidar extrinsic estimation is not supported in the stable V1 frontend");
        return false;
    }
    if (config.min_lidars < 1 || config.min_lidars > static_cast<int>(config.lidars.size())) {
        SetError(error, "multi_lidar.min_lidars is outside the configured lidar count");
        return false;
    }
    const auto& adaptive = config.adaptive_load;
    const int lidar_count = static_cast<int>(config.lidars.size());
    if (adaptive.enabled) {
        if (adaptive.tracking_min_lidars < 1 || adaptive.tracking_min_lidars > lidar_count ||
            adaptive.relocalization_min_lidars < adaptive.tracking_min_lidars ||
            adaptive.relocalization_min_lidars > lidar_count ||
            adaptive.cloud_publish_min_lidars < 1 || adaptive.cloud_publish_min_lidars > lidar_count) {
            SetError(error, "invalid multi_lidar.adaptive_load lidar-count thresholds");
            return false;
        }
        if (!std::isfinite(adaptive.target_latency_sec) ||
            !std::isfinite(adaptive.hard_latency_sec) || adaptive.target_latency_sec <= 0.0 ||
            adaptive.hard_latency_sec <= adaptive.target_latency_sec ||
            adaptive.degrade_processing_ratio <= 0.0 || adaptive.degrade_processing_ratio > 1.0 ||
            adaptive.recover_processing_ratio <= 0.0 ||
            adaptive.recover_processing_ratio >= adaptive.degrade_processing_ratio ||
            adaptive.degrade_consecutive_frames <= 0 || adaptive.recover_consecutive_frames <= 0 ||
            adaptive.point_strides.empty()) {
            SetError(error, "invalid multi_lidar.adaptive_load timing parameters");
            return false;
        }
        int previous_stride = 0;
        for (const int stride : adaptive.point_strides) {
            if (stride <= 0 || stride < previous_stride) {
                SetError(error, "multi_lidar.adaptive_load.point_strides must be positive and ordered");
                return false;
            }
            previous_stride = stride;
        }
    }
    std::set<int> ids;
    bool has_primary = false;
    for (const auto& sensor : config.lidars) {
        if (sensor.id < 0 || sensor.id > 255 || sensor.lidar_topic.empty() || !ids.insert(sensor.id).second) {
            SetError(error, "invalid or duplicate multi-lidar sensor entry");
            return false;
        }
        has_primary = has_primary || sensor.id == config.primary_lidar_id;
    }
    if (!has_primary) {
        SetError(error, "primary lidar id is not configured");
        return false;
    }
    return true;
}

}  // namespace

const MultiLidarSensorConfig* MultiLidarConfig::FindLidar(int id) const {
    for (const auto& lidar : lidars) {
        if (lidar.id == id) {
            return &lidar;
        }
    }
    return nullptr;
}

bool LoadMultiLidarConfig(const YAML::Node& root, MultiLidarConfig& config, std::string* error) {
    config = MultiLidarConfig();
    const YAML::Node multi = root["multi_lidar"];
    if (!multi) {
        return true;
    }
    try {
        config.enabled = multi["enabled"] ? multi["enabled"].as<bool>() : true;
        if (!config.enabled) {
            return true;
        }
        if (multi["primary_lidar_id"]) config.primary_lidar_id = multi["primary_lidar_id"].as<int>();
        if (multi["frame_period"]) config.frame_period = multi["frame_period"].as<double>();
        if (multi["match_tolerance"]) config.match_tolerance = multi["match_tolerance"].as<double>();
        if (multi["reorder_window"]) config.reorder_window = multi["reorder_window"].as<double>();
        if (multi["min_lidars"]) config.min_lidars = multi["min_lidars"].as<int>();
        if (multi["online_extrinsic_estimation"]) {
            config.online_extrinsic_estimation = multi["online_extrinsic_estimation"].as<bool>();
        }
        const YAML::Node adaptive = multi["adaptive_load"];
        if (adaptive) {
            auto& load = config.adaptive_load;
            load.enabled = adaptive["enabled"] ? adaptive["enabled"].as<bool>() : true;
            if (adaptive["tracking_min_lidars"]) {
                load.tracking_min_lidars = adaptive["tracking_min_lidars"].as<int>();
            }
            if (adaptive["relocalization_min_lidars"]) {
                load.relocalization_min_lidars = adaptive["relocalization_min_lidars"].as<int>();
            }
            if (adaptive["cloud_publish_min_lidars"]) {
                load.cloud_publish_min_lidars = adaptive["cloud_publish_min_lidars"].as<int>();
            }
            if (adaptive["cloud_publish_require_primary"]) load.cloud_publish_require_primary = adaptive["cloud_publish_require_primary"].as<bool>();
            if (adaptive["target_latency_sec"]) load.target_latency_sec = adaptive["target_latency_sec"].as<double>();
            if (adaptive["hard_latency_sec"]) load.hard_latency_sec = adaptive["hard_latency_sec"].as<double>();
            if (adaptive["degrade_processing_ratio"]) load.degrade_processing_ratio = adaptive["degrade_processing_ratio"].as<double>();
            if (adaptive["recover_processing_ratio"]) load.recover_processing_ratio = adaptive["recover_processing_ratio"].as<double>();
            if (adaptive["degrade_consecutive_frames"]) load.degrade_consecutive_frames = adaptive["degrade_consecutive_frames"].as<int>();
            if (adaptive["recover_consecutive_frames"]) load.recover_consecutive_frames = adaptive["recover_consecutive_frames"].as<int>();
            if (adaptive["point_strides"]) load.point_strides = adaptive["point_strides"].as<std::vector<int>>();
        }

        const YAML::Node topics = multi["topics"];
        const YAML::Node extrinsics = multi["extrinsics"];
        if (!topics || !extrinsics) {
            SetError(error, "multi_lidar.topics and multi_lidar.extrinsics are required");
            return false;
        }
        std::map<int, MultiLidarSensorConfig> sensors;
        for (const auto& item : topics) {
            const std::string key = item.first.as<std::string>();
            if (key.rfind("lidar", 0) != 0) {
                continue;
            }
            int id = 0;
            if (!ParseId(key, id)) {
                continue;
            }
            MultiLidarSensorConfig sensor;
            sensor.id = id;
            sensor.lidar_topic = item.second.as<std::string>();
            const std::string imu_key = "imu_" + std::to_string(id);
            sensor.imu_topic = topics[imu_key] ? topics[imu_key].as<std::string>() : InferImuTopic(sensor.lidar_topic);
            YAML::Node extrinsic = extrinsics["lidar" + std::to_string(id)];
            if (!extrinsic) extrinsic = extrinsics["lidar_" + std::to_string(id)];
            if (!extrinsic || !ReadTransform(extrinsic["T"], sensor.R_lidar_to_primary,
                                             sensor.t_lidar_to_primary, error)) {
                return false;
            }
            sensors.emplace(id, std::move(sensor));
        }
        for (auto& [id, sensor] : sensors) {
            (void)id;
            config.lidars.push_back(std::move(sensor));
        }
    } catch (const YAML::Exception& e) {
        SetError(error, e.what());
        return false;
    }
    return ValidateConfig(config, error);
}

bool LoadSelfPointFilterConfig(const YAML::Node& root, SelfPointFilterConfig& config,
                               std::string* error) {
    config = SelfPointFilterConfig();
    const YAML::Node filter = root["self_point_filter"];
    if (!filter) return true;

    try {
        config.enabled = filter["enabled"] ? filter["enabled"].as<bool>() : true;
        if (!config.enabled) return true;

        const double padding = filter["padding"] ? filter["padding"].as<double>() : 0.0;
        if (!std::isfinite(padding) || padding < 0.0) {
            SetError(error, "self_point_filter.padding must be finite and non-negative");
            return false;
        }
        if (filter["min_body"] || filter["max_body"]) {
            if (!filter["min_body"] || !filter["max_body"]) {
                SetError(error, "self_point_filter.min_body and max_body must be specified together");
                return false;
            }
            const auto minimum = filter["min_body"].as<std::vector<double>>();
            const auto maximum = filter["max_body"].as<std::vector<double>>();
            if (minimum.size() != 3 || maximum.size() != 3) {
                SetError(error, "self_point_filter.min_body and max_body must contain three values");
                return false;
            }
            config.min_body = Vec3d(minimum[0], minimum[1], minimum[2]);
            config.max_body = Vec3d(maximum[0], maximum[1], maximum[2]);
            if (!config.min_body.allFinite() || !config.max_body.allFinite() ||
                (config.min_body.array() >= config.max_body.array()).any()) {
                SetError(error, "self_point_filter explicit body bounds are invalid");
                return false;
            }
            config.min_body.array() -= padding;
            config.max_body.array() += padding;
            return true;
        }

        const char* required[] = {"front", "back", "left", "right", "bottom", "top"};
        for (const char* key : required) {
            if (!filter[key]) {
                SetError(error, std::string("self_point_filter.") + key + " is required");
                return false;
            }
        }
        const double front = filter["front"].as<double>();
        const double back = filter["back"].as<double>();
        const double left = filter["left"].as<double>();
        const double right = filter["right"].as<double>();
        const double bottom = filter["bottom"].as<double>();
        const double top = filter["top"].as<double>();
        const double values[] = {front, back, left, right, bottom, top, padding};
        if (!std::all_of(std::begin(values), std::end(values),
                         [](double value) { return std::isfinite(value) && value >= 0.0; })) {
            SetError(error, "self_point_filter extents and padding must be finite and non-negative");
            return false;
        }
        config.min_body = Vec3d(-back - padding, -right - padding, -bottom - padding);
        config.max_body = Vec3d(front + padding, left + padding, top + padding);
    } catch (const YAML::Exception& e) {
        SetError(error, e.what());
        return false;
    }
    return true;
}

std::size_t FilterSelfPoints(PointCloudType& cloud, const SelfPointFilterConfig& config,
                             const Mat3d& R_primary_to_body) {
    if (!config.enabled || cloud.empty()) return 0;
    const auto old_size = cloud.size();
    cloud.erase(std::remove_if(cloud.begin(), cloud.end(), [&](const PointType& point) {
                    const Vec3d p_body = R_primary_to_body * point.getVector3fMap().cast<double>();
                    return (p_body.array() >= config.min_body.array()).all() &&
                           (p_body.array() <= config.max_body.array()).all();
                }),
                cloud.end());
    cloud.width = cloud.size();
    cloud.height = 1;
    cloud.is_dense = false;
    return old_size - cloud.size();
}

CloudPtr DownsamplePreservingSource(const CloudPtr& cloud, double leaf_size) {
    CloudPtr filtered(new PointCloudType);
    if (!cloud || cloud->empty() || leaf_size <= 0.0) {
        if (cloud) *filtered = *cloud;
        return filtered;
    }
    using Key = std::tuple<long long, long long, long long>;
    struct Representative {
        PointType point;
        double center_distance = std::numeric_limits<double>::max();
    };
    std::map<Key, Representative> representatives;
    for (const auto& point : cloud->points) {
        const Vec3d p = point.getVector3fMap().cast<double>();
        const long long ix = static_cast<long long>(std::floor(p.x() / leaf_size));
        const long long iy = static_cast<long long>(std::floor(p.y() / leaf_size));
        const long long iz = static_cast<long long>(std::floor(p.z() / leaf_size));
        const Vec3d center((ix + 0.5) * leaf_size, (iy + 0.5) * leaf_size, (iz + 0.5) * leaf_size);
        const double distance = (p - center).squaredNorm();
        const Key key{ix, iy, iz};
        auto it = representatives.find(key);
        if (it == representatives.end() || distance < it->second.center_distance) {
            representatives[key] = Representative{point, distance};
        }
    }
    filtered->reserve(representatives.size());
    for (const auto& [key, representative] : representatives) {
        (void)key;
        filtered->push_back(representative.point);
    }
    filtered->header = cloud->header;
    filtered->width = filtered->size();
    filtered->height = 1;
    filtered->is_dense = false;
    return filtered;
}

void AdaptiveLidarLoadController::Reset(const MultiLidarConfig& config) {
    config_ = config;
    localization_good_ = false;
    degradation_step_ = 0;
    overload_frames_ = 0;
    recovery_frames_ = 0;
}

int AdaptiveLidarLoadController::MaximumDegradationStep() const {
    if (!config_.adaptive_load.enabled) return 0;
    const int point_steps = std::max(0, static_cast<int>(config_.adaptive_load.point_strides.size()) - 1);
    const int lidar_steps = std::max(0, static_cast<int>(config_.lidars.size()) -
                                           config_.adaptive_load.tracking_min_lidars);
    return point_steps + lidar_steps;
}

int AdaptiveLidarLoadController::PointStride() const {
    if (!config_.adaptive_load.enabled || config_.adaptive_load.point_strides.empty()) return 1;
    const int index = std::min(degradation_step_,
                               static_cast<int>(config_.adaptive_load.point_strides.size()) - 1);
    return config_.adaptive_load.point_strides[index];
}

int AdaptiveLidarLoadController::TargetLidarCount() const {
    const int total = static_cast<int>(config_.lidars.size());
    if (!config_.adaptive_load.enabled) return total;
    const int point_steps = std::max(0, static_cast<int>(config_.adaptive_load.point_strides.size()) - 1);
    const int lidar_reduction = std::max(0, degradation_step_ - point_steps);
    const int minimum = localization_good_.load() ? config_.adaptive_load.tracking_min_lidars
                                           : config_.adaptive_load.relocalization_min_lidars;
    return std::max(minimum, total - lidar_reduction);
}

void AdaptiveLidarLoadController::Observe(double processing_sec, double latency_sec,
                                           bool tracking_healthy) {
    if (!config_.adaptive_load.enabled) return;
    const auto& load = config_.adaptive_load;
    const bool overloaded = latency_sec >= load.target_latency_sec ||
                            processing_sec >= load.target_latency_sec * load.degrade_processing_ratio;
    const bool recovered = latency_sec <= load.target_latency_sec * load.recover_processing_ratio &&
                           processing_sec <= load.target_latency_sec * load.recover_processing_ratio;
    if (overloaded) {
        recovery_frames_ = 0;
        if (!tracking_healthy) {
            overload_frames_ = 0;
            return;
        }
        if (++overload_frames_ >= load.degrade_consecutive_frames) {
            degradation_step_ = std::min(MaximumDegradationStep(), degradation_step_ + 1);
            overload_frames_ = 0;
        }
        return;
    }
    overload_frames_ = 0;
    const int required_recovery_frames = tracking_healthy ? load.recover_consecutive_frames
                                                          : load.degrade_consecutive_frames;
    if ((!tracking_healthy || recovered) && degradation_step_ > 0) {
        if (++recovery_frames_ >= required_recovery_frames) {
            --degradation_step_;
            recovery_frames_ = 0;
        }
    } else {
        recovery_frames_ = 0;
    }
}

AdaptiveLidarSelection AdaptiveLidarLoadController::Select(const MultiLidarFrameStats& stats) const {
    AdaptiveLidarSelection selection;
    selection.point_stride = PointStride();
    selection.degradation_step = degradation_step_;
    const auto primary_it = std::find(stats.present_lidar_ids.begin(), stats.present_lidar_ids.end(),
                                      config_.primary_lidar_id);
    const int minimum = localization_good_.load() ? config_.adaptive_load.tracking_min_lidars
                                           : config_.adaptive_load.relocalization_min_lidars;
    if (primary_it == stats.present_lidar_ids.end() ||
        static_cast<int>(stats.present_lidar_ids.size()) < minimum) {
        return selection;
    }
    selection.lidar_ids.push_back(config_.primary_lidar_id);
    std::vector<int> secondary;
    for (const int id : stats.present_lidar_ids) {
        if (id != config_.primary_lidar_id) secondary.push_back(id);
    }
    std::sort(secondary.begin(), secondary.end(), [&](int lhs, int rhs) {
        const std::size_t lhs_points = stats.points_by_lidar.count(lhs) ? stats.points_by_lidar.at(lhs) : 0;
        const std::size_t rhs_points = stats.points_by_lidar.count(rhs) ? stats.points_by_lidar.at(rhs) : 0;
        return lhs_points != rhs_points ? lhs_points > rhs_points : lhs < rhs;
    });
    const int target = std::min(TargetLidarCount(), static_cast<int>(stats.present_lidar_ids.size()));
    for (const int id : secondary) {
        if (static_cast<int>(selection.lidar_ids.size()) >= target) break;
        selection.lidar_ids.push_back(id);
    }
    return selection;
}

bool AdaptiveLidarLoadController::CanPublishCloud(const MultiLidarFrameStats& stats) const {
    if (!config_.adaptive_load.enabled) return true;
    const auto& load = config_.adaptive_load;
    if (static_cast<int>(stats.present_lidar_ids.size()) < load.cloud_publish_min_lidars) return false;
    return !load.cloud_publish_require_primary ||
           std::find(stats.present_lidar_ids.begin(), stats.present_lidar_ids.end(),
                     config_.primary_lidar_id) != stats.present_lidar_ids.end();
}

bool AdaptiveLidarLoadController::IsHardStale(double latency_sec) const {
    return config_.adaptive_load.enabled && latency_sec >= config_.adaptive_load.hard_latency_sec;
}

CloudPtr SelectLidarPoints(const CloudPtr& cloud, const AdaptiveLidarSelection& selection) {
    CloudPtr selected(new PointCloudType);
    if (!cloud || selection.lidar_ids.empty()) return selected;
    const std::set<int> allowed(selection.lidar_ids.begin(), selection.lidar_ids.end());
    std::map<int, std::size_t> source_counters;
    const int stride = std::max(1, selection.point_stride);
    selected->reserve(cloud->size());
    for (const auto& point : cloud->points) {
        const int id = static_cast<int>(point.lidar_id);
        if (allowed.count(id) == 0) continue;
        const std::size_t index = source_counters[id]++;
        if (index % static_cast<std::size_t>(stride) == 0) selected->push_back(point);
    }
    selected->header = cloud->header;
    selected->width = selected->size();
    selected->height = 1;
    selected->is_dense = false;
    return selected;
}

void MultiLidarFrameAssembler::Reset(MultiLidarConfig config) {
    config_ = std::move(config);
    expected_lidar_ids_.clear();
    for (const auto& sensor : config_.lidars) expected_lidar_ids_.insert(sensor.id);
    frames_.clear();
    ready_frames_.clear();
    anchor_initialized_ = false;
    anchor_time_ = 0.0;
    max_seen_time_ = -std::numeric_limits<double>::infinity();
    last_emitted_bucket_ = std::numeric_limits<long long>::min();
    late_drop_count_ = 0;
    duplicate_drop_count_ = 0;
    tolerance_drop_count_ = 0;
    invalid_drop_count_ = 0;
    emitted_frame_count_ = 0;
    insufficient_lidar_drop_count_ = 0;
}

long long MultiLidarFrameAssembler::BucketFor(double timestamp) const {
    return static_cast<long long>(std::llround((timestamp - anchor_time_) / config_.frame_period));
}

double MultiLidarFrameAssembler::BucketTime(long long bucket) const {
    return anchor_time_ + static_cast<double>(bucket) * config_.frame_period;
}

bool MultiLidarFrameAssembler::IsComplete(const PartialFrame& frame) const {
    return frame.clouds.size() == expected_lidar_ids_.size();
}

bool MultiLidarFrameAssembler::IsExpired(long long bucket) const {
    return max_seen_time_ - BucketTime(bucket) >= config_.reorder_window;
}

bool MultiLidarFrameAssembler::AddCloud(int lidar_id, double timestamp, CloudPtr cloud) {
    if (!cloud || !std::isfinite(timestamp) || expected_lidar_ids_.count(lidar_id) == 0) {
        ++invalid_drop_count_;
        return false;
    }
    if (!anchor_initialized_) {
        anchor_time_ = timestamp;
        anchor_initialized_ = true;
    }
    const long long bucket = BucketFor(timestamp);
    if (bucket <= last_emitted_bucket_) {
        ++late_drop_count_;
        return false;
    }
    auto& frame = frames_[bucket];
    if (frame.clouds.count(lidar_id) != 0) {
        ++duplicate_drop_count_;
        return false;
    }
    if (!frame.timestamps.empty()) {
        double min_timestamp = timestamp;
        double max_timestamp = timestamp;
        for (const auto& [id, existing_timestamp] : frame.timestamps) {
            (void)id;
            min_timestamp = std::min(min_timestamp, existing_timestamp);
            max_timestamp = std::max(max_timestamp, existing_timestamp);
        }
        if (max_timestamp - min_timestamp > config_.match_tolerance) {
            ++tolerance_drop_count_;
            return false;
        }
    }
    frame.clouds.emplace(lidar_id, std::move(cloud));
    frame.timestamps.emplace(lidar_id, timestamp);
    max_seen_time_ = std::max(max_seen_time_, timestamp);
    PromoteReady(false);
    return true;
}

void MultiLidarFrameAssembler::PromoteReady(bool force) {
    while (!frames_.empty()) {
        auto it = frames_.begin();
        if (!IsComplete(it->second) && !force && !IsExpired(it->first)) {
            break;
        }
        const int assembly_min_lidars = config_.adaptive_load.enabled
                                            ? config_.adaptive_load.tracking_min_lidars
                                            : config_.min_lidars;
        if (static_cast<int>(it->second.clouds.size()) >= assembly_min_lidars) {
            ready_frames_.push_back(Assemble(it->second));
            ++emitted_frame_count_;
        } else {
            ++insufficient_lidar_drop_count_;
        }
        last_emitted_bucket_ = it->first;
        frames_.erase(it);
    }
}

bool MultiLidarFrameAssembler::PopReady(FusedLidarFrame& frame) {
    if (ready_frames_.empty()) return false;
    frame = std::move(ready_frames_.front());
    ready_frames_.pop_front();
    return true;
}

void MultiLidarFrameAssembler::Flush() { PromoteReady(true); }

FusedLidarFrame MultiLidarFrameAssembler::Assemble(const PartialFrame& frame) const {
    FusedLidarFrame fused;
    fused.cloud.reset(new PointCloudType);
    auto& stats = fused.stats;
    stats.begin_time = std::numeric_limits<double>::infinity();
    std::size_t total_points = 0;
    for (const auto& [id, cloud] : frame.clouds) {
        stats.begin_time = std::min(stats.begin_time, frame.timestamps.at(id));
        stats.present_lidar_ids.push_back(id);
        stats.points_by_lidar[id] = cloud ? cloud->size() : 0;
        total_points += stats.points_by_lidar[id];
    }
    for (const int id : expected_lidar_ids_) {
        if (frame.clouds.count(id) == 0) {
            stats.missing_lidar_ids.push_back(id);
            stats.points_by_lidar[id] = 0;
        }
    }
    stats.partial = !stats.missing_lidar_ids.empty();
    fused.cloud->reserve(total_points);
    double max_relative_ms = 0.0;
    for (const auto& [id, cloud] : frame.clouds) {
        if (!cloud) continue;
        const auto* sensor = config_.FindLidar(id);
        const double header_offset_ms = (frame.timestamps.at(id) - stats.begin_time) * 1e3;
        for (const auto& point : cloud->points) {
            PointType transformed = point;
            const Vec3d p_primary = sensor->R_lidar_to_primary * point.getVector3fMap().cast<double>() +
                                    sensor->t_lidar_to_primary;
            transformed.x = p_primary.x();
            transformed.y = p_primary.y();
            transformed.z = p_primary.z();
            transformed.time = std::max(0.0, point.time + header_offset_ms);
            transformed.lidar_id = static_cast<std::uint8_t>(id);
            max_relative_ms = std::max(max_relative_ms, transformed.time);
            fused.cloud->push_back(transformed);
        }
    }
    std::sort(fused.cloud->points.begin(), fused.cloud->points.end(),
              [](const PointType& lhs, const PointType& rhs) { return lhs.time < rhs.time; });
    fused.cloud->width = fused.cloud->size();
    fused.cloud->height = 1;
    fused.cloud->is_dense = false;
    fused.cloud->header.stamp = static_cast<std::uint64_t>(std::llround(stats.begin_time * 1e9));
    stats.merged_points = fused.cloud->size();
    stats.end_time = stats.begin_time + max_relative_ms * 1e-3;
    return fused;
}

}  // namespace lightning

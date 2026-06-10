//
// Created by user on 2026/3/18.
//

#ifndef LIGHTNING_IMU_FILTER_H
#define LIGHTNING_IMU_FILTER_H

#include "common/eigen_types.h"
#include "common/imu.h"

#include <deque>

namespace lightning {

/**
 * @brief IMU角速度滤波器。
 *
 * 当前实现只滤波 IMU::angular_velocity，不修改 linear_acceleration 和 timestamp。
 * Filter() 内部按顺序执行：
 * 1. 保存各轴角速度历史；
 * 2. 基于历史窗口检测并替换明显毛刺；
 * 3. 对角速度做中值滤波；
 * 4. 对角速度做移动平均；
 * 5. 根据上一帧滤波结果限制角速度变化率。
 *
 * 这个类带有历史缓存和统计量，应按时间顺序连续调用。若数据流重启，建议重新构造对象
 * 或增加显式Reset逻辑，避免旧统计量影响新数据。
 */
class IMUFilter {
   private:
    /// 滤波参数。窗口越大越平滑，但会引入更明显的滞后。
    struct Config {
        int median_window_size = 5;    // 中值滤波窗口大小，必须为不小于3的奇数
        int moving_avg_window = 3;     // 移动平均窗口，取最近N个角速度求平均
        double rate_limit = 3.0;       // 角速度变化率限制，单位 rad/s^2
        double spike_threshold = 3.0;  // 毛刺检测阈值，按角速度标准差倍数设置
        bool enable_adaptive = true;   // 是否根据在线统计的标准差调整毛刺阈值
    } config_;

    // 角速度历史缓存。滤波按轴独立处理，因此这里只保存三轴标量历史。
    std::deque<IMU> buffer_;             // 预留的整条IMU缓存，当前实现未实际写入
    std::deque<double> gyro_x_history_;  // x轴角速度历史
    std::deque<double> gyro_y_history_;  // y轴角速度历史
    std::deque<double> gyro_z_history_;  // z轴角速度历史

    IMU prev_filtered_;  // 上一条滤波后的IMU，用于计算dt并限制角速度突变

    // 在线统计信息，用于自适应毛刺检测。这里的std是绝对偏差的指数滑动估计，并非严格方差开根号。
    double gyro_mean_[3] = {0};  // 三轴角速度指数滑动均值
    double gyro_std_[3] = {0};   // 三轴角速度波动尺度估计
    int sample_count_ = 0;       // 已处理样本数，前期统计未稳定时不做毛刺替换

   public:
    /// 构造时把上一帧时间戳置为无效值，第一条数据不会进行变化率限制。
    IMUFilter() { prev_filtered_.timestamp = -1; }

    /// 设置中值滤波窗口。只接受不小于3的奇数，避免中位数定义不稳定。
    void SetMedianWindowSize(int size) {
        if (size >= 3 && size % 2 == 1) {
            config_.median_window_size = size;
        }
    }

    /// 设置角速度变化率上限。传入负数时取绝对值。
    void SetRateLimit(double limit) { config_.rate_limit = std::abs(limit); }

    /// 设置毛刺检测阈值倍数。值越小越容易把突变判定为毛刺。
    void SetSpikeThreshold(double threshold) { config_.spike_threshold = threshold; }

    /**
     * @brief 对单条IMU量测做角速度滤波。
     *
     * 返回值以raw_data为基础拷贝，因此加速度、时间戳等字段保持原始值；
     * 函数只覆盖angular_velocity三轴。
     *
     * @param raw_data 原始IMU量测。
     * @return 滤波后的IMU量测。
     */
    IMU Filter(const IMU &raw_data) {
        IMU filtered = raw_data;

        // 先写入当前原始角速度历史，后续滤波窗口会包含当前样本。
        updateBuffer(raw_data);

        // 三个角速度轴独立处理，避免某一轴异常影响其他轴。
        filtered.angular_velocity.x() = processAxis(raw_data.angular_velocity.x(), gyro_x_history_, 0);
        filtered.angular_velocity.y() = processAxis(raw_data.angular_velocity.y(), gyro_y_history_, 1);
        filtered.angular_velocity.z() = processAxis(raw_data.angular_velocity.z(), gyro_z_history_, 2);

        // 变化率限制使用滤波后的上一帧作为参考，只在dt合理时启用。
        if (prev_filtered_.timestamp > 0) {
            double dt = raw_data.timestamp - prev_filtered_.timestamp;
            if (dt > 0 && dt < 0.1) {  // 合理的dt范围
                filtered.angular_velocity.x() =
                    rateLimit(filtered.angular_velocity.x(), prev_filtered_.angular_velocity.x(), dt);
                filtered.angular_velocity.y() =
                    rateLimit(filtered.angular_velocity.y(), prev_filtered_.angular_velocity.y(), dt);
                filtered.angular_velocity.z() =
                    rateLimit(filtered.angular_velocity.z(), prev_filtered_.angular_velocity.z(), dt);
            }
        }

        // 用最终输出更新统计量，使毛刺阈值跟随滤波后的正常角速度分布。
        updateStatistics(filtered);

        prev_filtered_ = filtered;
        return filtered;
    }

   private:
    /// 对单个角速度轴执行毛刺替换、中值滤波和移动平均。
    double processAxis(double raw_value, std::deque<double> &history, int axis_idx) {
        double filtered = raw_value;

        // 统计量需要一定样本数才能稳定，启动初期只做窗口滤波，不做毛刺替换。
        if (history.size() >= config_.median_window_size && sample_count_ > 100) {
            filtered = detectAndRemoveSpike(raw_value, history, axis_idx);
        }

        // 中值滤波优先去掉孤立离群点。
        filtered = medianFilter(filtered, history);

        // 移动平均进一步平滑高频抖动，但会带来少量滞后。
        filtered = movingAverage(filtered, history);

        return filtered;
    }

    /// 写入原始角速度历史，并把历史长度裁剪到所有滤波步骤需要的最大窗口。
    void updateBuffer(const IMU &data) {
        gyro_x_history_.push_back(data.angular_velocity.x());
        gyro_y_history_.push_back(data.angular_velocity.y());
        gyro_z_history_.push_back(data.angular_velocity.z());

        // 至少保留10个样本，给统计和调试留一点缓冲；滤波只使用各自窗口大小。
        int max_history = std::max({config_.median_window_size, config_.moving_avg_window, 10});
        while (gyro_x_history_.size() > max_history) {
            gyro_x_history_.pop_front();
            gyro_y_history_.pop_front();
            gyro_z_history_.pop_front();
        }
    }

    /// 如果当前值相对近期中位数偏离过大，则认为是毛刺并替换为中位数。
    double detectAndRemoveSpike(double value, std::deque<double> &history, int axis_idx) {
        if (history.size() < config_.median_window_size) {
            return value;
        }

        // 使用最近median_window_size个原始角速度作为局部正常水平。
        std::vector<double> window(history.end() - config_.median_window_size, history.end());
        std::nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
        double median = window[window.size() / 2];

        // 与局部中位数比较，比直接和均值比较更不容易被单个异常值拖偏。
        double diff = std::abs(value - median);

        // 自适应阈值随角速度波动尺度变化；波动越大，允许的瞬时偏差越大。
        double threshold = config_.spike_threshold * gyro_std_[axis_idx];
        if (config_.enable_adaptive && gyro_std_[axis_idx] > 0) {
            threshold = std::max(threshold, config_.spike_threshold * 0.5);
        }

        // 判断为毛刺时用局部中位数替换，避免单点跳变进入后续积分。
        if (diff > threshold) {
            LOG(INFO) << "find imu spike: " << diff << ", " << threshold;
            return median;
        }

        return value;
    }

    /// 返回最近窗口内的中位数。窗口不足时直接返回输入值。
    double medianFilter(double value, std::deque<double> &history) {
        if (history.size() < config_.median_window_size) {
            return value;
        }

        std::vector<double> window(history.end() - config_.median_window_size, history.end());
        std::nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
        return window[window.size() / 2];
    }

    /// 返回最近窗口内的均值。窗口不足时直接返回输入值。
    double movingAverage(double value, std::deque<double> &history) {
        if (history.size() < config_.moving_avg_window) {
            return value;
        }

        double sum = 0;
        auto it = history.end() - config_.moving_avg_window;
        for (; it != history.end(); ++it) {
            sum += *it;
        }
        return sum / config_.moving_avg_window;
    }

    /// 限制相邻两帧滤波角速度的变化量，避免滤波输出出现过快跳变。
    double rateLimit(double current, double previous, double dt) {
        double max_change = config_.rate_limit * dt;
        double diff = current - previous;

        if (std::abs(diff) > max_change) {
            return previous + (diff > 0 ? max_change : -max_change);
        }
        return current;
    }

    /// 用指数滑动平均更新角速度均值和波动尺度，供后续毛刺检测使用。
    void updateStatistics(const IMU &data) {
        const double alpha = 0.01;  // 指数移动平均系数

        if (sample_count_ == 0) {
            gyro_mean_[0] = data.angular_velocity[0];
            gyro_mean_[1] = data.angular_velocity[1];
            gyro_mean_[2] = data.angular_velocity[2];
            gyro_std_[0] = 0.1;  // 初始值
            gyro_std_[1] = 0.1;
            gyro_std_[2] = 0.1;
        } else {
            // 更新均值
            double prev_mean[3] = {gyro_mean_[0], gyro_mean_[1], gyro_mean_[2]};
            gyro_mean_[0] = (1 - alpha) * gyro_mean_[0] + alpha * data.angular_velocity[0];
            gyro_mean_[1] = (1 - alpha) * gyro_mean_[1] + alpha * data.angular_velocity[1];
            gyro_mean_[2] = (1 - alpha) * gyro_mean_[2] + alpha * data.angular_velocity[2];

            // 更新标准差
            gyro_std_[0] = (1 - alpha) * gyro_std_[0] + alpha * std::abs(data.angular_velocity[0] - gyro_mean_[0]);
            gyro_std_[1] = (1 - alpha) * gyro_std_[1] + alpha * std::abs(data.angular_velocity[1] - gyro_mean_[1]);
            gyro_std_[2] = (1 - alpha) * gyro_std_[2] + alpha * std::abs(data.angular_velocity[2] - gyro_mean_[2]);
        }

        sample_count_++;
    }
};
}  // namespace lightning

#endif  // LIGHTNING_IMU_FILTER_H

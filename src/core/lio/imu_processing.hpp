#pragma once

#ifndef FASTER_LIO_IMU_PROCESSING_H
#define FASTER_LIO_IMU_PROCESSING_H

#include <glog/logging.h>
#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "common/eigen_types.h"
#include "common/measure_group.h"
#include "common/options.h"
#include "common/point_def.h"
#include "core/lightning_math.hpp"
#include "core/lio/eskf.hpp"
#include "core/lio/imu_filter.h"
#include "core/lio/pose6d.h"
#include "utils/timer.h"

namespace lightning {

/**
 * @brief IMU预处理和点云运动畸变补偿模块。
 *
 * ImuProcess位于Lidar量测和ESKF之间，主要负责三件事：
 * 1. 在系统启动阶段统计静止IMU均值和噪声，用于初始化重力、陀螺零偏和过程噪声；
 * 2. 对每帧Lidar时间段内的IMU量测做前向积分，驱动ESKF预测；
 * 3. 根据积分出的IMU位姿序列，把一帧点云内不同时刻采集的点补偿到扫描结束时刻。
 *
 * 注意：类中保留了上一帧末尾IMU、上一帧Lidar结束时刻等跨帧状态，因此同一个对象应按时间顺序处理连续帧。
 */
class ImuProcess {
   public:
    struct InitializationOptions {
        double min_duration = 0.05;
        int min_samples = 20;
        double max_mean_gyro_norm = std::numeric_limits<double>::infinity();
        double max_gyro_std = std::numeric_limits<double>::infinity();
        double max_acc_std = std::numeric_limits<double>::infinity();
        double min_mean_acc_norm = 0.0;
        double max_mean_acc_norm = std::numeric_limits<double>::infinity();
        double initial_yaw_deg = 0.0;
    };

    /// Eigen固定大小矩阵成员需要对齐分配，避免在容器或new对象时出现内存对齐问题。
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ImuProcess();
    ~ImuProcess();

    /// 重置IMU初始化状态和跨帧缓存，通常在重定位或重新开始处理数据时调用。
    void Reset();
    /// 设置Lidar到IMU的外参，transl和rot用于点云去畸变时在Lidar系与IMU系之间变换。
    void SetExtrinsic(const Vec3d &transl, const Mat3d &rot);
    /// 设置陀螺仪噪声缩放参数，初始化阶段会结合量测统计得到实际过程噪声。
    void SetGyrCov(const Vec3d &scaler);
    /// 设置加速度计噪声缩放参数，初始化阶段会结合量测统计得到实际过程噪声。
    void SetAccCov(const Vec3d &scaler);
    /// 设置陀螺仪零偏随机游走噪声。
    void SetGyrBiasCov(const Vec3d &b_g);
    /// 设置加速度计零偏随机游走噪声。
    void SetAccBiasCov(const Vec3d &b_a);
    void SetInitializationOptions(const InitializationOptions &options) { init_options_ = options; }

    /**
     * @brief 处理一组同步后的Lidar和IMU量测。
     *
     * 如果IMU尚未初始化，Process()会先累积若干帧IMU完成初始化；
     * 初始化完成后，会调用UndistortPcl()做ESKF预测和点云去畸变。
     *
     * @param meas 当前Lidar帧及其时间范围内的IMU量测。
     * @param kf_state ESKF状态，函数内部会修改其状态和协方差。
     * @param scan 输出的去畸变点云。
     */
    void Process(const MeasureGroup &meas, ESKF &kf_state, CloudPtr &scan);

    /// 返回IMU初始化是否已经完成。imu_need_init_为false表示可进入正常预测/去畸变流程。
    bool IsIMUInited() const { return imu_need_init_ == false; }
    /// 开关IMU量测滤波。开启时，去畸变前会对当前帧IMU拷贝做滤波。
    void SetUseIMUFilter(bool b) { use_imu_filter_ = b; }

    /// 获取初始化阶段估计出的加速度均值模长，可用于检查静止初始化是否接近重力加速度。
    double GetMeanAccNorm() const { return mean_acc_.norm(); }
    std::size_t GetInitializationSampleCount() const { return init_samples_.size(); }
    const SO3 &GetInitialRotation() const { return initial_rotation_; }
    Vec3d ScaleAccelerationForPrediction(const Vec3d &acc) const { return acc * acc_scale_factor_; }

    // 这些噪声参数会在IMU初始化和ESKF预测中使用，保持public是为了兼容原框架的配置方式。
    Eigen::Matrix<double, 12, 12> Q_;  // ESKF过程噪声协方差矩阵，顺序需与NavState::df_dw()的噪声定义一致
    Vec3d cov_acc_;                    // 加速度计噪声协方差（从IMU初始化中估计）
    Vec3d cov_gyr_;                    // 陀螺仪噪声协方差（从IMU初始化中估计）
    Vec3d cov_acc_scale_;              // 加速度计噪声协方差缩放因子（配置参数）
    Vec3d cov_gyr_scale_;              // 陀螺仪噪声协方差缩放因子（配置参数）
    Vec3d cov_bias_gyr_;               // 陀螺仪偏置不稳定性的协方差（配置参数）
    Vec3d cov_bias_acc_;               // 加速度计偏置不稳定性的协方差（配置参数）

   private:
    /// 静止初始化：估计加速度/角速度均值、噪声方差、初始重力方向和陀螺零偏。
    void IMUInit(const MeasureGroup &meas, ESKF &kf_state, int &N);
    /// 正常运行阶段：用IMU积分预测状态，并把点云补偿到扫描结束时刻。
    void UndistortPcl(const MeasureGroup &meas, ESKF &kf_state, CloudPtr &pcl_out);

    bool InitializationReady() const;

    PointCloudType::Ptr cur_pcl_un_ = nullptr;  // 当前帧去畸变点云缓存，Reset时重新分配
    lightning::IMUPtr last_imu_ = nullptr;      // 上一帧最后一条IMU，用于和当前帧第一条IMU形成连续积分区间
    std::deque<lightning::IMUPtr> imu_queue_;   // IMU队列缓存，当前实现主要在Reset中维护

    std::vector<Pose6D> imu_pose_;            // 当前Lidar帧内的IMU位姿序列，用于按时间段补偿点云
    Mat3d R_lidar_imu_ = Mat3d ::Identity();  // Lidar到IMU的旋转外参，点云先转到IMU系再做运动补偿
    Vec3d t_lidar_mu_ = Vec3d ::Zero();       // Lidar到IMU的平移外参，变量名保留原拼写
    Vec3d mean_acc_ = Vec3d::Zero();          // 初始化阶段累积的加速度计均值，用于估计重力方向
    Vec3d mean_gyr_ = Vec3d::Zero();          // 初始化阶段累积的陀螺仪均值，用作初始陀螺零偏
    Vec3d angvel_last_ = Vec3d ::Zero();      // 上一积分区间去零偏后的角速度，保存给下一帧起点使用
    Vec3d acc_s_last_ = Vec3d ::Zero();       // 上一积分区间转到世界系并去重力后的加速度
    double acc_scale_factor_ = 1.0;           // IMU加速度计缩放因子，初始化后用于修正加速度量测尺度

    double last_lidar_end_time_ = 0;  // 上一帧Lidar扫描结束时刻，用于避免跨帧IMU重复积分
    int init_iter_num_ = 1;           // IMU初始化累计次数，达到max_init_count_后结束初始化
    bool b_first_frame_ = true;       // 标记是否为第一帧，第一帧通常只用于建立初始时间和缓存
    bool imu_need_init_ = true;       // 是否仍处于IMU初始化阶段

    bool use_imu_filter_ = true;  // 是否对参与积分的IMU量测做滤波
    IMUFilter filter_;            // IMU滤波器实例，仅处理当前帧IMU拷贝
    InitializationOptions init_options_;
    std::deque<lightning::IMUPtr> init_samples_;
    SO3 initial_rotation_;
};

inline ImuProcess::ImuProcess() : b_first_frame_(true), imu_need_init_(true) {
    init_iter_num_ = 1;
    Q_.setZero();
    Q_.diagonal() << 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5, 0.0, 0.0, 0.0;
    cov_acc_ = Vec3d(0.1, 0.1, 0.1);
    cov_gyr_ = Vec3d(0.1, 0.1, 0.1);
    cov_bias_gyr_ = Vec3d(0.0001, 0.0001, 0.0001);
    cov_bias_acc_ = Vec3d(0.0001, 0.0001, 0.0001);
    mean_acc_ = Vec3d(0, 0, -1.0);
    mean_gyr_ = Vec3d(0, 0, 0);
    last_imu_.reset(new lightning::IMU());
}

inline ImuProcess::~ImuProcess() {}

inline void ImuProcess::Reset() {
    mean_acc_ = Vec3d(0, 0, -1.0);
    mean_gyr_ = Vec3d(0, 0, 0);
    angvel_last_.setZero();
    acc_scale_factor_ = 1.0;  // 重置IMU缩放因子

    imu_need_init_ = true;
    init_iter_num_ = 1;
    imu_queue_.clear();
    init_samples_.clear();
    imu_pose_.clear();
    initial_rotation_ = SO3();
    last_imu_.reset(new lightning::IMU());
    cur_pcl_un_.reset(new PointCloudType());
}

inline void ImuProcess::SetExtrinsic(const Vec3d &transl, const Mat3d &rot) {
    t_lidar_mu_ = transl;
    R_lidar_imu_ = rot;
}

inline void ImuProcess::SetGyrCov(const Vec3d &scaler) { cov_gyr_scale_ = scaler; }

inline void ImuProcess::SetAccCov(const Vec3d &scaler) { cov_acc_scale_ = scaler; }

inline void ImuProcess::SetGyrBiasCov(const Vec3d &b_g) { cov_bias_gyr_ = b_g; }

inline void ImuProcess::SetAccBiasCov(const Vec3d &b_a) { cov_bias_acc_ = b_a; }

/**
 * @brief 使用静止阶段IMU量测初始化ESKF状态和IMU噪声统计。
 *
 * @param meas 当前Lidar帧对应的IMU量测集合。
 * @param kf_state 待初始化的ESKF状态，函数会写入重力方向、陀螺零偏和初始协方差。
 * @param N 初始化阶段累计使用的IMU样本数，函数内会随量测递增。
 *
 * @details 初始化阶段假设载体基本静止：加速度均值主要反映重力方向，
 *          角速度均值作为初始陀螺零偏。函数同时用增量形式估计加速度计
 *          和陀螺仪噪声方差，避免保存全部历史IMU样本。
 */
inline void ImuProcess::IMUInit(const MeasureGroup &meas, ESKF &kf_state, int &N) {
    if (b_first_frame_) {
        Reset();
        b_first_frame_ = false;
    }
    for (const auto &imu : meas.imu_) init_samples_.push_back(imu);
    if (init_samples_.empty()) return;
    const double retention_duration = std::max(1.0, 5.0 * init_options_.min_duration);
    const double oldest_allowed = init_samples_.back()->timestamp - retention_duration;
    while (init_samples_.size() > 1 && init_samples_[1]->timestamp <= oldest_allowed) init_samples_.pop_front();

    mean_acc_.setZero();
    mean_gyr_.setZero();
    for (const auto &imu : init_samples_) {
        mean_acc_ += imu->linear_acceleration;
        mean_gyr_ += imu->angular_velocity;
    }
    mean_acc_ /= static_cast<double>(init_samples_.size());
    mean_gyr_ /= static_cast<double>(init_samples_.size());
    cov_acc_.setZero();
    cov_gyr_.setZero();
    for (const auto &imu : init_samples_) {
        const Vec3d da = imu->linear_acceleration - mean_acc_;
        const Vec3d dg = imu->angular_velocity - mean_gyr_;
        cov_acc_ += da.cwiseProduct(da);
        cov_gyr_ += dg.cwiseProduct(dg);
    }
    const double denominator = std::max<std::size_t>(1, init_samples_.size() - 1);
    cov_acc_ /= denominator;
    cov_gyr_ /= denominator;
    N = static_cast<int>(init_samples_.size());

    // 将静止统计结果写入ESKF：用加速度均值对齐初始姿态，重力固定在世界系-z，角速度均值为陀螺零偏。
    auto init_state = kf_state.GetX();
    init_state.timestamp_ = meas.imu_.back()->timestamp;
    const double mean_acc_norm = mean_acc_.norm();
    if (mean_acc_norm > 1e-6) {
        const Vec3d acc_dir = mean_acc_ / mean_acc_norm;
        const SO3 gravity_alignment(Quatd::FromTwoVectors(acc_dir, Vec3d::UnitZ()).normalized());
        const double initial_yaw_rad = init_options_.initial_yaw_deg * M_PI / 180.0;
        const SO3 heading_alignment = SO3::exp(Vec3d(0.0, 0.0, initial_yaw_rad));
        init_state.rot_ = heading_alignment * gravity_alignment;
    }
    initial_rotation_ = init_state.rot_;

    init_state.grav_ = Vec3d(0.0, 0.0, -G_m_s2);
    init_state.bg_ = mean_gyr_;
    init_state.ba_ = Vec3d::Zero();
    kf_state.ChangeX(init_state);

    // 计算并缓存IMU缩放因子
    // meas_acc_scale_ = G_m_s2 / mean_acc_.norm();  // TODO，解决冲突

    // Keep the initial covariance tight; the full ESKF lets lidar pose residuals update
    // velocity/bias/gravity through cross-covariance, so a unit prior can inject huge velocity.
    auto init_P = kf_state.GetP();
    init_P.setIdentity();
    init_P *= 1e-4;
    init_P.block<NavState::kBlockDim, NavState::kBlockDim>(NavState::kBgIdx, NavState::kBgIdx) =
        0.0001 * Mat3d::Identity();
    init_P.block<NavState::kBlockDim, NavState::kBlockDim>(NavState::kBaIdx, NavState::kBaIdx) =
        0.001 * Mat3d::Identity();
    init_P.block<NavState::kBlockDim, NavState::kBlockDim>(NavState::kGravIdx, NavState::kGravIdx) =
        1e-8 * Mat3d::Identity();
    kf_state.ChangeP(init_P);

    // LOG(INFO) << "P diag: " << init_P.diagonal().transpose();

    // 缓存当前批次最后一条IMU，供后续点云去畸变形成跨帧积分区间。
    last_imu_ = meas.imu_.back();
}

inline bool ImuProcess::InitializationReady() const {
    if (init_samples_.size() < static_cast<std::size_t>(std::max(1, init_options_.min_samples))) return false;
    if (init_samples_.back()->timestamp - init_samples_.front()->timestamp < init_options_.min_duration) return false;
    const double gyro_std = std::sqrt(std::max(0.0, cov_gyr_.maxCoeff()));
    const double acc_std = std::sqrt(std::max(0.0, cov_acc_.maxCoeff()));
    const double mean_acc_norm = mean_acc_.norm();
    return mean_gyr_.norm() <= init_options_.max_mean_gyro_norm && gyro_std <= init_options_.max_gyro_std &&
           acc_std <= init_options_.max_acc_std && mean_acc_norm >= init_options_.min_mean_acc_norm &&
           mean_acc_norm <= init_options_.max_mean_acc_norm;
}

/**
 * @brief 使用IMU前向积分对当前Lidar点云做运动畸变补偿。
 *
 * 本函数做两件事：
 * 1. 用IMU量测把ESKF状态从上一帧末尾积分到当前Lidar扫描结束时刻；
 * 2. 根据积分过程中保存的IMU位姿，把每个点补偿到扫描结束时刻，减小运动畸变。
 *
 * @param meas 当前帧Lidar点云和对应时间段内的IMU量测。
 * @param kf_state ESKF状态，会在函数内预测到当前Lidar扫描结束时刻。
 * @param pcl_out 输出的去畸变点云，所有点被补偿到扫描结束时刻的Lidar坐标系。
 */
inline void ImuProcess::UndistortPcl(const MeasureGroup &meas, ESKF &kf_state, CloudPtr &pcl_out) {
    /*** add the imu_ of the last frame-tail to the of current frame-head ***/
    // 这里把上一帧最后一个IMU插到当前帧IMU序列头部，是为了形成跨帧的第一个积分区间。
    // 否则当前帧起点到第一条IMU之间会缺少约束，点云开头附近的运动补偿容易不连续。
    auto v_imu = meas.imu_;
    v_imu.push_front(last_imu_);
    const double &imu_end_time = v_imu.back()->timestamp;

    // 点云时间戳使用“本帧扫描开始时刻”为零点，后面保存IMU位姿时也统一转成该相对时间。
    const double &pcl_beg_time = meas.lidar_begin_time_;
    const double &pcl_end_time = meas.lidar_end_time_;

    /*** Initialize IMU pose ***/
    // imu_pose_保存一串按时间排列的IMU状态，后面会用它把点云按时间分段补偿。
    // 第一项是扫描开始处的参考状态：offset_time=0，状态来自当前ESKF。
    auto imu_state = kf_state.GetX();
    imu_pose_.clear();
    imu_pose_.emplace_back(0.0, acc_s_last_, angvel_last_, imu_state.vel_, imu_state.pos_, imu_state.rot_.matrix());

    /*** forward propagation at each imu_ point ***/
    Vec3d angvel_avr, acc_avr, acc_imu, vel_imu, pos_imu;
    Mat3d R_imu;

    double dt = 0;
    Vec3d acc = Vec3d::Zero();
    Vec3d gyro = Vec3d::Zero();

    if (use_imu_filter_) {
        // 可选滤波只作用于本次用于积分的IMU拷贝，避免改写meas里的原始量测。
        for (auto &imu : v_imu) {
            auto imu_f = filter_.Filter(*imu);
            *imu = imu_f;
        }
    }

    for (auto it_imu = v_imu.begin(); it_imu < (v_imu.end() - 1); it_imu++) {
        // 每次取相邻两条IMU量测，构成一个积分区间 [head, tail]。
        // 后面用两端量测的平均值作为这个区间内的近似输入。
        auto &&head = *(it_imu);
        auto &&tail = *(it_imu + 1);

        if (tail->timestamp < last_lidar_end_time_) {
            // 这个区间完全落在上一帧Lidar结束之前，对当前帧补偿没有贡献。
            continue;
        }

        angvel_avr = .5 * (head->angular_velocity + tail->angular_velocity);
        acc_avr = .5 * (head->linear_acceleration + tail->linear_acceleration);

        // 先把加速度计量测缩放到标定尺度，再交给ESKF预测；角速度这里没有做额外缩放。
        acc_avr = ScaleAccelerationForPrediction(acc_avr);  // 使用缓存的缩放因子进行加速度计标定
        // 如果head早于上一帧Lidar结束时刻，只积分 [last_lidar_end_time_, tail] 这一段。
        // 这样可以避免把上一帧已经积分过的IMU时间重复计入当前帧。
        if (head->timestamp < last_lidar_end_time_) {
            dt = tail->timestamp - last_lidar_end_time_;
        } else {
            dt = tail->timestamp - head->timestamp;
        }

        acc = acc_avr;
        gyro = angvel_avr;

        if (dt <= 0.0) {
            continue;
        }

        const double warn_dt = std::max(0.1, 1.5 * static_cast<double>(lo::lidar_time_interval));
        if (dt > warn_dt) {
            LOG(WARNING) << "propagate over long imu interval: " << dt << ", lidar: " << pcl_beg_time << " -> "
                         << pcl_end_time;
        }
        // Q_是ESKF预测用的过程噪声。这里把初始化阶段估计出的陀螺仪、加速度计和零偏噪声写进去。
        // TODO.在IMU初始化完成后，其实这个方差就不会变
        Q_.block<3, 3>(0, 0).diagonal() = cov_gyr_;
        Q_.block<3, 3>(3, 3).diagonal() = cov_acc_;
        Q_.block<3, 3>(6, 6).diagonal() = cov_bias_gyr_;
        Q_.block<3, 3>(9, 9).diagonal() = cov_bias_acc_;
        kf_state.Predict(dt, Q_, gyro, acc);

        // LOG(INFO) << "gyro: " << gyro.transpose() << ", dt: " << dt;

        // LOG(INFO) << "acc: " << acc.transpose() << " grav: " << kf_state.GetX().grav_.norm()
        //           << ", vel: " << kf_state.GetX().vel_.transpose() << ", dt: " << dt;

        /* save the poses at each IMU measurements */
        imu_state = kf_state.GetX();
        // 保存去零偏后的角速度，以及转到世界系并扣除重力后的加速度。
        // 这些量后面会用于把区间内的点外推到其真实采集时刻。
        angvel_last_ = angvel_avr - imu_state.bg_;
        acc_s_last_ = imu_state.rot_ * (acc_avr - imu_state.ba_);
        for (int i = 0; i < 3; i++) {
            acc_s_last_[i] += imu_state.grav_[i];  // 去除重力向量
        }
        // 思考：为什么用的平均加速度/角速度，而时间戳不是用的中间时间戳，而是尾端时间戳
        // 1.ESKF预测使用的是前向欧拉积分：kf_state.Predict(dt, Q_, gyro, acc)
        // 2.在 [head_timestamp, tail_timestamp] 区间内，使用平均的输入值进行积分
        // 3.但积分结果对应的是区间末端的状态
        double &&offs_t = tail->timestamp - pcl_beg_time;
        imu_pose_.emplace_back(
            Pose6D(offs_t, acc_s_last_, angvel_last_, imu_state.vel_, imu_state.pos_, imu_state.rot_.matrix()));
    }

    /*** 计算帧结束时刻的位姿预测，确保IMU积分覆盖整个点云扫描周期 ***/
    // 前面的循环只积分到了最后一条IMU的时间 imu_end_time。
    // 点云补偿的目标帧是当前Lidar扫描结束时刻 pcl_end_time，所以还要把状态对齐到这个时刻。
    // note用于兼容两种情况：
    // - IMU早于点云结束：继续预测到点云结束；
    // - IMU晚于点云结束：按当前写法取两者时间差的非负值，沿用原框架的预测逻辑。
    double note = pcl_end_time > imu_end_time ? 1.0 : -1.0;
    dt = note * (pcl_end_time - imu_end_time);  // 正向或反向预测到点云结束时刻
    kf_state.Predict(dt, Q_, gyro, acc);

    imu_state = kf_state.GetX();  // 获取最终的位姿状态（作为运动补偿的参考帧）
    last_imu_ = meas.imu_.back();
    last_lidar_end_time_ = pcl_end_time;

    /*** 按时间戳对点云进行排序，便于后续按时间区间进行运动补偿 ***/
    pcl_out = meas.scan_;
    std::sort(pcl_out->points.begin(), pcl_out->points.end(),
              [](const PointType &p1, const PointType &p2) { return p1.time < p2.time; });

    /*** 开始点云运动畸变补偿（从后向前传播）***/
    // 点云已经按时间升序排序。这里从最后一个点往前处理，和“补偿到扫描结束时刻”的方向一致。
    // 外层循环按相邻IMU位姿划分时间区间，内层循环处理落在该区间内的点。
    if (pcl_out->empty()) {
        return;
    }
    // 从后向前处理点云，将所有点补偿到扫描结束时刻
    auto it_pcl = pcl_out->points.end() - 1;
    for (auto it_kp = imu_pose_.end() - 1; it_kp != imu_pose_.begin(); it_kp--) {
        auto head = it_kp - 1;  // 当前时间区间起始位姿
        auto tail = it_kp;      // 当前时间区间结束位姿
        // 对落在 [head, tail] 内的点，用区间起点状态加上该区间的速度/加速度外推点的采集位姿。
        // 获取区间起始时刻的运动状态
        R_imu = (head->rot);
        vel_imu = (head->vel);
        pos_imu = (head->pos);
        acc_imu = (tail->acc);     // 区间内的加速度（平均）
        angvel_avr = (tail->gyr);  // 区间内的角速度（平均）
        // 处理当前时间区间内的所有点云点
        for (; it_pcl->time / double(1000) > head->offset_time && it_pcl != pcl_out->points.begin(); it_pcl--) {
            // 点的time通常是毫秒，这里转成秒，并得到它相对当前IMU区间起点head的时间差。
            dt = it_pcl->time / double(1000) - head->offset_time;

            /// dt 有时候存在非法数据
            if (dt < 0 || dt > lo::lidar_time_interval) {
                // LOG(WARNING) << "find abnormal dt in cloud: " << dt;
                continue;
            }

            /* Transform to the 'end' frame, using only the rotation
             * Note: Compensation direction is INVERSE of Frame's moving direction
             * So if we want to compensate a point at timestamp-i to the frame-e
             * p_compensate = R_imu_e ^ T * (R_i * P_i + T_ei) where T_ei is represented in global frame */
            // [gj-2025-11-26] 修正：使用 0.5*dt 以匹配 exp 的 2*scale*|vec| 定义
            // R_i表示该点采集时刻的IMU姿态：从区间起点姿态R_imu按角速度积分dt得到。
            Mat3d R_i(R_imu * math::exp(angvel_avr, 0.5 * dt).matrix());  // 计算点采集时刻的IMU旋转矩阵
            // Mat3d R_i(R_imu * math::exp(angvel_avr, dt).matrix());

            Vec3d P_i(it_pcl->x, it_pcl->y, it_pcl->z);  // 点云原始位置
            // T_ei是“点采集时刻IMU位置”相对“扫描结束时刻IMU位置”的位移，表达在世界系。
            // 随后的p_compensate分三步：
            // 1. 用外参把原始Lidar点P_i转到IMU系；
            // 2. 用点采集时刻的IMU位姿变到世界系，并平移到扫描结束参考位置；
            // 3. 用扫描结束时刻的IMU位姿和外参逆变换，回到扫描结束时刻的Lidar系。
            Vec3d T_ei(pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt -
                       imu_state.pos_);  // 从点采集时刻到扫描结束时刻的平移向量
            Vec3d p_compensate = R_lidar_imu_.transpose() *
                                 (imu_state.rot_.inverse() * (R_i * (R_lidar_imu_ * P_i + t_lidar_mu_) + T_ei) -
                                  t_lidar_mu_);  // 执行运动补偿变换（将点补偿到扫描结束时刻）

            // 更新点云坐标
            it_pcl->x = p_compensate(0);
            it_pcl->y = p_compensate(1);
            it_pcl->z = p_compensate(2);

            // if (it_pcl == pcl_out->points.begin()) {
            //     break;
            // }
        }
    }
}

/**
 * @brief ImuProcess对外的主入口。
 *
 * 处理流程分两阶段：
 * 1. IMU初始化阶段：累积若干帧静止IMU，估计重力方向、陀螺零偏和噪声参数；
 * 2. 正常运行阶段：调用UndistortPcl()，用IMU预测ESKF并对当前点云去畸变。
 *
 * 初始化阶段不会输出去畸变点云，因为此时重力、零偏和噪声统计还不稳定。
 *
 * @param meas 当前Lidar帧和该帧时间范围内的IMU量测。
 * @param kf_state ESKF状态，初始化和正常预测都会修改该状态。
 * @param scan 输出点云，只有初始化完成后才会由UndistortPcl()写入。
 */
inline void ImuProcess::Process(const MeasureGroup &meas, ESKF &kf_state, CloudPtr &scan) {
    if (meas.imu_.empty()) {
        // 没有IMU时无法初始化或去畸变，直接跳过当前帧。
        return;
    }

    if (imu_need_init_) {
        // 初始化阶段每来一帧Lidar，就把该帧对应的IMU继续累积到均值/方差估计里。
        IMUInit(meas, kf_state, init_iter_num_);

        // IMUInit内部会递增init_iter_num_，但是否结束初始化由下面的阈值判断统一决定。
        imu_need_init_ = true;

        // 缓存当前帧最后一条IMU，初始化结束后的第一帧去畸变仍需要跨帧连续积分。
        last_imu_ = meas.imu_.back();

        auto imu_state = kf_state.GetX();
        if (InitializationReady()) {
            // cov_acc_ *= pow(meas_acc_scale_, 2);  // 使用缓存的缩放因子，方差则需要平方
            // 初始化累计足够多帧后，切换到正常预测/去畸变流程。
            imu_need_init_ = false;

            // 最终预测噪声使用外部配置的尺度参数，避免静止初始化统计值过小导致滤波器过度自信。
            // 这发生在 IMU 初始化累计次数超过 max_init_count_ 后。
            // 前面 IMUInit() 里在线估计出来的 cov_acc_ / cov_gyr_，到了这里会被外部配置的 cov_acc_scale_ / cov_gyr_scale_ 覆盖掉。
            cov_acc_ = cov_acc_scale_;
            cov_gyr_ = cov_gyr_scale_;
            const double mean_acc_norm = mean_acc_.norm();

            // 根据静止时加速度均值的模长推断原始加速度单位：
            // - 接近1：认为输入单位是g，预测前需要乘9.81转成m/s^2；
            // - 接近9.81：认为输入已经是m/s^2，不再缩放；
            // - 其他范围：单位不可信，保持不缩放并给出告警。
            if (mean_acc_norm > 0.5 && mean_acc_norm < 1.5) {
                acc_scale_factor_ = G_m_s2;
            } else if (mean_acc_norm > 7.0 && mean_acc_norm < 12.0) {
                acc_scale_factor_ = 1.0;
            } else {
                acc_scale_factor_ = 1.0;
                LOG(WARNING) << "imu init mean acc norm is abnormal for unit inference: " << mean_acc_norm
                             << ", keep accelerometer scale unchanged";
            }

            LOG(INFO) << "imu init done, bg: " << imu_state.bg_.transpose() << ", grav: " << imu_state.grav_.transpose()
                      << ", acc scale: " << acc_scale_factor_ << ", cov_acc: " << cov_acc_.transpose()
                      << ", cov_gyr: " << cov_gyr_.transpose() << ", mean: " << mean_acc_.transpose() << ", "
                      << mean_gyr_.transpose();
        } else {
            // 初始化未完成时继续等待更多IMU样本；当前帧不做点云去畸变输出。
            const double duration = init_samples_.empty()
                                        ? 0.0
                                        : init_samples_.back()->timestamp - init_samples_.front()->timestamp;
            LOG(INFO) << "waiting for stationary imu init, samples=" << init_samples_.size()
                      << ", duration=" << duration << ", gyro_std="
                      << std::sqrt(std::max(0.0, cov_gyr_.maxCoeff())) << ", acc_std="
                      << std::sqrt(std::max(0.0, cov_acc_.maxCoeff()));
        }

        return;
    }

    // 正常运行阶段：用当前帧IMU推进ESKF，并把点云补偿到扫描结束时刻。
    Timer::Evaluate([&, this]() { UndistortPcl(meas, kf_state, scan); }, "Undistort Pcl");
}
}  // namespace lightning

#endif

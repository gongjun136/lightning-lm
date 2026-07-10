//
// Created by xiang on 2022/2/15.
//

#ifndef FUSION_ESKF_HPP
#define FUSION_ESKF_HPP

#include "common/eigen_types.h"
#include "common/nav_state.h"
#include "core/lio/anderson_acceleration.h"

namespace lightning {

/**
 * @brief LIO使用的误差状态卡尔曼滤波器。
 *
 * 这个类维护两套量：
 * 1. 名义状态 x_：真实参与积分和输出的导航状态，类型为NavState；
 * 2. 误差状态协方差 P_：描述名义状态附近小扰动的不确定性。
 *
 * Predict() 用IMU进行前向预测，Update() 用Lidar、轮速、GPS等观测修正状态。
 *
 * 相比原MTK写法，这里不再用宏拼接任意状态块，而是固定使用NavState。
 * 这样状态块索引和流形操作都集中在NavState里，便于调试和阅读。
 *
 * 当前主要状态块包括：pos / rot / vel / bg，以及NavState中定义的其他固定块。
 */
class ESKF {
   public:
    static constexpr int process_noise_dim_ = 12;     // IMU预测噪声维度，通常包含陀螺仪、加速度计及零偏噪声
    static constexpr int pose_obs_dim_ = 6;           // 位姿观测维度，前三维平移，后三维旋转
    static constexpr int state_dim_ = NavState::dim;  // 误差状态维度，由NavState统一定义
    using StateVecType = NavState::VectState;         // 误差状态向量类型
    using CovType = Eigen::Matrix<double, state_dim_, state_dim_>;  // 误差状态协方差矩阵
    using ProcessNoiseType = Eigen::Matrix<double, process_noise_dim_, process_noise_dim_>;

    /**
     * @brief Update()支持的观测来源。
     *
     * 枚举值只负责选择观测函数。观测函数本身通过Options传入，并把线性化后的信息写入
     * CustomObservationModel。
     */
    enum class ObsType {
        LIDAR,                  // 开源版本只有Lidar
        WHEEL_SPEED,            // 单独的轮速观测
        WHEEL_SPEED_AND_LIDAR,  // 轮速+Lidar
        ACC_AS_GRAVITY,         // 重力作为加计观测量
        GPS,                    // GPS/RTk 六自由度位姿
        BIAS,
    };

   public:
    /**
     * @brief 构造ESKF。
     *
     * @param x 初始名义状态。
     * @param P 初始误差状态协方差。
     * @param use_aa 是否启用Anderson Acceleration加速迭代更新。
     */
    explicit ESKF(const NavState& x = NavState(), const CovType& P = CovType::Identity(), bool use_aa = true)
        : x_(x), P_(P), use_aa_(use_aa) {}

    ~ESKF() {}

    /**
     * @brief 一次观测线性化后的结果。
     *
     * ESKF::Update()本身不关心观测来自Lidar、轮速还是GPS，只要求外部观测函数把当前
     * 线性化点处的观测信息填到这个结构体里。以Lidar点云为例，LaserMapping::ObsModel()
     * 会完成数据关联、残差计算和雅可比计算，然后把所有有效点约束累加为信息形式
     * H^T H和H^T r。
     *
     * 这样做的原因是点云残差维度随匹配点数变化，直接传完整H会很大；传H^T H和H^T r
     * 可以先把大量点约束压缩成6维位姿信息，
     * ESKF::Update()再把它和先验协方差组合求解。
     */
    struct CustomObservationModel {
        bool valid_ = true;     // 观测模型是否有效；例如有效匹配点太少时会置false并放弃本次更新
        bool converge_ = true;  // 当前迭代是否收敛，由ESKF::Update()根据dx阈值更新

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> R_;  // 预留的观测噪声矩阵，当前Update主流程未直接使用

        /// 信息矩阵近似：累加每个残差项的H_j^T H_j，当前只约束6维位姿。
        Eigen::Matrix<double, pose_obs_dim_, pose_obs_dim_> HTH_;
        /// 信息向量：累加每个残差项的H_j^T r_j，符号约定由观测函数中的残差定义决定。
        Eigen::Matrix<double, pose_obs_dim_, 1> HTr_;

        double lidar_residual_mean_ = 0;  // Lidar残差统计值，当前用平方残差中位数衡量本轮匹配质量
        double lidar_residual_max_ = 0;   // Lidar最大平方残差，用于调试和异常观测分析
    };

    /// 用户定义的观测函数：输入当前线性化状态，输出CustomObservationModel中的矩阵和残差信息。
    using CustomObsFunction = std::function<void(NavState& s, CustomObservationModel& obs)>;

    /**
     * @brief ESKF运行参数和各类观测函数入口。
     *
     * Options由上层模块配置后通过Init()写入滤波器。观测函数可以为空，但调用Update()
     * 时对应ObsType的函数必须已经配置。
     */
    struct Options {
        CustomObsFunction lidar_obs_func_;           // 雷达观测函数
        CustomObsFunction wheelspeed_obs_func_;      // 轮速观测函数
        CustomObsFunction acc_as_gravity_obs_func_;  // 加计观测函数
        CustomObsFunction gps_obs_func_;             // GPS/RTK观测函数
        CustomObsFunction bias_obs_func_;            // 零偏观测函数
        int max_iterations_ = 4;                     // 单次Update最多迭代次数
        StateVecType epsi_;                          // 每个误差状态维度的收敛阈值
        bool use_aa_ = false;                        // 是否使用Anderson Acceleration加速收敛

        /// 速度增量裁剪参数，目前实现里相关逻辑处于注释状态。
        double vel_clip_norm_ = 1.0;  // 速度更新量范数上限
        double dv_ratio_ = 0.5;       // 速度更新量缩放比例

        double predict_cov_inflation_ = 1.01;        // 预测后协方差膨胀系数，避免滤波器过度自信
        double min_cov_diag_ = 1e-9;                 // 协方差对角线下限，防止数值退化到非正定
        double degeneracy_threshold_ratio_ = 1e-3;   // 退化方向判定阈值，相对最大特征值设置
        bool propagate_velocity_ = false;            // 是否在名义状态中传播速度；SANY MID360需打开以匹配Voxel-SLAM。
        double degeneracy_cov_inflation_ = 1.02;     // 发现退化方向后，对位姿协方差的膨胀系数
        double max_update_translation_step_ = 0.5;   // 单次迭代允许的最大平移修正量
        double max_update_rotation_step_deg_ = 5.0;  // 单次迭代允许的最大旋转修正量，单位deg
        double max_update_velocity_step_ = 2.0;      // <= 0 disables velocity-step rejection for full-state ESKF.
        double max_update_gyro_bias_step_ = 0.05;    // rad/s, <= 0 disables lidar inertial-step fallback.
        double max_update_acc_bias_step_ = 0.5;      // m/s^2, <= 0 disables lidar inertial-step fallback.
        double max_update_gravity_step_ = 0.05;      // m/s^2, <= 0 disables lidar inertial-step fallback.
        bool lidar_update_pose_only_ = false;        // Lidar pose observations do not directly update vel/bias/gravity.
        bool lidar_update_inertial_states_ = true;   // Lidar pose observations may update bg/ba/gravity via covariance.
    };

    /**
     * @brief 初始化滤波器运行选项。
     *
     * Init()不会重置x_和P_，只更新观测函数、迭代次数、收敛阈值和数值保护参数。
     *
     * @param options 上层配置好的ESKF选项。
     */
    void Init(Options options) {
        lidar_obs_func_ = options.lidar_obs_func_;
        wheelspeed_obs_func_ = options.wheelspeed_obs_func_;
        acc_as_gravity_obs_func_ = options.acc_as_gravity_obs_func_;
        gps_obs_func_ = options.gps_obs_func_;
        bias_obs_func_ = options.bias_obs_func_;
        maximum_iter_ = options.max_iterations_;
        limit_ = options.epsi_;
        use_aa_ = options.use_aa_;

        options_ = options;
    }

    /**
     * @brief 使用一段IMU输入预测名义状态和误差状态协方差。
     *
     * 该函数会调用NavState中的运动模型，把名义状态x_积分到下一时刻；
     * 同时根据状态雅可比和噪声雅可比传播P_。
     *
     * @param dt 积分时间间隔，单位秒。
     * @param Q IMU过程噪声协方差，维度为process_noise_dim_。
     * @param gyro 当前区间使用的角速度，通常为相邻IMU量测平均值。
     * @param acce 当前区间使用的加速度，通常为相邻IMU量测平均值。
     */
    void Predict(const double& dt, const ProcessNoiseType& Q, const Vec3d& gyro, const Vec3d& acce);

    /**
     * @brief 使用指定观测类型对状态进行迭代更新。
     *
     * Update()会根据obs选择对应的CustomObsFunction，反复线性化观测模型并求解误差状态增量。
     * 收敛或达到最大迭代次数后，最终更新名义状态x_和协方差P_。
     *
     * @param obs 本次融合的观测类型。
     * @param R 观测噪声缩放因子。数值越小，观测相对先验的权重越高。
     */
    void Update(ObsType obs, const double& R);

    // accessors
    /// 获取当前名义状态。
    const NavState& GetX() const { return x_; }
    /// 获取当前误差状态协方差。
    const CovType& GetP() const { return P_; }
    /// 获取滤波器时间戳缓存。当前主要时间戳也保存在x_.timestamp_中。
    const double& GetStamp() const { return stamp_; }

    /// 直接替换名义状态，常用于初始化或外部重定位。
    void ChangeX(const NavState& state) { x_ = state; }
    /// 直接替换协方差，常用于初始化或测试。
    void ChangeP(const CovType& P) { P_ = P; }
    /// 修改滤波器时间戳缓存。
    void ChangeStamp(const double& stamp) { stamp_ = stamp; }

    /// 运行时开关Anderson Acceleration。
    void SetUseAA(bool use_aa) { use_aa_ = use_aa; }
    /// 只修改名义状态中的时间戳。
    void SetTime(double timestamp) { x_.timestamp_ = timestamp; }

    /// 迭代次数
    int GetIterations() const { return iterations_; }
    /// 最终平均观测误差
    double GetFinalRes() const { return final_res_; }
    bool LastUpdateAccepted() const { return last_update_accepted_; }

   private:
    double stamp_ = 0.0;

    NavState x_;                       // 名义状态，包含位置、姿态、速度、零偏等
    CovType P_ = CovType::Identity();  // 误差状态协方差矩阵
    CovType F_x1_ = CovType::Identity();  // Predict阶段使用的离散误差状态转移矩阵
    CovType L_ = CovType ::Identity();    // Update末尾做流形协方差映射时的临时矩阵

    CustomObservationModel custom_obs_model_;
    CustomObsFunction lidar_obs_func_, wheelspeed_obs_func_, acc_as_gravity_obs_func_, gps_obs_func_, bias_obs_func_;

    int maximum_iter_ = 0;  // 最大迭代次数
    StateVecType limit_;    // 收敛阈值，每一维误差增量都小于对应阈值时认为收敛

    int iterations_ = 0;        // 最近一次Update实际使用的迭代次数
    double final_res_ = 0.0;    // 最近一次Update的相对残差
    bool last_update_accepted_ = false;

    /// 是否使用Anderson Acceleration加速迭代收敛。
    bool use_aa_ = false;
    AndersonAcceleration<double, state_dim_, 10> aa_;

    Options options_;  // 保存完整选项，供Predict/Update中的数值保护逻辑使用
};

}  // namespace lightning

#endif  // FUSION_ESKF_HPP

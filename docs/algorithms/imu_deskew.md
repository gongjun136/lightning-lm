@page imu_deskew IMU 初始化、时间积分与点云去畸变

# 外部契约

@ref lightning::ImuProcess "ImuProcess" 接收一帧扫描和覆盖其时间的 IMU。初始化未完成时不输出去畸变扫描；正常阶段先推进 ESKF 到扫描末时刻，再将各点变换到该时刻的主 LiDAR 系。系统上下游见 @ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算"。

本页按三个问题读：@ref imu_init_reading "什么时候允许结束初始化" → @ref imu_time_reading "每次究竟积分到哪个时间" → @ref imu_point_reading "一个点如何变到扫描末坐标系"。先用片段和手算例子确认行为，再进入完整函数；例子中的阈值与运动假设不代表车辆现场配置或实测结果。

@anchor imu_init_reading
# 1. 样本数够了，为什么初始化仍不能结束？

`IMUInit` 在保留窗口中计算加速度、角速度均值和样本方差；`InitializationReady` 同时要求样本数、持续时间、平均角速度、两类标准差及加速度模长通过。仅“收集够 N 条”不充分。

@snippet{lineno} core/lio/imu_processing.hpp imu-init-ready

先看两个提前 `return false`，再看最后串联的 `&&`。`gyro_std`、`acc_std` 是三轴样本方差中最大者的平方根，表示各轴标准差的最大值，不是三维向量模长的标准差。平均角速度模长的单位为 rad/s；加计统计来自原始输入，初始化后的单位缩放尚未发生。

<b>手算例子：</b>假设 `min_samples=20`、`min_duration=0.05` s，以严格均匀的 400 Hz 采样。20 个样本只跨越 19 个间隔，即 0.0475 s，所以数目达标仍会在第二个 `if` 返回 false；21 个样本才覆盖 0.05 s。即使时间也够，若配置平均角速度上限 0.1 rad/s、实际均值模长为 0.2 rad/s，最终 `&&` 仍拒绝结束初始化。结构体中的部分默认上限为无穷大，实际静止约束强弱取决于传入的 `InitializationOptions`。

姿态和偏置又怎样从均值得到？在 @ref lightning::ImuProcess::IMUInit() "IMUInit()" 中只看下面两段：

@snippet{lineno} core/lio/imu_processing.hpp imu-init-orientation

@snippet{lineno} core/lio/imu_processing.hpp imu-init-gravity-bias

`FromTwoVectors(acc_dir, UnitZ)` 将静止比力方向对齐到世界 +z，再左乘配置 yaw；随后设世界重力为 `(0,0,-9.81)`，陀螺偏置为角速度均值，加计偏置为零。例如均值为 `(0,0,9.81)`、配置 yaw=30°，倾斜对齐是单位旋转，最终姿态的 yaw 来自配置的 30°，而非由静止加速度测出来。`IMUInit()` 每次都暂写这些估计，是否准许进入正常处理仍由上面的就绪判断决定。

初始化完成后，预测噪声由配置尺度覆盖静止统计；加速度均值模长接近 1 时乘 9.81，接近 9.81 时不缩放，其他范围保持尺度并告警。不要把输入的 g 单位当作 m/s²。

这些缩放发生在通过就绪检查之后，不会自动换算初始化的加计标准差或模长阈值。还应区分“本次确认初始化完成”和“本帧已经去畸变”。先看入口 `if (imu_need_init_)` 与通过就绪检查后的赋值：

@snippet{lineno} core/lio/imu_processing.hpp imu-init-ready-transition

再看同一初始化分支的末尾。`imu_need_init_ = false` 不会跳出已经进入的 `if`；本次仍先 `return`，后续调用才走分支外的 `UndistortPcl()`：

@snippet{lineno} core/lio/imu_processing.hpp imu-init-return-before-deskew

<b>停止位置：</b>能指出样本、时间与统计门槛，并解释初始化完成这一帧为什么仍没有去畸变输出，就可以继续下一节。完整条件与配置类型见 @ref lightning::ImuProcess::InitializationReady() "InitializationReady()"、@ref lightning::ImuProcess::InitializationOptions "InitializationOptions"；完成后的控制流见 @ref lightning::ImuProcess::Process() "Process()"。

@anchor imu_time_reading
# 2. 跨帧 IMU 怎样衔接，为什么不能把负时间差取绝对值？

`UndistortPcl` 将上一帧末尾的 `last_imu_` 拼到序列头，保留跨帧区间。目标端点取 `min(tail.timestamp, scan_end)`，步长由当前状态时间计算；非正步跳过，最后仅补齐尚未到达扫描末的正时间差。长间隔会告警，不代表已得到可靠外推。

@snippet{lineno} core/lio/imu_processing.hpp imu-cross-frame-bridge

`last_imu_` 用于与当前首条 IMU 配成相邻量测；它不是把 ESKF 时间重设回上一条原始消息。循环先跳过完全落在上一扫描结束之前的区间，再进入下面的均值输入与目标端点计算：

@snippet{lineno} core/lio/imu_processing.hpp imu-interval-end

只追 `interval_end` 和 `kf_state.GetX().timestamp_`：两者都是传感器时间，单位 s；`dt` 从<b>当前状态已积分到的时间</b>计算，而非无条件使用 `tail-head`。裁剪端点也没有把两条 IMU 量测重新按裁剪比例插值，输入仍取两端均值。

<b>手算例子：</b>当前 ESKF 已到 10.02 s，下一个 IMU 区间尾端为 10.08 s，扫描末为 10.05 s。代码取 `interval_end=10.05`，只积分 0.03 s。若该扫描末已经覆盖，接下来的 `PredictTo(10.05, ...)` 会返回 false：

@snippet{lineno} core/lio/eskf.cc eskf-predict-time-gate

这里先比较目标时间与状态时间，合法后才调用 `Predict(target-current, ...)`。若状态已经为 10.06 s，它同样拒绝 10.05 s；不会用 `abs(-0.01)` 又向前推进，也不会偷偷把状态时间退回去。

末尾用同一个包装函数补齐剩余区间，并且只在预测成功时触发可选回调：

@snippet{lineno} core/lio/imu_processing.hpp imu-predict-scan-end

<b>读到这里：</b>能用当前状态时间和 `min(tail, scan_end)` 算出本次 dt 即可。名义状态的积分公式见 @ref eskf_prediction_reading "ESKF 一步预测"；更完整的缓存、长间隔告警和回调顺序见 @ref lightning::ImuProcess::UndistortPcl() "UndistortPcl()"。`PredictTo()` 通过时间检查只证明执行了该段预测，不证明 IMU 数据或长时间外推在物理上可靠。

@anchor imu_point_reading
# 3. 一个点怎样变换到扫描末时刻的 LiDAR 系？

设点 p 在时刻 t 的 LiDAR 系，外参为 \f$T_{il}\f$，扫描末时刻 e。对应同一世界点：
\f[
p_w=R_{wi}(t)(R_{il}p_l+t_{il})+p_{wi}(t).
\f]
变回末时刻的 LiDAR：
\f[
p_{l,e}=R_{il}^T\left[R_{wi}(e)^T\left(R_{wi}(t)(R_{il}p_l+t_{il})+
p_{wi}(t)-p_{wi}(e)\right)-t_{il}\right].
\f]
`p_compensate` 正是这个顺序。区间内使用
\f$R(t)=R_h\operatorname{Exp}(\omega\Delta t)\f$、
\f$p(t)=p_h+v_h\Delta t+\frac12a\Delta t^2\f$，并不重新运行滤波器。
加速度表达在世界系，包含重力项；`math::exp(angvel_avr,0.5*dt)` 的半角约定见 @ref geometry "坐标系、右扰动与李群雅可比"。

核对公式前先认清本段的变量，避免把 `head` 与扫描末混为同一状态：

| 源码量 | 公式中的含义 | 单位/坐标系 |
|---|---|---|
| `R_lidar_imu_`, `t_lidar_mu_` | 外参 \f$R_{il},t_{il}\f$ | LiDAR → IMU；平移 m |
| `R_imu`, `pos_imu`, `vel_imu` | 区间起点 h 的姿态、位置与速度 | 姿态 IMU → 世界；位置 m、速度 m/s，均在世界系 |
| `angvel_avr`, `acc_imu` | 本区间去偏角速度与世界净加速度 | rad/s、m/s² |
| `imu_state.rot_`, `imu_state.pos_` | 扫描末 e 的状态 | IMU → 世界、世界位置 m |
| `P_i` | 点在其采样时刻的 LiDAR 坐标 | m |
| `dt` | 点采样时刻减当前 IMU 区间起点 | s；不等于整段扫描时长 |

点的 `time` 先从相对扫描开始的毫秒转换为秒，再减去区间起点的相对秒数。进入补偿的时间条件在这里：

@snippet{lineno} core/lio/imu_processing.hpp imu-point-time-gate

然后对这个点构造姿态，并按上述公式变换坐标：

@snippet{lineno} core/lio/imu_processing.hpp imu-point-rotation

@snippet{lineno} core/lio/imu_processing.hpp imu-point-compensation

`T_ei` 在实际赋值中是 \f$p_{wi}(t)-p_{wi}(e)\f$，表达在世界系。按括号从内向外读：先用外参把 `P_i` 转到 IMU；用 `R_i` 和 `T_ei` 表达它相对扫描末 IMU 原点的世界向量；最后应用末时刻旋转和外参的逆，得到末时刻 LiDAR 坐标。`0.5 * acc_imu * dt * dt` 仅用于这里的区间内点位置近似，不能反推 ESKF 的名义位置也使用同一积分式。

<b>手算例子：</b>假设外参为单位变换、姿态不变、匀速沿 x 前进，某个进入循环的点采样时 IMU 世界位置为 0，扫描末为 1 m，点的原始坐标为 `(10,0,0)` m。则 `T_ei=(-1,0,0)`，补偿后点为 `(9,0,0)` m：同一个静止世界点，在前移 1 m 的末时刻传感器中距离减小 1 m。若传感器完全不动，`T_ei=0` 且前后旋转相同，外参的正、逆变换抵消，坐标保持不变。

<b>停止位置与边界：</b>能从 `P_i → T_ei → p_compensate` 手算上述结果即可。源码最后只把补偿结果写回 x/y/z；这里没有再次对点运行滤波更新。完整实现见 @ref lightning::ImuProcess::UndistortPcl() "UndistortPcl()"，上游同步与无输出条件见 @ref lio_drain_reading "LIO 排空导读"。

# 按需核对的实现边界

点按相对毫秒时间排序，从尾到头补偿。当前内层循环条件排除了 `points.begin()`，因此不能写成“每一个点都保证执行补偿”；这是源码现状，本次文档整理不改变算法。可选 IMU 滤波分支通过共享指针改写 IMU 对象，不能仅根据局部容器拷贝认定原始消息完全不变。

<b>源码入口：</b>@ref imu_processing.hpp "imu_processing.hpp"（`IMUInit`、`InitializationReady`、`UndistortPcl`、`Process`），@ref lightning::ESKF::PredictTo "PredictTo"。

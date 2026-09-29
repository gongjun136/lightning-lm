@page imu_deskew IMU 初始化、时间积分与点云去畸变

# 外部契约

@ref lightning::ImuProcess "ImuProcess" 接收一帧扫描和覆盖其时间的 IMU。初始化未完成时不输出去畸变扫描；正常阶段先推进 ESKF 到扫描末时刻，再将各点变换到该时刻的主 LiDAR 系。系统上下游见 @ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算"。

# 静止初始化

`IMUInit` 在保留窗口中计算加速度、角速度均值和样本方差；`InitializationReady` 同时要求样本数、持续时间、平均角速度、两类标准差及加速度模长通过。仅“收集够 N 条”不充分。
用 `FromTwoVectors(acc_dir, UnitZ)` 求倾斜对齐，再左乘配置初始 yaw；设世界重力为 (0,0,-9.81)，陀螺偏置为角速度均值，加计偏置初始化为零。静止加速度只能确定倾斜，不能观测绝对航向。

初始化完成后，预测噪声由配置尺度覆盖静止统计；加速度均值模长接近 1 时乘 9.81，接近 9.81 时不缩放，其他范围保持尺度并告警。不要把输入的 g 单位当作 m/s²。

# 时间边界

`UndistortPcl` 将上一帧末尾的 `last_imu_` 拼到序列头，保留跨帧区间。目标端点取 `min(tail.timestamp, scan_end)`，步长由当前状态时间计算；非正步跳过，最后仅补齐尚未到达扫描末的正时间差。长间隔会告警，不代表已得到可靠外推。

# 逐点变换推导

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

# 实现边界

点按相对毫秒时间排序，从尾到头补偿。当前内层循环条件排除了 `points.begin()`，因此不能写成“每一个点都保证执行补偿”；这是源码现状，本次文档整理不改变算法。可选 IMU 滤波分支通过共享指针改写 IMU 对象，不能仅根据局部容器拷贝认定原始消息完全不变。

**源码入口：** @ref imu_processing.hpp "imu_processing.hpp"（`IMUInit`、`InitializationReady`、`UndistortPcl`、`Process`），@ref lightning::ESKF::PredictTo "PredictTo"。

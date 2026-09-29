@page timing_contracts 时间一致性检查

# 按边界核对

1. 消息 header 与点内时间：秒与毫秒转换只做一次，扫描末时间与 IMU 覆盖一致。
2. `PredictTo`：状态时间是积分前沿，跨扫描边界截断区间，不对负 dt 取绝对值。
3. 跨线程：sensor 队列保持到达次序，不自动对任意乱序输入按时间排序。
4. PGO：重定位清理旧相对运动队列，防止跨恢复时间段插值。
5. 输出：传感器时间用于位姿关联，steady clock 用于判断真实等待；两种年龄都必须可解释。

**源码依据：** @ref imu_processing.hpp "imu_processing.hpp"、@ref eskf.cc "eskf.cc"、@ref localization.cpp "localization.cpp"、@ref pgo_impl.cc "pgo_impl.cc"、@ref sany_localization_output.cc "sany_localization_output.cc"。
具体积分推导见 @ref imu_deskew "IMU 初始化、时间积分与点云去畸变"；图状态见 @ref pose_graph_theory "位姿图、增量求解与高频平滑"；输出失效门控见 @ref output_contracts "定位状态、车体参考点与 ROS 输出"。

# 回归方法

覆盖重复/回退时间、跨扫描末 IMU、无后续 IMU、长间隔、重定位后旧数据与输出时间单调性。单元测试覆盖局部边界，bag 回放观察各队列丢弃、帧龄和有效输出。WSL 回放通过只能证明该输入和平台条件，不代表域控全栈时延合格。

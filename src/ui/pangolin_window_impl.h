#pragma once

#include <pangolin/pangolin.h>

#include <pcl/filters/voxel_grid.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/impl/pcl_base.hpp>

#include <atomic>
#include <deque>
#include <mutex>
#include <string>
#include <thread>

#include "common/keyframe.h"
#include "common/loop_candidate.h"

#include "ui/pangolin_window.h"
#include "ui/ui_car.h"
#include "ui/ui_cloud.h"
#include "ui/ui_trajectory.h"

namespace lightning::ui {

/**
 */
class PangolinWindowImpl {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PangolinWindowImpl() = default;
    ~PangolinWindowImpl() = default;

    PangolinWindowImpl(const PangolinWindowImpl &) = delete;
    PangolinWindowImpl &operator=(const PangolinWindowImpl &) = delete;
    PangolinWindowImpl(PangolinWindowImpl &&) = delete;
    PangolinWindowImpl &operator=(PangolinWindowImpl &&) = delete;

    /// 初始化，创建可用于渲染的点云和轨迹
    bool Init();

    void Reset(const std::vector<Keyframe::Ptr> &keyframes);

    /// 注销
    bool DeInit();

    /// 渲染所有信息
    void Render();

    /// 获取窗口名称
    std::string GetWindowName() const;

   public:
    /// 后台渲染线程
    std::thread render_thread_;

    /// 一些辅助的锁和原子变量
    std::mutex mtx_map_cloud_;
    std::mutex mtx_current_scan_;  ///< 保护待显示的当前扫描点云。
    std::mutex mtx_nav_state_;  ///< 保护待显示的导航状态与位姿。
    std::mutex mtx_gps_pose_;  ///< 保护 GPS 位姿显示数据。
    std::mutex mtx_loop_info_;  ///< 保护回环连接与候选显示数据。

    std::mutex mtx_reset_;  ///< 串行保护显示状态的重置。

    std::atomic<bool> exit_flag_;  ///< 请求渲染线程结束的标志。
    std::atomic<bool> should_quit_ = true;  ///< 窗口或上层是否请求退出。
    std::atomic<bool> render_ready_ = false;  ///< 渲染资源是否已初始化完成。
    std::atomic<bool> render_finished_ = false;  ///< 渲染循环是否已结束。

    std::atomic<bool> cloud_global_need_update_;   ///< 全局点云是否需要更新
    std::atomic<bool> cloud_dynamic_need_update_;  ///< 动态点云是否需要更新
    std::atomic<bool> kf_result_need_update_;      ///< 卡尔曼滤波结果
    std::atomic<bool> current_scan_need_update_;   ///< 更新当前扫描
    std::atomic<bool> lidarloc_need_update_;       ///< 雷达位置？

    pcl::PointCloud<PointType>::Ptr current_scan_ = nullptr;  ///< 当前scan
    SE3 newest_frontend_pose_;                                ///< 最新pose
    SE3 predicted_pose_;  ///< 最近的高频预测位姿，用于显示。
    SE3 newest_backend_pose_;  ///< 最新pose
    SE3 current_scan_pose_;    ///< 当前scan对应的pose or Twb/Twi
    std::deque<std::pair<int, int>> loop_info_;  ///< 算法侧等待更新到 UI 的回环连接。
    std::vector<LoopCandidate> new_loop_candidate_;  ///< 待显示的新增回环候选。

    // 地图点云
    std::map<int, CloudPtr> cloud_global_map_;  ///< 按地图块 ID 保存的静态地图显示点云。
    std::map<int, CloudPtr> cloud_dynamic_map_;  ///< 按地图块 ID 保存的动态地图显示点云。

    /// 滤波器状态
    Sophus::SE3d pose_;
    double confidence_;  ///< 当前显示的定位置信度。
    Vec3d vel_;  ///< 当前显示的导航速度。
    Vec3d bias_acc_;  ///< 当前显示的加速度计偏置。
    Vec3d bias_gyr_;  ///< 当前显示的陀螺仪偏置。
    Vec3d grav_;  ///< 当前显示的重力向量。

    Sophus::SE3d T_imu_lidar_;  ///< 显示点云与 IMU 状态时使用的固定外参。
    int max_size_of_current_scan_ = 200;  ///< 当前扫描数据保留多少个
    std::vector<std::shared_ptr<Keyframe>> all_keyframes_;  ///< 供全局轨迹与点云显示使用的关键帧集合。

    //////////////////////////////// 以下和render相关 ///////////////////////////
   private:
    /// 创建OpenGL Buffers
    void AllocateBuffer();
    void ReleaseBuffer();

    void CreateDisplayLayout();

    void DrawAll();  // 作图：画定位窗口

    /// 渲染点云，调用各种Update函数
    void RenderClouds();
    bool UpdateGlobalMap();
    bool UpdateDynamicMap();
    bool UpdateState();
    bool UpdateCurrentScan();

    void RenderLabels();

   private:
    /// 窗口layout相关
    int win_width_ = 1920;
    int win_height_ = 1080;  ///< 渲染窗口高度，单位像素。
    static constexpr float cam_focus_ = 5000;  ///< 主视图相机焦距配置。
    static constexpr float cam_z_near_ = 1.0;  ///< 主视图相机近裁剪平面距离。
    static constexpr float cam_z_far_ = 1e10;  ///< 主视图相机远裁剪平面距离。
    static constexpr int menu_width_ = 210;  ///< 控制菜单宽度，单位像素。
    const std::string win_name_ = "UI";  ///< 主渲染窗口名称。
    const std::string dis_main_name_ = "main";  ///< 主显示容器名称。
    const std::string dis_3d_name_ = "Cam 3D";  ///< 三维显示视图名称。
    const std::string dis_3d_main_name_ = "Cam 3D Main";  ///< main
    const std::string dis_plot_name_ = "Plot";  ///< 统计曲线视图名称。
    const std::string dis_imgs_name = "Images";  ///< 图像视图名称。

    bool following_loc_ = true;       ///< 相机是否追踪定位结果
    bool draw_frontend_traj_ = true;  ///< 可视化前端轨迹
    bool draw_backend_traj_ = true;   ///< 可视化后端轨迹

    // text
    pangolin::GlText gltext_label_global_;  ///< 全局地图视图的文字标签。
    pangolin::GlText gltext_label_state_;  ///< 状态信息视图的文字标签。

    // camera
    pangolin::OpenGlRenderState s_cam_main_;  ///< 主三维视图的相机投影与观察状态。

    /// cloud rendering
    ui::UiCar backend_car_{Vec3f(0.2, 0.2, 0.8)};   // 白色车
    ui::UiCar frontend_car_{Vec3f(0.2, 0.2, 0.8)};  ///< 显示前端估计位姿的车辆图形。
    ui::UiCar pred_car_{Vec3f(0.8, 0.5, 0.8)};  ///< 显示高频预测位姿的车辆图形。

    std::map<int, std::shared_ptr<ui::UiCloud>> cloud_map_ui_;  ///< 用来渲染的点云地图
    std::map<int, std::shared_ptr<ui::UiCloud>> cloud_dyn_ui_;  ///< 用来渲染的点云地图
    std::shared_ptr<ui::UiCloud> current_scan_ui_;              ///< current scan
    std::deque<std::shared_ptr<ui::UiCloud>> scans_;            ///< current scan 保留的队列
    std::deque<std::pair<int, int>> loop_info_ui_;  ///< 渲染线程使用的回环连接缓存。

    // trajectory
    std::shared_ptr<ui::UiTrajectory> traj_scans_ = nullptr;         ///< 激光扫描的轨迹
    std::shared_ptr<ui::UiTrajectory> traj_newest_state_ = nullptr;  ///< 最新state的轨迹

    // 滤波器状态相关 Data logger object
    pangolin::DataLog log_vel_;           ///< odom frame下的速度
    pangolin::DataLog log_vel_baselink_;  ///< baselink frame下的速度
    pangolin::DataLog log_bias_acc_;      ///< accelerometer bias
    pangolin::DataLog log_confidence_;    ///< confidence
    pangolin::DataLog log_error_;  ///< 显示定位误差的曲线数据日志。

    std::unique_ptr<pangolin::Plotter> plotter_vel_ = nullptr;  ///< 世界坐标速度曲线绘制器。
    std::unique_ptr<pangolin::Plotter> plotter_vel_baselink_ = nullptr;  ///< 车体坐标速度曲线绘制器。
    std::unique_ptr<pangolin::Plotter> plotter_bias_acc_ = nullptr;  ///< 加速度计偏置曲线绘制器。
    std::unique_ptr<pangolin::Plotter> plotter_confidence_ = nullptr;  ///< 定位置信度曲线绘制器。
    std::unique_ptr<pangolin::Plotter> plotter_err_ = nullptr;  ///< 定位误差曲线绘制器。
    std::unique_ptr<pangolin::Plotter> plotter_err_eval_ = nullptr;  ///< 轨迹评估误差曲线绘制器。
};

}  // namespace lightning::ui

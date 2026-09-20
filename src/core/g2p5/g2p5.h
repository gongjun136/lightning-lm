//
// Created by xiang on 25-6-23.
//

#ifndef LIGHTNING_G2P5_H
#define LIGHTNING_G2P5_H

#include "common/eigen_types.h"
#include "common/keyframe.h"
#include "core/g2p5/g2p5_map.h"
#include "core/system/async_message_process.h"

#include <thread>

namespace lightning::g2p5 {

/**
 * @brief 2.5D栅格地图渲染与三维点云到二维scan的转换模块。
 *
 * @details
 * g2p5 是一个2.5D地图模块，它计算从雷达中心点出发，打到障碍物部分的射线，然后进行occupancy map投影。
 * 可以很方便地将g2p5转换成opencv或者ros格式的地图进行发布和存储
 *
 * g2p5假设雷达是水平安装的。在雷达水平面以下，存在一个地面的平面（floor），地面的高度应该低于雷达安装高度
 * 这个地面参数可以动态估计，也可以事先指定
 * 比如，我们默认雷达水平面为z=0，雷达安装在1m的高度上，那么地面默认为z=-1.0
 * 在地面上方一定区间内的点云才会参与栅格地图的计算，默认取 (min_floor, max_floor)
 * 之间的点云，比如地面以上的0.5至1.0米的障碍物 那么取 min_floor = 0.5, max_floor = 1.0
 *
 * 如果想过滤掉天花板，应该将max_floor适当设小一点
 *
 * 如果是这个区间的障碍物，g2p5会模拟一条从雷达高度的点射向障碍物点，并计算沿途的射线高度
 * 如果栅格中的物体高度低于此射线，则不会被刷新；如果高于此射线，则会被刷白
 *
 * 在lightning-lm中，栅格地图仅用于显示和存储（输出给导航），并不用于定位，所以默认不需要过高的分辨率
 *
 * 地图渲染主要有两个步骤组成：
 * 1. 前端在计算完关键帧时，会立即更新到最新的地图上，此时应该发布最新的地图
 * 2. 后端在发现回环时，触发一次地图重绘过程。由于重绘时间较长，需要在重绘完成后再放到1,
 * 此时1仍然是活跃的，仍然会发布地图
 *
 * g2p5的3转2目前角分辨率精度较低，但速度很快，正常情况一帧点云应该在10ms以内
 */
class G2P5 {
   public:
    /**
     * @brief G2P5渲染与地面过滤参数。
     *
     * 参数可由构造函数传入，也可在Init()中从YAML配置覆盖部分字段。
     */
    struct Options {
        Options() {}

        bool online_mode_ = true;  ///< 是否为在线模式；在线模式下前端渲染由异步线程处理。
        bool esti_floor_ = false;  ///< 是否需要从关键帧点云中动态估计地面平面。

        double lidar_height_ = 0.0;  ///< 雷达的安装高度。

        double default_floor_height_ = -1.0;  ///< 默认地面高度，通常位于雷达下方。
        double min_th_floor_ = 0.5;           ///< 距离地面超过该高度的点才可能作为障碍物写入栅格地图。
        double max_th_floor_ = 1.2;           ///< 距离地面低于该高度的障碍物点才会参与栅格地图计算。
        float usable_scan_range_ = 50.0;      ///< 用于计算栅格地图的最远障碍物距离。

        double grid_map_resolution_ = 0.1;  ///< 栅格地图分辨率；室外建议0.1，室内可使用0.05。

        bool verbose_ = true;  ///< 是否输出调试日志。
    };

    /**
     * @brief 构造G2P5地图渲染器。
     * @param options 初始渲染参数。
     */
    explicit G2P5(Options options = Options()) : options_(options) {}

    /**
     * @brief 析构渲染器并停止后台线程。
     */
    ~G2P5();

    /**
     * @brief 地图更新回调类型。
     *
     * 回调参数为最新栅格地图指针，调用方应按需复制或只读使用。
     */
    using MapUpdateCallback = std::function<void(G2P5MapPtr map)>;

    /**
     * @brief 从YAML文件读取配置信息并启动必要的渲染线程。
     * @param yaml_path YAML配置文件路径。
     */
    void Init(std::string yaml_path);

    /**
     * @brief 设置地图更新回调。
     * @param func 地图发生变化时调用的回调函数。
     */
    void SetMapUpdateCallback(MapUpdateCallback func) { map_update_cb_ = func; }

    /**
     * @brief 请求退出渲染器并等待后台渲染线程结束。
     */
    void Quit();

    /**
     * @brief 增加一个关键帧并触发前端地图更新。
     * @param kf 待加入栅格地图的关键帧。
     */
    void PushKeyframe(Keyframe::Ptr kf);

    /**
     * @brief 触发一次全局地图重绘。
     *
     * 如果后台已经在重绘，再次触发会使当前重绘过程停止，并在下一轮使用最新关键帧重新绘制。
     */
    void RedrawGlobalMap();

    /**
     * @brief 获取最新地图的深拷贝。
     * @return 最新地图副本；尚未生成地图时返回nullptr。
     */
    G2P5MapPtr GetNewestMap();

    /**
     * @brief 设置是否启用并行渲染模式。
     * @param enable_parallel 为true时请求启用并行渲染。
     * @return 当前接口固定返回true。
     */
    bool SetParallelRendering(bool enable_parallel = false) {
        parallel_render_ = enable_parallel;
        return true;
    }

    /**
     * @brief 查询后台是否仍有未完成的重绘任务。
     * @return 后台重绘忙碌时返回true，否则返回false。
     */
    bool IsBusy() { return is_busy_; }

   private:
    /**
     * @brief 将一组关键帧增量写入已有地图。
     * @param[in] kfs 待写入地图的关键帧集合。
     * @param[in,out] map 待更新的栅格地图。
     * @return 写入流程正常结束返回true。
     */
    bool AddKfToMap(const std::vector<Keyframe::Ptr> &kfs, G2P5MapPtr &map);

    /**
     * @brief 根据关键帧轨迹和点云范围扩展地图边界。
     * @param[in] kfs 用于计算边界的关键帧集合。
     * @param[in,out] map 待缩放的栅格地图。
     * @return 地图边界有效时返回true，否则返回false。
     */
    bool ResizeMap(const std::vector<Keyframe::Ptr> &kfs, G2P5MapPtr &map);

    /**
     * @brief 将单个关键帧渲染到前端地图。
     * @param kf 待渲染的关键帧。
     */
    void RenderFront(Keyframe::Ptr kf);

    /**
     * @brief 后台循环重绘全局地图。
     */
    void RenderBack();

    /**
     * @brief 将关键帧中的三维点云转换为二维scan并更新栅格地图。
     * @param kf 待转换的关键帧。
     * @param[in,out] map 待更新的栅格地图。
     */
    void Convert3DTo2DScan(Keyframe::Ptr kf, G2P5MapPtr &map);

    /**
     * @brief 从关键帧点云中检测地面平面参数。
     * @param kf 用于地面估计的关键帧。
     * @return 成功检测到近似水平地面时返回true，否则返回false。
     */
    bool DetectPlaneCoeffs(Keyframe::Ptr kf);

    /**
     * @brief 沿二维scan射线将可通行区域写入地图。
     * @param[in] ang_distance_height 每个角度对应的距离和高度信息，期望大小为360。
     * @param kf 提供位姿的关键帧。
     * @param[in,out] map 待更新的栅格地图。
     */
    void SetWhitePoints(const std::vector<Vec2d> &ang_distance_height, Keyframe::Ptr kf, G2P5MapPtr &map);

    Options options_;  ///< 渲染与地面过滤参数。

    std::atomic_bool parallel_render_ = false;  ///< 是否需要并行渲染。
    std::atomic_bool is_busy_ = false;          ///< 是否仍有未渲染完成的后台任务。
    MapUpdateCallback map_update_cb_;           ///< 地图更新回调函数。
    std::atomic_bool quit_flag_ = false;        ///< 退出标记位。

    std::mutex kf_mutex_;                       ///< 保护全部关键帧列表的互斥锁。
    std::vector<Keyframe::Ptr> all_keyframes_;  ///< 全部关键帧。

    std::mutex newest_map_mutex_;       ///< 保护最新地图指针的互斥锁。
    G2P5MapPtr newest_map_ = nullptr;   ///< 当前可对外读取的最新地图。

    /// @name 前端渲染状态
    /// @{
    std::mutex frontend_mutex_;                                         ///< 前端状态互斥锁。
    G2P5MapPtr frontend_map_ = nullptr;                                 ///< 前端最新绘制的地图。
    Keyframe::Ptr frontend_current_ = nullptr;                          ///< 前端正在绘制的关键帧。
    sys::AsyncMessageProcess<Keyframe::Ptr> draw_frontend_map_thread_;  ///< 前端绘制地图的异步处理线程。
    /// @}

    /// @name 后端重绘状态
    /// @{
    std::thread draw_backend_map_thread_;           ///< 后端重绘线程。
    std::atomic_bool backend_redraw_flag_ = false;  ///< 后端重绘请求标记。
    G2P5MapPtr backend_map_ = nullptr;              ///< 后端正在绘制的地图。
    /// @}

    Vec4d floor_coeffs_ = Vec4d(0, 0, 1.0, 1.0);  ///< 地面方程参数。
};

}  // namespace lightning::g2p5

#endif  // LIGHTNING_G2P5_H

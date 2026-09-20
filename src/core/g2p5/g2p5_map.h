//
// Created by xiang on 25-6-23.
//

#ifndef LIGHTNING_G2P5_MAP_H
#define LIGHTNING_G2P5_MAP_H

#include "common/eigen_types.h"
#include "common/std_types.h"
#include "core/g2p5/g2p5_subgrid.h"

#include <bitset>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <opencv2/core.hpp>

namespace lightning::g2p5 {

/**
 * @brief G2P5自定义2.5D栅格地图数据结构。
 *
 * @details
 * 主要用于对外的绘制，并没有匹配部分
 * 数据存储在grids_中，是一个二维数组，内部还有4层subgrids
 * x为行指针，y为列指针，y优先增长
 *
 * 可以通过ToCV 或者 ToROS 转换成OpenCV格式或者ROS格式进行显示和存储
 */
class G2P5Map {
   public:
    /**
     * @brief 栅格地图参数。
     */
    struct Options {
        float resolution_ = 0.05;      ///< 最细栅格分辨率。
        float occupancy_ratio_ = 0.3;  ///< 判定占用时使用的占用比例阈值。
    };

    /**
     * @brief 构造栅格地图对象。
     * @param options 地图分辨率和占用阈值参数。
     */
    G2P5Map(Options options) : options_(options) {
        grid_reso_ = options_.resolution_ * sub_grid_width_;
        grids_ = nullptr;
    }

    /**
     * @brief 析构地图并释放内部SubGrid数组。
     */
    ~G2P5Map();

    /// 子网格索引位宽。
    static inline constexpr int SUB_GRID_SIZE = SubGrid::SUB_GRID_SIZE;

    /**
     * @brief 由当前地图内容创建深拷贝。
     * @return 新地图指针，内部SubGrid数据与当前对象独立。
     */
    std::shared_ptr<G2P5Map> MakeDeepCopy();

    /**
     * @brief 转换为ROS OccupancyGrid消息。
     * @return 可用于发布或保存的ROS占用栅格。
     */
    nav_msgs::msg::OccupancyGrid ToROS();

    /**
     * @brief 转换为OpenCV图像。
     * @return 用于显示或保存的cv::Mat。
     */
    cv::Mat ToCV();

    /**
     * @brief 清空内部地图数据并重置边界。
     */
    void ReleaseResources();

    /**
     * @brief 初始化地图边界并分配SubGrid数组。
     * @param temp_min_x 地图最小x坐标。
     * @param temp_min_y 地图最小y坐标。
     * @param temp_max_x 地图最大x坐标。
     * @param temp_max_y 地图最大y坐标。
     * @return 初始化成功返回true；边界无效时返回false。
     */
    bool Init(const float &temp_min_x, const float &temp_min_y, const float &temp_max_x, const float &temp_max_y);

    /**
     * @brief 调整地图边界并保留旧地图中仍然落在新边界内的数据。
     * @param temp_min_x 新地图最小x坐标。
     * @param temp_min_y 新地图最小y坐标。
     * @param temp_max_x 新地图最大x坐标。
     * @param temp_max_y 新地图最大y坐标。
     * @return 调整完成返回true。
     */
    bool Resize(const float &temp_min_x, const float &temp_min_y, const float &temp_max_x, const float &temp_max_y);

    /**
     * @brief 按地图坐标写入一个命中点。
     * @param px 命中点x坐标。
     * @param py 命中点y坐标。
     * @param if_hit 是否按占用点写入；通常true表示障碍命中，false表示空闲。
     * @param height 命中点相对地面高度。
     */
    void SetHitPoint(const float &px, const float &py, const bool &if_hit, float height);

    /**
     * @brief 沿雷达到目标点的射线写入miss点。
     *
     * 白色区域的直线填充算法。除了2D填充以外，还需要给出雷达高度和目标位置高度，
     * 用于在射线上插值高度并避免错误清空高于射线的障碍物。
     *
     * @param point_x 射线终点x坐标。
     * @param point_y 射线终点y坐标。
     * @param laser_origin_x 雷达原点x坐标。
     * @param laser_origin_y 雷达原点y坐标。
     * @param height 射线终点高度。
     * @param lidar_height 雷达高度。
     */
    void SetMissPoint(const float &point_x, const float &point_y, const float &laser_origin_x,
                      const float &laser_origin_y, float height, float lidar_height);

    /**
     * @brief 判断细栅格索引是否落在当前地图范围内。
     * @param point 细栅格索引。
     * @return 索引位于已分配地图范围内返回true，否则返回false。
     */
    bool IsObstacle(const Vec2i &point) {
        // 当前实现仅判断点是否落在已分配地图范围内。
        int xi = point.x() >> SUB_GRID_SIZE;
        int yi = point.y() >> SUB_GRID_SIZE;
        if (xi < 0 || xi >= grid_size_x_ || yi < 0 || yi >= grid_size_y_) {
            return false;
        }
        return true;
    }

    /**
     * @brief 按细栅格索引更新单个栅格状态。
     * @param point_index 以最细分辨率表示的栅格索引。
     * @param if_hit 是否写入命中状态。
     * @param height 当前点或射线位置的高度。
     */
    void UpdateCell(const Vec2i &point_index, const bool &if_hit, float height);

    /**
     * @brief 将地图坐标转换为最细栅格索引。
     * @param x 地图坐标x。
     * @param y 地图坐标y。
     * @param[out] x_index 输出的x方向细栅格索引。
     * @param[out] y_index 输出的y方向细栅格索引。
     * @return 坐标落在地图范围内返回true，否则返回false。
     */
    bool GetDataIndex(const float x, const float y, int &x_index, int &y_index);

    /**
     * @brief 查询地图是否尚未分配内部数据。
     * @return 内部SubGrid数组为空返回true，否则返回false。
     */
    bool IsEmpty() { return (grids_ == nullptr); }

    /**
     * @brief 获取地图边界。
     * @param[out] min_x 地图最小x坐标。
     * @param[out] min_y 地图最小y坐标。
     * @param[out] max_x 地图最大x坐标。
     * @param[out] max_y 地图最大y坐标。
     */
    inline void GetMinAndMax(float &min_x, float &min_y, float &max_x, float &max_y) {
        min_x = min_x_;
        min_y = min_y_;
        max_x = max_x_;
        max_y = max_y_;
    }

    /**
     * @brief 设置地图边界。
     * @param min_x 地图最小x坐标。
     * @param min_y 地图最小y坐标。
     * @param max_x 地图最大x坐标。
     * @param max_y 地图最大y坐标。
     */
    inline void SetMinAndMax(float &min_x, float &min_y, float &max_x, float &max_y) {
        min_x_ = min_x;
        min_y_ = min_y;
        max_x_ = max_x;
        max_y_ = max_y;
    }

    /**
     * @brief 获取最细栅格分辨率。
     * @return 单个细栅格的物理尺寸。
     */
    float GetGridResolution() const { return options_.resolution_; }

   private:
    /**
     * @brief 将二维SubGrid索引转换为一维数组索引。
     * @param sx x方向尺寸。
     * @param x x方向索引。
     * @param y y方向索引。
     * @return 一维数组索引。
     */
    inline int MapIdx(int sx, int x, int y) { return (sx) * (y) + (x); }

   private:
    float grid_reso_ = 0.0;  ///< 单个SubGrid覆盖的地图边长。
    float min_x_ = 0;        ///< 地图最小x坐标。
    float min_y_ = 0;        ///< 地图最小y坐标。
    float max_x_ = 0;        ///< 地图最大x坐标。
    float max_y_ = 0;        ///< 地图最大y坐标。
    int grid_size_x_ = 0;    ///< SubGrid数组在x方向的尺寸。
    int grid_size_y_ = 0;    ///< SubGrid数组在y方向的尺寸。

    inline static const int sub_grid_width_ = (1 << SUB_GRID_SIZE);  ///< 单个SubGrid包含的细栅格宽度。

    Options options_;  ///< 地图参数。

    SubGrid **grids_ = nullptr;  ///< 实际存储SubGrid数据的位置。
};

/// G2P5Map共享指针类型。
typedef std::shared_ptr<G2P5Map> G2P5MapPtr;

}  // namespace lightning::g2p5

#endif  // LIGHTNING_G2P5_MAP_H

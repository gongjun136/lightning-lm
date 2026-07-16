//
// Created by xiang on 23-2-7.
//

#ifndef LIGHNING_TILED_MAP_H
#define LIGHNING_TILED_MAP_H

#include <opencv2/core/core.hpp>
#include <set>
#include <map>
#include <unordered_map>

#include "common/functional_points.h"
#include "common/point_def.h"
#include "common/std_types.h"
#include "core/lightning_math.hpp"
#include "core/maps/tiled_map_chunk.h"

namespace lightning {

/**
 * @brief 同时管理静态地图和动态地图的分块点云地图。
 *
 * 静态地图从文件目录中读取，动态地图实时从Lidar odom中构建。
 * 地图按二维chunk切分，用于按位姿加载附近地图并卸载远处地图。
 *
 * TODO:
 * - **done** 动态区域也需要分块，不然太大
 * - 动态区域应该有生命周期和最大采样点（目前在option中设置）
 * - 动态图层的保存与读取
 * - thread safety
 */
class TiledMap {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**
     * @brief 动态点云图层的生命周期策略。
     */
    enum class DynamicCloudPolicy {
        /// 短时间保留，车辆移出该区域就会清空
        SHORT = 0,

        /// 长时间保留在内存，车辆移出区域后仍然保留，但不会存盘并在下次开机时使用
        LONG,

        /// 更长时间的保留，会存盘，下次运行也会读取，超过一定时间后才会清除
        PERSISTENT,
    };

    /**
     * @brief 分块地图运行参数。
     */
    struct Options {
        Options() {}
        std::string map_path_ = "./data/maps/";                       ///< 地图存储路径。
        float chunk_size_ = 100.0;                                    ///< 单个chunk边长，单位m。
        float inv_chunk_size_ = 1.0 / 100.0;                          ///< chunk边长倒数，用于坐标到网格索引转换。
        float voxel_size_in_chunk_ = 0.1;                             ///< chunk内部点云体素滤波分辨率。
        const int load_nearby_chunks_ = 8;                            ///< 载入邻近范围内的点云数量上限。
        bool save_dynamic_layer_ = false;                             ///< 是否将动态图层也存入地图文件夹。
        int load_map_size_ = 2;                                       ///< 当前位置周围需要载入的网格半径。
        int unload_map_size_ = 3;                                     ///< 超出该网格半径的chunk会被卸载。
        bool enable_dynamic_polygon_ = false;                         ///< 是否启用动态区域多边形过滤。
        size_t max_pts_in_dyn_chunk_ = 50000;                         ///< 每个动态chunk中最多保留的点数。
        DynamicCloudPolicy policy_ = DynamicCloudPolicy::PERSISTENT;  ///< 动态区域处理策略。

        bool delete_when_unload_ = false;   ///< 卸载区域时是否清空点云。
        bool load_dyn_cloud_ = true;        ///< 初始化时是否加载动态点云。
        bool save_dyn_when_quit_ = true;    ///< 退出时是否保存动态点云。
        bool save_dyn_when_unload_ = true;  ///< 卸载时是否保存动态点云。
    };

    /**
     * @brief 动态区域多边形定义。
     */
    struct DynamicPolygon {
        DynamicPolygon() {}
        DynamicPolygon(int id, const Vec2d& pt) : id_(id) { polygon_.emplace_back(cv::Point2f(pt[0], pt[1])); }
        int id_ = 0;                         ///< 动态区域ID。
        std::vector<cv::Point2f> polygon_;   ///< 动态区域顶点，使用地图平面坐标。
    };

    /**
     * @brief NDT目标地图构建时选取的邻域范围。
     */
    enum class NearbyType {
        CENTER,   ///< 只使用当前中心体素。
        NEARBY6,  ///< 使用当前体素及其6邻域。
    };

    /**
     * @brief NDT栅格地图参数。
     */
    struct NdtMapOptions {
        NdtMapOptions() {}
        float static_cloud_weight_ = 1.0;              ///< 静态点云权重。
        float dynamic_cloud_weight_ = 0.1;             ///< 动态点云权重。
        NearbyType nearby_type_ = NearbyType::CENTER;  ///< 邻近类型。
        float voxel_size_ = 1.0;                       ///< voxel边长，单位m。
        int max_iteration_ = 5;                        ///< 最大迭代次数。
        double eps_ = 1e-2;                            ///< 收敛判定阈值。
        double res_outlier_th_ = 50.0;                 ///< 残差异常值拒绝阈值。
        float inv_voxel_size_ = 1.0;                   ///< voxel边长倒数。
        int min_pts_in_voxel_ = 3;                     ///< 每个voxel中估计NDT所需的最小点数。
        int min_effective_pts_ = 500;                  ///< 配准所需的最小有效点数。
        size_t capacity_ = 100000;                     ///< 缓存的体素数量上限。
        int max_pts_in_voxel_ = 50;                    ///< 每个voxel中参与估计的最大点数。

        bool verbose_ = false;  ///< 是否输出NDT调试信息。
    };

    using KeyType = Eigen::Matrix<int, 3, 1>;  ///< 体素三维索引。

    /**
     * @brief NDT体素内的点和统计量。
     */
    struct VoxelData {
        VoxelData() {}
        VoxelData(size_t id) { idx_.emplace_back(id); }
        VoxelData(const Vec3f& pt) {
            pts_.emplace_back(pt);
            num_pts_ = 1;
        }

        void AddPoint(const Vec3f& pt) {
            pts_.emplace_back(pt);
            if (!ndt_estimated_) {
                num_pts_++;
            }
        }

        std::vector<Vec3f> pts_;   ///< 内部点，多于一定数量之后再估计均值和协方差。
        std::vector<size_t> idx_;  ///< 点云中点的索引。
        Vec3f mu_ = Vec3f::Zero();         ///< 均值。
        Mat3f sigma_ = Mat3f::Identity();  ///< 协方差。
        Mat3f info_ = Mat3f::Zero();       ///< 协方差逆矩阵。

        bool ndt_estimated_ = false;  ///< NDT统计量是否已经估计。
        int num_pts_ = 0;             ///< 总点数，用于增量更新估计。
    };

    /// NDT体素哈希表，key为三维体素索引。
    using GridHashMap = std::unordered_map<KeyType, VoxelData, math::hash_vec<3>>;

    /// 分块地图哈希表，key为二维chunk网格索引。
    using ChunkHashMap = std::unordered_map<Vec2i, std::shared_ptr<MapChunk>, math::hash_vec<2>>;

    /**
     * @brief 构造分块地图。
     * @param options 分块地图运行参数。
     */
    TiledMap(Options options = Options()) : options_(options) {}

    /**
     * @brief 从单个完整PCD地图转换为分块地图。
     * @param map 输入的完整点云地图。
     * @param start_pose 地图起点位姿，用作分块地图原点参考。
     * @param map_path 分块地图保存目录。
     * @return 转换成功返回true，否则返回false。
     */
    bool ConvertFromFullPCD(CloudPtr map, const SE3& start_pose, const std::string& map_path = "./data/maps/");

    /// 载入地图索引文件。
    bool LoadMapIndex();

    /**
     * @brief 按当前位姿载入附近chunk，并卸载远处chunk。
     * @param pose 当前车辆或传感器位姿。
     */
    void LoadOnPose(const SE3& pose);

    // /// 载入某个位置附近的点云
    // void LoadOnPose(const SE3& pose, bool wndt = true);

    /// 获取当前载入的静态点云，按chunk ID返回深拷贝。
    std::map<int, CloudPtr> GetStaticCloud();

    /// 获取当前载入的动态点云，按chunk ID返回深拷贝。
    std::map<int, CloudPtr> GetDynamicCloud();

    /// 获取当前载入chunk拼接后的完整点云地图。
    CloudPtr GetAllMap();

    bool GetGlobalStaticBounds(Vec3d& min_point, Vec3d& max_point, std::size_t& point_count);

    /**
     * @brief 使用最近配准之后的scan更新动态图层。
     * @param cloud_world 世界坐标系下点云。
     * @param remove_old_pts 是否移除对应动态区域中的历史点。
     */
    void UpdateDynamicCloud(CloudPtr cloud_world, bool remove_old_pts = true);

    /// 重置动态图层。
    void ResetDynamicCloud();

    /// 获取地图原点。
    Vec3d GetOrigin() const { return origin_; }

    /// 添加地图功能点。
    void AddFP(const FunctionalPoint& fp) { func_points_.emplace_back(fp); }

    /// 获取所有功能点。
    std::vector<FunctionalPoint> GetAllFP() const { return func_points_; }

    /// 获取当前已载入chunk数量。
    int NumActiveChunks() const { return loaded_chunks_.size(); }

    /// 查询静态地图是否有更新。
    bool MapUpdated() const { return map_updated_; }

    /// 查询动态图层是否有更新。
    bool DynamicMapUpdated() const { return dynamic_map_updated_; }

    /// 清除静态地图和动态图层的更新标记。
    void CleanMapUpdate() {
        map_updated_ = false;
        dynamic_map_updated_ = false;
    }

    /// 更新静态点云地图。
    void AddStaticCloud(CloudPtr cloud);

    /// 更新动态点云地图。
    void AddDynamicCloud(CloudPtr cloud);

    /// 更新体素内部数据，根据新加入的pts和历史估计状态更新NDT统计量。
    void UpdateVoxel(VoxelData& v, bool flag_first_scan = true);

    /// 获取静态NDT栅格地图。
    GridHashMap GetStaticGridMap();

    /// 获取动态NDT栅格地图。
    GridHashMap GetDynamicGridMap();

    GridHashMap static_grids_;   ///< 静态地图栅格。
    GridHashMap dynamic_grids_;  ///< 动态地图栅格。

    /// 获取动态区域多边形集合，key为动态区域ID。
    std::map<int, DynamicPolygon> GetDynamicPolygons() const { return dynamic_polygon_; }

    /**
     * @brief 将地图chunk存储到二进制PCD文件。
     * @param only_dynamic true时只保存动态图层，false时保存静态地图和动态图层。
     */
    void SaveToBin(bool only_dynamic = false);

    /**
     * @brief 将当前载入chunk设置为NDT配准目标。
     * @tparam T NDT对象指针或智能指针类型，需要提供AddTarget/ComputeTargetGrids/initCompute接口。
     * @param ndt 待设置目标点云的NDT对象。
     */
    template <class T>
    void SetNewTargetForNDT(T ndt) {
        bool has_cloud = false;
        {
            UL lock(static_data_mutex_);
            /// 静态chunk始终按当前loaded_chunks_集合加入目标地图。
            for (auto& idx : loaded_chunks_) {
                if (!static_chunks_[idx]->cloud_->empty()) {
                    has_cloud = true;
                }
                ndt->AddTarget(static_chunks_[idx]->cloud_);
            }
        }

        {
            UL lock(dynamic_data_mutex_);
            /// 动态chunk可能不存在或未加载，只有有效点云才加入目标地图。
            for (auto& idx : loaded_chunks_) {
                if (dynamic_chunks_.find(idx) != dynamic_chunks_.end()) {
                    if (dynamic_chunks_[idx]->cloud_ != nullptr) {
                        has_cloud = true;
                        ndt->AddTarget(dynamic_chunks_[idx]->cloud_);
                    }
                }
            }
        }

        /// 先计算目标栅格；只有存在有效点云时才初始化配准计算。
        ndt->ComputeTargetGrids();
        if (has_cloud) {
            ndt->initCompute();
        }
    }

    /// 清理地图。
    void ClearMap();

   private:
    /**
     * @brief 测试某个点是否落在动态区域。
     * @param pt 给定点，世界坐标系。
     * @return 动态区域ID；未落入动态区域时返回实现定义的无效ID。
     */
    int FallsInDynamicArea(const Vec3f& pt);

   private:
    /**
     * @brief 将float二维平面坐标转换为chunk网格索引。
     * @param pt 世界坐标系下的二维点。
     * @return 所在chunk的二维网格索引。
     */
    inline Vec2i Pos2Grid(const Vec2f& pt) const {
        Vec2d p = ((pt.cast<double>() * options_.inv_chunk_size_) + Vec2d(0.5, 0.5));
        return Vec2i(floor(p[0]), floor(p[1]));
    }

    /**
     * @brief 将double二维平面坐标转换为chunk网格索引。
     * @param pt 世界坐标系下的二维点。
     * @return 所在chunk的二维网格索引。
     */
    inline Vec2i Pos2Grid(const Vec2d& pt) const {
        Vec2d p = ((pt * options_.inv_chunk_size_) + Vec2d(0.5, 0.5));
        return Vec2i(floor(p[0]), floor(p[1]));
    }

    Options options_;                ///< 分块地图运行参数。
    NdtMapOptions ndt_map_options_;  ///< NDT栅格地图参数。

    std::mutex static_data_mutex_;   ///< 静态数据锁。
    std::mutex dynamic_data_mutex_;  ///< 动态数据锁。

    Vec3d origin_ = Vec3d::Zero();     ///< 地图原点。
    ChunkHashMap static_chunks_;       ///< 静态地图chunk集合。
    ChunkHashMap dynamic_chunks_;      ///< 动态地图chunk集合。
    std::map<int, Vec2i> id_to_grid_;  ///< chunk ID到二维网格索引的映射。
    int chunk_id_ = 0;                 ///< chunk生成时分配的递增ID。

    Vec2i last_load_grid_ = Vec2i::Zero();  ///< 上次触发加载时所在的网格。
    bool last_load_grid_set_ = false;       ///< 是否已有上次加载网格记录。

    std::set<Vec2i, math::less_vec<2>> loaded_chunks_;  ///< 已经载入的chunk，以网格索引记录。
    bool map_updated_ = false;                          ///< 静态地图是否更新。
    bool dynamic_map_updated_ = false;                  ///< 动态地图是否更新。

    std::map<int, DynamicPolygon> dynamic_polygon_;  ///< 从txt读取的动态区域多边形。

    std::vector<FunctionalPoint> func_points_;  ///< 功能点集合。

    bool flag_first_dynamic_scan_ = true;  ///< 首帧动态点云特殊处理标记。
    bool flag_first_static_scan_ = true;   ///< 首帧静态点云特殊处理标记。
};

}  // namespace lightning

#endif

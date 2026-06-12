//
// Created by xiang on 25-6-23.
//

#include "core/g2p5/g2p5.h"
#include "common/constant.h"

#include <pcl/ModelCoefficients.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <map>
#include <pcl/filters/impl/voxel_grid.hpp>
#include <pcl/segmentation/impl/sac_segmentation.hpp>

#include "utils/timer.h"
#include "yaml-cpp/yaml.h"

namespace lightning::g2p5 {

G2P5::~G2P5() { Quit(); }

void G2P5::Quit() {
    quit_flag_ = true;

    if (options_.online_mode_) {
        draw_frontend_map_thread_.Quit();
    }

    if (draw_backend_map_thread_.joinable()) {
        draw_backend_map_thread_.join();
    }
}

/**
 * @brief 加入关键帧并触发前端地图渲染。
 * @param kf 待加入地图的关键帧。
 */
void G2P5::PushKeyframe(Keyframe::Ptr kf) {
    /// 保存完整关键帧序列，供后端全局重绘时重新构建地图。
    UL lock(kf_mutex_);
    all_keyframes_.emplace_back(kf);

    /// 在线模式下交给异步前端线程；离线模式下直接同步渲染，便于顺序回放和调试。
    if (options_.online_mode_) {
        draw_frontend_map_thread_.AddMessage(kf);
    } else {
        RenderFront(kf);
    }
}

/**
 * @brief 将单个关键帧增量渲染到前端地图，并发布最新地图。
 * @param kf 待渲染的关键帧。
 */
void G2P5::RenderFront(Keyframe::Ptr kf) {
    {
        /// 记录前端当前处理到的关键帧，后端重绘完成后会用它追齐前端进度。
        UL lock(frontend_mutex_);
        frontend_current_ = kf;
    }

    {
        /// 更新前端地图和对外可见的最新地图时使用同一把锁，避免外部读到中间状态。
        UL lock{newest_map_mutex_};

        lightning::Timer::Evaluate([&]() { AddKfToMap({kf}, frontend_map_); }, "G2P5 Occupancy Mapping", true);
        newest_map_ = frontend_map_;
    }

    /// 向外回调
    if (map_update_cb_) {
        map_update_cb_(newest_map_);
    }
}

void G2P5::RedrawGlobalMap() { backend_redraw_flag_ = true; }

void G2P5::RenderBack() {
    while (!quit_flag_) {
        while (!backend_redraw_flag_ && !quit_flag_) {
            sleep(1);
        }

        if (quit_flag_) {
            break;
        }

        is_busy_ = true;

        /// 后端重绘被触发
        backend_redraw_flag_ = false;

        /// 重新绘制整张地图，如果中间过程又被重绘了，则退出
        std::vector<Keyframe::Ptr> all_keyframes;
        {
            UL lock(kf_mutex_);
            all_keyframes = all_keyframes_;
        }

        if (all_keyframes.empty()) {
            is_busy_ = false;
            continue;
        }

        G2P5Map::Options opt;
        opt.resolution_ = options_.grid_map_resolution_;
        backend_map_ = std::make_shared<G2P5Map>(opt);
        auto cur_kf = all_keyframes.begin();
        bool abort = false;

        for (; cur_kf != all_keyframes.end(); ++cur_kf) {
            AddKfToMap({*cur_kf}, backend_map_);
            if (backend_redraw_flag_) {
                LOG(INFO) << "backend redraw triggered in process, abort";
                abort = true;
                break;
            }

            if (quit_flag_) {
                abort = true;
                break;
            }
        }

        if (abort) {
            /// 继续重绘
            is_busy_ = false;
            continue;
        }

        /// 绘制过程中前端可能发生了更新，要保证后端绘制和前端的一致性
        int cur_idx = all_keyframes.back()->GetID();
        while (true) {
            Keyframe::Ptr frontend_kf = nullptr;
            {
                UL lock(frontend_mutex_);
                frontend_kf = frontend_current_;
            }

            if (cur_idx == frontend_kf->GetID()) {
                break;
            }

            int frontend_idx = frontend_kf->GetID();
            std::vector<Keyframe::Ptr> kfs;

            {
                UL lock(kf_mutex_);
                for (int i = cur_idx + 1; i <= frontend_idx; ++i) {
                    kfs.emplace_back(all_keyframes_[i]);
                }
            }

            AddKfToMap(kfs, backend_map_);
            cur_idx = frontend_idx;
        }

        {
            /// 同步前后端地图，替换newest map
            UL lock{newest_map_mutex_};
            frontend_map_ = backend_map_;
            newest_map_ = frontend_map_;
        }

        /// 向外回调
        if (map_update_cb_) {
            map_update_cb_(newest_map_);
        }

        is_busy_ = false;
    }

    LOG(INFO) << "backend render quit";
}

/**
 * @brief 根据关键帧位姿和有效点云范围扩展栅格地图边界。
 * @param[in] kfs 用于估计地图边界的关键帧集合。
 * @param[in,out] map 待扩展的栅格地图。
 * @return 地图边界有效时返回true，否则返回false。
 */
bool G2P5::ResizeMap(const std::vector<Keyframe::Ptr> &kfs, G2P5MapPtr &map) {
    /// 以当前地图边界为初值，只在新关键帧超出范围时扩展。
    float init_min_x, init_min_y, init_max_x, init_max_y;
    float min_x, min_y, max_x, max_y;
    map->GetMinAndMax(init_min_x, init_min_y, init_max_x, init_max_y);
    min_x = init_min_x;
    min_y = init_min_y;
    max_x = init_max_x;
    max_y = init_max_y;

    /// 遍历待加入的关键帧，将关键帧位置和采样后的点云位置都纳入边界估计。
    for (auto it = kfs.begin(); it != kfs.end(); ++it) {
        if (quit_flag_) {
            return true;
        }

        auto kf = *it;
        SE3 pose = kf->GetOptPose();
        auto cloud = kf->GetCloud();

        if (pose.translation().x() < min_x) {
            min_x = pose.translation().x();
        }

        if (pose.translation().y() < min_y) {
            min_y = pose.translation().y();
        }

        if (pose.translation().x() > max_x) {
            max_x = pose.translation().x();
        }

        if (pose.translation().y() > max_y) {
            max_y = pose.translation().y();
        }

        if (min_x > max_x || min_y > max_y) {
            return false;
        }

        /// 点云只做稀疏采样来估计边界，避免扩图检查成为前端渲染瓶颈。
        for (size_t i = 0; i < cloud->points.size(); i += 10) {
            float range = cloud->points[i].getVector3fMap().norm();

            if (range > options_.usable_scan_range_ || range <= 0.01 || std::isnan(range)) {
                continue;
            }

            Vec3d point = pose * cloud->points[i].getVector3fMap().cast<double>();

            if ((point.x() - 1) < min_x) {
                min_x = point.x() - 1;
            }
            if ((point.y() - 1) < min_y) {
                min_y = point.y() - 1;
            }
            if ((point.x() + 1) > max_x) {
                max_x = point.x() + 1;
            }
            if ((point.y() + 1) > max_y) {
                max_y = point.y() + 1;
            }
        }
    }

    if (min_x > max_x || min_y > max_y) {
        return false;
    }

    /// 如果估计边界超出现有地图，则按栅格分辨率对齐后扩展地图。
    if (min_x < init_min_x || min_y < init_min_y || max_x > init_max_x || max_y > init_max_y) {
        min_x = (min_x > init_min_x) ? init_min_x : min_x;
        min_y = (min_y > init_min_y) ? init_min_y : min_y;
        max_x = (max_x < init_max_x) ? init_max_x : max_x;
        max_y = (max_y < init_max_y) ? init_max_y : max_y;

        float r = map->GetGridResolution();
        min_x = static_cast<int>((floor)(min_x / r)) * r;
        min_y = static_cast<int>((floor)(min_y / r)) * r;
        max_x = static_cast<int>((ceil)(max_x / r)) * r;
        max_y = static_cast<int>((ceil)(max_y / r)) * r;
        map->Resize(min_x, min_y, max_x, max_y);

        // LOG(INFO) << "map resized to " << min_x << ", " << min_y << ", " << max_x << ", " << max_y;
    }

    return true;
}

/**
 * @brief 将关键帧集合增量写入栅格地图。
 * @param[in] kfs 待写入地图的关键帧集合。
 * @param[in,out] map 待更新的栅格地图。
 * @return 当前实现固定返回true。
 */
bool G2P5::AddKfToMap(const std::vector<Keyframe::Ptr> &kfs, G2P5MapPtr &map) {
    /// 如果需要，更新地图大小
    ResizeMap(kfs, map);

    /// 地图边界准备好后，逐帧执行三维点云到二维scan的投影与栅格更新。
    for (const auto &kf : kfs) {
        Convert3DTo2DScan(kf, map);
    }

    return true;
}

G2P5MapPtr G2P5::GetNewestMap() {
    UL lock{newest_map_mutex_};
    LOG(INFO) << "getting newest map";
    if (newest_map_ == nullptr) {
        return nullptr;
    }

    return newest_map_->MakeDeepCopy();
}

/**
 * @brief 将关键帧点云按角度聚合为二维scan，并据此更新占用栅格。
 * @param kf 待转换的关键帧。
 * @param[in,out] map 待更新的栅格地图。
 */
void G2P5::Convert3DTo2DScan(Keyframe::Ptr kf, G2P5MapPtr &map) {
    /// 根据配置选择动态估计地面，或使用默认水平地面模型。
    if (options_.esti_floor_) {
        if (!DetectPlaneCoeffs(kf)) {
            /// 如果动态检测失败，就回退到默认水平地面参数。
            floor_coeffs_ = Vec4d(0, 0, 1, -options_.default_floor_height_);
        } else {
            if (options_.verbose_) {
                LOG(INFO) << "floor coeffs: " << floor_coeffs_.transpose();
            }
        }
    } else {
        floor_coeffs_ = Vec4d(0, 0, 1, -options_.default_floor_height_);
    }

    /// Step 1: 统计360个整数角度方向上的射线高度分布。
    /// NOTE: 角度被round到整数度，存在角分辨率损失，但能换取较快的渲染速度。
    std::vector<std::map<double, double>> rays(360);               ///< map键值：距离-相对高度，以距离排序。
    std::vector<Vec2d> angle_distance_height(360, Vec2d::Zero());  ///< 每个角度最终选中的距离-高度值。
    std::vector<Vec3d> pts_3d;                                     ///< 有效障碍物点，雷达坐标系下。

    SE3 Twb = kf->GetOptPose();

    double min_th = options_.min_th_floor_;
    double max_th = options_.max_th_floor_;

    /// 关键帧点云在雷达坐标系下，写栅格时通过Twb转换到地图坐标系。
    auto cloud = kf->GetCloud();

    CloudPtr obstacle_cloud(new PointCloudType);

    /// 障碍物高度带内的点会写为hit点，同时也参与各角度射线末端估计。
    int cnt_valid = 0;
    for (size_t i = 0; i < cloud->points.size(); ++i) {
        const auto &pt = cloud->points[i];
        if (quit_flag_) {
            return;
        }

        Vec3d pc = Vec3d(pt.x, pt.y, pt.z);
        Vec4d pn = Vec4d(pt.x, pt.y, pt.z, 1);

        /// 计算激光点所在角度、水平距离以及到地面的有符号距离。
        Vec2d p = pc.head<2>();
        double dis = p.norm();

        if (dis > options_.usable_scan_range_) {
            continue;
        }

        double dis_floor = pn.dot(floor_coeffs_);  /// 该点到地面的距离
        double dangle = atan2(p[1], p[0]) * constant::kRAD2DEG;
        int angle = int(round(dangle) + 360) % 360;

        if (dis_floor > min_th) {
            if (dis_floor < max_th) {
                /// 特别矮的和特别高的都不计入障碍物。
                pts_3d.emplace_back(pc);

                rays[angle].insert({dis, dis_floor});

                /// 障碍物点直接写为hit点。
                Vec3d p_world = Twb * pc;
                map->SetHitPoint(p_world[0], p_world[1], true, dis_floor);

                cnt_valid++;

                obstacle_cloud->points.push_back(pt);
            }
        } else if (dis_floor > -min_th) {
            /// 地面附近的点不写hit，但保留在射线分布中，用于判断可通行段。
            rays[angle].insert({dis, dis_floor});
            cnt_valid++;
        }
    }

    // if (options_.verbose_) {
    //     obstacle_cloud->is_dense = false;
    //     obstacle_cloud->height = 1;
    //     obstacle_cloud->width = obstacle_cloud->size();

    //     pcl::io::savePCDFile("./data/obs.pcd", *obstacle_cloud);
    // }

    // LOG(INFO) << "valid obs: " << cnt_valid << ", total: " << cloud->size();

    std::vector<double> floor_esti_data;  // 地面高度估计值

    /// Step 2: 考察每个方向上的高度分布曲线。
    /// 正常场景中，每个方向由较低高度开始（地面），转到较高高度（物体）；
    /// 若测量点不足，则认为该方向可能被遮挡或进入盲区。
    constexpr double default_ray_distance = -1;
    const double floor_rh = floor_coeffs_[3];

    for (int i = 0; i < 360; ++i) {
        if (quit_flag_) {
            return;
        }

        if (rays[i].size() < 2) {
            /// 该方向测量数据很少，认为有严重遮挡或失效，使用默认最小距离。
            angle_distance_height[i] = Vec2d(default_ray_distance, floor_rh);
            continue;
        }

        /// 从远到近寻找地面到障碍物的跳变位置，作为该方向的scan端点。
        for (auto iter = rays[i].rbegin(); iter != rays[i].rend(); ++iter) {
            if (iter->second < options_.min_th_floor_) {
                angle_distance_height[i] = Vec2d(iter->first, iter->second);
                continue;
            }

            auto next_iter = iter;
            next_iter++;

            if (next_iter != rays[i].rend()) {
                if (iter->second > options_.min_th_floor_ && next_iter->second < options_.min_th_floor_) {
                    /// 当前点是障碍、下一个更近点不是障碍，说明这里是可通行段尽头。
                    angle_distance_height[i] = Vec2d(iter->first, iter->second);
                    break;
                }
            } else {
                angle_distance_height[i] = Vec2d(iter->first, iter->second);
            }
        }
    }

    /// Step 3: 以二维scan方式沿射线写入miss点，刷新可通行区域。
    SetWhitePoints(angle_distance_height, kf, map);
}

/**
 * @brief 沿二维scan射线将可通行区域写入地图。
 * @param[in] pt2d 每个整数角度方向的射线端点，元素格式为距离-高度。
 * @param[in] kf 提供当前雷达位姿的关键帧。
 * @param[in,out] map 待写入miss点的栅格地图。
 *
 * 根据每个角度方向选出的scan端点，将雷达原点到端点之间的区域标记为miss。
 * 无有效测量或超出可用范围的方向会被跳过；距离很近的默认端点仍用于清理车体附近区域。
 */
void G2P5::SetWhitePoints(const std::vector<Vec2d> &pt2d, Keyframe::Ptr kf, G2P5MapPtr &map) {
    assert(pt2d.size() == 360);

    SE3 pose = kf->GetOptPose();
    Vec3d orig = pose.translation();

    for (int i = 0; i < 360; ++i) {
        if (quit_flag_) {
            return;
        }

        double angle = float(i) * constant::kDEG2RAD;
        float r = pt2d[i][0];
        float h = pt2d[i][1];

        Vec3d p_local(r * cos(angle), r * sin(angle), h);
        Vec3d p_world = pose * p_local;

        /// 无测量值或超过可用scan范围时，认为该方向无效。
        if (r <= 0 || r > options_.usable_scan_range_) {
            /// 默认近距离端点仍用于涂白车体附近区域。
            if (r < 0.1) {
                map->SetMissPoint(p_world[0], p_world[1], orig[0], orig[1], h, options_.lidar_height_);
            }
            continue;
        }

        map->SetMissPoint(p_world[0], p_world[1], orig[0], orig[1], h, options_.lidar_height_);
    }
}

/**
 * @brief 从关键帧点云中用RANSAC估计地面平面系数。
 * @param kf 用于估计地面的关键帧。
 * @return 成功估计到近似水平且内点数量足够的地面时返回true，否则返回false。
 */
bool G2P5::DetectPlaneCoeffs(Keyframe::Ptr kf) {
    /// 使用PCL RANSAC平面模型拟合候选地面点。
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
    pcl::SACSegmentation<PointType> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_PLANE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setDistanceThreshold(0.25);

    /// 只取默认地面高度附近、低于雷达安装高度的点作为地面候选点。
    CloudPtr cloud(new PointCloudType);
    for (auto &pt : kf->GetCloud()->points) {
        if (pt.z < options_.lidar_height_ + options_.default_floor_height_) {
            cloud->points.push_back(pt);
        }
    }

    cloud->width = cloud->points.size();
    cloud->height = 1;
    cloud->is_dense = false;

    /// 候选点过少时不做拟合，避免得到不稳定的地面参数。
    if (cloud->size() < 200) {
        // LOG(ERROR) << "not enough points cloud->size(): " << cloud->size();
        return false;
    }

    // if (options_.verbose_) {
    //     pcl::io::savePCDFile("./data/floor_candi.pcd", *cloud);
    // }

    seg.setInputCloud(cloud);
    seg.segment(*inliers, *coefficients);

    /// 当前G2P5假设雷达水平安装，因此只接受法向接近z轴的地面平面。
    if (coefficients->values[2] < 0.99) {
        LOG(ERROR) << "floor is not horizontal. ";
        return false;
    }

    /// 内点数量过少说明平面证据不足，继续沿用默认或上一轮地面参数。
    if (inliers->indices.size() < 100) {
        LOG(ERROR) << "cannot get enough points on floor: " << inliers->indices.size();
        return false;
    }

    /// 保存平面方程系数，后续用点到平面的距离判断障碍物高度。
    for (int i = 0; i < 4; ++i) {
        floor_coeffs_[i] = coefficients->values[i];
    }
    cloud->clear();

    return true;
}

void G2P5::Init(std::string yaml_path) {
    auto yaml = YAML::LoadFile(yaml_path);
    options_.esti_floor_ = yaml["g2p5"]["esti_floor"].as<bool>();
    options_.min_th_floor_ = yaml["g2p5"]["min_th_floor"].as<float>();
    options_.max_th_floor_ = yaml["g2p5"]["max_th_floor"].as<float>();
    options_.lidar_height_ = yaml["g2p5"]["lidar_height"].as<float>();
    options_.grid_map_resolution_ = yaml["g2p5"]["grid_map_resolution"].as<float>();
    options_.default_floor_height_ = yaml["g2p5"]["floor_height"].as<float>();

    G2P5Map::Options opt;
    opt.resolution_ = options_.grid_map_resolution_;
    frontend_map_ = std::make_shared<G2P5Map>(opt);

    if (options_.online_mode_) {
        draw_frontend_map_thread_.SetProcFunc([this](Keyframe::Ptr kf) { RenderFront(kf); });
        draw_frontend_map_thread_.Start();
    }

    draw_backend_map_thread_ = std::thread([this]() { RenderBack(); });
}

}  // namespace lightning::g2p5

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

Lightning-LM 是一个极速激光雷达建图与定位系统（Lightning-Speed Lidar SLAM），基于AA-FasterLIO的快速LIO前端，支持实时回环检测、高精度定位、3D到2D地图转换等功能。该系统纯CPU环境下即可运行，资源占用极低。

## 构建和依赖管理

### 环境要求
- Ubuntu 22.04+ (Ubuntu 20.04 未测试但应该可行)
- ROS2 Humble+

### 依赖安装
```bash
# 安装所有必需依赖
./scripts/install_dep.sh
```

主要依赖包括：
- ros2 humble及以上
- Pangolin (用于可视化)
- OpenCV, PCL, yaml-cpp, glog, gflags, pcl_conversions

### 编译命令
```bash
# 编译项目
colcon build

# 设置环境
source install/setup.bash
```

### 编译输出
编译后会生成以下可执行文件：
- `run_slam_online`: 在线SLAM建图
- `run_slam_offline`: 离线SLAM建图
- `run_loc_online`: 在线定位
- `run_loc_offline`: 离线定位
- `run_frontend_offline`: 前端测试
- `run_loop_offline`: 回环检测测试
- `test_ui`: UI测试

## 系统架构

### 核心模块 (`src/core/`)
- **lio/**: 激光雷达惯性里程计前端，基于AA-FasterLIO
  - `eskf.hpp`: 扩展卡尔曼滤波器
  - `laser_mapping.cc`: 激光建图主模块
  - `pointcloud_preprocess.cc`: 点云预处理和运动畸变校正
- **localization/**: 定位模块
  - `localization.cpp`: 定位系统主接口
  - `lidar_loc/`: 基于NDT的激光定位算法
  - `pose_graph/`: 位姿图优化(PGO)
- **loop_closing/**: 基于NDT的回环检测模块
- **maps/**: 地图管理，支持分块动态加载
  - `tiled_map.cc`: 分块地图管理
- **g2p5/**: 3D到2D栅格地图转换模块
- **ivox3d/**: 高效点云索引结构
- **miao/**: 轻量图优化库(类似g2o，支持增量优化)
- **system/**: 系统集成接口
  - `slam.cc`: SLAM系统主接口
  - `loc_system.cc`: 定位系统接口

### 数据流程
1. **SLAM建图**: IMU+激光雷达 → 点云预处理 → LIO前端 → 关键帧检测 → 地图构建 → 回环检测 → 后端优化 → 2D地图生成
2. **定位**: 地图加载 → 初始位姿估计 → NDT精确定位 → 动态物体处理 → 位姿图优化 → 高频位姿输出

## 配置系统

### 配置文件位置 (`config/`)
- `default.yaml`: 默认配置
- `default_nclt.yaml`: NCLT数据集配置
- `default_vbr.yaml`: VBR数据集配置
- `default_utbm.yaml`: UTBM数据集配置

### 关键配置参数
```yaml
# 系统开关
system:
  with_loop_closing: true    # 回环检测
  with_ui: true              # 3D可视化
  with_2dui: false           # 2D可视化
  with_g2p5: true            # 2D栅格地图
  map_path: "new_map"        # 地图存储路径

# LIO参数
fasterlio:
  lidar_type: 2              # 1:Livox, 2:Velodyne, 3:Ouster
  point_filter_num: 10       # 点云采样数
  use_aa: true               # Anderson加速
  extrinsic_T: [0, 0, 0.28]  # IMU-LiDAR外参

# 回环检测
loop_closing:
  loop_kf_gap: 20            # 关键帧检查间隔
  with_height: true          # 高度约束(防Z轴飘移)

# 定位参数
lidar_loc:
  init_with_fp: true         # 特征初始化
  update_dynamic_cloud: true # 动态点云更新
```

## 运行方式

### 建图
```bash
# 在线建图
ros2 run lightning run_slam_online --config ./config/default_nclt.yaml

# 离线建图
ros2 run lightning run_slam_offline --input_bag bag_file --config ./config/default_nclt.yaml

# 保存地图 (服务调用)
ros2 service call /lightning/save_map lightning/srv/SaveMap "{map_id: new_map}"
```

### 定位
```bash
# 在线定位
ros2 run lightning run_loc_online --config ./config/default_nclt.yaml

# 离线定位
ros2 run lightning run_loc_offline --input_bag bag_file --config ./config/default_nclt.yaml
```

### 地图查看
```bash
# 查看完整点云地图
pcl_viewer ./data/new_map/global.pcd

# 查看栅格地图
# 地图保存在 data/new_map/map.pgm
```

## 支持的激光雷达类型
- **Livox系列** (`lidar_type: 1`)
- **Velodyne** (`lidar_type: 2`)
- **Ouster** (`lidar_type: 3`)

## 新设备调试要点

1. **设置激光雷达类型**：在配置文件中设置正确的 `fasterlio.lidar_type`
2. **配置话题名称**：修改 `common.lidar_topic` 和 `common.imu_topic`
3. **时间戳处理**：关注 `fasterlio.time_scale` 参数，确保点云时间戳正确
4. **外参标定**：IMU和LiDAR外参默认为零即可，系统对此不敏感
5. **建议调试流程**：先录包 → 调通离线模式 → 再调试在线模式

## 性能特点
- **低资源占用**: 在线定位0.8核，建图1.2核（32线激光雷达，无UI）
- **高频输出**: 100Hz IMU频率的位姿输出
- **大规模场景**: 支持地图分区动态加载
- **实时性**: 支持多线程并发和动态跳帧

## 开发指南
- **新功能开发**: 主要在`src/core/`下添加新模块
- **算法调参**: 通过配置文件调整，避免重新编译
- **调试**: 使用离线模式进行断点调试
- **数据支持**: 通过`wrapper/bag_io.h`添加新数据格式支持
- **可视化**: 通过`ui/`模块扩展可视化功能
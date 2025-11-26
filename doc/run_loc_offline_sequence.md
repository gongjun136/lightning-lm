# run_loc_offline.cc 流程时序图

```mermaid
sequenceDiagram
    participant Main as main()
    participant Config as 配置解析
    participant Localization as 定位系统
    participant RosbagIO as ROSbag读取器
    participant IMU as IMU数据处理
    participant Lidar as 激光数据处理
    participant Timer as 计时器

    Note over Main,Timer: 程序启动阶段
    Main->>Main: 初始化Google日志
    Main->>Main: 设置日志参数
    Main->>Main: 解析命令行参数
    Main->>Main: 检查输入数据包参数

    alt 输入数据包为空
        Main->>Main: 输出错误日志
        Main->>Main: 返回-1退出程序
    end

    Note over Main,Localization: 初始化阶段
    Main->>Localization: 创建RosbagIO对象
    Main->>Localization: 设置定位选项(offline模式)
    Main->>Localization: 创建Localization实例
    Main->>Localization: 调用Init(config, map_path)

    Note over Localization,Config: 读取配置文件
    Localization->>Config: 读取YAML配置文件
    Config-->>Localization: 获取lidar_topic和imu_topic

    Note over Main,RosbagIO: 注册回调函数阶段
    Main->>RosbagIO: 添加IMU数据处理回调
    Main->>RosbagIO: 添加PointCloud2处理回调
    Main->>RosbagIO: 添加Livox激光数据处理回调

    Note over RosbagIO,Main: 数据处理循环
    RosbagIO->>RosbagIO: 开始处理ROS bag数据

    loop 遍历所有消息
        alt 消息类型为IMU
            RosbagIO->>IMU: 调用IMU回调函数
            IMU->>Localization: ProcessIMUMsg(imu)
            IMU->>IMU: usleep(1000)
            IMU-->>RosbagIO: 返回true继续处理
        else 消息类型为PointCloud2
            RosbagIO->>Lidar: 调用激光点云回调函数
            Lidar->>Localization: ProcessLidarMsg(cloud)
            Lidar->>Lidar: usleep(1000)
            Lidar-->>RosbagIO: 返回true继续处理
        else 消息类型为Livox
            RosbagIO->>Lidar: 调用Livox激光回调函数
            Lidar->>Localization: ProcessLivoxLidarMsg(cloud)
            Lidar->>Lidar: usleep(1000)
            Lidar-->>RosbagIO: 返回true继续处理
        end
    end

    Note over Main,Timer: 结束阶段
    RosbagIO-->>Main: 数据处理完成
    Main->>Timer: 打印所有计时统计
    Main->>Localization: 调用Finish()清理资源
    Main->>Main: 输出"done"日志
    Main->>Main: 返回0正常退出
```

## 关键组件说明

### 1. 初始化阶段
- **日志初始化**: 设置Google Logging参数，启用彩色日志输出
- **参数解析**: 解析命令行参数，包括输入bag路径、配置文件路径、地图路径
- **参数验证**: 检查必需的输入bag参数是否存在

### 2. 定位系统初始化
- **RosbagIO创建**: 用于读取ROS bag数据的接口
- **Localization创建**: 核心定位系统实例，设置为离线模式
- **系统初始化**: 读取配置文件和地图数据

### 3. 回调函数注册
- **IMU回调**: `ProcessIMUMsg()` - 处理IMU数据
- **PointCloud2回调**: `ProcessLidarMsg()` - 处理标准激光点云数据
- **Livox回调**: `ProcessLivoxLidarMsg()` - 处理Livox激光雷达数据

### 4. 数据处理循环
- 按时间顺序遍历ROS bag中的所有消息
- 根据消息类型调用相应的回调函数
- 每次处理后休眠1ms以控制处理速度

### 5. 程序结束
- 打印性能统计信息
- 清理定位系统资源
- 正常退出程序

## 数据流向

```
ROS Bag → RosbagIO → 回调函数 → Localization → 定位结果
```

每个传感器消息都会按照时间顺序被处理，IMU数据和激光数据在Localization系统中进行融合，最终输出定位结果。
# SLAM 建图定位方案技术文档

## 1. 概述

本文面向一个以 LiDAR-IMU 融合为核心的 SLAM 建图定位系统，给出建图、回环、地图表达、重定位和在线定位的整体技术方案。系统以激光雷达点云提供几何约束，以 IMU 提供高频运动预测和点云运动畸变补偿，通过误差状态迭代卡尔曼滤波实现前端里程计，通过关键帧、回环检测和位姿图优化保证全局地图一致性，并在定位阶段利用分块点云地图、NDT 匹配和滑窗位姿图融合输出高频稳定位姿。

系统的核心目标如下：

1. **实时建图**：在 LiDAR 与 IMU 输入下估计连续位姿，生成全局一致的三维点云地图。
2. **全局一致性**：通过关键帧回环检测与位姿图优化抑制前端累积漂移。
3. **地图可用性**：同时输出全局 PCD、分块点云地图和可选 2D 栅格地图。
4. **在线定位**：加载已有地图，结合 LIO 短时预测与全局点云匹配得到稳定定位结果。
5. **大场景适应**：通过地图分块动态加载、动态层更新和高频外推降低内存与延迟压力。

## 2. 坐标系与状态定义

本文使用以下坐标系：

| 坐标系 | 含义 |
| --- | --- |
| $\mathcal{F}_w$ | 世界坐标系，也是建图和定位输出的全局参考系 |
| $\mathcal{F}_i$ | IMU 坐标系 |
| $\mathcal{F}_l$ | LiDAR 坐标系 |
| $\mathcal{F}_{l_j}$ | 当前扫描中第 $j$ 个点采样时刻的 LiDAR 坐标系 |
| $\mathcal{F}_{l_e}$ | 当前扫描结束时刻的 LiDAR 坐标系 |

旋转矩阵 $\mathbf{R}_{wi}$ 表示从 $\mathcal{F}_i$ 到 $\mathcal{F}_w$ 的旋转；平移 ${}^{w}\mathbf{t}_{wi}$ 表示从 $o_w$ 指向 $o_i$ 的向量，在 $\mathcal{F}_w$ 下表达。LiDAR 到 IMU 的外参写作 $\mathbf{T}_{il}=(\mathbf{R}_{il},{}^{i}\mathbf{t}_{il})$。

前端名义状态采用惯导导航状态：

$$
\hat{\mathbf{x}} =
\left(
\hat{\mathbf{R}}_{wi},
{}^{w}\hat{\mathbf{p}},
{}^{w}\hat{\mathbf{v}},
{}^{i}\hat{\mathbf{b}}_g,
{}^{w}\hat{\mathbf{g}}
\right).
$$

其中 ${}^{i}\hat{\mathbf{b}}_g$ 为陀螺零偏，${}^{w}\hat{\mathbf{g}}$ 为世界系重力向量。当前框架的在线前端以位置、姿态、速度、陀螺零偏为主要可修正误差块，重力由静态初始化给定并在运行中保持稳定；LiDAR-IMU 外参和加速度计零偏采用配置或离线标定结果，不作为默认在线估计量。这样的状态设计牺牲一部分自标定自由度，换取更小的状态维度和更稳定的实时估计。

误差状态采用切空间小扰动：

$$
\tilde{\mathbf{x}} =
\begin{bmatrix}
\delta{}^{w}\mathbf{p}^{\mathsf{T}} &
\delta\boldsymbol{\theta}^{\mathsf{T}} &
\delta{}^{w}\mathbf{v}^{\mathsf{T}} &
\delta{}^{i}\mathbf{b}_g^{\mathsf{T}}
\end{bmatrix}^{\mathsf{T}},
$$

并通过广义加法注入名义状态：

$$
\mathbf{x} = \hat{\mathbf{x}} \boxplus \tilde{\mathbf{x}},
\quad
\mathbf{R}_{wi} =
\hat{\mathbf{R}}_{wi}\operatorname{Exp}(\delta\boldsymbol{\theta}).
$$

## 3. 系统总体架构

系统分为建图链路和定位链路。建图链路侧重生成高质量地图，定位链路侧重利用已有地图输出稳定低延迟位姿。两者共享 LiDAR-IMU 前端、点云预处理、地图表达和位姿图优化等基础能力。

```mermaid
flowchart TB
    classDef input fill:#EEF6FF,stroke:#4C7EDB,color:#172B4D;
    classDef front fill:#FFF8E8,stroke:#D99A25,color:#3D2A00;
    classDef mapping fill:#EEF9F1,stroke:#36A26B,color:#113D26;
    classDef loc fill:#F6F0FF,stroke:#8B65D9,color:#2C1E4A;
    classDef output fill:#FFF1F1,stroke:#DC6767,color:#4A1F1F;

    subgraph Input["传感器输入与同步"]
        direction LR
        Lidar["LiDAR<br/>点云"]:::input
        IMU["IMU<br/>角速度/加速度"]:::input
        Sync["预处理<br/>时间同步"]:::input

        Lidar --> Sync
        IMU --> Sync
    end

    subgraph Frontend["LIO 前端：实时状态估计"]
        direction LR
        Init["IMU<br/>静态初始化"]:::front
        Deskew["IMU 预测<br/>点云去畸变"]:::front
        Match["scan-to-map<br/>几何约束"]:::front
        IEKF["迭代 ESKF<br/>状态更新"]:::front
        IVox["iVOX<br/>局部地图"]:::front
        KF["关键帧<br/>选择"]:::front

        Init --> Deskew --> Match --> IEKF --> KF
        IEKF --> IVox
        IVox --> Match
    end

    subgraph Mapping["建图后端：全局一致性与地图生成"]
        direction LR
        LC["回环候选<br/>检测"]:::mapping
        NDTLoop["多分辨率<br/>NDT 验证"]:::mapping
        PGO["位姿图<br/>优化"]:::mapping
        OptKF["优化关键帧<br/>位姿"]:::mapping
        Map3D["全局/分块<br/>3D 地图"]:::output
        Grid["2D<br/>栅格地图"]:::output

        LC --> NDTLoop --> PGO --> OptKF
        OptKF --> Map3D
        OptKF --> Grid
    end

    subgraph Localization["在线定位：地图匹配与高频输出"]
        direction LR
        Tiled["分块地图<br/>动态加载"]:::loc
        LocNDT["NDT_OMP<br/>全局匹配"]:::loc
        LocPGO["定位滑窗<br/>PGO"]:::loc
        HF["高频外推<br/>平滑输出"]:::output

        Tiled --> LocNDT --> LocPGO --> HF
    end

    Sync --> Init
    KF --> LC
    Map3D --> Tiled
    IEKF -.->|短时预测| LocNDT
```

建图阶段的主流程为：

1. 对 LiDAR 和 IMU 数据进行格式统一、时间检查和缓存。
2. 按一帧 LiDAR 扫描的起止时间同步 IMU 序列。
3. 启动阶段统计 IMU 均值，初始化重力方向和陀螺零偏。
4. 正常运行阶段用 IMU 对当前扫描周期做前向传播，补偿点云运动畸变。
5. 将去畸变点云与 iVOX 局部地图匹配，构造点到面几何残差。
6. 使用迭代 ESKF 更新当前位姿，并增量维护局部地图。
7. 当运动超过距离或角度阈值时生成关键帧。
8. 对关键帧序列进行回环检测、NDT 验证和位姿图优化。
9. 根据优化后关键帧位姿输出全局点云地图、分块地图和 2D 栅格地图。

定位阶段的主流程为：

1. 加载建图阶段保存的分块地图索引和起点功能点。
2. 根据外部初值或功能点初始化全局位姿。
3. LIO 前端持续输出短时相对运动估计，作为全局匹配初值。
4. 根据预测位姿动态加载附近地图块，并构造 NDT 目标地图。
5. 将当前扫描与局部全局地图做 NDT 匹配，得到绝对定位观测。
6. 使用定位滑窗 PGO 融合 LidarLoc 绝对约束、LidarOdom 相对约束和 DR/IMU 递推约束。
7. 通过外推和平滑器输出与 IMU 接近同频的定位结果。

## 4. 传感器预处理与时间同步

### 4.1 点云统一表达

系统支持多种 LiDAR 数据源，预处理层需要将不同消息格式统一为带时间戳的点云：

$$
{}^{l}\mathbf{p}_{l_jj} =
\begin{bmatrix}
x_j & y_j & z_j
\end{bmatrix}^{\mathsf{T}},
\quad
\tau_j \in [0,T_s],
$$

其中 $\tau_j$ 为点在当前扫描帧内的相对采样时间，$T_s$ 为扫描周期。逐点时间是后续去畸变的关键。如果 $\tau_j$ 的单位、起点或范围错误，点云补偿会沿错误轨迹执行，点到面残差将出现系统性偏差。

预处理还包含以下质量控制：

| 步骤 | 目的 |
| --- | --- |
| 盲区过滤 | 去除距离过近、几何不稳定或传感器盲区内的点 |
| ROI 高度过滤 | 去除明显无关的高空点或地面以下异常点 |
| 点云抽样 | 控制当前扫描规模，降低 IEKF 观测构造成本 |
| 强度/标签过滤 | 对部分雷达剔除低质量回波、多路径或异常回波 |
| 时间戳检查 | 检测 LiDAR/IMU 回退和断流，避免缓存污染 |

### 4.2 LiDAR-IMU 时间同步

对第 $k$ 帧点云，设扫描起始时刻为 $t_k^b$，结束时刻为 $t_k^e$：

$$
t_k^e = t_k^b + \max_j \tau_j.
$$

系统从 IMU 缓存中取出满足

$$
t_k^b \leq t_m \leq t_k^e
$$

的 IMU 序列，与当前点云组成一个测量组。只有当最新 IMU 时间覆盖到 $t_k^e$ 之后，才允许进入本帧处理；否则等待后续 IMU 数据。这一策略保证了去畸变和滤波预测拥有完整的扫描周期运动信息。

```mermaid
sequenceDiagram
    participant L as LiDAR 缓存
    participant I as IMU 缓存
    participant S as 同步器
    participant F as LIO 前端

    L->>S: 取最早点云帧 t_b
    S->>S: 根据逐点时间估计 t_e
    I->>S: 查询 IMU 是否覆盖 t_e
    alt IMU 不足
        S-->>F: 等待更多 IMU
    else IMU 足够
        S->>I: 弹出不晚于 t_e 的 IMU
        S-->>F: 输出 MeasureGroup
    end
```

## 5. LIO 前端里程计

LIO 前端是系统实时性的核心。它将 IMU 预测、点云去畸变、scan-to-map 几何匹配和局部地图更新放在同一个闭环内，形成“预测-配准-更新-建图”的迭代过程。

### 5.1 IMU 静态初始化

启动阶段假设载体近似静止。IMU 测量模型为：

$$
\begin{aligned}
{}^{i}\overline{\boldsymbol{\omega}}
&= {}^{i}\boldsymbol{\omega}
 + {}^{i}\mathbf{b}_g
 + \mathbf{n}_g, \\
{}^{i}\overline{\mathbf{a}}
&= \mathbf{R}_{iw}
\left({}^{w}\mathbf{a}-{}^{w}\mathbf{g}\right)
 + {}^{i}\mathbf{b}_a
 + \mathbf{n}_a .
\end{aligned}
$$

静止时 ${}^{i}\boldsymbol{\omega}\approx\mathbf{0}$，${}^{w}\mathbf{a}\approx\mathbf{0}$。因此可用 IMU 统计均值估计陀螺零偏和重力方向：

$$
{}^{i}\hat{\mathbf{b}}_g
=
\frac{1}{M}\sum_{m=1}^{M}{}^{i}\overline{\boldsymbol{\omega}}_m,
\quad
{}^{w}\hat{\mathbf{g}}
=
-\frac{\bar{\mathbf{a}}}{\|\bar{\mathbf{a}}\|}g.
$$

这里 $\bar{\mathbf{a}}$ 为加速度计均值，$g\approx 9.81\,\mathrm{m/s^2}$。如果加速度均值模长接近 $1$，说明输入可能以 $g$ 为单位，需进行尺度修正；如果模长接近 $9.81$，则通常已经是 $\mathrm{m/s^2}$。

静态初始化的可靠性直接影响后续姿态水平度和点云去畸变质量。实际使用中应尽量保证启动前若干帧 IMU 无明显加减速和大角速度。

### 5.2 IMU 预测模型

设去零偏角速度和加速度为：

$$
\begin{aligned}
{}^{i}\hat{\boldsymbol{\omega}}_m
&= {}^{i}\overline{\boldsymbol{\omega}}_m
 - {}^{i}\hat{\mathbf{b}}_{g,m},\\
{}^{i}\hat{\mathbf{a}}_m
&= {}^{i}\overline{\mathbf{a}}_m.
\end{aligned}
$$

当前默认方案不在线估计 ${}^{i}\mathbf{b}_a$，因此加速度计零偏误差主要被过程噪声、LiDAR 观测和工程参数吸收。完整惯导预测模型可写为：

$$
\begin{aligned}
\hat{\mathbf{R}}_{wi,m+1}^{-}
&=
\hat{\mathbf{R}}_{wi,m}^{+}
\operatorname{Exp}\!\left({}^{i}\hat{\boldsymbol{\omega}}_m\Delta t\right),\\
{}^{w}\hat{\mathbf{v}}_{m+1}^{-}
&=
{}^{w}\hat{\mathbf{v}}_m^{+}
+
\left(
\hat{\mathbf{R}}_{wi,m}^{+}{}^{i}\hat{\mathbf{a}}_m
+{}^{w}\hat{\mathbf{g}}
\right)\Delta t,\\
{}^{w}\hat{\mathbf{p}}_{m+1}^{-}
&=
{}^{w}\hat{\mathbf{p}}_m^{+}
+
{}^{w}\hat{\mathbf{v}}_m^{+}\Delta t
+
\frac{1}{2}
\left(
\hat{\mathbf{R}}_{wi,m}^{+}{}^{i}\hat{\mathbf{a}}_m
+{}^{w}\hat{\mathbf{g}}
\right)\Delta t^2,\\
{}^{i}\hat{\mathbf{b}}_{g,m+1}^{-}
&=
{}^{i}\hat{\mathbf{b}}_{g,m}^{+}.
\end{aligned}
$$

在工程实现中，为抑制加速度零偏和时间同步误差导致的速度漂移，可对速度积分和协方差膨胀采用更保守的策略。无论采用完整积分还是保守外推，核心原则都是：IMU 给出扫描周期内的连续运动先验，LiDAR 观测负责校正长期漂移。

误差协方差预测采用一阶离散化：

$$
\tilde{\mathbf{P}}_{m+1}^{-}
=
\boldsymbol{\Phi}_m
\tilde{\mathbf{P}}_m^{+}
\boldsymbol{\Phi}_m^{\mathsf{T}}
+
\mathbf{G}_m\mathbf{Q}_m\mathbf{G}_m^{\mathsf{T}},
$$

其中 $\boldsymbol{\Phi}_m$ 为误差状态转移矩阵，$\mathbf{G}_m$ 为过程噪声雅可比，$\mathbf{Q}_m$ 为 IMU 噪声协方差。为了避免滤波器过度自信，协方差预测后可进行轻微膨胀，并强制保持对称和正定下限。

### 5.3 点云运动畸变补偿

一帧机械式或固态 LiDAR 点云并非在同一时刻采集。若载体在扫描周期内运动，原始点云会出现扭曲。系统通过 IMU 在扫描周期内积分得到离散位姿序列，并将每个点补偿到扫描结束时刻。

设第 $j$ 个点采样时刻为 $t_j=t_k^b+\tau_j$，扫描结束时刻为 $t_e$。点在采样时刻 LiDAR 坐标系下为 ${}^{l_j}\mathbf{p}_{l_jj}$。补偿链路为：

$$
{}^{i_j}\mathbf{p}_{i_jj}
=
\mathbf{R}_{il}\,{}^{l_j}\mathbf{p}_{l_jj}
+
{}^{i}\mathbf{t}_{il},
$$

$$
{}^{w}\mathbf{p}_{wj}
=
\mathbf{R}_{wi_j}\,{}^{i_j}\mathbf{p}_{i_jj}
+
{}^{w}\mathbf{t}_{wi_j},
$$

$$
{}^{l_e}\mathbf{p}_{l_ej}
=
\mathbf{R}_{li}
\left[
\mathbf{R}_{i_ew}
\left(
{}^{w}\mathbf{p}_{wj}
-
{}^{w}\mathbf{t}_{wi_e}
\right)
-
{}^{i}\mathbf{t}_{il}
\right].
$$

补偿后的点云全部表达在 $\mathcal{F}_{l_e}$ 中，作为当前帧 scan-to-map 匹配输入。

```mermaid
flowchart TD
    A[原始点: 采样时刻 LiDAR 系] --> B[LiDAR -> IMU 外参]
    B --> C[采样时刻 IMU -> 世界]
    C --> D[世界 -> 扫描结束 IMU]
    D --> E[IMU -> 扫描结束 LiDAR]
    E --> F[去畸变点云]
```

### 5.4 Scan-to-Map 观测模型

当前帧去畸变点云降采样后，与 iVOX 局部地图建立对应关系。对第 $j$ 个点，先变换到世界系：

$$
{}^{w}\hat{\mathbf{p}}_{wj}
=
\hat{\mathbf{R}}_{wi}
\left(
\mathbf{R}_{il}\,{}^{l}\mathbf{p}_{lj}
+
{}^{i}\mathbf{t}_{il}
\right)
+
{}^{w}\hat{\mathbf{p}}.
$$

在局部地图中查询近邻点并拟合局部平面：

$$
{}^{w}\mathbf{n}_{j}^{\mathsf{T}}{}^{w}\mathbf{p}
+ d_j = 0,
\quad
\|{}^{w}\mathbf{n}_{j}\|=1.
$$

点到面残差定义为：

$$
r_j
=
{}^{w}\mathbf{n}_{j}^{\mathsf{T}}
{}^{w}\hat{\mathbf{p}}_{wj}
+ d_j.
$$

令

$$
{}^{i}\mathbf{p}_{ij}
=
\mathbf{R}_{il}\,{}^{l}\mathbf{p}_{lj}
+
{}^{i}\mathbf{t}_{il},
$$

对位置和右扰动姿态的一阶雅可比可写为：

$$
\mathbf{J}_j
=
\begin{bmatrix}
{}^{w}\mathbf{n}_{j}^{\mathsf{T}} &
-
{}^{w}\mathbf{n}_{j}^{\mathsf{T}}
\hat{\mathbf{R}}_{wi}
\left[{}^{i}\mathbf{p}_{ij}\right]_{\times}
\end{bmatrix}.
$$

点到面约束利用局部平面结构，通常比点到点 ICP 收敛更快；但在单一平面、走廊、长直道路等场景中存在可观性退化。系统会通过有效点数量、残差分布和信息矩阵特征值判断观测质量。当 $\mathbf{H}^{\mathsf{T}}\mathbf{H}$ 出现低秩或病态时，可对退化方向降低更新权重或膨胀协方差，避免错误约束把状态拉偏。

可选的点到点 ICP 项用于增强局部几何不足时的约束。设最近邻目标点为 ${}^{w}\mathbf{q}_{wj}$，点到点残差为：

$$
\mathbf{r}^{\mathrm{pt}}_j
=
{}^{w}\hat{\mathbf{p}}_{wj}
-
{}^{w}\mathbf{q}_{wj}.
$$

实际更新中，点到面和点到点约束会按权重累积为 6 维位姿信息矩阵和信息向量：

$$
\mathbf{H}^{\mathsf{T}}\mathbf{H}
=
\sum_j \mathbf{J}_j^{\mathsf{T}}\mathbf{J}_j,
\quad
\mathbf{H}^{\mathsf{T}}\mathbf{r}
=
\sum_j \mathbf{J}_j^{\mathsf{T}}r_j.
$$

这种形式避免直接存储大规模观测矩阵，使计算复杂度主要与状态维度相关，而不是与点数线性放大到矩阵求逆维度。

### 5.5 迭代 ESKF 更新

LiDAR 观测更新采用迭代误差状态滤波。每次迭代在当前名义状态 $\hat{\mathbf{x}}^\kappa$ 处重新建立对应关系、拟合平面并线性化残差：

$$
\mathbf{r}^{\kappa}
\approx
\mathbf{r}(\hat{\mathbf{x}}^\kappa)
-
\mathbf{H}^{\kappa}\tilde{\mathbf{x}}^\kappa.
$$

更新可从最大后验角度理解为：

$$
\min_{\tilde{\mathbf{x}}^\kappa}
\left\|
\tilde{\mathbf{x}}^\kappa
-
\tilde{\boldsymbol{\mu}}^{\kappa-}
\right\|_{\left(\tilde{\mathbf{P}}^{\kappa-}\right)^{-1}}^2
+
\left\|
\mathbf{r}^{\kappa}
-
\mathbf{H}^{\kappa}\tilde{\mathbf{x}}^\kappa
\right\|_{\mathbf{N}^{-1}}^2.
$$

由于点云残差数量很大，系统采用信息形式更新。观测信息被累加到位姿 6 维块中，再与先验协方差融合：

$$
\left[
\left(\tilde{\mathbf{P}}^{\kappa-}\right)^{-1}
+
\mathbf{H}^{\mathsf{T}}\mathbf{N}^{-1}\mathbf{H}
\right]
\Delta\tilde{\mathbf{x}}
=
\mathbf{H}^{\mathsf{T}}\mathbf{N}^{-1}\mathbf{r}
-
\left(\tilde{\mathbf{P}}^{\kappa-}\right)^{-1}
\tilde{\boldsymbol{\mu}}^{\kappa-}.
$$

求得增量后通过 $\boxplus$ 注入名义状态。若增量小于阈值或达到最大迭代次数，则结束当前帧更新。为了提高收敛稳定性，可启用 Anderson Acceleration；若加速后残差反而变大，则回退到上一可靠迭代状态。

```mermaid
flowchart TD
    P[IMU 预测先验] --> M[当前状态投影点云到地图]
    M --> N[iVOX 近邻搜索]
    N --> F[局部平面拟合/对应筛选]
    F --> R[构造残差和信息矩阵]
    R --> U[信息形式 IEKF 更新]
    U --> C{收敛?}
    C -- 否 --> M
    C -- 是 --> O[输出后验位姿]
```

## 6. 局部地图与关键帧

### 6.1 iVOX 局部地图

局部地图采用增量式体素索引结构。三维空间被划分为分辨率为 $r_v$ 的体素，每个体素只保存有限数量的代表点。查询某个点的邻域时，可在中心体素及其 6/18/26 邻域内查找候选点，并选择距离最近的若干点用于平面拟合。

iVOX 的作用有三点：

1. **提供低延迟近邻查询**：避免每帧对全局点云做大规模 KD-Tree 重建。
2. **控制局部地图密度**：通过体素中心距离和已有近邻判断是否加入新点。
3. **服务实时观测模型**：为当前扫描的点到面残差提供局部几何面片。

增量更新时，当前帧点云经过后验位姿变换到世界系：

$$
{}^{w}\mathbf{p}_{wj}
=
\hat{\mathbf{R}}_{wi}
\left(
\mathbf{R}_{il}\,{}^{l}\mathbf{p}_{lj}
+
{}^{i}\mathbf{t}_{il}
\right)
+
{}^{w}\hat{\mathbf{p}}.
$$

如果当前点所在体素中已有更靠近体素中心的代表点，则可拒绝加入；如果该区域稀疏或缺少近邻，则直接加入。这种自适应策略使地图在边缘和稀疏区域保留更多信息，在平坦密集区域抑制冗余点。

### 6.2 关键帧策略

关键帧是全局地图、回环检测和 2D 栅格地图的基本单元。当当前位姿相对上一关键帧满足任一条件时创建新关键帧：

$$
\left\|
{}^{w}\mathbf{p}_{k}
-
{}^{w}\mathbf{p}_{k_{\mathrm{last}}}
\right\|
>
d_{\mathrm{kf}},
$$

或

$$
\left\|
\operatorname{Log}
\left(
\mathbf{R}_{wi,k_{\mathrm{last}}}^{\mathsf{T}}
\mathbf{R}_{wi,k}
\right)
\right\|
>
\theta_{\mathrm{kf}}.
$$

关键帧保存当前去畸变点云、LIO 位姿、优化位姿和时间戳。LIO 位姿表示前端原始轨迹，优化位姿表示经过回环或后端修正后的全局一致轨迹。新关键帧的初始优化位姿由上一关键帧优化位姿递推得到：

$$
\mathbf{T}_{w k}^{\mathrm{opt}}
=
\mathbf{T}_{w,k-1}^{\mathrm{opt}}
\left(
\mathbf{T}_{w,k-1}^{\mathrm{lio}}
\right)^{-1}
\mathbf{T}_{w k}^{\mathrm{lio}}.
$$

这样即使前端 LIO 轨迹后续被全局修正，关键帧序列仍能保持连续一致。

## 7. 回环检测与全局优化

前端 LIO 不可避免存在累积漂移。回环检测通过发现当前关键帧与历史关键帧的重复观测，构造长距离约束，并通过位姿图优化把漂移分摊到整条轨迹中。

### 7.1 候选帧筛选

候选帧筛选同时使用时间和空间约束：

| 约束 | 作用 |
| --- | --- |
| 关键帧间隔 | 避免每帧检测，降低计算成本 |
| ID 最小间隔 | 排除时间上相邻的局部帧 |
| 最近历史阈值 | 防止把短期重叠误认为回环 |
| 平面距离阈值 | 只保留优化位姿下空间接近的历史帧 |

设当前关键帧为 $k$，历史关键帧为 $i$。若

$$
|k-i| > N_{\mathrm{id}},
\quad
\left\|
\left(
{}^{w}\mathbf{p}_{k}^{\mathrm{opt}}
-
{}^{w}\mathbf{p}_{i}^{\mathrm{opt}}
\right)_{xy}
\right\|
<
d_{\mathrm{loop}},
$$

则 $i$ 可作为回环候选。

### 7.2 多分辨率 NDT 验证

候选帧通过几何配准验证。系统围绕历史候选关键帧构建局部子图，将若干邻近关键帧点云按优化位姿拼接为目标子地图；当前关键帧或其局部点云作为源点云，使用多分辨率 NDT 从粗到细优化。

NDT 将目标地图按体素划分，并在每个体素中估计均值和协方差：

$$
{}^{w}\boldsymbol{\mu}_v
=
\frac{1}{M_v}\sum_{m=1}^{M_v}{}^{w}\mathbf{q}_{m},
\quad
{}^{w}\boldsymbol{\Sigma}_v
=
\frac{1}{M_v-1}
\sum_{m=1}^{M_v}
\left({}^{w}\mathbf{q}_{m}-{}^{w}\boldsymbol{\mu}_v\right)
\left({}^{w}\mathbf{q}_{m}-{}^{w}\boldsymbol{\mu}_v\right)^{\mathsf{T}}.
$$

源点变换后落入体素 $v(j)$，负对数似然中的主要残差为：

$$
e_j
=
\left(
{}^{w}\hat{\mathbf{p}}_{wj}
-
{}^{w}\boldsymbol{\mu}_{v(j)}
\right)^{\mathsf{T}}
{}^{w}\boldsymbol{\Sigma}_{v(j)}^{-1}
\left(
{}^{w}\hat{\mathbf{p}}_{wj}
-
{}^{w}\boldsymbol{\mu}_{v(j)}
\right).
$$

多分辨率策略先用大体素扩大收敛域，再用小体素细化位姿。若最终 NDT 概率分数超过阈值，则认为回环有效，并生成相对位姿约束。

### 7.3 位姿图优化

位姿图以关键帧位姿为顶点：

$$
\mathcal{V}=\{\mathbf{T}_{wk}\}.
$$

边包括相邻关键帧运动约束、回环约束和可选高度约束：

$$
\min_{\{\mathbf{T}_{wk}\}}
\sum_{(i,j)\in\mathcal{E}_{\mathrm{odom}}}
\left\|
\operatorname{Log}
\left(
\mathbf{Z}_{ij}^{-1}
\mathbf{T}_{wi}^{-1}
\mathbf{T}_{wj}
\right)
\right\|_{\boldsymbol{\Omega}_{ij}}^2
+
\sum_{(i,j)\in\mathcal{E}_{\mathrm{loop}}}
\rho
\left(
\left\|
\operatorname{Log}
\left(
\mathbf{Z}_{ij}^{-1}
\mathbf{T}_{wi}^{-1}
\mathbf{T}_{wj}
\right)
\right\|_{\boldsymbol{\Omega}_{ij}}^2
\right).
$$

其中 $\mathbf{Z}_{ij}$ 是相对位姿观测，$\boldsymbol{\Omega}_{ij}$ 是信息矩阵，$\rho(\cdot)$ 是鲁棒核函数。运动约束通常权重更高，回环约束使用 Cauchy 等鲁棒核降低误匹配影响。

```mermaid
flowchart LR
    K[新关键帧] --> C[时空候选筛选]
    C --> S[构建历史子地图]
    S --> N[多分辨率 NDT]
    N --> V{分数通过?}
    V -- 否 --> X[丢弃候选]
    V -- 是 --> E[生成回环边]
    E --> G[增量位姿图优化]
    G --> M[更新关键帧优化位姿]
    M --> R[触发地图重绘]
```

当启用高度约束时，可加入

$$
r_{h,k} = z_k - z_0
$$

作为先验约束，用于室外平面或单层场景抑制 Z 轴漂移。但在多层室内、坡道、立体结构场景中，高度约束可能压制真实三维运动，应谨慎开启。

## 8. 地图构建与地图表达

### 8.1 全局三维地图

建图结束时，将所有关键帧点云按位姿拼接为全局地图：

$$
\mathcal{M}_{3D}
=
\bigcup_k
\left\{
\mathbf{T}_{wk}^{\mathrm{map}}
{}^{k}\mathbf{p}_{j}
\right\},
$$

其中 $\mathbf{T}_{wk}^{\mathrm{map}}$ 根据是否启用回环，选择优化位姿或 LIO 位姿。全局地图通常进行体素滤波，平衡精度、体积和加载速度。

系统输出两类三维地图：

| 地图 | 用途 |
| --- | --- |
| 全局 PCD | 便于整体检查、可视化和离线评估 |
| 分块点云地图 | 供定位阶段按需加载，支持大场景和动态层 |

### 8.2 分块地图

分块地图将全局点云按二维平面网格切分。设分块边长为 $s_c$，点 ${}^{w}\mathbf{p}$ 对应的块索引为：

$$
\mathbf{g}
=
\left\lfloor
\frac{
\begin{bmatrix}
{}^{w}p_x & {}^{w}p_y
\end{bmatrix}^{\mathsf{T}}
}{s_c}
+
\begin{bmatrix}
0.5 \\ 0.5
\end{bmatrix}
\right\rfloor .
$$

每个地图块保存为独立 PCD，并在索引文件中记录：

1. 地图原点。
2. 块 ID 与二维网格坐标。
3. 每个块对应的点云文件。
4. 功能点，例如建图起点、恢复点或人工标注定位点。

定位时，根据当前预测位姿只加载附近若干块：

$$
\|\mathbf{g}-\mathbf{g}_{\mathrm{cur}}\|_1
\leq
N_{\mathrm{load}}.
$$

远离当前位置的地图块会被卸载。这样可使定位计算量与局部环境大小相关，而不是与全局地图大小相关。

### 8.3 动静态图层

定位阶段除了静态地图，还可维护动态点云层。静态地图来自建图结果，动态图层来自定位过程中可靠匹配后的在线扫描。动态图层的作用是适应临时障碍、场景布置变化和局部结构更新。

动态层可采用三种策略：

| 策略 | 行为 | 适用场景 |
| --- | --- | --- |
| 短期 | 离开区域后清空 | 临时障碍频繁变化 |
| 长期 | 留在内存但不落盘 | 单次任务内的环境变化 |
| 持久 | 保存到磁盘，下次启动加载 | 稳定变化的长期场景 |

更新动态图层需要满足匹配成功、定位置信度足够高、与上次更新距离或时间超过阈值等条件。为避免把车辆自身或低矮地面噪声写入地图，更新前会进行高度过滤。

### 8.4 3D 到 2D 栅格地图

可选的 g2p5 模块将三维关键帧点云投影为二维占据栅格。该模块假设 LiDAR 近似水平安装，或可通过地面估计获得地面平面：

$$
\pi_f:\quad
\mathbf{n}_f^{\mathsf{T}}\mathbf{p}+d_f=0.
$$

对每个 LiDAR 点计算其到地面的高度：

$$
h_j =
\mathbf{n}_f^{\mathsf{T}}{}^{l}\mathbf{p}_{j}+d_f.
$$

若

$$
h_{\min} < h_j < h_{\max},
$$

则该点被视为可投影障碍物，对应世界坐标下的栅格标记为占用。与此同时，从 LiDAR 原点到障碍点方向上的可见区域通过射线更新为自由空间。

```mermaid
flowchart TD
    A[关键帧点云] --> B[地面平面/默认地面高度]
    B --> C[计算点到地面高度]
    C --> D{高度在障碍区间?}
    D -- 是 --> O[写入占用栅格]
    D -- 否 --> G[地面/无效点处理]
    O --> R[射线释放自由空间]
    G --> R
    R --> M[ROS 兼容 OccupancyGrid]
```

输出栅格可转换为 ROS 导航兼容格式：

| 栅格值 | 含义 | 保存图像 |
| --- | --- | --- |
| $0$ | 自由空间 | 白色 |
| $100$ | 占用空间 | 黑色 |
| $-1$ | 未知空间 | 灰色 |

2D 栅格主要服务导航和显示，不参与默认三维定位。回环发生后，由于关键帧优化位姿改变，需要触发全局重绘，以保持 2D 地图与优化后的三维轨迹一致。

## 9. 在线定位方案

定位阶段的输入是在线 LiDAR/IMU 数据和建图阶段保存的分块地图。定位输出是世界系下的高频位姿 $\mathbf{T}_{wb}$ 或 $\mathbf{T}_{wl}$。

### 9.1 定位总体流程

```mermaid
flowchart TD
    A[加载地图索引和功能点] --> B[设置外部初值/功能点初始化]
    B --> C[LIO 输出短时相对运动]
    C --> D[按预测位姿加载地图块]
    D --> E[构建 NDT 目标地图]
    E --> F[当前扫描 NDT 匹配]
    F --> G[得到 LidarLoc 绝对观测]
    C --> H[PGO 相对约束]
    G --> I[定位滑窗 PGO]
    H --> I
    I --> J[DR/IMU 外推]
    J --> K[平滑与 TF 输出]
```

定位由两类信息共同决定：

1. **相对运动**：LIO/IMU 给出短时间连续运动，频率高、局部平滑，但会漂移。
2. **绝对匹配**：当前点云与全局地图 NDT 匹配，频率较低、依赖初值和地图质量，但可消除漂移。

系统不直接用全局匹配替代里程计输出，而是通过滑窗 PGO 和高频外推融合二者，使输出既平滑又不长期漂移。

### 9.2 初始位姿

定位需要一个足够接近真实位置的初值。系统支持两类初值：

1. **外部初值**：由人工、上位机、GNSS 或其他系统给定 $\mathbf{T}_{w0}$。
2. **功能点初值**：地图索引中保存的功能点，如建图起点或恢复点。

给定初值后，系统加载该位置附近地图块，并执行 NDT 匹配。如果置信度超过初始化阈值，则定位进入正常跟踪状态。若初值存在较大 yaw 不确定性，可采用 yaw 网格搜索：固定位置、roll、pitch，在一定角度范围内采样多个 yaw 初值，选择 NDT 分数最高者，再进入精配准。

### 9.3 NDT 全局匹配

定位匹配使用 NDT_OMP。地图块被组合为当前目标点云，并预计算体素高斯分布。当前扫描以 LIO 递推位姿为初值：

$$
\mathbf{T}_{wk}^{\mathrm{guess}}
=
\mathbf{T}_{w,k-1}^{\mathrm{abs}}
\left(
\mathbf{T}_{o,k-1}^{\mathrm{lio}}
\right)^{-1}
\mathbf{T}_{o,k}^{\mathrm{lio}},
$$

其中 $\mathbf{T}_{w,k-1}^{\mathrm{abs}}$ 是上一帧绝对定位结果，$\mathbf{T}_{o,k}^{\mathrm{lio}}$ 是 LIO 在自身里程计坐标系中的位姿。NDT 在该初值附近求解：

$$
\mathbf{T}_{wk}^{\mathrm{ndt}}
=
\arg\min_{\mathbf{T}}
\sum_j
\left(
{}^{w}\mathbf{p}_{wj}(\mathbf{T})
-
{}^{w}\boldsymbol{\mu}_{v(j)}
\right)^{\mathsf{T}}
{}^{w}\boldsymbol{\Sigma}_{v(j)}^{-1}
\left(
{}^{w}\mathbf{p}_{wj}(\mathbf{T})
-
{}^{w}\boldsymbol{\mu}_{v(j)}
\right).
$$

为了避免全局匹配在退化场景中产生跳变，可对 NDT 相对初值的修正量做比例融合：

$$
\mathbf{T}_{wk}^{\mathrm{loc}}
=
\mathbf{T}_{wk}^{\mathrm{guess}}
\operatorname{Exp}
\left(
\alpha
\operatorname{Log}
\left[
\left(\mathbf{T}_{wk}^{\mathrm{guess}}\right)^{-1}
\mathbf{T}_{wk}^{\mathrm{ndt}}
\right]
\right),
\quad
0<\alpha\leq1.
$$

当启用强制 2D 定位时，输出会将 roll、pitch 和 z 分量投影到约束平面：

$$
\mathrm{roll}=0,\quad
\mathrm{pitch}=0,\quad
z=0.
$$

这适合平面车辆定位，但不适合多层、坡道或明显三维运动场景。

### 9.4 定位滑窗 PGO

定位 PGO 将低频绝对匹配和高频相对运动融合。每个定位帧为一个 SE(3) 顶点，约束包括：

| 约束 | 类型 | 作用 |
| --- | --- | --- |
| LidarLoc | 绝对先验 | 将当前帧拉回全局地图 |
| LidarOdom | 相对边 | 维持短时连续性和平滑性 |
| DR/IMU | 相对边 | 在 LiDAR 匹配延迟或失败时外推 |
| Prior | 边缘化先验 | 滑窗移除旧帧后保留历史信息 |

优化目标为：

$$
\min_{\{\mathbf{T}_{wk}\}}
\sum_k
\left\|
\operatorname{Log}
\left[
\left(\mathbf{Z}_{wk}^{\mathrm{loc}}\right)^{-1}
\mathbf{T}_{wk}
\right]
\right\|_{\boldsymbol{\Omega}_{\mathrm{loc}}}^2
+
\sum_{(i,j)}
\left\|
\operatorname{Log}
\left[
\left(\mathbf{Z}_{ij}^{\mathrm{rel}}\right)^{-1}
\mathbf{T}_{wi}^{-1}\mathbf{T}_{wj}
\right]
\right\|_{\boldsymbol{\Omega}_{\mathrm{rel}}}^2.
$$

LidarLoc 分数较低时，绝对约束权重降低；LidarOdom 被判定异常时，其相对约束可被短时间降权或跳过。滑窗中早期帧收敛后可边缘化为先验约束，以限制计算规模。

### 9.5 高频输出与平滑

全局匹配通常低于 IMU 频率，且 NDT 匹配存在计算延迟。系统使用最新 PGO 结果作为低频基准，再用 DR/IMU 或 LIO 队列外推到最新时刻：

$$
\mathbf{T}_{w,t}^{\mathrm{out}}
=
\mathbf{T}_{w,k}^{\mathrm{pgo}}
\left(
\mathbf{T}_{r,k}^{-1}
\mathbf{T}_{r,t}
\right),
\quad
t\geq t_k.
$$

其中 $\mathbf{T}_{r}$ 是相对运动源的位姿。最后通过平滑器抑制小幅抖动，输出稳定的 TF 或定位消息。若车辆处于静止状态，可直接保持上一定位结果，避免静止时匹配噪声导致位姿跳动。

## 10. 关键质量控制

### 10.1 前端退化检测

点云几何可能无法完整约束 6 自由度。例如：

| 场景 | 弱约束方向 |
| --- | --- |
| 单一地面 | 平面内平移、绕法向旋转 |
| 长走廊 | 沿走廊方向平移 |
| 稀疏开阔区域 | 多个方向均弱 |
| 动态物体占比高 | 对应关系不稳定 |

系统通过 $\mathbf{H}^{\mathsf{T}}\mathbf{H}$ 的特征值判断观测退化。若最小特征值相对最大特征值过小，说明存在弱观测方向。处理方式包括：

1. 对退化方向降低更新强度。
2. 膨胀位姿协方差，避免过度自信。
3. 引入 IMU、轮速、GNSS 或高度先验辅助约束。
4. 延迟关键帧或地图更新，等待更好的几何视角。

### 10.2 对应关系筛选

错误对应通常比优化器选择更容易导致失败。有效对应需要满足：

1. 近邻数量足够。
2. 局部平面拟合残差低于阈值。
3. 点到平面距离不过大。
4. 激光束方向与平面关系不过于退化。
5. 当前帧有效匹配数量超过下限。

若有效点数过少，应跳过 LiDAR 更新或降低观测权重，避免少量错误面片主导滤波结果。

### 10.3 回环可靠性

回环误匹配会破坏全局地图。可靠回环应满足：

1. 候选帧与当前帧有足够时间间隔。
2. 优化位姿下空间接近。
3. 多分辨率 NDT 能稳定收敛。
4. 配准分数超过阈值。
5. 回环边经过鲁棒核和卡方误差检测。

对多层建筑、重复走廊、货架区等强别名场景，应提高回环阈值，或加入高度、语义、强度、GNSS 等额外判别信息。

### 10.4 定位健康度

定位阶段需要同时监测绝对匹配和相对运动：

| 指标 | 异常含义 |
| --- | --- |
| NDT 置信度连续偏低 | 当前扫描与地图不一致或初值偏差过大 |
| LidarLoc 与 LidarOdom 相对运动差异过大 | LIO 漂移、地图匹配跳变或动态物体干扰 |
| IMU/DR 断流 | 高频外推不可用 |
| 地图块加载为空 | 初值偏离地图或分块索引异常 |
| 静止状态位姿抖动 | 匹配噪声主导输出，需要静止保持 |

若绝对匹配连续失败，可进入跟随 DR/LIO 状态；若失败持续时间过长，则应触发重新初始化或请求外部初值。

## 11. 参数设计建议

### 11.1 前端参数

| 参数类别 | 建议 |
| --- | --- |
| 点云降采样 | 应保证有效点数充足；过大分辨率会导致平面约束不足 |
| iVOX 分辨率 | 与场景结构尺度匹配；室内可小，室外大场景可适当增大 |
| 平面拟合阈值 | 过小会丢弃有效面，过大会引入曲面和动态噪声 |
| ESKF 最大迭代次数 | 通常 3-5 次即可，过多迭代收益有限且增加延迟 |
| IMU 噪声 | 过小会过信 IMU，过大会使预测不稳定；应结合设备标定 |
| 关键帧阈值 | 距离过小地图冗余，过大回环和栅格地图稀疏 |

### 11.2 回环参数

| 参数类别 | 建议 |
| --- | --- |
| 回环检测间隔 | 根据关键帧密度设置，避免高频无效检测 |
| 空间搜索半径 | 过小漏检回环，过大增加误匹配 |
| NDT 分数阈值 | 重复结构场景应更严格 |
| 高度约束 | 单层平面场景可开启，多层三维场景应关闭 |
| 鲁棒核阈值 | 应允许小漂移闭环，但拒绝明显错误边 |

### 11.3 定位参数

| 参数类别 | 建议 |
| --- | --- |
| 地图加载范围 | 应覆盖一次匹配可能收敛的空间范围 |
| NDT 分辨率 | 粗分辨率扩大收敛域，细分辨率提高精度 |
| 初始化阈值 | 应高于正常跟踪阈值，避免错误初始化 |
| 强制 2D | 仅用于平面车辆或固定高度场景 |
| 动态层策略 | 长期变化选 persistent，短时障碍选 short/long |
| 平滑因子 | 过大跟随快但抖动，过小平滑但延迟 |

## 12. 方案特点与适用边界

### 12.1 技术特点

1. **紧耦合 LIO 前端**：IMU 预测和 LiDAR 几何观测在同一 ESKF 框架内融合。
2. **直接点云配准**：不依赖手工边缘/平面特征提取，适配 Livox、Velodyne、Ouster、RoboSense 等多种雷达。
3. **高效局部地图**：iVOX 支持快速邻域查询和增量更新，适合实时 scan-to-map。
4. **回环全局一致**：NDT 回环验证与位姿图优化可显著降低长距离漂移。
5. **地图工程化表达**：同时支持全局 PCD、分块地图、功能点和 2D 栅格地图。
6. **定位高频输出**：全局匹配低频修正，相对运动高频外推，兼顾精度和实时性。
7. **动态场景适应**：通过动态图层在线更新缓解地图长期变化问题。

### 12.2 适用场景

该方案适用于室外园区、道路、厂区、仓储、机器人巡检、室内外混合通道等具备稳定几何结构的场景。若环境中存在足够墙面、地面、立柱、建筑轮廓或其他几何结构，点到面 LIO 和 NDT 定位通常能获得较高稳定性。

### 12.3 局限性

1. **强动态环境**：大量移动物体会污染局部地图和全局匹配。
2. **几何退化场景**：长直隧道、开阔平面、玻璃墙面等会削弱点云约束。
3. **初值依赖**：NDT/ICP 属于局部优化方法，定位初始化偏差过大时可能收敛到错误位置。
4. **时间同步敏感**：逐点时间和 IMU-LiDAR 时间覆盖错误会直接影响去畸变。
5. **多层场景回环困难**：仅依赖平面距离筛选候选时，多层结构可能产生错误候选或漏检。
6. **外参依赖**：默认不在线估计外参，外参误差会体现为系统性配准残差。

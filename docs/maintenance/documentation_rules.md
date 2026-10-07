@page documentation_rules 基于代码 diff 的文档更新规则

# 目标

文档审查不以“改了多少行”为准，而以 diff 是否改变读者依赖的契约为准。代码作者在提交前查看：

```bash
git diff -- CMakeLists.txt cmake package.xml src scripts config
git diff --name-only -- CMakeLists.txt cmake package.xml src scripts config
```

# 路径到文档的映射

| diff 命中 | 必查文档 | 必查内容 |
|---|---|---|
| 顶层/`src/**/CMakeLists.txt`, `package.xml`, `cmake/**` | @ref build_and_run "构建页" | 目标、依赖、安装位置、构建命令、开关 |
| `src/app/**` | @ref index "首页"、@ref architecture "架构页"、对应 flow | 程序入口、flags、返回码、顶层所有权和调用顺序 |
| `scripts/**/*.sh` | @ref script_contracts "脚本契约"、对应 flow/build | 参数、环境、cwd、子调用、副作用、输出/退出契约 |
| `src/core/lio/**` | @ref laser_mapping_module "前端模块页"、@ref offline_slam_flow "离线流程" | 数据同步、状态更新、关键帧、锁和不变量 |
| `src/core/backend/**`, `loop_closing/**` | @ref backend_module "后端模块页"、@ref offline_slam_flow "离线流程" | 后端 mode、线程、优化/提交顺序、输出 |
| `src/core/localization/**`, `system/loc_system.*`, `app/run_loc_online.cc` | @ref online_localization_flow "在线定位流程"、@ref localization_module "定位模块页"、@ref online_systems "在线编排页" | 输入顺序、重定位策略、map/odom、发布门限、线程与退出 |
| `src/core/system/**` | @ref online_systems "在线编排页"、受影响 flow | 所有权、消息队列、ROS 回调、启动/停止顺序 |
| `src/core/maps/**`, `g2p5/**` | @ref map_and_support_modules "地图与支撑模块页"、导出 flow | 坐标系、格式、原点、覆盖/原子性 |
| `src/wrapper/**`, `src/io/**`, message dependencies | @ref map_and_support_modules "地图与支撑模块页"、所有受影响 flow | topic type、回调线程、序列化/时间语义 |
| `src/ui/**`, `src/utils/**`, `src/core/miao/**` | @ref map_and_support_modules "地图与支撑模块页" | UI/文件副作用、诊断边界、优化器所有权与求解约束 |
| `config/**` | 启动页和受影响 module/flow | 默认值、单位、启用条件、互斥/前置关系 |
| 公共头文件 API | 对应模块页 + Doxygen | 签名、所有权、线程安全、失败语义 |

# 语义触发器

即使文件路径没有新增模块，只要 diff 出现下列变化，也必须更新文档：

- 新增/删除/重命名 `add_executable`、`add_library`、install target 或 ROS dependency。
- 新增 flag、YAML key、环境变量、topic/service、输出文件或退出码。
- `unique_ptr/shared_ptr`、栈对象或回调捕获改变，导致所有权/生命周期变化。
- `std::thread/async/OpenMP`、mutex、condition variable、atomic 或队列策略变化。
- 输入排序、时间同步、坐标变换、frame ID、单位、阈值或容错变化。
- pipeline 中新增 early return、drop、fallback、重试、覆盖或删除行为。

# 结论证据规则

1. 每个关键结论至少关联“真实文件 + 符号或稳定行段”；优先写 Doxygen `@ref`，行号只用于帮助首次定位。
2. API 语义必须从声明与调用者交叉确认；仅凭名称或注释不能写成事实。
3. 性能、线程安全、执行顺序若未由实现保证，前缀写 **推断：** 并说明依据。
4. 不复制大段源码；写不变量、状态转换和失败边界，并链接 API/源码。
5. 行号可能漂移，评审 diff 时应验证符号链接仍成立；重命名后运行 docs target 让断链变成构建失败。
6. `index.md` → `guide/` 章节页 → 专题页使用 `@subpage` 建立唯一父节点的导航树；跨章节与专题间跳转使用 `@ref`。所有页面必须能从首页到达。双向 `@subpage` 会形成页面父子环，并可能让 Doxygen 1.9.1 长时间停在符号解析阶段；生成前先运行源文件检查。

# PR/提交检查

```text
[ ] 我检查了目标/入口/脚本/API/配置的 diff 分类。
[ ] 数据流、生命周期、线程或关键约束变化已同步到对应页。
[ ] shell 的参数、环境、cwd、调用链和副作用已同步。
[ ] 新结论有 @ref 或文件+符号证据；推断已标记。
[ ] 若涉及业务代码或依赖变化，普通构建和链接通过。
[ ] `BUILD_DOCS=ON` 的 `docs` target 通过，且 `doxygen-warnings.log` 中没有 `docs` 告警。
```

# 分层维护与校验

流程放在 `flows/`，模块职责放在 `modules/`，公式和参数化放在 `algorithms/`；`guide/` 将相关专题组织到同一个学习章节。ESKF 状态或噪声变化同时检查算法页与章节中的代码对照表；残差变更核对正负号、扰动方向、单位与信息权重。首页保持章节入口和学习路线，细节留在专题。

```bash
python3 docs/maintenance/check_links.py
cmake -S docs -B build_docs
cmake --build build_docs --target docs
python3 docs/maintenance/check_links.py --html build_docs/docs/html
```

文档构建会拒绝手写页面的 Doxygen 告警；历史源码注释告警仍留在日志，不将其误报为本轮零告警。链接检查覆盖全部手册页、API 目标、锚点和图像资源；公式与图布局还需浏览器查看。

Shell 源码通过 `shell_sources.dox` 的 `@example{lineno}` 直接引用真实文件，手册用 `@ref` 跳转。新增、删除或重命名脚本时同步注册与调用说明；检查器会拒绝注册遗漏、重复或不存在的脚本。生成的脚本源码页也参加 HTML 链接检查。

源码注释变化检查非注释 token 不变。改动文档 CMake 时核验独立 docs 和 `BUILD_DOCS=ON` 两个入口；涉及业务代码或依赖变化时，再使用同一依赖环境验证业务构建。文档迁移之后检查整个仓库与工作区入口的路径，不只检查 Markdown。

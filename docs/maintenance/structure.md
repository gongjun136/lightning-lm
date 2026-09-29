@page documentation_structure 文档结构与保留标准

# 目录职责

| 路径 | 放什么 | 不放什么 |
|---|---|---|
| `docs/index.md`, `architecture.md` | 主入口与总架构 | 长推导和逐次实验汇报 |
| `getting_started/` | 构建、工作区、脚本契约 | 与运行无关的个人笔记 |
| `flows/` | 从入口到结果的完整路径 | 复制各模块实现细节 |
| `modules/` | 职责、生命周期、并发、失败边界 | 所有函数逐条复述 |
| `algorithms/` | 核心公式、参数化、实现差异、源码映射 | 无实现关联的数学教材 |
| `reference/` | 配置与可重现评测协议 | 无版本上下文的性能承诺 |
| `operations/` | 当前诊断、对齐、回归与必要故障证据 | 过时部署补丁指令 |
| `maintenance/` | 更新规则、结构约束、校验方式 | 临时工作记录 |
| `assets/` | 当前页面/README 实际引用的资源 | 未引用图片、整批实验产物 |

只使用一个 `docs/`。生成站点不提交到源目录，Doxygen 配置在 `docs/Doxyfile.in`。这里的层次代表阅读深度，不用 `internal/` 暗示只有内部用户才能访问算法页。

# 设计参考与取舍

[Ceres 文档入口](https://github.com/ceres-solver/ceres-solver/blob/master/docs/source/index.rst) 将安装、教程、导数、建模与求解原理分开；本手册借鉴阅读层次，不引入它的 Sphinx 工具链。
[GTSAM API](https://gtsam.org/doxygen/) 按逻辑模块导航并保留类/文件索引；本项目继续使用 Doxygen + Graphviz + Markdown，将主题页和本仓库 API 连通。

# 本轮清理

保留既有流程与模块页，并纠正在线 SLAM 队列等与源码不符的描述。将有用的几何、ESKF、统计推导吸收到直接关联实现的算法页；移除 `gj/` 的独立通用笔记、重复 DOCX/XLSX、旧架构图与汇总文档。
移除 2026 年 7 月的一次性评测汇报、重复报告、manifest 和专属图表，保留评测方法与运行脚本。它们记录特定旧版本，不作为当前代码说明。保留 9 月定位故障证据及配套 JSON，因为它们解释当前失效门控与回归边界；不把历史现场结果描述成当前版本已验证性能。

README 仍引用的演示图归入 `assets/examples/`。保留材料必须能回答当前维护问题；不能仅因“曾经生成过”就持续放入仓库。

# 导航约束

只有 `index.md` 使用 `@subpage`，每页唯一一次。模块与原理、原理与 API 之间使用 `@ref`，防止父子页环。新增主题先确定稳定 page id，再补首页、调用方链接和源码反向链接。

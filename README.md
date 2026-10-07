# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

**最新源码与交付入口（2026-10-07）**：本仓库是版本来源；最小构建包由
[`tools/package_minimal.py`](tools/package_minimal.py) 从指定 Git 提交生成，使用
`WCNS_BUILD_TOOLS=OFF` 只编译求解器，不附带网格和计算结果。
见 [最小包构建说明](docs/minimal-build-readme.md)。已有案例的新增大网格以无损 XZ 归档随仓库保存，
新克隆后先运行 `python tools/restore_case_grids.py` 恢复，再使用案例配置。
[CBC 九组合脚本与操作步骤](campaigns/cbc_9cases/README.md) 继续保留，已有上传包可继续使用。
本地计算结果和临时验证数据保留在原目录，新生成的结果与构建产物由 `.gitignore` 排除。

当前内部版本为 **WCNS v2.6**。新增均匀各向同性自由衰减和强迫湍流，包含 MPI 分布 FFT/FFTW 初始化、低波数强迫、53项统计、能谱、相关函数及跨分区续算。详见 [HIT 详细报告](docs/hit-v2.6/HIT算例与v2.6实现报告.md)、[CBC 算例](cases/manual/case11_hit_decay/README.md) 和 [JHTDB 算例](cases/manual/case12_hit_forced/README.md)。配套16³至128³八块网格与原始参考数据已准备，本机仅进行配置、解析与三步测试，未做长计算。

2026-10-07 初始化更新：新增 CBC **预演化后逐模态匹配谱**及所给 C++ / Samtaney IC4 的**解析谱随机相位初始化**，通过 `hit.initialization` 切换。准备阶段独立输出和重启，正式时钟与统计清零；旧配置行为保留。见 [详细说明与参数](docs/hit-v2.6/HIT初始化与预演化说明.md)、[CBC 准备配置](cases/manual/case11_hit_decay/prepared_production.wcns)、[解析谱配置与网格](cases/manual/case13_hit_analytic/README.md)、[本次验证](docs/hit-v2.6/初始化更新验证.md)。本次仅做初始化检查和短步测试，未进行实际长计算。

v2.5 历史说明：新增 SD7003 和压缩拐角多块结构网格 ILES，见 [第五章调研报告](docs/chapter5-v2.5-research.md)、[v2.5 使用说明](docs/chapter5-v2.5-implementation.md) 和 [验证记录](docs/v2.5-validation.md)。本机只做短步测试；物理统计尚未收敛验证。

v2.4 历史说明：在 v2.3 基础上新增 Re_h=1400 周期山 ILES 初场、动态流量驱动、时间／展向统计和重分区续算；配套 101376 单元网格及 smoke/production 配置已准备。设置依据和补充假设见 [周期山报告](docs/periodic-hill-v2.4.md)，测试见 [v2.4 验收](docs/v2.4-validation.md)，变更见 [发行说明](docs/release-notes-2.4.md)。本机仅进行短步验证，未启动统计长算。v2.4 源码快照使用 `tools/package_v2_4.py` 打包，按 SHA-256 清单核对实际内容。以下 v2.3 段落保留为历史背景。

`v2.3` 仍属于物理验证受限的内部版本：没有在本机运行 Case06 大网格或目标三维翼型；本次
优化只完成公式等价、单元、小网格 SSPRK3、串行/MPI 回归和短时性能卡口，不宣称长期低马赫湍流
或复杂激波精度已经通过验证。全速度 Roe 见 [`docs/all-speed-roe.md`](docs/all-speed-roe.md)，
求解器覆盖审查见 [`docs/riemann-solvers-v2.2.md`](docs/riemann-solvers-v2.2.md)，v2.3 优化与
自动验收见 [`docs/performance-optimization-v2.3.md`](docs/performance-optimization-v2.3.md) 和
[`docs/v2.3-validation.md`](docs/v2.3-validation.md)。v2.0_pre 的范围和卡口见
[`docs/v2.0-pre-development-plan.md`](docs/v2.0-pre-development-plan.md)，逐项支持见
[`docs/v2.0-pre-capability-matrix.md`](docs/v2.0-pre-capability-matrix.md)，数学定义见
[`算法补充.md`](算法补充.md)，使用方法和限制见
[`docs/user-manual.md`](docs/user-manual.md) 与
[`docs/known-limitations.md`](docs/known-limitations.md)。代码复审和数值一致性证据见
[`docs/v2.0-pre-code-review.md`](docs/v2.0-pre-code-review.md)，Linux 服务器迁移与运行步骤见
[`docs/linux-server-guide.md`](docs/linux-server-guide.md)。

v1.1.0 的详细范围、P--U 阶段、自动卡口、人工判断及
Git 闭环见 [`docs/v1.1.0-development-plan.md`](docs/v1.1.0-development-plan.md)；
物理容许性、局部通量降阶、壁面载荷/热流、输运和性能公式见
[`算法补充.md`](算法补充.md) 第 11 节；各阶段实现状态以对应设计和验收报告为准。

V、W、X、AA、Y、Z、AB0 和 AB1 的既有设计与验收记录保留在开发仓库。原 AC--AE 中依赖
服务器和目标算例的物理验证继续属于正式 `v2.0.0` 后续工作，不阻塞 `v2.0_pre` 精简源码
目录；历史计划见 [`docs/v2.0.0-development-plan.md`](docs/v2.0.0-development-plan.md)。

本开发仓库保留阶段设计、自动测试、人工算例及验收证据。v2.3 的确定性 Linux 源码包由
`tools/package_v2_3.py` 从版本提交直接生成，只附带完整配置模板和 t=250 分叉配置模板，不包含网格、结果或
完整算例。历史上的独立 `wcns_v3_release` 精简仓库不再作为版本来源真值。

当前程序具备 CGNS 结构多块网格读取、两套独立高阶几何 profile、确定性运行时剖分与 MPI halo、
六种界面重构、Rusanov/HLL/HLLC/Roe/Li--Gu 与 Rieper all-speed Roe/旋转 Roe、层流 Navier--Stokes、SSPRK3、定常/
非定常 LU-SGS、Roe/Weiss--Smith、SA-neg、SST-2003m、实验级标准 k-epsilon 及阶段 Z 的五种
LES 候选实现。边界输出
覆盖壁面量、压力/黏性分载荷和三维展向分箱；接受步时间统计覆盖 mean/RMS/covariance/Favre，
并随可改变 rank 数的 CGNS checkpoint 连续恢复。v1.1 的 SSPRK 稳健化路径继续保留且默认关闭。
检查点默认使用严格数值签名续算；显式 `restart.mode=algorithm_change` 可在同一网格上导入守恒
状态并切换核心算法，同时重置多步历史、收敛基准和时间统计。
逐步使用说明见 [`docs/user-manual.md`](docs/user-manual.md)，源码二次开发见
[`docs/developer-guide.md`](docs/developer-guide.md)，可复制的完整配置见
[`examples/full_case_template_v2.3.wcns`](examples/full_case_template_v2.3.wcns)；简明运行速查见
[`docs/runtime-guide.md`](docs/runtime-guide.md)。

新增的 `turbulent_channel` 初场、y-z 截面监测和专用槽道壁摩擦/`Re_tau` 统计已用于
[`case05`](cases/manual/case05_3d_turbulent_channel/README.md) 的 4-rank、5 步 Linux 迁移前可行性卡口；
该稀疏网格结果不是湍流统计或 DNS 验收。

四块圆柱 O 网格、Re=20/40/100/200 层流 Navier--Stokes 结果以及 Mach 5 Euler 钝体绕流的
完整配置、复现方法、图像和精度边界见
[`case07`](cases/manual/case07_2d_cylinder/README.md)。

## 构建与测试

```powershell
cmake -S . -B build -G "MinGW Makefiles" -DWCNS_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix build\install
```

启用 MPI（Windows/MinGW 下使用 Intel MPI）：

```powershell
cmake -S . -B build-mpi -G "MinGW Makefiles" -DWCNS_ENABLE_MPI=ON -DWCNS_BUILD_TESTS=ON
cmake --build build-mpi
ctest --test-dir build-mpi --output-on-failure
```

CGNS 4.4.0 源码归档随仓库提供，CMake 会以静态 ADF 后端构建，不需要联网下载或单独安装 HDF5。
安装树的 `bin` 包含正式求解器、开发仓库中的网格生成器、独立验证器和上游 CGNS 工具，`share/wcns` 包含配置模板、
算法/运行文档和第三方通知。MinGW 运行库不会自动复制；具体环境与安装检查见
[`docs/runtime-guide.md`](docs/runtime-guide.md)。

## 设计约定

- 内部索引从零开始。
- 物理单元范围为 `[0, n)`，ghost 索引允许为负数。
- 二维网格仍使用三维索引和五分量 Euler 状态。
- 全局块编号、MPI 所属进程和进程内数组下标相互独立。
- Euler 状态采用五分量守恒量；重构空间和 Riemann 求解器由严格配置选择。
- WCNS 求解需要至少三层 cell-centered ghost，块连接通信守恒量，接收后转换原始量。

# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

当前内部版本为 **WCNS v2.3**。本版本保持 v2.2 的数值功能和配置语义，重点清理求解器
热路径中的重复状态转换、重复参数校验和通用 Roe 特征矩阵乘法；标准 Roe、Li--Gu 与
Rieper all-speed Roe 改为等价的 Roe--Pike 闭式耗散，预处理 Roe 仍保留其专用矩阵。
RANS/LES、Weiss--Smith+LU-SGS、三维载荷和可重启统计继续保留。本项目采用
本机独立开发，不启用外部 CI；质量卡口由本机构建、CTest、算法规格和发布契约检查承担。
对外发布和 WCNS 自有代码许可证决定均暂缓。

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

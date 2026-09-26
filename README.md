# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

当前稳定版本为 **WCNS v1.1.0**；正在整理的 **v2.0_pre** 是功能完整、物理验证受限的内部
前瞻源码版本。本项目采用本机独立开发，不启用外部 CI；质量卡口由本机串行/MPI CMake、
CTest、算法规格和发布契约检查承担。对外发布和 WCNS 自有代码许可证决定均暂缓。

`v2.0_pre` 包含原 v2 路线的 RANS/LES、低 Mach 预处理、定常/非定常 LU-SGS、三维载荷与
可重启统计能力，但不运行 Case06 大网格、长期湍流或尚未提供的目标三维翼型。因此它证明
通用算法路径可构建、可执行并与冻结公式一致，不宣称工程精度已经通过验证。范围和卡口见
[`docs/v2.0-pre-development-plan.md`](docs/v2.0-pre-development-plan.md)，逐项支持见
[`docs/v2.0-pre-capability-matrix.md`](docs/v2.0-pre-capability-matrix.md)，数学定义见
[`算法补充.md`](算法补充.md)，使用方法和限制见
[`docs/user-manual.md`](docs/user-manual.md) 与
[`docs/known-limitations.md`](docs/known-limitations.md)。

v1.1.0 的详细范围、P--U 阶段、自动卡口、人工判断及
Git 闭环见 [`docs/v1.1.0-development-plan.md`](docs/v1.1.0-development-plan.md)；
物理容许性、局部通量降阶、壁面载荷/热流、输运和性能公式见
[`算法补充.md`](算法补充.md) 第 11 节；各阶段实现状态以对应设计和验收报告为准。

V、W、X、AA、Y、Z、AB0 和 AB1 的既有设计与验收记录保留在开发仓库。原 AC--AE 中依赖
服务器和目标算例的物理验证继续属于正式 `v2.0.0` 后续工作，不阻塞 `v2.0_pre` 精简源码
目录；历史计划见 [`docs/v2.0.0-development-plan.md`](docs/v2.0.0-development-plan.md)。

本开发仓库保留阶段设计、自动测试、人工算例及验收证据。v1.1.0 的确定性内部源码包由
`tools/package_release.py` 从版本提交直接生成，并包含版本验收所需的 Case07；NACA0012 不在
本版本范围内。历史上的独立 `wcns_v3_release` 精简仓库不再作为版本来源真值。

当前程序具备 CGNS 结构多块网格读取、两套独立高阶几何 profile、确定性运行时剖分与 MPI halo、
六种界面重构、Rusanov/HLLC/Roe、层流 Navier--Stokes、SSPRK3、定常/非定常 LU-SGS、Roe/
Weiss--Smith、SA-neg、SST-2003m、实验级标准 k-epsilon 及阶段 Z 的五种 LES 候选实现。边界输出
覆盖壁面量、压力/黏性分载荷和三维展向分箱；接受步时间统计覆盖 mean/RMS/covariance/Favre，
并随可改变 rank 数的 CGNS checkpoint 连续恢复。v1.1 的 SSPRK 稳健化路径继续保留且默认关闭。
逐步使用说明见 [`docs/user-manual.md`](docs/user-manual.md)，源码二次开发见
[`docs/developer-guide.md`](docs/developer-guide.md)，可复制的完整配置见
[`examples/full_case_template.wcns`](examples/full_case_template.wcns)；简明运行速查见
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

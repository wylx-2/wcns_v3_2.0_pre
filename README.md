# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

当前生产版本为 **WCNS v1.1.0**。`v1.1.0-rc.1` 已通过人工核验，并在一次不改变程序功能的
源码清理和完整本地回归后进入主干。本项目当前采用本机独立开发，不启用外部 CI；GitHub
Actions 工作流已在发布后移除，质量卡口继续由本机串行/MPI CMake、CTest 和验收脚本承担。
对外发布和代码许可证决定均暂缓。程序在本地发布矩阵基础上，已经过二维
Riemann、三维泊肃叶流、扭曲网格等熵涡、双马赫反射、三维槽道流迁移/长算和二维圆柱
低速/高超声速绕流等检查。数学与算法约定见 [`算法补充.md`](算法补充.md)，完整使用方法见
[`docs/user-manual.md`](docs/user-manual.md)，源码扩展方法见
[`docs/developer-guide.md`](docs/developer-guide.md)，v1.1.0 变化和迁移见
[`docs/release-notes-1.1.0.md`](docs/release-notes-1.1.0.md)。当前能力边界与许可状态分别见
[`docs/known-limitations.md`](docs/known-limitations.md) 和 [`LICENSE.md`](LICENSE.md)。

v1.1.0 的详细范围、P--U 阶段、自动卡口、人工判断及
Git 闭环见 [`docs/v1.1.0-development-plan.md`](docs/v1.1.0-development-plan.md)；
物理容许性、局部通量降阶、壁面载荷/热流、输运和性能公式见
[`算法补充.md`](算法补充.md) 第 11 节；各阶段实现状态以对应设计和验收报告为准。

v2.0.0 已进入分阶段开发：V、W、X、AA、Y、Z 已通过相应自动/人工卡口并合入
`release/v2.0.0`。Z 已实现五种三维 LES、SGS 诊断、可重启接受步时间统计、展向载荷分箱，
以及 BDF2 双时间 LU-SGS 耦合。完整 HIT、槽道、能谱与长时间统计仍属于服务器 S 级
`pending`，短 smoke 不作为物理湍流验证。当前进入阶段 AB，首先处理 Case06 O 网格同一块面上
壁面段与尾迹自连接段并存的守恒权重支持；在工况冻结和服务器预算批准前不启动大型计算。
后续依次进入三维通用能力和用户提供的三维翼型/机翼目标算例。
详细阶段 V--AE、
自动卡口、人工放行和 Git 规则见
[`docs/v2.0.0-development-plan.md`](docs/v2.0.0-development-plan.md)，计划算法见
[`算法补充.md`](算法补充.md) 第 12 节，Case06 只读登记见
[`docs/v2.0.0/case06-intake.md`](docs/v2.0.0/case06-intake.md)，AA 前置人工审查单见
[`docs/v2.0.0/stage-aa-reorder-review.md`](docs/v2.0.0/stage-aa-reorder-review.md)。

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
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix build\install
```

启用 MPI（Windows/MinGW 下使用 Intel MPI）：

```powershell
cmake -S . -B build-mpi -G "MinGW Makefiles" -DWCNS_ENABLE_MPI=ON
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

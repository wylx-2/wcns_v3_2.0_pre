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

v2.0.0 已进入分阶段开发：V、W 已通过自动与人工验收，当前 `stage/v2.0.0-x`
正在实现和验证首个湍流模型 SA-neg。公式、输运、边界、流场/历史输出和模型
checkpoint 已接入候选实现，但 TMR 平板与 NACA0012 定量收敛卡口尚未关闭，
因此 SA-neg 还不得标记为生产验证通过，也尚无 X 候选标签。后续计划仍包括
k--omega SST、标准 k--epsilon、原始 Smagorinsky、尺度相似、混合/动态 Smagorinsky、
WALE、定常/非定常 LU-SGS 和 Roe/Weiss--Smith 低 Mach 预处理，最终依次完成
二维 Case06 NACA0012 及用户后续提供的三维翼型/机翼算例。详细阶段 V--AE、
自动卡口、人工放行和 Git 规则见
[`docs/v2.0.0-development-plan.md`](docs/v2.0.0-development-plan.md)，计划算法见
[`算法补充.md`](算法补充.md) 第 12 节，Case06 只读登记见
[`docs/v2.0.0/case06-intake.md`](docs/v2.0.0/case06-intake.md)。

本开发仓库保留阶段设计、自动测试、人工算例及验收证据。v1.1.0 的确定性内部源码包由
`tools/package_release.py` 从版本提交直接生成，并包含版本验收所需的 Case07；NACA0012 不在
本版本范围内。历史上的独立 `wcns_v3_release` 精简仓库不再作为版本来源真值。

当前程序具备 CGNS 结构多块网格读取、两套独立高阶几何 profile、单 zone 受约束二次剖分、同 rank/MPI 非阻塞 halo 交换、六种界面重构（含保持六点调用契约的 `zero_order`）、Rusanov/HLLC/Roe、WCNS-Euler 空间离散、层流 Navier--Stokes 黏性通量、显式源项和 SSPRK3 推进。v1.1 新增默认关闭的 SSPRK 候选态物理容许性检查、troubled-cell 离散支持传播、逐面受控降阶和整步缩步重试；关闭时保持 v1.0 数值路径。严格配置现完整支持常黏度/Sutherland/Prandtl，边界面 `p_w/T_w/mu_w/Cp/Cf/q_wall/traction` 和积分力、力矩、`Cd/Cl/Cm`。阶段 X 分支另已接入待完成物理验收的 SA-neg 输运和湍流黏性闭合，并支持 `nu_tilde`、`mu_t_over_mu`、源项诊断与壁距离输出。正式入口还支持可配置 MDCD 色散/耗散系数、定常/非定常停止、MPI 全局残差、精确时间事件、CGNS/Tecplot 流场、TXT/Tecplot 历史与统计、多截面 x-z/y-z 监测、manifest，以及可改变 rank 数和叶块划分的 CGNS 检查点重启；二维经典双马赫反射已有专用初场和时变边界。逐步使用说明见 [`docs/user-manual.md`](docs/user-manual.md)，源码二次开发见 [`docs/developer-guide.md`](docs/developer-guide.md)，可复制的完整配置见 [`examples/full_case_template.wcns`](examples/full_case_template.wcns)；简明运行速查仍见 [`docs/runtime-guide.md`](docs/runtime-guide.md)，发布算例的生成、独立重读和矩阵入口见 [`docs/release-validation.md`](docs/release-validation.md)。

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

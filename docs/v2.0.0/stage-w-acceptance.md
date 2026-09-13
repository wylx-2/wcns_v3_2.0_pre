# v2.0.0 阶段 W 自动验收

结论：**通过（2026-09-13）**。阶段 W 已建立通用 RANS/LES 基础设施，但没有注册任何具体
RANS/LES 闭合；当前可运行的 schema 2 组合仍只有 `none + ssprk3 + none`。本报告不把
SA/SST/LES、LU-SGS、低 Mach 或 Case06 物理计算写成已实现。

## 1. 范围与 Git

- 验收实现提交：`dcd0a2a277bc6ce7d35d1ffd42b81060a9785848`；阶段分支：
  `stage/v2.0.0-w`；基线是 V 已合入的 `5abdfdb`。
- 实现包括独立模型场/descriptor、模型 registry、通用标量面通量与边界/容许性、MPI
  `all-gather`、壁距 BVH、唯一黏性闭合入口、源 Jacobian、输出量注册、重启 payload、最小
  schema 2，以及 Case06 平面二维和同 zone 自连接兼容。
- `TurbulenceModelRegistry::create_builtin()` 只注册 `none`；所有后续模型名、`lu_sgs` 和
  `weiss_smith` 在启动前 fail closed。
- 用户目录 `cases/manual/case06_2d_naca0012/` 与
  `cases/manual/case05_3d_turbulent_channel/results/lowmach-longrun-segment01/` 未暂存、未修改。

## 2. 独立 Release 构建与完整回归

两套构建均从不存在的新目录配置，并在实现提交上完成构建与安装：

| 项目 | 结果 |
|---|---|
| OS / CPU | Microsoft Windows NT 10.0.26200.0 / Intel64 Family 6 Model 151，12 logical CPU |
| CMake / generator | 3.28.0 / MinGW Makefiles |
| C/C++ | MinGW-w64 GCC/G++ 8.1.0，Release |
| MPI | Intel MPI，`C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe` |
| 串行 | `build-v2-w-serial`；configure/build/install 成功；CTest **60/60**，20.57 s |
| MPI | `build-v2-w-mpi`；configure/build/install 成功；CTest **110/110**，36.89 s |
| 算法规格 | `python tools/verify_algorithm_spec.py`，**6/6** |
| whitespace | `git diff --check` 通过；仅提示既有 CRLF→LF 规范化，不存在错误 |

MPI 相比阶段 V 的 108 项增加 `wcns.halo_exchange.4` 和 `wcns.mpi_runtime.4`；完整矩阵仍覆盖
release、光滑/间断/强激波、黏性、源项、Case07、输出、checkpoint/restart、失败注入及
1/2/4/8-rank 原回归。

## 3. 通用模型场与输运卡口

1. 平均流 `flow.conservative.components()` 在分配 1/2/3 个模型场前后均严格为 5；空
   `TurbulenceFieldSet` 不持有数据分配。
2. descriptor 名称/唯一性/量纲/下界校验、确定顺序 pack/unpack、签名不匹配拒绝和非法模型值
   拒绝通过；恢复失败不会产生部分合法状态。
3. 周期制造问题对 1/2/3 个标量分别检查二阶面中心对流、中心扩散和解析源组合，观测阶断言
   均为 `p > 1.95`；线性耗散源 SSPRK3 时间加密断言为 `p > 2.9`。
4. 标量面通量逐项验证
   `F_net = mass_flux * phi_upwind - D * grad(phi)·n`，扩散系数为负、尺寸不同或非有限值均拒绝；
   Dirichlet 镜像和外推 ghost 逐项通过。
5. 模型场沿现有多块轴置换 halo 路径执行，1/2/3 分量在 1/2/4 rank 下逐 ghost 与 donor
   interior 完全相等；标量只变换索引，不旋转数值。

## 4. 壁距、闭合与隐式源卡口

- 二维 segment、三维 triangle 最近点、BVH 叶/内部节点、平板/槽道、等距规范键、非有限输入
  和物理 wall patch 提取均通过；16→32 段圆的误差比断言 `>3.9`，8→16 纬向带球面三角化
  的误差比断言 `>3.7`，对应二阶几何收敛。
- MPI 将每个几何元编码为固定十个 `Real`，全收集后规范排序去重；1/2/4 rank 的全局元表、
  最近距离和等距 primitive 选择一致。
- 显式零湍流贡献与原层流黏性通量的 `x/y/z` 五分量和黏度逐位相等；非零对称应力、能量热流
  与 `mu_t` 的解析增量通过，二维非法 z 通量一致拒绝。
- 两变量源项对角块逐项验证为 `alpha*I-J`；源/Jacobian 尺寸、有限性和正时间对角均检查。

## 5. `model=none`、配置、输出和重启

- 阶段 V 与 W 的二维/三维代表性 release freestream 最终 CGNS 用独立验证器以容差 0 比较：
  分别为 640 和 5120 个样本，`max_abs=0`。
- 层流热路径分配探针保持 `first_allocations=6, first_bytes=208`，第二次仍为 `6/208`；新增的
  空模型场没有增加推进期分配。
- schema 1 原完整矩阵全过且不接受 v2 键；schema 2 必填模型/积分器/预处理三个键，当前仅
  `none/ssprk3/none` 可运行。未知名称、模型无关参数和未实现组合全部有拒绝测试。
- descriptor 驱动的模型场可注册到 field quantity registry，并具有量纲缩放；无活动模型时不
  暴露占位字段。模型场 restart payload 校验 descriptor 与真实区长度；schema 2
  `none/ssprk3/none` 明确映射到 schema 1 层流 checkpoint 身份。生产 CGNS 中出现首个模型场
  时由阶段 X 复用该 descriptor/payload 契约写入数组，不能改变当前顺序或兼容规则。

## 6. Case06 只读复探

| 项目 | 结果 |
|---|---|
| 文件 | `naca0012_897x257_str.cgns`，5,545,984 byte |
| SHA-256（验收前/后） | `ADBEC980FB1C1CC2956DCFDFB918ED25F7159F2287C958DB5056C926D725351E` / 相同 |
| `cgnscheck -v` | exit 0；0 error，6 warning，与 intake 一致 |
| schema 2 串行 dry-run | 越过平面二维与自连接，精确退出于部分面张量权重缺口，exit 1 |
| schema 2、2-rank dry-run | 两 rank 同一精确退出，`mpiexec` exit 1，无挂起 |

本轮证明阶段 V 的两个基础兼容缺口已解决，但没有证明 Case06 自由流或物理计算通过。新暴露的
`a partially connected block face does not admit tensor-product conservation weights` 已在 intake
中解释：下 `j` 面同时含物理壁段和自连接尾迹段，需要非张量子区间守恒权重。它是 AB 正式
Case06 验收前置项，不允许通过忽略连接、扩大壁面或改写用户网格绕过。

## 7. 阶段结论

W0--W6 的公共接口和 fail-closed 边界满足阶段目标，独立串行/MPI 自动卡口全过，且层流基线
无数值或分配回归。依据项目负责人 2026-09-13 对 V→W 连续执行且中间不设人工审批的授权，
阶段 W 可创建 `v2.0.0-w-candidate.1` 并以 `--no-ff` 合入 `release/v2.0.0`。

连续授权到 W 为止；阶段 X（SA-neg）仍须等待新的人工批准。本报告也不免除 X 对生产 CGNS
模型场 checkpoint、SA 边界/降阶及高阶 MMS 的阶段专属卡口。

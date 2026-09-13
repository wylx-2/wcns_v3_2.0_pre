# v2.0.0 阶段 V 设计：规格、输入与基线冻结

状态：**2026-09-13 设计冻结。项目负责人已批准 v2.0.0 计划，并授权 V→W 连续推进，中间不设
人工等待；本文件把该授权作为阶段 V 的 G0/G9 书面依据。**

基线：`release/v2.0.0` 提交 `f72df39`。阶段分支：`stage/v2.0.0-v`。

## 1. 阶段边界

V 只冻结需求、公式选择、配置草案、用户输入身份、测试矩阵和可重复基线，不修改生产 C++、
当前 schema 1 parser 或数值轨迹。W 及以后才能实现接口和代码。

必做交付物：

1. 本设计和正式验收报告；
2. [`config-schema-2-draft.md`](config-schema-2-draft.md)；
3. [`target-wing-intake-template.md`](target-wing-intake-template.md)；
4. [`case06-v-smoke.wcns`](case06-v-smoke.wcns) 及 Case06 只读检查证据；
5. 当前串行/MPI Release、规格检查、Case06 dry-run 和环境基线；
6. V--AE 测试映射与人工判断清单。

## 2. 冻结的模型选择

### 2.1 RANS

| 配置名 | 冻结变体 | v2.0.0 目标支持级别 |
|---|---|---|
| `none` | v1.1.0 层流路径 | 生产，兼容基线 |
| `sa_neg` | NASA TMR fully turbulent SA-noft2-neg；trip/压缩修正关闭 | 生产 |
| `k_omega_sst` | SST-2003m；明确 $F_1/F_2$、cross diffusion 和 production limiter | 生产 |
| `k_epsilon` | 标准高 Reynolds 数模型，只与显式 wall function 组合 | 实验目标；阶段 Y 决定是否升级 |

RANS 总能不含 $k$；各向同性 $-2\rho k/3$ 显式进入应力。候选 $Pr_t=0.9$，最终键和边界值在
W/Y 设计冻结。模型常数以《算法补充》12.4--12.6 和阶段 X/Y 的逐式表为准。

### 2.2 LES

LES 只允许三维非定常 BDF2/双时间计算。必做模型为：

| 配置名 | 冻结基线 | 支持目标 |
|---|---|---|
| `smagorinsky` | $C_s=0.17$；默认无壁面阻尼，`van_driest` 必须显式选择 | 生产 |
| `scale_similarity` | Bardina/Leonard 结构项，$C_B=1$，允许 backscatter | 受限生产 |
| `mixed_smagorinsky_similarity` | 上述相似项加 Smagorinsky 耗散 | 生产 |
| `dynamic_smagorinsky` | Germano--Lilly，测试滤波比 2，局部确定性平均 | 生产 |
| `wale` | $C_w=0.325$ | 生产 |

Vreman 等更多 SGS 模型为条件增强，不预留可被误用的配置名。动态系数的具体 stencil、分母保护、
负值/上限和物理边界滤波必须在 Z0 冻结并写入 restart signature；V 的常数只是版本级初始选择，
不授权不经 Z0 设计直接实现。

离散滤波基线为三点张量积核 `[1/4, 1/2, 1/4]`，严格保持常量；$\Delta=V_c^{1/3}$，
$\widehat\Delta/\Delta=2$。高度各向异性和边界闭合必须以滤波响应及 MPI 等价证据修订，不能
按 rank 隐式改变。

## 3. 冻结的隐式与低 Mach 路线

1. 唯一新增隐式算法为 scalar-diagonal 起步、允许模型源分块对角的 Yoon--Jameson LU-SGS；
   不把 Jacobi 或 point-implicit 对外命名为 LU-SGS。
2. 定常：本地伪时间、LU-SGS 前扫/后扫、残差与载荷窗口联合停止。
3. 非定常：BDF1 首步，之后 BDF2；每个物理层通过双时间 LU-SGS 收敛，伪迭代不推进物理时间、
   不输出物理快照、不累计统计。
4. 未预处理 LU-SGS 覆盖 Rusanov/HLLC/Roe 权威空间残差；近似 Jacobian 可以用一致谱半径，
   但 manifest 必须同时记录 residual flux 与 Jacobian 类型。
5. Weiss--Smith 首版只和 Roe 组合，同时进入伪导数、耗散、谱半径、远场和回退。定常收敛到
   $R(Q)=0$；非定常收敛到守恒 BDF 方程，预处理不得进入真实物理时间导数。
6. MPI 首版采用 rank 间 block-Jacobi、rank 内按冻结块/单元序执行 LU-SGS；每个完整前后扫后
   交换增量/状态 halo。必须验证随 rank 改变只影响收敛路径、不影响收敛解。

## 4. 配置与兼容决策

- schema 2 由 W 开始实现；schema 1 严格迁移为 `model=none, preconditioner=none,
  integrator=ssprk3`。
- 新增数值键均进入 summary、manifest 和 restart signature。纯输出选择不进入数值签名。
- LES 模型在二维、steady 或缺少物理时间配置时初始化前拒绝。
- k-epsilon 与 `resolved` 壁面、低 Mach 与非 Roe、BDF2 缺少历史且未声明 BDF1 启动等组合拒绝。
- 本阶段草案不是生产 parser；W 必须从合法/非法表生成解析测试。

## 5. 用户文件和算例输入

| 路径 | 状态 | 阶段处理 |
|---|---|---|
| `cases/manual/case06_2d_naca0012/` | 未跟踪用户网格；SHA-256 已登记 | 只读 CGNS 检查和 dry-run，不暂存 |
| `cases/manual/case05_3d_turbulent_channel/results/lowmach-longrun-segment01/` | 未跟踪长算结果，16 文件、118,067,120 byte | 完全保留，不读取数值内容、不暂存 |

Case06 的 CGNS 结构检查成功，但当前运行路径的兼容性探针稳定暴露两个既有缺口：二维网格以
`PhysicalDimension=3` 保存且 `Z=0`，以及 O 型网格尾迹切口使用同 zone 自连接。它们分别触发
`2D metrics currently require two-dimensional physical space` 和
`self-connectivity is not supported by stage D`。V 将这两个精确失败签名作为 intake 结果，不把
它们解释成网格损坏，也不声称自由流已经启动；W 增加基础兼容与回归，AB 再做正式 Case06
自由流/工况验证。没有用户工况时不得给物理结论。最终三维目标使用 intake 模板，原文件永不覆盖。

## 6. V--AE 测试映射

| 能力 | 最早公式/单元 | 数值验证 | 系统算例 |
|---|---|---|---|
| 动态附加场/halo/restart | W | 标量 MMS | W/AE |
| 壁距 | W | 平板/槽道/圆/球 | X/Y/AB/AD |
| SA-neg | X | MMS、平板、TMR NACA | AB/AE |
| SST/k-epsilon | Y | MMS、槽道、逆压梯度 | AB/AE |
| LES 滤波/SGS | Z | Fourier、TGV、HIT、槽道 | AC/AD/AE |
| 壁面/时间/截面统计 | Z | 解析信号、加和、restart | AB/AC/AD |
| LU-SGS steady | AA | MMS、Couette、圆柱 | AB/AD/AE |
| LU-SGS BDF2/双时间 | AA | 时间制造解、涡/TGV | AB(URANS)/AC/AD |
| Weiss--Smith | AA | $M=10^{-1..-3}$、$\beta=1$ | AB/AD/AE |

## 7. 阶段 V 自动卡口

1. `git diff --check`、本地文档链接和《算法补充》规格检查通过；
2. 从干净、独立的串行/MPI Release 目录配置、构建和安装；完整 CTest 全过；
3. `case06-v-smoke.wcns` 由正式 `wcns_run --dry-run` 串行和至少 2-rank 执行兼容性探针；
   `cgnscheck -v` 无 error、已登记 warning；探针必须成功或精确命中本设计第 5 节登记的两个
   既有缺口，其他错误一律失败；原网格哈希前后不变；
4. schema 2 草案的合法/非法组合、迁移、摘要和 restart 影响无矛盾；
5. 验收报告记录 OS/CPU/CMake/C++/MPI、提交、命令、测试数、耗时、工作树和用户文件；
6. 所有生产源码与 `f72df39` 一致。

自动全过后创建 `v2.0.0-v-candidate.1`。依据 2026-09-13 连续授权，不等待 `v-approved` 人工
标签，以 `--no-ff` 合入 `release/v2.0.0` 并立即建立阶段 W；任一卡口出现未登记失败则停止连续
流程。登记的兼容缺口不被视为功能已经通过，必须保留到对应修复回归。

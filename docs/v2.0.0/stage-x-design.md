# v2.0.0 阶段 X 设计：SA-neg

状态：**2026-09-13 设计已冻结；实现与非长算自动卡口已通过，定量物理验收因
定常加速阶段依赖未关闭。** 阶段 W 已经人工验收；本阶段只实现首个 RANS 模型
`sa_neg`，不提前实现 SST、LES、LU-SGS 或低 Mach 预处理。

## 1. 模型身份

`sa_neg` 严格采用 NASA Turbulence Modeling Resource（TMR）的 SA-neg：

- 正分支保留标准 $f_{t2}=c_{t3}\exp(-c_{t4}\chi^2)$；fully turbulent 只关闭
  $f_{t1}\Delta U^2$ trip 源，因此不得再称 `SA-noft2-neg`；
- 采用 modified-vorticity 方法 1(c)，负分支采用 $f_n$、负分支生产和正号破坏贡献；
- 压缩、旋转/曲率、粗糙壁、转捩和 QCR 修正全部关闭；
- $\widetilde\nu<0$ 是合法求解状态，不能用零裁剪替代 SA-neg；负分支令 $\nu_t=0$；
- 常数、逐式定义、守恒形式和边界语义以《算法补充》12.4 为唯一代码规格。

推进量为 $\rho\widetilde\nu$，场存储的 `nu_tilde` 为运动黏度尺度工作变量。平均流湍流应力采用
Boussinesq 偏应力，RANS 总能不另含 $k$；湍流热流使用独立 $Pr_t=0.9$。模型对流复用平均流
质量通量，按 TMR 数值验证建议使用二阶上风；扩散和源项逐项守恒组合，不能把源项并入未记录
的人工耗散。

## 2. 软件契约

### 2.1 配置

schema 2 新增且只在 `turbulence.model=sa_neg` 时合法：

```text
turbulence.prandtl = 0.9
turbulence.wall_treatment = resolved
turbulence.sa.farfield_nu_tilde_ratio = 3
turbulence.sa.source_treatment = explicit | local_implicit
```

远场比必须在 `[3,5]`。SA-neg 要求黏性方程和至少一个 resolved no-slip wall；阶段 X 不支持
wall function。`local_implicit` 只对 SA 局部源使用解析 Jacobian，外层仍是 SSPRK3；它不是
阶段 AA 的 LU-SGS。SA-neg 暂不与整步重试/局部降阶事务组合，非法组合启动前拒绝。

### 2.2 场、边界和输出

descriptor 顺序固定为：

1. `nu_tilde`：唯一推进场；
2. `mu_t_over_mu`：非负诊断；
3. `sa_production`：带符号方程贡献；
4. `sa_destruction`：带符号方程贡献，正分支通常为负；
5. `wall_distance`：全局最近 no-slip wall 距离；
6. `sa_negative_branch`：0/1 分支标志。

连接面先复用阶段 W 的标量 halo。无滑移壁以镜像 ghost 强制 $\widetilde\nu_w=0$；farfield 和
inflow 使用配置远场值；outflow、symmetry 和其余非壁物理面零阶外推。边/角 ghost 不属于输运
模板契约。黏性闭合和模型输运使用独立 MPI tag，不能交叉匹配。

checkpoint 版本在活动 SA 模型时升级为 2，保存模型 descriptor signature 和真实区
`NuTilde`，ghost、壁距及诊断在恢复后重建。`none` 继续逐位使用 v1 文件格式；模型、常数、
源处理或 descriptor 不匹配必须在推进前失败。

## 3. 离散与稳定性

- `nu_tilde` 梯度使用当前 profile 的守恒高阶度量梯度；扩散面系数为左右
  $\rho(\nu+\widetilde\nu f_n)/\sigma$ 的算术平均。
- 模型对流使用同一计算面的平均流质量通量和二阶上风值。模型面通量按 owner 交换后再形成
  高阶通量差分，以维持跨块唯一共享面。
- 显式稳定步同时包括分子/涡黏扩散谱半径和正的局部源 Jacobian。局部隐式更新为
  $\Delta q=\Delta t R/(1-\Delta t\,\partial S/\partial\widetilde\nu)$；分母非正即失败。
- 不设置 `nu_tilde` floor，不隐式限制负分支占比，不用大范围裁剪掩盖源项或边界错误。

## 4. 分层验收算例

验收输入、散列和阈值的机器可读副本位于
[`cases/validation/sa_tmr/manifest.json`](../../cases/validation/sa_tmr/manifest.json)。层级不可互相替代：

1. **公式/MMS（代码自有）**：正负分支逐点标量、源 Jacobian 有限差分、SA 专属连续制造解的
   $L_1/L_2/L_\infty$ 和观测阶；这是代码验证，不是物理模型验证。
2. **TMR finite flat plate（独立物理验证）**：$M=0.2$、$Re=5\times10^6$、$Pr=0.72$、
   $Pr_t=0.9$、Sutherland、绝热壁、$\widetilde\nu_\infty/\nu_\infty=3$。只使用 TMR 生成器
   产生的嵌套网格；level 8 仅 smoke，至少 level 6/5/4 构成趋势，比较 $C_f(x=0.8697742)$、
   总阻力、速度/工作变量剖面和 $y^+$。
3. **TMR NACA0012 validation family（独立物理验证）**：原始公开结构 C-grid，$M=0.15$、
   $Re_c=6\times10^6$、$\alpha=0,10,15^\circ$、远场约 $500c$、绝热壁、完全湍流。
   113x33 只做导入/有限性 smoke；225x65、449x129、897x257 检查系统网格趋势，并只和同一
   工况的 TMR `CL/CD/Cp/Cf` 数据比较。
4. **TMR NACA0012 numerical family II（离散敏感性）**：$\alpha=10^\circ$，使用修正后的闭合
   尖尾缘 Family II 和明确的“有/无 point-vortex”远场分支，检查 `CM` 及积分量渐近区间。
   Family II 与 validation family 几何不同，禁止混用结果或用一种网格的参考值验收另一种。

用户 `cases/manual/case06_2d_naca0012` 明确排除在 X 的模型验收之外：它属于阶段 AB 的项目算例，
既不能替代 TMR 独立网格，也不会在 X 修改。结构多块曲面 smoke 使用仓库生成网格，只验证接口、
halo 与 rank 一致性，同样不能替代上述两个物理基准。

## 5. 自动卡口

X 候选必须同时满足：

- 串行与 MPI Release 构建、完整 CTest、算法规格检查和 `git diff --check` 全过；
- 公式常数和正/负分支与冻结字面量在 `5e-13` 相对/绝对容差内，解析源 Jacobian与中心差分在
  `2e-6` 相对容差内；SA MMS 三范数观测阶均不低于 1.90；
- 1/2/4 rank 的模型场、残差、负分支计数及 restart 连续结果在冻结舍入容差内一致；
- TMR 平板至少三个非 smoke 网格呈系统收敛，level 4 的 $C_f$ 与三代码公开包络一致，积分阻力、
  剖面和 $y^+$ 通过 manifest 阈值；
- NACA validation family 的三层非 smoke 网格完成，897x257 的 `CL/CD/Cp/Cf` 通过阈值；
  numerical Family II 的 `CM` 与选择的 point-vortex 分支渐近区间一致；
- 均匀远场、壁面零值、显式/局部隐式对照、负分支和 checkpoint/restart 都无非有限值、隐藏
  floor、未登记回退或部分输出；
- 自动报告必须列出命令、输入散列、rank、运行时间、迭代/残差、积分量、曲线误差和失败原因。

任一 TMR 网格没有被当前结构 CGNS/PLOT3D 入口忠实读取、工况/边界未冻结、或计算尚未收敛，
阶段状态只能是“进行中/未通过”，不能用较宽容差、最粗网格或 Case06 生成候选标签。

## 6. Git 与人工判断

全部 X 工作留在 `stage/v2.0.0-x`。通过自动卡口后才允许创建不可移动的
`v2.0.0-x-candidate.1`；随后停下，由人工查看源项/破坏项分布、负分支占比、平板对数层与
摩擦、NACA `Cp/Cf` 曲线和网格外推。人工接受后才以 `--no-ff` 合入 `release/v2.0.0`，再开始 Y。

权威输入：NASA TMR SA equations、Finite Flat Plate Numerical Analysis、NACA0012 Validation、
NACA0012 Numerical Analysis 及其各自 grid/results 页面。抓取脚本固定 URL 与 SHA-256；网页内容
更新时必须人工审查并更新 manifest，不能静默接受新文件。

本阶段直接核对的官方页面为 [SA 方程](https://tmbwg.github.io/turbmodels/spalart.html)、
[finite flat plate 数值验证](https://tmbwg.github.io/turbmodels/finiteflatplatenumerics_val.html)、
[NACA0012 validation](https://tmbwg.github.io/turbmodels/naca0012_val.html) 和
[NACA0012 numerical analysis](https://tmbwg.github.io/turbmodels/naca0012numerics_val.html)。

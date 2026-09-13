# v2.0.0 阶段 AA 设计：LU-SGS、双时间与低 Mach

状态：**2026-09-13 AA0 设计冻结；按已批准的 `AA → X-B → Y` 连续授权执行。**

基线为 `stage/v2.0.0-aa` 建立点 `332f23d`，其生产代码锚点仍为 X-A 的 `d32010a`。
阶段只修改层流/SA 与模型无关 implicit contract；SST/k-epsilon/LES 专属耦合分别在 Y/Z 验收。

## 1. 数学与离散身份

权威公式为《算法补充》12.7、12.8、12.15。空间残差仍由用户选择的 WCNS profile、重构和
Rusanov/HLLC/Roe 形成，隐式推进不得另建低阶物理残差。一次定常伪迭代求解

$$
\left(\frac{I}{\Delta\tau_i}+\frac{\partial R_i}{\partial Q}\right)\Delta Q_i=R_i,
\qquad
\Delta\tau_i=CFL_\tau/\sigma_i,
$$

其中程序现有 `residual` 的符号是 $dQ/dt=R$，$\sigma_i$ 为单位体积无粘、黏性和模型扩散谱半径
之和。采用 scalar-spectral LU-SGS：$D_i=(1/\Delta\tau_i+\sigma_i)I-J_{src,i}$，本地词典序
`(block,k,j,i)` 做完整前扫和逆序后扫；本 rank 已更新邻居进入 $L/U$，远程和块连接增量在一次
sweep 内固定为零，相当于 block-Jacobi/块内 SGS。每个外层迭代前刷新物理 halo，1/2/4 rank
必须收敛到同一稳态解；并行路径可以有不同迭代数，但不得有不同终态或停止语义。

平均流五分量共享同一标量谱对角；SA 的 $\rho\widetilde\nu$ 另加解析源 Jacobian。更新先写候选
缓冲，逐单元检查 $\rho,p,T,e$ 与模型合法域；任一分量非法时按 $1,1/2,\ldots,1/128$ 对全体增量
统一回溯。仍非法则整次伪迭代失败并恢复输入，不做 clipping 或部分提交。

## 2. 非定常 BDF2 双时间

非定常 `lu_sgs` 的每个物理步固定求解

$$
G(Q^{n+1})=R(Q^{n+1})-
\frac{3Q^{n+1}-4Q^n+Q^{n-1}}{2\Delta t}=0.
$$

首个没有合法历史层的物理步使用 BDF1。内层使用与定常完全相同的 LU-SGS 扫序，在对角增加
$1/\Delta t$（BDF1）或 $3/(2\Delta t)$（BDF2）。只有 `L2 <= abs_tol` 或同时满足
`L2/L2_initial <= rel_tol` 且守恒缺陷通过，物理层才整体接受；达到最大内迭代时恢复 $Q^n$ 并
返回数值失败。输出、统计和正式 step/time 只在接受后更新，最终物理步由 driver 截短以精确命中
`run.t_end`。

BDF 历史保存在每个块的可选 `ImplicitTimeHistory`，包括上一物理层平均流和全部 transported
模型变量。非定常 LU-SGS checkpoint 版本升级并写入历史层、有效标志、物理步长、积分器和
预处理签名；缺历史时只允许显式记录的 BDF1 重启，不能静默伪造 BDF2。

## 3. Weiss--Smith 预处理

配置只允许 `algorithm.riemann=roe + time.integrator=lu_sgs + preconditioner.type=weiss_smith`。
参考速度和预处理参数严格采用 12.7：

$$
U_r=\min[a,\max(|\boldsymbol u|,M_{cut}a,C_\nu\nu/\ell)],\qquad
\beta=(U_r/a)^2.
$$

$\ell=2V/(\sum_f |S_f|)$；无黏计算令黏性项为零。面上使用左右单元限制的较大 $\beta$，避免跨面
欠耗散。实现提供可独立测试的 $\Gamma$、$\Gamma^{-1}$、五个预处理特征速度和 Roe 耗散变换；
`beta=1` 必须在 `5e-13` 内恢复普通 Roe 通量和谱半径。低 Mach Roe 失败时只允许使用相同
$\beta$ 的最低阶预处理耗散，不能回退到未预处理 HLLC/Rusanov。预处理只进入数值耗散、谱半径
和伪时间对角；BDF 物理导数仍对守恒量 $Q$ 离散。

首版默认 `M_cut=1e-3`、`C_nu=1.0`，合法范围分别为 `(0,1]` 和 `[0,10]`。局部
$M\ge1$ 时严格取 `beta=1`。配置摘要、manifest、checkpoint signature 同时记录 Riemann 残差
和 approximate Jacobian 类型。

## 4. 配置与失败语义

schema 2 新增：

```text
time.integrator = ssprk3 | lu_sgs
time.physical.scheme = bdf2
time.physical.step = <positive>                 # 仅 unsteady lu_sgs 必填
time.dual_time.max_iterations = 100
time.dual_time.absolute_tolerance = 1e-10
time.dual_time.relative_tolerance = 1e-8
time.dual_time.cfl = 5
lu_sgs.sweeps = 1
lu_sgs.jacobian = scalar_spectral
lu_sgs.relaxation = 1
preconditioner.type = none | weiss_smith
preconditioner.mach_cutoff = 1e-3
preconditioner.viscous_cutoff = 1
```

定常 `lu_sgs` 使用 `run.cfl` 作为本地伪时间 CFL；非定常物理步使用 `time.physical.step`，内层使用
`time.dual_time.cfl`。首版 `lu_sgs.sweeps` 必须为 1，`jacobian` 必须为 `scalar_spectral`，避免
接受尚未实现的配置。`ssprk3` 不接受任何 `time.physical.*`、`time.dual_time.*` 或 `lu_sgs.*`；
`preconditioner=none` 不接受 cutoff 键。SA 的显式和 local-implicit 源配置均可运行，但 LU-SGS
始终把同一解析源 Jacobian 纳入模型对角。

## 5. 自动验收阈值

1. 公式：$\Gamma\Gamma^{-1}$ 最大误差 `<=5e-12`；特征速度与独立标量式 `<=5e-13`；
   `beta=1` Roe 通量/谱半径 `<=5e-13`；法向旋转和二维退化 `<=2e-12`。
2. LU-SGS：常量场增量 `<=1e-14`；线性 scalar 系统与稠密解 `<=2e-11`；前/后扫顺序、块边界
   固定增量和回溯事务由哨兵测试覆盖。
3. 稳态：制造源/Couette/Poiseuille 的 LU-SGS 终态与解析或收紧 SSPRK3 的三范数差
   `<=2e-8`；达到同一残差阈值所需外迭代和墙钟均优于全局 SSPRK3，具体改善率在运行前写入
   `stage-aa-acceptance.md`，未达即失败而不是修改阈值。
4. 非定常：平移涡或线性衰减问题用三层物理步长，BDF2 $L_2$ 观测阶 `>=1.90`；每层最终内残差
   `<=max(1e-10,1e-8*initial)`；未收敛层不推进 time/step。
5. 低 Mach：$M=10^{-1},10^{-2},10^{-3}$ 的压力扰动除以 $M^2$ 保持在冻结包络内，速度误差不随
   $1/M$ 增长；`M>=1` 与无预处理 Roe 一致。
6. SA：TMR 平板 level 6 和 NACA0012 225 完成启动、有限性和同停止标准加速对照；完整物理趋势
   属于 X-B，不用最粗网格或 Case06 替代。
7. 工程：串行/MPI Release 干净构建与安装、全部 CTest、1/2/4 rank、同/异 rank checkpoint
   恢复通过；schema 1 和 `ssprk3+none` 的冻结摘要、restart 与数值路径不变。

## 6. 提交与候选

按“设计/配置与公式核—LU-SGS—双时间/restart—低 Mach—测试和文档”拆分提交。所有 AA 自动
卡口通过后创建 `v2.0.0-aa-candidate.N`，按连续授权非快进回合到 X；任一卡口失败则保留报告并
停止。用户 Case06 和 case05 long-run segment01 始终不暂存、不修改。


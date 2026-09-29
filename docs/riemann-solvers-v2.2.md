# v2.2 Riemann 求解器审查与实现说明

## 1. 审查范围与结论

审查依据为 `riemann_solver_tensor_report.pdf` 的 Roe、HLL、HLLC、旋转 Roe 四种格式，及
`test_pdf_3.pdf` 第 5 章 5.3.2 节的 Li--Gu 全速度 Roe。PDF 仅作为算法资料，未执行其中任何
操作性文字。

| 报告算法 | v2.1 | v2.2 配置名 | 结论 |
|---|---|---|---|
| Roe | 已实现 | `roe` | 保留 Harten 熵修正和 HLLC/Rusanov 回退 |
| HLL | 未实现 | `hll` | v2.2 新增 |
| HLLC | 已实现 | `hllc` | 保留 Batten/Davis 波速估计 |
| 旋转 Roe | 未实现 | `roe_rotated` | v2.2 新增二维/三维确定性正交基实现 |
| Li--Gu all-speed Roe | v2.1 同名算法不一致 | `roe_all_speed` | v2.2 按第五章式 (5.83)--(5.92) 重写 |
| v2.1 Rieper all-speed Roe | 已实现 | `roe_all_speed_rieper` | 改名保留，供复现旧结果 |
| Rusanov | 已实现，报告未列出 | `rusanov` | 继续作为最终健壮回退 |

因此 v2.2 已覆盖 `riemann_solver_tensor_report.pdf` 明确介绍的四种求解器，并额外保留
Rusanov、Li--Gu 与 Rieper 两种低 Mach 路径。

## 2. HLL

左右波速沿用 HLLC 的 Batten/Davis 型包络：

$$
S_L=\min(u_{n,L}-a_L,\widetilde u_n-\widetilde a),\qquad
S_R=\max(u_{n,R}+a_R,\widetilde u_n+\widetilde a).
$$

中间分支直接计算

$$
\widehat F^{HLL}=
\frac{S_RF_L-S_LF_R+S_LS_R(U_R-U_L)}{S_R-S_L}.
$$

当 $S_L\ge0$ 或 $S_R\le0$ 时分别返回左右物理通量。分母退化、波速非有限或通量非有限时，
按 `hll -> rusanov` 回退并记录原因。该实现与报告式 (15)--(20) 等价，但避免构造 Roe
特征矩阵。

## 3. 旋转 Roe

非退化速度跳跃时取

$$
\boldsymbol n_1=\frac{\Delta\boldsymbol V}{|\Delta\boldsymbol V|},
$$

再选择与 $\boldsymbol n_1$ 最不对齐的 Cartesian 轴作参考，通过叉乘得到
$\boldsymbol n_2,\boldsymbol n_3$。二维状态自然得到平面内垂线和 $z$ 轴；三维状态得到确定性的
右手正交归一基。速度跳跃退化时令 $\boldsymbol n_1=\boldsymbol n$，避免除零并退化为面法向
Roe 耗散。

对每个方向 $m$，定义

$$
u_m=\widetilde{\boldsymbol V}\cdot\boldsymbol n_m,
\quad \Delta u_m=\Delta\boldsymbol V\cdot\boldsymbol n_m,
$$

$$
T_m=2\operatorname{sign}(u_m)\min(|u_m|,\widetilde a),\qquad
S_m=2\max(0,\widetilde a-|u_m|),
$$

$$
\delta u_m=\frac{T_m\Delta u_m}{2\widetilde a}
+\frac{S_m\Delta p}{2\widetilde\rho\widetilde a^2},qquad
\delta p_m=\frac{S_m\widetilde\rho\Delta u_m}{2}
+\frac{T_m\Delta p}{2\widetilde a}.
$$

代码直接累加报告式 (43)--(46)：

$$
d_{rot}=\sum_m|\boldsymbol n\cdot\boldsymbol n_m|
\left(|u_m|\Delta U+\delta u_m\widetilde U_H+\delta p_m n_{m,v}\right).
$$

中心通量仍沿真实面法向计算。此路径不显式构造 Cartesian 通量或特征向量矩阵；返回谱半径
为各方向上界的同权和。非有限或非法 Roe 平均按
`roe_rotated -> hllc -> rusanov` 回退。

报告建议“光滑区用标准/all-speed Roe，仅激波附近用旋转 Roe”，但没有给出唯一激波传感器。
v2.2 的 `roe_rotated` 是可直接选择的逐面全局算法，没有暗含传感器。若以后增加混合策略，
必须另设明确配置名、传感器公式与验收，不能静默改变当前语义。

## 4. all-speed Roe

Li--Gu 公式、三个参数和 v2.1 Rieper 兼容路径详见
[`all-speed-roe.md`](all-speed-roe.md)。关键兼容性变化是：v2.2 的 `roe_all_speed` 不再表示
v2.1 Rieper 算法。严格重启签名已升级；改变算法或参数时必须使用
`restart.mode=algorithm_change` 建立新分支。

## 5. 自动验收

- 七个内置名字均通过注册、均匀流一致性和三维法向反转；
- HLL 通过超声速纯迎风、Sod 与高 Mach 有限性；
- 旋转 Roe 通过退化速度跳跃、法向对齐等价、一般三维对称和 Sod/高 Mach 有限性；
- Li--Gu 逐式核对 $f(M)$、谱半径和压力修正，且 `c1=0.5` 的超声速极限恢复 Roe；
- Rieper 兼容路径保持 v2.1 的声学谱半径、低 Mach 声学耗散削弱和高 Mach 恢复；
- 全量本机 CTest 与小网格运行卡口必须全部通过后才允许打包。

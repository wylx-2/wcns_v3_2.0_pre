# WCNS schema 2 配置草案（阶段 V 冻结输入）

状态：**AA 已实现 LU-SGS、BDF2 双时间和 Weiss--Smith；阶段 Y 候选已实现 SST-2003m、
实验性标准 k-epsilon 与两方程隐式耦合，正在执行定量算例卡口。** 当前可运行选择为
`none|sa_neg|k_omega_sst|k_epsilon`；LES 键仍由 Z 实现，尚未实现的键必须明确拒绝而不是
静默降级。每次扩展必须同步合法/非法配置测试、summary、manifest 和 restart signature。

## 1. 兼容迁移

schema 1 未出现 v2 键时严格迁移为：

```text
turbulence.model = none
time.integrator = ssprk3
preconditioner.type = none
```

旧 `run.mode`、`run.cfl`、`run.t_end` 和 steady 停止键继续解释。schema 2 拒绝未知键、重复键、
非有限数、模型无关参数和不支持组合。

层流最小合法配置片段为：

```text
schema_version = 2
turbulence.model = none
time.integrator = ssprk3
preconditioner.type = none
```

schema 1 不接受上述 v2 键；其缺省迁移发生在内部语义层，不改变 schema 1 的文本摘要和旧重启
签名。`none` 还拒绝 `turbulence.prandtl` 与 `turbulence.wall_treatment`，防止无效参数污染身份。

## 2. 键空间

| 键 | 候选值/类型 | 数值签名 | 约束 |
|---|---|---|---|
| `turbulence.model` | `none|sa_neg|k_omega_sst|k_epsilon|smagorinsky|scale_similarity|mixed_smagorinsky_similarity|dynamic_smagorinsky|wale` | 是 | LES 仅 3D unsteady |
| `turbulence.experimental` | bool | 是 | 当前只允许 k-epsilon，且必须为 true |
| `turbulence.prandtl` | positive real | 是 | RANS 模型启用时 |
| `turbulence.wall_treatment` | `resolved|wall_function` | 是 | k-epsilon 必须 wall_function |
| `turbulence.freestream.intensity` | real in `(0,1]` | 是 | SST/k-epsilon |
| `turbulence.freestream.length_scale` | positive real | 是 | SST/k-epsilon |
| `turbulence.model_floor` | positive real | 是 | SST/k-epsilon；运算保护及有计数的接受态正值投影 |
| `turbulence.two_equation.source_treatment` | `explicit|local_implicit` | 是 | SST/k-epsilon |
| `turbulence.wall_function.y_plus_min/max` | positive real, max > min | 是 | wall_function 适用区间与诊断 |
| `turbulence.sa.farfield_nu_tilde_ratio` | real in `[3,5]` | 是 | 仅 `sa_neg`；基准 3 |
| `turbulence.sa.source_treatment` | `explicit|local_implicit` | 是 | 仅 `sa_neg`；只处理局部源 Jacobian |
| `les.filter.type` | `box3_tensor` | 是 | LES 必需 |
| `les.filter.width_ratio` | positive real | 是 | 基线 1 |
| `les.test_filter.ratio` | real > 1 | 是 | 动态/相似模型必需，基线 2 |
| `les.smagorinsky.cs` | `[0,0.3]` | 是 | Smag/混合，基线 0.17 |
| `les.similarity.cb` | finite real | 是 | 相似/混合，基线 1 |
| `les.dynamic.*` | average/clipping/floor | 是 | 仅动态模型 |
| `les.wale.cw` | positive real | 是 | WALE，基线 0.325 |
| `les.sgs_prandtl` | positive real | 是 | LES 必需 |
| `time.integrator` | `ssprk3|lu_sgs` | 是 | v2 隐式选 lu_sgs |
| `time.physical.scheme` | `bdf2` | 是 | unsteady lu_sgs；无历史首步自动 BDF1 |
| `time.physical.step` | positive real | 是 | unsteady lu_sgs 必填；最终层可截短命中 t_end |
| `time.dual_time.max_iterations` | positive integer | 是 | 默认 100 |
| `time.dual_time.absolute_tolerance` | positive real | 是 | 默认 `1e-10` |
| `time.dual_time.relative_tolerance` | positive real | 是 | 默认 `1e-8` |
| `time.dual_time.cfl` | positive real | 是 | 默认 5 |
| `lu_sgs.jacobian` | `scalar_spectral` | 是 | 首版唯一值；共享标量对角，平均流面块保留 Euler $A_n$ |
| `lu_sgs.sweeps` | integer in `[1,4]` | 是 | 默认 1；额外扫为同一线性系统的缺陷修正 |
| `lu_sgs.relaxation` | real in `(0,1]` | 是 | 默认 1；非法候选另做全局回溯 |
| `preconditioner.type` | `none|weiss_smith` | 是 | weiss_smith 只配 Roe+LU-SGS |
| `preconditioner.mach_cutoff` | real in `(0,1]` | 是 | Weiss--Smith 默认 `1e-3` |
| `preconditioner.viscous_cutoff` | real in `[0,10]` | 是 | Weiss--Smith 默认 1 |
| `boundary.point_vortex.enabled` | bool | 是 | 默认 false；仅二维 farfield，NACA numerical-analysis 对照用 |
| `boundary.point_vortex.lift_coefficient` | finite real | 是 | 给定 $C_{L,pv}$，不从迭代载荷反馈 |
| `boundary.point_vortex.center_x/y` | finite real | 是 | 默认 `(0.25,0)` |
| `boundary.point_vortex.chord` | positive real | 是 | 默认 1；进入点涡环量 |
| `geometry.metric.fallback` | `strict|phenglei_finite_volume` | 是 | 默认 strict；后者仅 PH profile 合法且必须报告计数 |
| `geometry.metric.maximum_reference_relative_difference` | non-negative real | 是 | 默认 0.20；Family II 尖尾缘冻结 0.35 |
| `statistics.time.*` | start/end/every/weight | 否；累加器身份单列 | 只累计接受物理步 |

## 3. 合法草案

阶段 X 的 SA-neg 显式基准：

```text
schema_version = 2
turbulence.model = sa_neg
turbulence.prandtl = 0.9
turbulence.wall_treatment = resolved
turbulence.sa.farfield_nu_tilde_ratio = 3
turbulence.sa.source_treatment = explicit
time.integrator = ssprk3
preconditioner.type = none
```

`local_implicit` 只把解析源 Jacobian 加入每个 SSPRK stage 的局部标量更新，不等同于 LU-SGS，
也不改变 `time.integrator=ssprk3`。SA-neg 当前要求黏性求解和至少一个 no-slip resolved wall；
阶段 X 不启用 wall function、trip 或压缩修正。

Family II 带 point-vortex 对照在同一 SA 配置上增加：

```text
boundary.point_vortex.enabled = true
boundary.point_vortex.lift_coefficient = 1.09125
boundary.point_vortex.center_x = 0.25
boundary.point_vortex.center_y = 0
boundary.point_vortex.chord = 1
```

其局部远场速度和拒绝条件见《算法补充》12.11；无 PV 主分支省略这些键。

定常低 Mach SST：

```text
schema_version = 2
turbulence.model = k_omega_sst
turbulence.prandtl = 0.9
turbulence.wall_treatment = resolved
turbulence.freestream.intensity = 0.0003872983346207417
turbulence.freestream.length_scale = 0.000006928203230275509
turbulence.model_floor = 1e-12
turbulence.two_equation.source_treatment = local_implicit
time.integrator = lu_sgs
lu_sgs.jacobian = scalar_spectral
lu_sgs.sweeps = 1
preconditioner.type = weiss_smith
algorithm.riemann = roe
run.mode = steady
```

实验性标准 k-epsilon：

```text
schema_version = 2
turbulence.model = k_epsilon
turbulence.experimental = true
turbulence.prandtl = 0.9
turbulence.wall_treatment = wall_function
turbulence.freestream.intensity = 0.01
turbulence.freestream.length_scale = 0.1
turbulence.model_floor = 1e-12
turbulence.two_equation.source_treatment = local_implicit
turbulence.wall_function.y_plus_min = 30
turbulence.wall_function.y_plus_max = 300
time.integrator = lu_sgs
lu_sgs.jacobian = scalar_spectral
lu_sgs.sweeps = 1
preconditioner.type = none
run.mode = steady
```

三维非定常 WALE：

```text
schema_version = 2
turbulence.model = wale
les.filter.type = box3_tensor
les.filter.width_ratio = 1.0
les.wale.cw = 0.325
les.sgs_prandtl = 0.9
time.integrator = lu_sgs
time.physical.scheme = bdf2
time.dual_time.enabled = true
time.dual_time.max_iterations = 100
time.dual_time.l2_relative = 1e-8
preconditioner.type = none
run.mode = unsteady
```

## 4. 必须拒绝的草案

| 组合 | 原因 |
|---|---|
| 2D + 任一 LES model | LES 只支持三维 |
| steady + 任一 LES model | LES 必须物理时间推进 |
| unsteady + LU-SGS + dual_time=false | 未定义隐式物理层求解 |
| BDF2 无 BDF1 启动或合法两历史层 | 时间历史不完整 |
| Weiss--Smith + HLLC/Rusanov | v2 首版只支持 Roe |
| Weiss--Smith + SSPRK3 | 预处理伪时间系统不完整 |
| k-epsilon + resolved | 标准高 Re 模型壁面条件不适用 |
| k-epsilon 未显式给出 `turbulence.experimental=true` | 当前物理证据不足以开放默认配置 |
| 非 k-epsilon + `turbulence.experimental` | 实验声明不能污染其他模型身份 |
| `none` + 任意 RANS/LES 专属键 | 模型无关参数污染签名 |
| dynamic model 缺 test filter | Germano 恒等式不完整 |
| 非动态模型出现 `les.dynamic.*` | 无效/误导参数 |
| 输出统计累计伪迭代 | 非物理采样 |

## 5. 输出和重启

字段候选包括 wall distance、`nu_tilde/k/omega/epsilon`、`mu_model/mu`、SGS 应力六分量、
`Pi_sgs`、动态系数和模型源/耗散。checkpoint 保存模型 descriptor、BDF 历史层、物理/伪迭代、
动态平均状态及统计累加器。任何缺字段、变体/常数/滤波/积分器签名不同必须在推进前拒绝。

RANS 无滑移壁的 `output.boundary.quantities` 允许
`wall_distance,friction_velocity,wall_y_plus,wall_y_plus_class`。解析壁面由权威切向黏性牵引
重构 $u_\tau$ 和 $y^+$；壁面函数路径输出模型实际使用的 $y^+$。分类字段固定为 0（$y^+\le5$）、
1（$5<y^+<30$）、2（$30\le y^+\le300$）、3（$y^+>300$）。无活动 RANS 模型、无黏性运行
或非无滑移 patch 请求上述量必须拒绝；分类只供后处理，不允许触发模型自动切换。

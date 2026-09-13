# WCNS schema 2 配置草案（阶段 V 冻结输入）

状态：**设计草案，不是当前 schema 1 parser 已支持功能。** 最终键在阶段 W/X/Y/Z/AA 逐步
实现，每次实现必须同步合法/非法配置测试、summary、manifest 和 restart signature。

## 1. 兼容迁移

schema 1 未出现 v2 键时严格迁移为：

```text
turbulence.model = none
time.integrator = ssprk3
preconditioner.type = none
```

旧 `run.mode`、`run.cfl`、`run.t_end` 和 steady 停止键继续解释。schema 2 拒绝未知键、重复键、
非有限数、模型无关参数和不支持组合。

## 2. 键空间

| 键 | 候选值/类型 | 数值签名 | 约束 |
|---|---|---|---|
| `turbulence.model` | `none|sa_neg|k_omega_sst|k_epsilon|smagorinsky|scale_similarity|mixed_smagorinsky_similarity|dynamic_smagorinsky|wale` | 是 | LES 仅 3D unsteady |
| `turbulence.prandtl` | positive real | 是 | RANS 模型启用时 |
| `turbulence.wall_treatment` | `resolved|wall_function` | 是 | k-epsilon 必须 wall_function |
| `turbulence.farfield.*` | 模型专属正量 | 是 | 只允许当前模型所需键 |
| `les.filter.type` | `box3_tensor` | 是 | LES 必需 |
| `les.filter.width_ratio` | positive real | 是 | 基线 1 |
| `les.test_filter.ratio` | real > 1 | 是 | 动态/相似模型必需，基线 2 |
| `les.smagorinsky.cs` | `[0,0.3]` | 是 | Smag/混合，基线 0.17 |
| `les.similarity.cb` | finite real | 是 | 相似/混合，基线 1 |
| `les.dynamic.*` | average/clipping/floor | 是 | 仅动态模型 |
| `les.wale.cw` | positive real | 是 | WALE，基线 0.325 |
| `les.sgs_prandtl` | positive real | 是 | LES 必需 |
| `time.integrator` | `ssprk3|lu_sgs` | 是 | v2 隐式选 lu_sgs |
| `time.physical.scheme` | `bdf1|bdf2` | 是 | unsteady lu_sgs；生产基线 bdf2 |
| `time.dual_time.enabled` | bool | 是 | unsteady lu_sgs 必须 true |
| `time.dual_time.*` | iteration/tolerance/cfl | 是 | 内迭代停止与失败策略 |
| `lu_sgs.jacobian` | `scalar_spectral|block_source` | 是 | AA0 冻结 |
| `lu_sgs.sweeps` | positive integer | 是 | 一次内迭代的完整前后扫数 |
| `preconditioner.type` | `none|weiss_smith` | 是 | weiss_smith 只配 Roe+LU-SGS |
| `preconditioner.*` | Mach/viscous cutoff | 是 | 低 Mach 必需 |
| `statistics.time.*` | start/end/every/weight | 否；累加器身份单列 | 只累计接受物理步 |

## 3. 合法草案

定常低 Mach SST：

```text
schema_version = 2
turbulence.model = k_omega_sst
turbulence.prandtl = 0.9
turbulence.wall_treatment = resolved
time.integrator = lu_sgs
lu_sgs.jacobian = block_source
lu_sgs.sweeps = 1
preconditioner.type = weiss_smith
algorithm.riemann = roe
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
| `none` + 任意 RANS/LES 专属键 | 模型无关参数污染签名 |
| dynamic model 缺 test filter | Germano 恒等式不完整 |
| 非动态模型出现 `les.dynamic.*` | 无效/误导参数 |
| 输出统计累计伪迭代 | 非物理采样 |

## 5. 输出和重启

字段候选包括 wall distance、`nu_tilde/k/omega/epsilon`、`mu_model/mu`、SGS 应力六分量、
`Pi_sgs`、动态系数和模型源/耗散。checkpoint 保存模型 descriptor、BDF 历史层、物理/伪迭代、
动态平均状态及统计累加器。任何缺字段、变体/常数/滤波/积分器签名不同必须在推进前拒绝。

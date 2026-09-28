# WCNS v2.1 发布说明

v2.1 是在 v2.0_pre 内部前瞻版上的小版本更新，核心新增项是可与 SSPRK3 组合的显式
all-speed Roe 通量。版本仍按内部源码发布；不启用外部 CI，不改变当前许可证和对外发布状态。

## 新增

- 新增 `algorithm.riemann=roe_all_speed`。它按 Rieper 型低马赫修正缩放两支声学波强度中的
  法向速度跳跃贡献，中央物理通量、Roe 特征速度、熵修正和谱半径保持不变。
- schema 1 可直接用 SSPRK3；schema 2 使用
  `time.integrator=ssprk3 + preconditioner.type=none`。高 Mach 时逐项退化为普通 Roe。
- 回退链为 `roe_all_speed -> hllc -> rusanov`，诊断和 restart signature 保留真实请求算法名。
- 提供 case05 从 t=250 checkpoint 分叉的 SCMM6/线性 MDCD/`diss=0.001` 配置模板，可用
  `restart.mode=algorithm_change` 导入现有守恒态。
- 新增显式 all-speed Roe 数学、配置、限制和 Linux 操作文档。

## 兼容性

- 旧配置中的 `rusanov`、`hllc`、`roe` 和 Weiss--Smith 行为不变。
- v2.0_pre checkpoint 可用原算法做 `strict` 续算；切换到 `roe_all_speed` 必须使用
  `algorithm_change`。该模式保留 step/time 和守恒态，但重置时间统计、隐式历史和收敛状态。
- `run.t_end` 始终是绝对终止时间。例如 t=250 后再推进 50，应写 `run.t_end=300`。

## 验证范围

本机卡口覆盖：内建注册、均匀流一致性、低 Mach 修正确实生效、物理谱半径不变、高 Mach 与
普通 Roe 一致、法向反转、Sod/高 Mach 有限性、schema 2 SSPRK3 配置、小网格实际推进，以及
从旧算法 checkpoint 改算法导入守恒态。完整结果见 [`v2.1-validation.md`](v2.1-validation.md)。

没有执行大型 NACA0012、长期槽道统计或目标三维翼型；新增格式也不解除显式声学 CFL 限制。
这些边界见 [`known-limitations.md`](known-limitations.md)。

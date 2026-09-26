# WCNS schema 2 配置参考

本文件是 `v2.0_pre` 的稳定配置入口。完整逐键说明和示例见《用户手册》；运行时严格拒绝未知
键、重复键、非有限值、模型无关参数和不支持组合。

## 模型与推进组合

| 类别 | 可选值 | 约束 |
|---|---|---|
| `turbulence.model` | `none`, `sa_neg`, `k_omega_sst`, `k_epsilon`, `smagorinsky`, `scale_similarity`, `mixed_smagorinsky_similarity`, `dynamic_smagorinsky`, `wale` | k--epsilon 必须显式 experimental；LES 仅三维黏性非定常 |
| `time.integrator` | `ssprk3`, `lu_sgs` | 定常/非定常均可使用 LU-SGS；非定常 LU-SGS 必须双时间 |
| `time.physical.scheme` | `bdf2` | 无历史首步自动使用 BDF1 |
| `preconditioner.type` | `none`, `weiss_smith` | Weiss--Smith 只允许 Roe+LU-SGS |
| `algorithm.riemann` | `roe`, `rusanov`, `hllc` | Rusanov/HLLC 没有低马赫预处理版本 |

SA-neg 使用 resolved wall；SST-2003m 使用 resolved wall；标准高 Reynolds 数 k--epsilon 只允许
wall function 且当前保持实验支持等级。LES 的测试滤波比固定为 2；动态平均固定为
`local_box3_tensor`。Smagorinsky/混合模型的生产 wall damping 当前只允许 `none`；请求
`van_driest` 会明确失败。

## 统计

瞬时空间统计由下列配置选择：

```text
output.statistics.enabled = true
output.statistics.format = txt
output.statistics.every_steps = 1
output.statistics.quantities = total_mass,total_momentum_x,total_energy
```

接受步时间统计在同一选择上增加：

```text
statistics.time.enabled = true
statistics.time.start = 0
statistics.time.end = 1
statistics.time.every_steps = 1
statistics.time.weight = accepted_dt
```

输出包含每个量的 Reynolds/Favre mean、RMS，以及所有唯一量对的 Reynolds/Favre covariance。
仅接受的物理步按实际重叠时间加权；定常伪时间不应被解释为物理时间统计。

三维规则面可启用 x-z 平面、y-z 截面和槽道壁面统计。三维物理边界可配置显式右手
`drag/lift/span` 方向和严格递增的 `span_bin_edges`，各箱压力/黏性/总载荷必须回收到同一输出
事件的整体载荷。面心分箱不切割跨箱面，生产箱边界应与网格截面对齐。

## 关键拒绝项

- 二维、定常、无黏或 SSPRK3 的 LES；
- Weiss--Smith 与 Rusanov/HLLC 或 SSPRK3；
- k--epsilon 未声明 `turbulence.experimental=true`，或使用 resolved wall；
- LES 测试滤波比不等于 2，或模型专属参数出现在其他模型；
- 未注册统计量、重复量、越界截面、非共面规则截面；
- 与 checkpoint 中模型、滤波、积分器、统计身份或网格不一致的重启。

运行前始终执行：

```text
wcns_run --config <case.wcns> --dry-run
```

dry-run 只证明读取、配置、网格、分区及工作区初始化成功，不是物理收敛证据。

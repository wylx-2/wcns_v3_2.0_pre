# 解析谱随机相位 HIT：Samtaney IC4 初始条件

使用无显式湍流模型、三维周期八块结构网格。配置与网格分别为：

| 配置 | 网格 | 用途 |
|---|---|---|
| `smoke.wcns` | `grids/hit16.cgns` | 最多 3 步程序测试；谱峰未解析，不用于物理验证 |
| `coarse32.wcns` | `grids/hit32.cgns` | 服务器粗网格比较 |
| `production.wcns` | `grids/hit64.cgns` | 服务器正式计算 |
| `refine128.wcns` | `grids/hit128.cgns` | 网格敏感性 |

`hit.initialization=analytic_random_phase`，A=0.00013，k₀=8。采用确定性随机相位、共轭对称、无散投影及一次全局 K₀ 归一化；初始密度与温度为 1。谱峰使用无量纲角波数。此方法不能同时指定外部谱文件或 `hit.initial_energy`。

参考参数为 Samtaney 等（2001）D4 的初始 Mt=0.3、连续谱 Reλ=72；K₀=0.500523533878318，τ=0.5424170415173236。全部连续积分见 [initial_reference.json](initial_reference.json)。离散随机场的梯度和 Reλ 不会由总动能归一化同时精确固定。

此配置保留当前 HIT 常黏度设置；Samtaney 的 μ∝T^0.76、李新亮等的 Sutherland 关系未在本次初始化更新中实现。它是同类初场的常黏度 ILES 变体，不能作为两文献完整复现的声明。

在源码根目录执行短测或仅检查：

```sh
mpiexec -n 4 build/wcns_run --config cases/manual/case13_hit_analytic/smoke.wcns
mpiexec -n 4 build/wcns_run --config cases/manual/case13_hit_analytic/production.wcns --dry-run
```

正式配置计算到 7τ，仅在服务器去掉 `--dry-run` 后执行。默认输出目录相对启动目录；重新运行使用新目录。退出码 2 表示达到设定的步数/时间资源限制时应查看日志，不应视作已达到最终物理时间。

统计继续写出 HIT 历史、三维壳谱、一维谱、相关函数、应力、压强/密度涨落、梯度偏斜度和检查点。元数据记录解析 K₀ 与 τ；对文献画图使用 t/τ 和 K/K₀。`lambda_isotropic` 及 `Mt` 的演化后定义差别见 [初始化与预演化详细说明](../../../docs/hit-v2.6/HIT初始化与预演化说明.md)。

# 显式 all-speed Roe：算法、配置与 t=250 分叉

WCNS v2.1 的 `roe_all_speed` 是 Rieper 型低马赫 Roe 耗散修正。它仍离散原始守恒 Euler/
Navier--Stokes 方程，只改变界面 Roe 耗散中法向速度跳跃对两支声学波的贡献。因此它可以直接
用于显式 SSPRK3，也可以作为未启用 Weiss--Smith 时的普通空间通量用于 LU-SGS。

它和 `roe + preconditioner.type=weiss_smith` 不是同一种算法。后者改变伪时间系统、特征速度和
LU-SGS 更新，只允许 LU-SGS；`roe_all_speed` 不使用预处理矩阵，也不改变物理时间导数。

## 1. 离散定义

单位法向为 $\boldsymbol n$ 时，数值通量仍写成

$$
\widehat{\boldsymbol F}_n=
\frac12\left(\boldsymbol F_n(\boldsymbol Q_L)+\boldsymbol F_n(\boldsymbol Q_R)\right)
-\frac12\widetilde{\boldsymbol R}|\widetilde{\boldsymbol\Lambda}|
\boldsymbol\alpha^{AS}.
$$

波速保持普通 Roe 的
$\widetilde u_n-\widetilde a,\widetilde u_n,\widetilde u_n,\widetilde u_n,
\widetilde u_n+\widetilde a$，并继续使用相同的 Harten 熵修正。令

$$
\Delta u_n=(\boldsymbol u_R-\boldsymbol u_L)\cdot\boldsymbol n,
\qquad
z=\min\left(1,\frac{\|\widetilde{\boldsymbol u}\|}{\widetilde a}\right),
$$

其中速度和声速均为 Roe 平均。$z$ 使用三维欧氏速度模，是原始二维低马赫 Roe 修正的旋转
不变三维扩展。普通 Roe 的两支声学波强度为

$$
\alpha_- = \frac{\Delta p}{2\widetilde a^2}
-\frac{\widetilde\rho\Delta u_n}{2\widetilde a},
\qquad
\alpha_+ = \frac{\Delta p}{2\widetilde a^2}
+\frac{\widetilde\rho\Delta u_n}{2\widetilde a}.
$$

all-speed 版本只作

$$
\alpha_-^{AS} = \frac{\Delta p}{2\widetilde a^2}
-z\frac{\widetilde\rho\Delta u_n}{2\widetilde a},
\qquad
\alpha_+^{AS} = \frac{\Delta p}{2\widetilde a^2}
+z\frac{\widetilde\rho\Delta u_n}{2\widetilde a}.
$$

实现中等价地对普通特征投影增加

$$
\delta\alpha_-=(1-z)\frac{\widetilde\rho\Delta u_n}{2\widetilde a},
\qquad
\delta\alpha_+=-\delta\alpha_-.
$$

其余三支波强度、中央物理通量、熵修正和回退判据不变。$z=1$ 时逐项恢复普通 Roe；异常 Roe
平均仍按 `roe_all_speed -> hllc -> rusanov` 回退并记录真实请求算法名。

## 2. 能做什么、不能做什么

- 低 Mach 区域削弱 Roe 中随声速放大的速度跳跃耗散，避免把低速涡结构过度抹平。
- 物理 Euler 特征速度及返回给时间步控制器的谱半径仍为
  $|\widetilde u_n|+\widetilde a$，所以显式 SSPRK3 的声学 CFL 限制仍然存在。
- 该修正不是不可压缩投影、压力修正或伪时间预处理，不会让低 Mach 显式计算获得与对流速度
  成比例的时间步。
- 高 Mach 区域严格退化为普通 Roe；跨声速处仍保留普通 Roe 的熵修正。
- 本机小网格卡口证明代码路径可执行和离散性质符合上述公式，不替代槽道长期统计或外流精度
  验证。

算法依据为 F. Rieper, “A low-Mach number fix for Roe's approximate Riemann solver”,
Journal of Computational Physics 230 (2011) 5263--5287,
<https://doi.org/10.1016/j.jcp.2011.03.025>。WCNS 的三维局部 Mach 定义采用上述旋转不变扩展，
因此复现实验时必须记录程序版本和数值签名，不能与论文中其他局部 Mach 定义混写。

## 3. 配置方法

schema 1 的显式路径只需选择新通量；时间推进本来就是 SSPRK3：

```text
schema_version = 1
algorithm.riemann = roe_all_speed
```

schema 2 应明确写出：

```text
schema_version = 2
algorithm.riemann = roe_all_speed
time.integrator = ssprk3
preconditioner.type = none
```

不要同时设置 `preconditioner.type=weiss_smith`。Weiss--Smith 仍只允许
`algorithm.riemann=roe + time.integrator=lu_sgs`。

## 4. 从 case05 的 t=250 checkpoint 建立新分支

源码包中的 `examples/channel_retau180_from_t250_all_speed_roe.wcns` 冻结了 SCMM6、
线性 MDCD、`diss=0.001`、all-speed Roe 和 SSPRK3。使用前只替换模板末尾的 checkpoint
路径，并按服务器作业时限修改 `run.max_wall_time`。核心差异必须是：

```text
algorithm.profile = scmm6_wcns
algorithm.reconstruction = mdcd_linear
algorithm.reconstruction_variables = primitive
algorithm.riemann = roe_all_speed
algorithm.mdcd.disp = 0.0463783
algorithm.mdcd.diss = 0.001

run.cfl = 0.3
run.t_end = 300.0

restart.path = REPLACE_WITH_T250_CHECKPOINT.cgns
restart.mode = algorithm_change
```

`run.t_end` 是绝对物理终止时间，不是“从 checkpoint 再计算多久”；t=250 起算而希望再推进
50，应写 300。因为 Riemann 算法已改变，第一次必须使用 `algorithm_change`，它保留 checkpoint
中的 step/time 和守恒状态，但重置时间统计及算法历史。必须使用新的 `case.name` 和
`output.directory`，不能覆盖普通 Roe 分支。

`mesh.path` 必须指向生成该 checkpoint 时使用的**原始 CGNS 网格文件**。网格签名包含 base/
zone 身份和逐点坐标；重新生成的同名、同尺寸网格也可能签名不同。现有 t=250 checkpoint 的
签名为 `18099232003167909757`。v2.1 已用这份 checkpoint 及其原始网格完成真实 dry-run，恢复
到 `step=1047500, time=250`；若日志报 `checkpoint mesh signature differs`，应找回原始网格，
不得关闭校验或强制导入。

服务器顺序为：

```bash
install/serial/bin/wcns_run --config case05-all-speed-from-t250.wcns --dry-run
mpirun -np 4 install/mpi/bin/wcns_run --config case05-all-speed-from-t250.wcns --dry-run
mpirun -np 4 install/mpi/bin/wcns_run --config case05-all-speed-from-t250-smoke.wcns
```

短测配置应把 `run.max_steps` 暂时改为 2--5，并使用独立输出目录。核对启动摘要同时出现
`riemann_solver=roe_all_speed` 和 `time(integrator=ssprk3`，首个 checkpoint 的 time 仍从
约 250 延续，且没有非有限值或异常回退。短测通过后再恢复正式步数/墙钟。新分支自产生第一个
checkpoint 后，后续不再改算法的续段改用 `restart.mode=strict`。

## 5. v2.1 本机卡口

- 单元：均匀流一致性、低 Mach 通量确实不同于普通 Roe、物理谱半径不变、高 Mach 与普通 Roe
  逐分量一致、法向反转对称、Sod/高 Mach 状态有限。
- 配置：`roe_all_speed + ssprk3 + preconditioner.none` 可解析。
- 运行：小型二维低 Mach 均匀流用 SSPRK3 推进并保持有限。
- 重启：从普通 HLLC checkpoint 以 `algorithm_change` 切到 `roe_all_speed`，逐值核对导入守恒场。
- 实际种子：case05 的 64-rank t=250 Roe checkpoint 配原始网格 dry-run 成功，恢复
  `step=1047500, time=250` 和 mesh signature `18099232003167909757`。

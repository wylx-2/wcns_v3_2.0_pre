# v2.2 显式 all-speed Roe：Li--Gu 实现、兼容模式与 t=250 分叉

## 1. v2.1 审查结论

v2.1 的 `roe_all_speed` 实际采用 Rieper 型修正：普通 Roe 特征值、中央通量和 CFL 谱半径
不变，只把两支声学波强度中的法向速度跳跃乘以局部 Mach 因子。它与 `test_pdf_3.pdf`
第 5 章 5.3.2 节的 Li--Gu 全速度 Roe **不一致**。后者还要求：

1. 用 $\widetilde a'=f(M)\widetilde a$ 修改两支声学特征值；
2. 加入抑制低 Mach 压力棋盘振荡的压力跳跃项；
3. 显式给出参考 Mach 数 $M_{ref}$、Roe 耗散系数 $c_1$ 和压力系数 $c_2$。

v2.2 将 `roe_all_speed` 改为下面的 Li--Gu 实现；v2.1 算法以
`roe_all_speed_rieper` 保留。旧 checkpoint 若从原 `roe_all_speed` 切换到 v2.2 的同名算法，
必须使用 `restart.mode=algorithm_change`，不能把二者视为严格相同的离散。

## 2. Li--Gu 离散定义

单位面法向为 $\boldsymbol n$，Roe 平均法向速度、声速和密度分别为
$\widetilde u_n,\widetilde a,\widetilde\rho$。定义

$$
M=\max\left(\frac{|\widetilde u_n|}{\widetilde a},M_{ref}\right),
$$

以及

$$
f(M)=\min\left[
1,
M\sqrt{\frac{4+(1-M^2)^2}{1+M^2}}
\right],
\qquad \widetilde a'=f(M)\widetilde a.
$$

当 $M\ge 1$ 时直接取 $f=1$，避免不必要的大数中间量。Roe--Pike 左右特征向量保持不变，
五个特征值改为

$$
\widetilde u_n-\widetilde a',\quad
\widetilde u_n,\quad\widetilde u_n,\quad\widetilde u_n,\quad
\widetilde u_n+\widetilde a'.
$$

令 $\widehat U=(U_L+U_R)/2$、$\widehat p=p(\widehat U)$，并定义

$$
\widehat Q=\widehat U+(0,0,0,0,\widehat p)^T.
$$

v2.2 的最终法向通量严格按第 5 章式 (5.92) 实现：

$$
\widehat F_n=
\frac{F_n(U_L)+F_n(U_R)}{2}
-c_1\widetilde R|\widetilde\Lambda^{AS}|\widetilde L(U_R-U_L)
-[1-f(M)]\frac{c_2(p_R-p_L)}
{M_{ref}\widetilde\rho\widetilde a}\widehat Q.
$$

默认值采用该章标定值

$$
M_{ref}=0.1,\qquad c_1=0.02,\qquad c_2=0.05.
$$

Li--Gu 路径按文档中的 $|\widetilde\Lambda^{AS}|$ 直接取绝对值，不额外施加普通 Roe 的
Harten 熵修正。返回给时间推进器的谱半径为

$$
\rho_A=|\widetilde u_n|+\widetilde a'.
$$

因此它可与 SSPRK3 直接组合，而且低 Mach 声学刚性会减小；但稳定步长仍必须由程序返回的
谱半径和实际网格共同确定，不能手工按流速任意放大。非法 Roe 平均按
`roe_all_speed -> hllc -> rusanov` 回退。

## 3. 配置

schema 2 推荐显式写出：

```text
algorithm.riemann = roe_all_speed
algorithm.roe_all_speed.reference_mach = 0.1
algorithm.roe_all_speed.dissipation_scale = 0.02
algorithm.roe_all_speed.pressure_coefficient = 0.05

time.integrator = ssprk3
preconditioner.type = none
```

三个参数分别对应 $M_{ref},c_1,c_2$，并进入配置摘要和 restart signature。它们只允许在
`algorithm.riemann=roe_all_speed` 时出现。取 `c1=0.5` 且局部 $M\ge1$ 时，压力修正消失，
通量恢复未施加熵修正的标准 Roe 形式。报告给出的 `c1=0.02` 是 LES 标定值，不应未经验证
直接用于强激波生产计算。

如需复现 v2.1 的 Rieper 路径：

```text
algorithm.riemann = roe_all_speed_rieper
time.integrator = ssprk3
preconditioner.type = none
```

Rieper 兼容路径仍保持普通 Roe 的声学谱半径与 Harten 熵修正。Weiss--Smith 则是另一套
`algorithm.riemann=roe + time.integrator=lu_sgs` 伪时间预处理，不能与上述两个显式通量混用。

## 4. 从 case05 的 t=250 checkpoint 分叉

`examples/channel_retau180_from_t250_all_speed_roe.wcns` 已配置 SCMM6、线性 MDCD、
`diss=0.001`、Li--Gu all-speed Roe 和 SSPRK3。替换 checkpoint 路径后，第一次启动必须使用：

```text
run.t_end = 300.0
restart.path = REPLACE_WITH_T250_CHECKPOINT.cgns
restart.mode = algorithm_change
```

`run.t_end` 是绝对物理终止时间；从 250 再推进 50 应写 300。新分支产生首个 checkpoint 后，
参数不再变化的后续续段改用 `strict`。每个不同 Riemann/参数组合必须使用独立 `case.name` 和
`output.directory`。

checkpoint 必须与生成它的原始 CGNS 网格配套。现有 t=250 种子的网格签名为
`18099232003167909757`；同尺寸的重建网格也可能签名不同，不得关闭校验。服务器上先依次执行
串行 dry-run、MPI dry-run 和 2--5 步 smoke，再恢复正式墙钟和步数。

## 5. v2.2 小规模卡口

- 逐式检查 $f(M)$、修正谱半径和压力跳跃项；
- 均匀流一致性、法向反转、低 Mach 有限性及 `c1=0.5` 的高 Mach Roe 恢复；
- Rieper 兼容路径保持 v2.1 的谱半径和声学波强度行为；
- schema 2 参数解析、非法组合拒绝、SSPRK3 小网格运行及算法变更重启；
- 不把这些小规模结果解释为长期槽道统计或大型翼型物理验证。

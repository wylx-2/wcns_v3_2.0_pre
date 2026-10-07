# HIT 两种初始化与预演化准备阶段

2026-10-07，WCNS v2.6 初始化增量更新。已有 v2.6 发布包保留；本次修改位于当前源码树。

## 1. 新增行为和入口

新增 `hit.initialization` 参数，在以下两条路线之间切换。均采用三维周期均匀多块网格、无显式湍流模型及 MPI 分布傅里叶变换。

| 参数值 | 速度场及归一化 | 用途 |
|---|---|---|
| `shell_spectrum`（默认） | 保留 v2.6 的随机无散场和逐壳层配谱；外部实验谱不补偿未解析动能 | CBC / JHTDB；CBC 可进一步开启预演化与 Kang 重匹配 |
| `analytic_random_phase` | 解析谱、独立随机相位、无散投影、一次全局动能归一化；不逐壳拟合 | 所给 C++ 程序的初始化方法，Samtaney IC4 / 李新亮等的均匀热力学初场 |

主要源码为 `include/wcns/runtime/{case_config,hit}.hpp`、`src/runtime/{case_config,hit}.cpp` 和 `src/app/wcns_run.cpp`。主求解器、通量格式及时间积分算法不变。未显式设置新参数的旧算例继续采用原行为和原 strict 检查点签名。

配套文件：

- [CBC 64³ 正式配置](../../cases/manual/case11_hit_decay/prepared_production.wcns)、[CBC 16³ 短测](../../cases/manual/case11_hit_decay/prepared_smoke.wcns)。另有 `prepared_coarse32.wcns`、`prepared_refine128.wcns`。
- [解析谱 64³ 正式配置](../../cases/manual/case13_hit_analytic/production.wcns)、[解析谱 16³ 短测](../../cases/manual/case13_hit_analytic/smoke.wcns)。另有 32³、128³配置及各自八块 CGNS 网格。
- [解析初始积分参考值](../../cases/manual/case13_hit_analytic/initial_reference.json)。这些是连续谱公式计算值，不是演化后的 DNS 数据。
- [新增验收脚本](../../tests/test_hit_initialization.py)，用 NumPy DFT 独立核验导出的实场，且只做 16³ 短步计算。

## 2. CBC：预演化后再匹配初始谱

### 2.1 原作者流程与本程序时间

[Bae 的 CBC 初场页面](https://web.stanford.edu/~hjbae/CBC) 指定：先生成随机初场，以 QR 模型从实验时间标记 t′=0 演化到 42，再用 `kang_ic_64_qr.m` / `kangrescale.m` 将此时的速度场重新匹配到初始实验谱，然后开始测站间衰减。其配套 `notes_nondim.pdf` 给出时间尺度。本文遵循这一准备流程，预演化使用用户要求的 WCNS ILES；不启用 QR 显式模型。

定义 t′=t U₀/M，程序时间 t*=t Uref/Lref。CBC 使用 M=0.0508 m、U₀=10 m/s、Lref=11M=0.5588 m、Uref=0.2718933614489327 m/s，因此：

\[
t^*_{prep}=42\frac{M}{U_0}\frac{U_{ref}}{L_{ref}}
=0.10381382891686522.
\]

**42 不能直接填入程序时间参数。** 完成准备后，正式时间 0 对应实验测站 t′=42；正式采样时间依旧为 0、0.13841843855582028、0.31885676024465742，分别比较 42、98、171 测站。网页示例把继续计算的终点写成 191；本算例沿用已交付、具有实验数据的 171 测站。

```ini
hit.initialization = shell_spectrum
hit.spectrum_file = reference/cbc_spectrum_42_nondimensional.dat
hit.preparation.time = 0.10381382891686522
hit.preparation.max_steps = 10000000
```

`hit.preparation.time=0` 关闭准备阶段。新旧配置分开保存：原 `production.wcns` 仍是直接随机初场，新入口为 `prepared_production.wcns`。

### 2.2 逐模态重匹配

原作者的 `kangrescale.m` 对每个傅里叶模态的向量幅值进行缩放，并非把一个整数壳内所有模态统一缩放。程序采用正变换 DFT/N³、逆变换直接求和。令整数波矢为 m，Δk=2π/L，κ=|m|Δk，则重匹配目标为：

\[
\left|\widehat{\boldsymbol u}_{target}(\boldsymbol m)\right|^2
=\frac{2E_{ref}(\kappa)\Delta k}{4\pi |\boldsymbol m|^2},
\qquad
\widehat{\boldsymbol u}_{new}
=\sqrt{\frac{2E_{ref}(\kappa)\Delta k}
{4\pi |\boldsymbol m|^2|\widehat{\boldsymbol u}_{s}|^2}}
\widehat{\boldsymbol u}_{s}.
\]

重匹配中采用 `energy.m` 的对数插值及两端对数外推；这与原 v2.6 随机初场在最低测量波数以下采用 κ⁴ 外推的处理明确区分。非零目标对应零幅值模态时直接报错，不用新随机相位填补。零模态与截断范围外的模态清零。离散球壳模态数不等于连续球面面积，因此输出壳谱应与离散目标求和比较，不能要求每个壳的 E 数值恰等于壳中心连续曲线。

WCNS 是可压缩求解器，准备阶段会产生小量纵向速度；故先对演化场进行一次横向投影，随后乘正实数保持横向模态相位与方向。原作者不需要此投影，因为输入已为不可压缩无散场。本次保留已有随机初场生成方式；相位分布、ILES 演化、初始球形 N/3 截断和谱导数投影均是适配，不声称逐点重现 MATLAB 的 QR/交错网格结果。谱无散也不等于非线性 WCNS 离散连续方程残差严格为零。

重匹配仅改变速度；保留准备阶段的密度和内能，通过精确更新动能项维持各点压力和温度。此次人为动能改变不计入正式衰减预算，也不作为一次 RK 时间步。

### 2.3 输出、停止与重启

准备阶段写入 `output.directory/preparation/`，使用自己的时间和步数。准备历史中的 `statistics_weight` 始终为零。完成后清空 HIT 积分、采样历史、时间和步数，再写正式 t=0 状态。

- `*_hit_history.csv`：准备阶段及正式阶段各有一份；现有能谱、相关函数、梯度偏斜度、密度/压力涨落等输出继续提供。
- `preparation/*_rematch.csv`：重匹配前壳谱、重匹配后壳谱及逐模态求和的离散目标，可直接检查谱调整。
- `*_hit_metadata.txt`：标记 `phase`、初始化方法、是否完成准备、准备时长与已完成准备步数；解析谱路线还记录 A、k₀、K₀ 和 τ。
- 准备阶段达到步数上限、墙钟上限或收到停止信号时，不开始正式计算；按配置保存检查点。两阶段共享本次运行墙钟额度，步数上限分别设置。
- strict 检查点保存阶段标记。准备中重启继续准备；准备结束但尚未重匹配的检查点执行一次重匹配；正式检查点不重复预演化或重匹配。已验证 MPI 重分区重启。
- 准备路线禁用 `restart.mode=algorithm_change` 和通用场时间平均模块，避免清除阶段标记或混入准备统计；使用 HIT 自带统计。

新示例启用 `output.checkpoint.write_initial=true`，正式 step=0 即保留重匹配后的初始场。续算采用 `tools/prepare_hit_restart.py` 生成配置；可增大准备步数上限，但不得改初始化方法、种子、谱内容、准备时长和物理/数值参数后仍作 strict 续算。新续算使用新输出目录，阶段历史按各自时钟拼接。

## 3. 解析谱随机相位法：C++ 与两篇文献

### 3.1 初始化公式

Samtaney 等（2001），第 II.C 节式 (2.9)、(2.12) 定义：

\[
E(\kappa)=A\kappa^4\exp[-2\kappa^2/\kappa_0^2],\qquad
K_0=\frac{3A}{64}\sqrt{2\pi}\kappa_0^5.
\]

李新亮、傅德薰、马延文（2002），第 1.1 节也采用该谱、初始无散速度与均匀密度/压力/温度。所给 C++ 中 k₀=8、A=0.00013；它以三个独立随机相位构造复向量，再投影到垂直于波矢的平面，逆变换后整体归一到 K₀。本次方法保留这一数学流程。

令三个分量的 φd 由 `(seed, canonical(m), d)` 的确定性哈希产生，构造：

\[
b_d=\sqrt{\frac23\frac{A\kappa^2e^{-2\kappa^2/\kappa_0^2}\Delta k^3}{4\pi}}
(\sin\phi_d+i\cos\phi_d),\qquad
\widehat{\boldsymbol u}=\left(I-\frac{\boldsymbol m\boldsymbol m^T}{|\boldsymbol m|^2}\right)\boldsymbol b.
\]

全场乘同一系数 √(K₀/Kraw)。归一化依据谱 Parseval 动能，与逆变换后求物理空间平均动能等价；已经用独立 DFT 检查。原 C++ 使用相反的投影整体符号，不改变统计分布。

相对所给 C++ 的明确调整：

1. 使用成对波矢及共轭对称，保证逆变换虚部只含舍入误差。原片段独立生成 ±m 再仅取实部，未显式保证此条件；本实现不逐点复制原随机场。
2. 增加可配置随机种子。分区或 MPI 进程数改变不改变初场。
3. 沿用现有球形截断：`hit.cutoff=0` 表示整数波数 N/3，正值要求小于 N/2，避免 Nyquist 共轭歧义。原片段遍历整个笛卡尔波数立方体。
4. A、k₀ 均定义在程序无量纲**角波数** κ 上；L=2π 时 κ=|m|。其他盒长自动使用 Δk=2π/L。

解析谱路线有意采用原 C++ 的全局 K₀ 归一化，意味着被截掉的能量会通过统一倍数反映到保留模态；这与 CBC 实验谱不补能的原则不同。网格太粗或 k₀ 太大时不能据此宣称初始谱已充分解析。16³、k₀=8 的配置仅验证程序流程；正式使用至少比较 32³、64³、128³及不同种子。

### 3.2 热力学初场与参数

```ini
hit.initialization = analytic_random_phase
hit.spectrum_amplitude = 0.00013
hit.peak_wave = 8
hit.seed = 20261003
hit.preparation.time = 0
# 不设置 hit.spectrum_file；不设置 hit.initial_energy
```

此路线对应 **IC4**：ρ*=T*=1，p*=1/(γ Ma_ref²)，平均速度为零，总能量按内能加动能构造。没有实现 IC1–IC3 的压力泊松方程、等熵热力学涨落或非零初始散度。解析初场的启动瞬态本身属于这些文献的问题定义，因此不默认套用 CBC 的预演化重匹配。

| 新参数 | 默认值 | 约束 |
|---|---|---|
| `hit.initialization` | `shell_spectrum` | 只接受表中两种名称 |
| `hit.spectrum_amplitude` | 0.00013 | 有限且为正；解析路线使用 |
| `hit.preparation.time` | 0 | 非负；正数不小于 1e-12；仅外部目标谱的衰减路线可用 |
| `hit.preparation.max_steps` | 10000000 | 正整数；准备阶段独立停止限制 |

原有 `hit.peak_wave` 默认仍为 4；所给片段需要显式填 8。`hit.initial_energy` 仅保留给旧解析壳谱路线，新方法拒绝显式指定它，避免与 A 推导的 K₀ 冲突。解析方法可与常功率/常低波段能量强迫结合，但拒绝 `jhtdb_shells`，因为其强制覆盖前两个谱壳会破坏该初始化定义。

## 4. 对照初始量及比较边界

A=0.00013、k₀=8 的连续积分为：K₀=0.500523533878318，Ω₀=40.04188271026545，u′₀=0.5776524525342888，λ₀=2/k₀=0.25，L₀=√(2π)/k₀=0.31332853432887503，τ=0.5424170415173236。使用 t/τ 作为文献中的时间坐标。

新增正式示例指定 Mt₀=0.3、连续谱 Reλ₀=72（Samtaney 表 I 的 D4 初始参数），取 ν*=u′₀λ₀/72=0.0020057376824107248，Ma_ref=Mt₀/√(2K₀)=0.2998430630684572。有限网格/单个随机实现的梯度和 Reλ 会有偏差，必须以输出初始实测值同时报告；总 K₀ 与 Mt₀ 则由归一化与均匀声速确定。

当前 HIT 模块仍要求常黏度。Samtaney 使用 μ∝T^0.76，李新亮等使用 Sutherland 关系；因此本次示例是**相同类型初始条件下的常黏度 ILES 变体**，并非两篇论文的完整物理复现。新方法不自动把 `reference.mach` 等价为 Mt；必须按上式设置参考温度/速度。若后续要求逐曲线复现，还需匹配输运定律、分辨率、初始 Reλ/Mt、统计定义和随机样本。

已有 `lambda_isotropic` / `Re_lambda_isotropic` 使用耗散关系估算；论文 λ²=u′²/⟨(∂u₁/∂x₁)²⟩ 对可压缩、有限样本流场不必与之相等。现有 `Mt` 使用 √⟨T⟩ 推算声速，论文用平均局部声速；初始均匀温度时相同，演化后需注意差别。`K_mass` 可用于密度加权动能对比。高波数谱、梯度偏斜度及耗散必须考虑粗网格滤波和数值耗散。

CBC 准备时长来自原作者流程，但该时长是否足以使当前 ILES 相位关联成熟，需要在服务器观察偏斜度、谱传递及衰减趋势，并检查种子/准备时长敏感性。此次短测不构成成熟湍流或实验/DNS 精度验证。

## 5. 使用与并行实现

在源码根目录构建后，可仅做配置与初始化检查：

```sh
mpiexec -n 4 build/wcns_run --config cases/manual/case11_hit_decay/prepared_production.wcns --dry-run
mpiexec -n 4 build/wcns_run --config cases/manual/case13_hit_analytic/production.wcns --dry-run
```

去掉 `--dry-run` 才执行两阶段/正式求解，生产配置留待服务器运行。本机应使用两个 `*smoke.wcns`。相对网格/谱路径基于配置所在目录解析；`output.directory` 则相对启动程序的工作目录。默认不允许覆盖已有输出。

继续使用 v2.6 的 MPI slab 变换与块—slab 的 Alltoallv 重分布，每进程保存约 O(N³/P) 的场数据。FFTW 开启时执行本地一维变换，跨进程转置由本项目 MPI 层完成；这不是直接链接 `fftw_mpi_plan_dft_3d`，但整个三维初始化和重匹配均为分布式并行，不在根进程生成完整初场。`WCNS_ENABLE_FFTW=OFF` 可使用已有 radix-2 后备实现。

## 6. 参考材料

1. [Bae：CBC 初始条件代码及流程](https://web.stanford.edu/~hjbae/CBC)，核对 `kangrescale.m`、`energy.m`、`makefft.m`、`notes_nondim.pdf`，访问日期 2026-10-07。
2. Samtaney, R., Pullin, D. I., Kosović, B. (2001). *Direct numerical simulation of decaying compressible turbulence and shocklet statistics*. Physics of Fluids 13, 1415–1430。以用户提供 PDF 的 II.C、II.D 和表 I 为初始化依据。
3. 李新亮、傅德薰、马延文（2002）.《可压缩均匀各向同性湍流的直接数值模拟》，中国科学(A辑)，32(8)。以用户提供 PDF 第 1.1 节、第 2 节表 1 为依据。
4. 用户提供 `iso_wcns_2/hit_turbulence.cpp` 的 `generate_full_turbulence_parallel`，核对随机相位、投影及全局动能归一化；原文件未修改。

附：本次 [验证记录](初始化更新验证.md) 和参考来源 SHA-256 记录保存在同一目录。

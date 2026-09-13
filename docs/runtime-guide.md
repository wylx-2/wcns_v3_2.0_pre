# WCNS 运行、配置、输出与重启指南

本文以 WCNS `1.1.0` 已发布能力为主，并记录 v2.0.0 阶段 X 分支已实现、
但尚未通过完整物理验收的 `schema_version = 2` 能力；正式入口仍为
`wcns_run`，定位为简明速查。
逐步用户手册见 [`user-manual.md`](user-manual.md)，源码扩展指南见
[`developer-guide.md`](developer-guide.md)，完整配置模板见
[`examples/full_case_template.wcns`](../examples/full_case_template.wcns)，算法数学约定见
[`算法补充.md`](../算法补充.md)。配置文件采用严格的 UTF-8 `key = value` 格式：空行和以
`#` 开头的行被忽略，键不可重复；未知键、缺失必填键、非法枚举、`NaN/Inf` 和空列表项均在
分配流场前失败。

## 1. 构建与运行

v1.1.0 本机已验证环境包括 Windows 11、CMake 3.28.0、MinGW-w64 GCC 8.1.0、本机 Python
工具链和 Intel MPI 2021.10。CMake 最低声明版本为 3.20；未列出的编译器、操作系统和 MPI
组合目前属于未验证环境，而不是已知不兼容。仓库不启用外部 CI，合并前应在本机执行适用的
串行/MPI 构建、CTest 和算法规格检查。

串行 Release 构建（运行示例前需把目标结构网格放到 `examples/freestream.cgns`，或修改模板
中的 `mesh.path`）：

```powershell
cmake -S . -B build -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_MPI=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
cmake --install build --prefix build\install
build\wcns_run.exe --config examples\freestream.wcns --dry-run
build\wcns_run.exe --config examples\freestream.wcns
```

MPI 构建及运行：

```powershell
cmake -S . -B build-mpi -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-mpi --parallel 4
mpiexec -n 4 build-mpi\wcns_run.exe --config examples\freestream.wcns
```

安装后正式入口位于 `<prefix>\bin`，配置模板和文档位于 `<prefix>\share\wcns`。安装 smoke
可先用 `wcns_generate_release_cgns` 生成 ADF 网格，再复制并修改 `freestream.wcns.in` 的
占位符后执行 `wcns_run --config`。MinGW 二进制需要与构建器相容的
`libgcc_s_sjlj-1.dll`、`libstdc++-6.dll`；MPI 构建另需 Intel MPI 的 `impi.dll` 和启动器
在 `PATH` 中。安装规则不复制这些编译器/MPI 运行库。

`--dry-run` 完成配置广播、CGNS 元数据/叶块读取、分区、度量、初场或检查点恢复及所有启动
校验，但不进入时间推进和输出状态机。`mesh.path` 与 `restart.path` 的相对路径以配置文件所在
目录解析；`output.directory` 当前以进程工作目录解析，发布算例建议写绝对路径或从仓库根运行。

程序在 rank 0 输出完整配置/分区摘要、网格签名以及由参考量导出的 `Re` 和 `Ma`。配置中
禁止直接输入 `Re`、`Ma`、`reference.reynolds` 或 `reference.mach`。

## 2. 必填配置组

### 2.1 算例、网格和算法

| 键 | 可选值或含义 |
|---|---|
| `schema_version` | `1`，或 v2 开发分支的 `2` |
| `case.name` | 文件名前缀；不安全字符在输出名中替换为 `_` |
| `mesh.path` | 结构多块 CGNS 网格 |
| `algorithm.profile` | `phenglei_wcns` 或 `scmm6_wcns`；两套度量/算子独立使用 |
| `algorithm.reconstruction` | `weno_js`、`weno_z`、`mdcd_linear`、`mdcd_hybrid` |
| `algorithm.reconstruction_variables` | `conservative`、`primitive`、`characteristic` |
| `algorithm.riemann` | `rusanov`、`hllc`、`roe` |
| `robustness.enabled` | 默认 `false`；开启 SSPRK 候选物理容许性检查和局部逐面降阶 |
| `robustness.max_local_recomputations` | 每个 RK 阶段最大局部重算轮数，默认 `3` |
| `robustness.max_step_retries` | 最大整步缩步重试数，默认 `4` |
| `robustness.time_step_reduction` | 缩步因子，必须在 `(0,1)`，默认 `0.5` |
| `robustness.minimum_time_step` | 允许的最小时间步，默认 `1e-12` |

schema 2 必须显式给出时间积分和预处理器。保持层流路径时使用：

```text
turbulence.model = none
time.integrator = ssprk3
preconditioner.type = none
```

阶段 X 分支可选 SA-neg：

```text
turbulence.model = sa_neg
turbulence.prandtl = 0.9
turbulence.wall_treatment = resolved
turbulence.sa.farfield_nu_tilde_ratio = 3
turbulence.sa.source_treatment = explicit | local_implicit
time.integrator = ssprk3
preconditioner.type = none
```

`sa_neg` 要求 `run.viscous=true`、至少一个 resolved no-slip wall，远场比限定在
`[3,5]`；当前不允许与整步稳健化事务同时开启。`local_implicit` 仅对 SA 局部源项作
对角更新，它不是局部伪时间步或 LU-SGS。SST、k--epsilon、LES、`lu_sgs` 和
`weiss_smith` 仍会在启动前明确拒绝。schema 1 不接受这些 v2 键，其摘要和旧检查点
身份保持不变。SA-neg 当前是阶段 X 验证中功能，不应用于对外生产结论。完整迁移边界见
[`v2.0.0/config-schema-2-draft.md`](v2.0.0/config-schema-2-draft.md)。

低 Mach 预处理尚未实现。无粘界面通量只走所选
Riemann 求解器；发生非法中间状态时按冻结的确定性回退链处理并计数。稳健化开启后按
“原方案、同重构 primitive、linear5/primitive、zero_order/conservative+Rusanov”的有效去重
梯子升级真实残差直接支持及一层转置支持保护面；候选始终在独立缓冲区验证，失败整步不推进时间。详细数学定义见
[`../算法补充.md`](../算法补充.md) 11.2.3。

### 2.2 气体和参考量

`gas.gamma` 必填。`gas.molar_mass` 与 `gas.specific_gas_constant` 二选一，不能同时输入。
以下五个参考量均必填并必须为正：

- `reference.velocity`：`U_ref`；
- `reference.density`：`rho_ref`；
- `reference.temperature`：`T_ref`；
- `reference.length`：`L_ref`；
- `reference.viscosity`：`mu_ref`。

压力按 `rho_ref U_ref^2` 无量纲化。程序由这些输入和气体模型计算并仅在日志中报告
`Re = rho_ref U_ref L_ref / mu_ref` 与
`Ma = U_ref / sqrt(gamma R T_ref)`。

### 2.3 层流输运

不写任何 `transport.*` 键时保持 v1.0 默认：常黏度 `mu/mu_ref=1`、`Pr=0.72`。显式常黏度为：

```text
transport.model = constant
transport.prandtl = 0.72
transport.constant.viscosity_ratio = 1.0
```

Sutherland 模型必须显式给出参考温度黏度比，并以 K 或无量纲比值二选一给常数：

```text
transport.model = sutherland
transport.prandtl = 0.72
transport.sutherland.reference_viscosity_ratio = 1.0
transport.sutherland.temperature = 110.4
# 或 transport.sutherland.temperature_ratio = 0.38313378448724628
```

模型专属键不可交叉使用；两个温度键不可并存。程序启动时报告 `mu(T_ref)/mu_ref`、Pr、
`S/T_ref`、真实单元温度范围及对应黏度范围。规范化后的全部输运参数进入 restart signature，
等价的 K/比值表示产生相同签名。

### 2.4 分区

| 键 | 含义 |
|---|---|
| `partition.mode` | `zones_only` 不切分；`auto_split` 与当前 `force_split` 均进入受约束切分流程；`force_split` 暂未提供额外的强制切分语义 |
| `partition.allow_idle_ranks` | 不可形成足够合法叶块时是否允许空闲 rank |
| `partition.max_load_ratio` | 最大/平均负载目标，必须不小于 1 |
| `partition.min_cells_per_active_direction` | 每个叶块活动方向的最少单元数 |

当 CGNS zone 少于 rank 数时，`auto_split` 对结构 zone 作确定性受约束二分，同时切片物理边界、
原块连接并生成兄弟连接。若模板宽度或最少单元数使目标 rank 数不可行，程序在流场分配前给出
诊断，不静默生成过窄叶块。`phenglei_wcns` 的严格下限为 4，`scmm6_wcns` 为 5；实际配置可
选择更大的安全下限。

## 3. 初场、边界和源项

`initial.type` 支持：

- `uniform`：`rho,u,v,w`，以及 `temperature` 或 `pressure`；
- `sod_x`：`x0` 和 `left_/right_` 前缀的 `rho,u,v,p`；
- `quadrant_riemann`：`x0,y0` 和 `ne_/nw_/sw_/se_` 前缀的 `rho,u,v,p`；
- `isentropic_vortex`：`x0,y0,beta,background_u,background_v,period_x,period_y`；
  `period_x/period_y` 缺省或为零时使用普通距离，正值时分别使用该周期内的最短有符号距离。
  双向周期扭曲网格必须显式给出周期，避免超出坐标包围盒的周期映像破坏初场连续性；
- `couette`：`y0,y1,lower_velocity,upper_velocity,lower_temperature,`
  `upper_temperature,velocity_curvature,temperature_curvature,pressure`；速度和温度分别为
  端点线性插值加 `velocity_curvature*eta*(1-eta)` 与
  `temperature_curvature*eta*(1-eta)`，密度由常压状态方程得到；曲率项仅用于从非平衡初场
  验证定常收敛，边界端点值不变；
- `poiseuille`：`y0,y1,centerline_velocity,temperature,temperature_curvature,pressure`；
  速度为 `4*centerline_velocity*eta*(1-eta)`。温度为壁温 `temperature` 加
  `temperature_curvature*(eta/6-eta^2/2+2*eta^3/3-eta^4/3)`；常黏度 Poiseuille
  能量平衡可取该系数为
  `(gamma-1)*Ma^2*Pr*(4*centerline_velocity)^2`，零值则是等温低 Mach 解析初场；
- `linear_conduction`：`y0,y1,lower_temperature,upper_temperature,pressure`；速度为零、
  温度按 y 线性变化，密度由常压状态方程得到；
- `manufactured_periodic`：`beta,background_u,background_v`。

未显式给出的初场参数使用相应初始化器的内建默认值。生产初始化只写真实单元；随后由统一
halo/物理边界路径填充所需 ghost。物理边界 ghost 只保证边界条件产生的物理量及其守恒量
有效，不得读取 ghost 坐标、度量或二级派生量。

`boundary.default` 以及按 CGNS 边界名覆盖的 `boundary.<patch>.type` 支持 `farfield`、
`inflow`、`outflow`、`slip_wall`、`no_slip_adiabatic_wall`、
`no_slip_isothermal_wall`、`symmetry`、`periodic`。按同一 CGNS 边界名还可设置：

- 移动壁速度 `boundary.<patch>.wall_velocity_x/y/z`；
- 等温壁温度 `boundary.<patch>.wall_temperature`；
- 入口、远场或可选出口目标态 `boundary.<patch>.rho/u/v/w`，以及
  `boundary.<patch>.temperature` 或 `boundary.<patch>.pressure` 二选一。

一旦提供目标态中的任一字段，必须提供正的 `rho` 和恰好一个正的 `temperature|pressure`；
未给出的速度分量为零。目标态不能配置到壁面，壁温只能配置到等温无滑移壁，壁速度只能
配置到壁面。未提供分边界数据时，farfield/inflow 状态仍取初场在原点的解析值，等温壁温度
仍回退到 `initial.temperature`（再回退为 1）。运行时切分后的边界 patch 保留原 CGNS 名，
因此自动继承相同的分边界数据。边界数据进入配置摘要和数值重启签名。

源项由 `source.enabled` 控制。关闭时不得设置非零源项参数；开启时 `source.models` 可列出：

- `uniform_conservative`：五个 `source.uniform.*` 守恒源；
- `body_force`：`source.body.ax/ay/az`；
- `pressure_gradient`：`source.pressure_gradient.x/y/z` 保存常量体积力
  $\boldsymbol G=-\nabla p$，动量源为 $\boldsymbol G$、总能量源为
  $\boldsymbol u\cdot\boldsymbol G$，适合周期 Poiseuille 通道；
- `manufactured`：五个 `source.manufactured.*` 幅值。

多个模型在每个 SSPRK 子步同址求值并相加。二维算例的 z 动量源必须为零。

## 4. 运行和停止

公共必填键为 `run.mode = steady | unsteady`、`run.viscous`、`run.cfl` 和
`run.max_steps`。`run.max_steps` 始终是硬上限；`run.max_wall_time = 0` 表示关闭墙钟上限，
正值要求启用检查点。

定常模式中的 `time` 是伪时间，使用下列残差键：

- `steady.min_steps`、`steady.check_interval_steps`、`steady.consecutive_checks`；
- `steady.reference_floor`；
- `steady.l2_absolute`、`steady.l2_relative`；
- `steady.linf_enabled`、`steady.linf_absolute`、`steady.linf_relative`。

五个守恒分量分别使用真实单元、正 Jacobian 和 profile 守恒权重形成 MPI 全局 `L2/Linf`。
首次检查冻结参考残差；所有分量连续满足规定次数后返回 `steady_converged`。

非定常模式必须给出正的 `run.t_end`。残差只监测、不触发收敛；时间步会裁剪到 `t_end` 和
下一个精确时间输出事件，因此正常结束为 `physical_time_reached`，不会越过目标时间。

停止优先级为数值失败、模式正常完成、用户信号/墙钟安全停止、最大步数。退出码为：

| 退出码 | 含义 |
|---:|---|
| 0 | `steady_converged` 或 `physical_time_reached` |
| 1 | 启动、配置、CGNS、重启签名或 I/O 异常 |
| 2 | `maximum_steps`、`wall_time_checkpoint`、`user_signal_checkpoint` |
| 3 | `numerical_failure` |

`SIGINT/SIGTERM` 和墙钟上限只在完整步边界生效。安全停止写最终安全检查点；数值失败保留历史
与 manifest，但不把已知非法状态写成正常流场或可重启检查点。

## 5. 输出配置

`output.directory`、`output.allow_existing` 和 `output.dimensional` 必填。默认推荐
`allow_existing = false`，目录已存在即失败；重启追加或覆盖同名 step 文件必须显式设为
`true`。各输出类别使用以下通用调度键：

```text
every_steps = 0
every_time = 0
explicit_times = 0.1,0.25,0.5
write_initial = false
write_final = true
```

步、时间、显式、初场和最终事件取并集，同一 `(step,time)` 只写一次。`every_time` 与
`explicit_times` 在非定常模式下通过裁剪时间步精确命中。

### 5.1 流场

`output.field.format = cgns | tecplot | both`。CGNS 按输入的原 zone 重组为 `CellCenter`
FlowSolution；Tecplot ASCII 按原 zone 写 cell-center ordered zone。支持字段：

| 名称 | 含义 |
|---|---|
| `rho,u,v,w,p,T` | 原始物理量 |
| `rho_u,rho_v,rho_w,rho_E` | 守恒量 |
| `sound_speed,mach,total_enthalpy,entropy_proxy` | 热力学派生量 |
| `viscosity` | 当前层流输运模型黏度 |
| `nu_tilde` | SA-neg 工作变量；仅模型活动时可用 |
| `mu_t_over_mu` | SA 涡黏度/分子黏度 |
| `sa_production,sa_destruction` | SA 产生项和破坏项诊断 |
| `wall_distance,sa_negative_branch` | 壁距离和 SA-neg 负分支标识 |
| `jacobian` | `partial(x,y,z)/partial(xi,eta,zeta)`；二维为面积尺度、三维为体积尺度 |

`output.dimensional = true` 时按参考量恢复量纲；否则输出内部无量纲值。需要 ghost 坐标/度量
或梯度的涡量、散度和 Q 判据尚未注册；壁面热流应通过 `output.boundary.quantities=q_wall`
请求。若在 `output.field.quantities` 请求未知量，会在首次输出前明确报
“unknown field quantity”，不会输出占位零值。

### 5.2 边界面工程量和载荷

`output.boundary.format = txt | tecplot`。启用时必须给出物理 patch、逐面量、来流参考状态、
参考面积/长度、力矩中心及单位 drag/lift/tangent 方向。内建量为 `p_w,T_w,mu_w,Cp,Cf,q_wall`、
压力/黏性/总牵引的 `x/y/z` 分量；无粘运行不能请求黏性量。

程序在已接受状态上复用当前 profile 的压力迹和黏性壁面强迹，使用全局守恒边界权重
`w_b|S|`。所有 rank 将固定宽度面记录汇总到 root，并按
`(patch,source_zone,k,j,i,axis,side)` 排序查重。每个事件生成
`<case>.boundary.r<ranks>.step....txt|dat`，载荷固定 schema 写入
`<case>.loads.r<ranks>.txt`。二维头记录 `force_per_unit_span=true`。这些配置进入摘要、digest
和 manifest，但不改变数值 restart signature。完整公式和符号见《算法补充》11.4。

### 5.3 残差历史、统计和 manifest

`output.history.format = txt | tecplot`。历史列固定，包含 step/time/dt/CFL/wall time、总残差、
五个平均流分量及活动模型场的 `L2/Linf`、冻结参考值、归一化值、连续通过次数、
重构/Riemann 回退、稳健化各级
owner 面数、troubled-cell/局部重算/整步 retry、proposed/accepted dt、候选最小 `rho/p/T/e`、
是否进行残差检查和停止原因。固定 schema 不接受 `output.history.quantities`；TXT 直接写停止原因字符串，
Tecplot 写数值 `stop_reason_code` 并在 `AUXDATA STOP_REASON_CODES` 中给出映射。

`output.statistics.format = txt | tecplot`。当前内建可选量为 `total_mass`、
`total_momentum_x/y/z`、`total_energy`，采用与残差一致的原 zone 守恒积分权重。

每次运行的 `<case>.manifest.r<ranks>.txt` 记录版本、Git 提交、编译器、构建类型、MPI 数、配置/分区
摘要、网格和重启签名、最终 step/time/dt/wall time、停止原因及成功提交的输出文件。流场、
历史、统计、检查点和 manifest 均先写同目录临时文件，成功关闭后再改名。

## 6. 检查点和不同 rank 重启

检查点只能是项目内部 CGNS 布局，始终无损保存原 zone 上五个无量纲守恒场；
SA-neg 活动时使用 checkpoint schema 2 并额外保存 `NuTilde`、模型描述符、模型残差参考和
连续计数，`none` 仍保持 schema 1。所有检查点另含格式版本、step/time/dt、网格签名和
数值重启签名。除带 step/time 的文件外，
还更新 `<case>.checkpoint.latest.cgns`。`latest` 是本次运行拥有的滚动别名，因此即使
`output.allow_existing=false`，后续检查点事件也会安全替换它；运行前既有输出目录以及带
step/time 的不可变文件仍按覆盖策略严格拒绝。

重启配置加入：

```text
restart.path = output/run-a/case.checkpoint.latest.cgns
```

程序先按原 zone 重读，再按当前 `StructuredPartitionPlan` 分发，因此可由 1 rank 检查点以
2/4/8 rank 和新的合法叶块划分续算。profile、重构/Riemann、气体、参考量、边界、源项、
黏性开关、网格或格式不兼容会在第一次推进前由所有 rank 一致失败；输出目录、输出调度、
rank 数和 `t_end` 可以改变。Tecplot 流场不是检查点，不能用于重启。

## 7. 自定义字段和统计量

自定义字段实现 `IFieldQuantity`，自定义统计量实现 `IStatisticQuantity`。两者的
`QuantityDescriptor` 必须给出唯一名称、`Cell` 位置、依赖、无量纲/有量纲单位和缩放类型；
字段求值只遍历真实单元。registry 拒绝空名、重名、未知依赖、循环依赖和非有限结果。

生产对象支持注入扩展后的 registry：

```cpp
auto fields = wcns::FieldQuantityRegistry::create_builtin();
fields.register_quantity(std::make_shared<MyFieldQuantity>());
wcns::ProductionFieldWriter writer(
    mpi, config, plan, local_blocks, metrics, quantity_context,
    mesh_path, std::move(fields));

auto statistics = wcns::StatisticRegistry::create_builtin();
statistics.register_quantity(std::make_shared<MyStatisticQuantity>());
wcns::RuntimeOutputManager output(
    mpi, config, plan, mesh_signature, &statistic_context,
    event_writer, std::move(statistics));
```

随后只需把新名称写入 `output.field.quantities` 或 `output.statistics.quantities`。标准
`wcns_run` 未加载动态插件；要让标准可执行程序识别项目自定义 C++ 量，应在其启动装配处创建
并注入 registry，无需修改 `SimulationDriver`、输出调度或求解器推进循环。

## 8. 常见故障诊断

- `unknown configuration key`：拼写错误或使用了当前 schema 未实现的键；不要绕过严格 parser。
- `partition ... infeasible`：rank 太多、块太窄或 `min_cells_per_active_direction` 太大；降低 rank、
  允许 idle rank，或使用满足模板要求的网格。
- `output directory already exists`：使用新目录；只有明确要续写时才设 `allow_existing=true`。
- `checkpoint ... signature differs`：网格或影响数值轨迹的配置改变；不能强行忽略。
- `unknown field/statistic quantity`：名称不在当前 registry；检查拼写或注入自定义实现。
- 退出码 2：算例没有正常完成，检查 manifest 的 `stop_reason`；它不是数值崩溃。
- 退出码 3：发生数值失败；检查最后历史记录、回退计数、初场正性、边界条件和 CFL。

完整的实现边界和发布限制集中列在 [`known-limitations.md`](known-limitations.md)。

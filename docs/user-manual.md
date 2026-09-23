# WCNS 用户手册

本文面向第一次接触本程序的算例使用者，以 WCNS `1.1.0` 已发布功能为主，并记录 v2.0.0
阶段 Z 分支已经实现、但仍受候选/服务器验证边界约束的 `schema_version = 2` 能力；生产入口仍为
`wcns_run`。按本文顺序操作，可以从源码构建程序、准备
CGNS 网格、填写配置、完成串行或 MPI 计算、识别停止状态、读取输出并从检查点续算。

本手册描述的是当前程序已经实现的行为。数学定义见[`算法补充.md`](../算法补充.md)，源码扩展见[`developer-guide.md`](developer-guide.md)，实现边界见[`known-limitations.md`](known-limitations.md)。可复制的完整配置见[`examples/full_case_template.wcns`](../examples/full_case_template.wcns)。

## 1. 开始前必须知道的约定

1. 程序只读取二维或三维**结构多块** CGNS 网格；当前随附 CGNS 4.4.0 只启用 ADF 后端。HDF5-CGNS、非结构网格、重叠网格、滑移接口、运动网格和 AMR 不在当前范围内。
2. 输入网格坐标、初场、边界数据、源项、计算时间全部使用程序的**内部无量纲量**。只有 `gas.*` 和 `reference.*` 是用于定义量纲的有量纲参考输入。
3. 二维状态仍保存五个 Euler 分量 `(rho,rho*u,rho*v,rho*w,rho*E)`，但二维的 `w` 及 z 动量必须为零。
4. `algorithm.profile = phenglei_wcns` 与 `scmm6_wcns` 是两套独立的度量、线性插值、通量导数和物理边界闭合组合；不能从两套 profile 中交叉抽取部件。
5. schema 1 时间推进固定为显式 SSPRK3；schema 2 还支持定常 LU-SGS 以及非定常 BDF2/BDF1
   双时间 LU-SGS。`run.max_steps` 对两者始终是硬上限，未收敛内迭代不能接受物理时间层。
6. 配置解析是严格的：键区分大小写，未知键、重复键、缺失必填键、空值、`NaN/Inf`、非法枚举和逗号列表空项都会在计算前失败。
7. 建议每次计算使用新的输出目录，并保留最终 manifest。退出码为 0 才表示正常达到定常收敛或非定常目标时间；最大步数、墙钟和信号停止返回 2，不能当作“计算成功收敛”。

## 2. 源码目录与可执行程序

在仓库的 `wcns` 目录运行后续命令。主要目录如下：

```text
wcns/
|-- include/wcns/          C++ 公共头文件
|-- src/                   程序实现和正式入口
|-- tests/                 单元、MPI、CGNS 与端到端测试
|-- examples/              可复制的普通配置
|-- cases/config/          发布矩阵使用的占位符模板
|-- cases/manual/          已执行的人工验收算例及结果
|-- tools/                 网格生成、结果验证和矩阵驱动
|-- third_party/cgns/      随仓库提供的 CGNS 4.4.0 源码归档
|-- docs/                  设计、验收和使用文档
|-- CMakeLists.txt         顶层构建入口
`-- 算法补充.md            数学及离散算法约定
```

启用 CGNS 后会生成五个面向用户的程序：

| 程序 | 用途 |
|---|---|
| `wcns_run` | 正式求解器；读取一个 `.wcns` 配置并运行 |
| `wcns_generate_release_cgns` | 生成均匀、扭曲、加密、周期和通道测试网格 |
| `wcns_validate_release_case` | 独立重读输出 CGNS，检查有限性、解析误差或两场差异 |
| `wcns_compare_metric_profiles` | 对同一网格独立计算两套 profile 的度量并报告差异 |
| `wcns_inspect_structured_mesh` | 只读报告结构网格拓扑、二维单元质量及指定壁面几何 |

`wcns_unit_tests` 等测试程序只在 `WCNS_BUILD_TESTS=ON` 时构建，不是生产求解入口。

## 3. 准备编译环境

### 3.1 已验证的 Windows 环境

v1.0 实测组合包括：Windows 10、CMake 3.28、MinGW-w64 GCC 8.1、Python 3.14.5 和
Intel MPI 2021.10。CMake 最低声明版本是 3.20。其他编译器或系统不一定不能使用，但属于
未进入发布矩阵的组合，首次使用必须完整运行测试。

确认工具可见：

```powershell
cmake --version
g++ --version
python --version
git --version
```

MPI 构建还应确认：

```powershell
mpiexec -help
Get-ChildItem Env:I_MPI_ROOT
```

若 `I_MPI_ROOT` 不存在，CMake 还会尝试
`C:\Program Files (x86)\Intel\oneAPI\mpi\latest`。MinGW 生成的可执行程序运行时需要相容的
`libgcc_s_sjlj-1.dll`、`libstdc++-6.dll`；MPI 版本还需要 `impi.dll` 和 Intel MPI 环境。

### 3.2 Linux 或其他工具链

非 Windows/MinGW 的 MPI 构建使用 CMake 的 `find_package(MPI COMPONENTS CXX)`。可使用
系统 GCC/Clang 和 OpenMPI/MPICH，但当前没有冻结的 Linux 发布基线。建议先建串行 Release，
再建 MPI Release，并分别执行完整 CTest。

### 3.3 不需要额外安装 CGNS/HDF5

`third_party/cgns/CGNS-4.4.0.zip` 已在仓库中。配置阶段会从本地归档构建静态 CGNS，
不访问网络。当前强制 `CGNS_ENABLE_HDF5=OFF`，所以输入文件必须能由 ADF 后端打开。

## 4. 编译、测试和安装

### 4.1 串行 Release 构建

在 `wcns` 根目录执行：

```powershell
cmake -S . -B build-user-serial -G "MinGW Makefiles" `
  -DWCNS_ENABLE_CGNS=ON `
  -DWCNS_ENABLE_MPI=OFF `
  -DWCNS_BUILD_TESTS=ON `
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-user-serial --parallel 4
ctest --test-dir build-user-serial --output-on-failure
```

各参数含义：

- `-S .` 指定源码根目录；`-B build-user-serial` 把生成文件隔离在单独目录。
- `WCNS_ENABLE_CGNS=ON` 才会构建正式求解器和 CGNS 工具。
- `WCNS_ENABLE_MPI=OFF` 生成串行后端；仍使用同一并行抽象，rank 数固定为 1。
- `WCNS_BUILD_TESTS=ON` 构建并注册测试；只想制作最小安装包时可设为 `OFF`。
- `CMAKE_BUILD_TYPE=Release` 启用优化。不要用 Debug 性能评价生产算例。
- `--parallel 4` 是编译并发数，可按机器调整，与求解 MPI rank 数无关。
- `ctest --output-on-failure` 只在失败时展开对应测试输出。

若之前用不同编译器或 MPI 选项配置过同一个构建目录，不要在其中混改缓存；使用新的`-B` 目录最安全。

### 4.2 MPI Release 构建

```powershell
cmake -S . -B build-user-mpi -G "MinGW Makefiles" `
  -DWCNS_ENABLE_CGNS=ON `
  -DWCNS_ENABLE_MPI=ON `
  -DWCNS_BUILD_TESTS=ON `
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-user-mpi --parallel 4
ctest --test-dir build-user-mpi --output-on-failure
```

Windows/MinGW 下该配置要求 Intel MPI 的 `mpi.h`、`impi.lib` 和 `mpiexec`。如果 CMake 报
`WCNS_ENABLE_MPI requires an Intel MPI installation`，先修复 MPI 安装或环境变量，不要把
MPI 头文件和库从不同版本混合链接。

Linux 下可由 CMake 查找系统 MPI，命令示例为：

```bash
cmake -S . -B build-user-linux-mpi \
  -DWCNS_ENABLE_CGNS=ON \
  -DWCNS_ENABLE_MPI=ON \
  -DWCNS_BUILD_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-user-linux-mpi --parallel 4
ctest --test-dir build-user-linux-mpi --output-on-failure
mpiexec -n 2 build-user-linux-mpi/wcns_run --config path/to/case.wcns --dry-run
```

这是未进入当前冻结发布矩阵的环境。若 CMake 找到多个 MPI，实现的头文件、链接库和`mpiexec` 必须来自同一套安装。

### 4.3 安装到独立目录

```powershell
cmake --install build-user-mpi --prefix build-user-install
Get-ChildItem build-user-install\bin
Get-ChildItem build-user-install\share\wcns
```

安装树中：

- `bin` 保存 WCNS 四个正式程序以及上游 CGNS 安装规则带入的工具；
- `share/wcns` 保存 README、算法文档、用户/开发手册、示例和配置模板；
- 工具链与 MPI 的运行 DLL 不会被复制，部署到另一台机器时仍需提供相容运行环境。

### 4.4 最小构建验收

至少完成以下三项后再投入长时间计算：

```powershell
ctest --test-dir build-user-serial --output-on-failure
ctest --test-dir build-user-mpi --output-on-failure
build-user-mpi\wcns_generate_release_cgns.exe periodic-square `
  build-user-mpi\manual-smoke.cgns 16 16 1.0
```

已有构建中的测试数可能随开发增加，不要把某个固定数量硬编码成成功条件；以 CTest 零失败和命令退出码 0 为准。

## 5. 十分钟完成第一个计算

以下步骤生成一个 32×16 的二维结构网格，复制配置并做串行短计算。

### 步骤 1：建立不与源码混杂的算例目录

```powershell
New-Item -ItemType Directory -Force work\quickstart
```

### 步骤 2：生成 ADF-CGNS 网格

```powershell
build-user-serial\wcns_generate_release_cgns.exe rectangle `
  work\quickstart\mesh.cgns 32 16 1 1.0 1.0 false
```

这里的参数依次是输出文件、x/y 总单元数、x 向原生 zone 数、物理长度 `Lx/Ly`、是否 x
周期。生成的非周期物理边界名为 `left`、`right`、`bottom`、`top`。

### 步骤 3：复制完整模板

```powershell
Copy-Item examples\full_case_template.wcns work\quickstart\quickstart.wcns
```

编辑下列值：

```text
case.name = quickstart
mesh.path = mesh.cgns
run.max_steps = 10
run.t_end = 0.01
output.directory = work/quickstart/output
output.checkpoint.enabled = false
```

由于 `mesh.path` 相对配置文件解析，网格与配置在同一目录时只写 `mesh.cgns`。由于从仓库根启动，`output.directory` 写成 `work/quickstart/output`。

### 步骤 4：只做启动检查

```powershell
build-user-serial\wcns_run.exe --config work\quickstart\quickstart.wcns --dry-run
if ($LASTEXITCODE -ne 0) { throw "WCNS dry-run failed" }
```

`--dry-run` 会实际完成配置解析、CGNS 读取、分区、块连接构建、度量计算、初场或检查点恢复和启动校验，但不推进、不创建输出目录。rank 0 应打印配置摘要、分区摘要、网格签名以及由参考量导出的 `Re`、`Ma`，最后打印 `WCNS dry-run completed`。

### 步骤 5：正式运行

```powershell
build-user-serial\wcns_run.exe --config work\quickstart\quickstart.wcns
$code = $LASTEXITCODE
Write-Host "wcns exit code = $code"
```

正常非定常结束应打印 `reason=physical_time_reached`，退出码为 0。若输出目录已存在且模板中`output.allow_existing=false`，程序会拒绝运行；改用新目录最安全。

### 步骤 6：检查结果

```powershell
Get-ChildItem work\quickstart\output
Get-Content work\quickstart\output\quickstart.manifest.r1.txt
Get-Content work\quickstart\output\quickstart.history.r1.txt | Select-Object -Last 3
```

还可独立检查最终 CGNS：

```powershell
$final = Get-ChildItem work\quickstart\output\quickstart.field.*.cgns |
  Sort-Object Name | Select-Object -Last 1
build-user-serial\wcns_validate_release_case.exe finite $final.FullName
```

## 6. CGNS 网格输入要求

### 6.1 基本结构

一个可用网格必须满足：

- 至少一个 `CGNSBase_t`，`CellDimension` 为 2 或 3，且物理维数不小于单元维数；
- zone 类型为 `Structured`；每个活动方向至少有两个顶点；
- 必须有 `CoordinateX`、`CoordinateY`；三维物理空间还必须有 `CoordinateZ`；
- 坐标为有限值，网格不能产生非正或退化 Jacobian 体积/面积；
- 物理边界使用顶点位置的 `BC_t/PointRange`，且范围必须描述完整块面，不可只描述边或角；
- 块间共形连接使用 `GridConnectivity1to1_t`，接收与供体范围、轴变换和层数必须一致；
- 周期连接也用 1-to-1 connectivity，并通过 CGNS 周期平移/旋转数据表达；二维旋转周期当前因 CGNS 表意歧义被拒绝。

CGNS 索引是 1-based；读入后程序转换为 0-based 半开区间并检查边界。用户不需要在配置中写索引范围，只需用准确的 `BC_t` 名称覆盖边界类型和物理数据。

### 6.2 网格坐标必须预先无量纲化

若物理坐标为 `x_dim`，写入求解网格的值应为 `x=x_dim/L_ref`。例如 2 m 长通道选`reference.length=2` 时，网格 x 范围应写 0 到 1；若网格仍写 0 到 2，程序会把它理解为 2 个参考长度。`output.dimensional=true` 时输出坐标再乘 `L_ref`。

### 6.3 物理边界名与连接的职责

- CGNS `BC_t` 的**名称和范围**决定配置如何找到边界；正式入口会用`boundary.default`/`boundary.<name>.type` 覆盖读入的 CGNS BC 类型。
- CGNS 1-to-1 connectivity 决定块间或周期 halo 通信，不能用两个物理 BC 假装成块连接。
- 配置中虽然能解析 `periodic` 边界字符串，但生产物理边界填充不接受“周期 BC”；周期必须已在 CGNS connectivity 中给出。不要写 `boundary.left.type=periodic` 来替代连接。
- 自动二次剖分产生的人造切面由程序建立兄弟连接，不需要写回 CGNS。

### 6.4 用随附工具生成网格

通用模式：

```text
wcns_generate_release_cgns output.cgns dimension cells_i cells_j cells_k zones_i warp periodic_x
```

二维必须 `dimension=2,cells_k=1`；`zones_i` 是 x 向原生 zone 数，`cells_i` 必须可整除；`periodic_x` 是 `true|false`。其他专用模式：

```text
periodic-square output.cgns cells_i cells_j length
warped-periodic-square output.cgns cells_i cells_j length Ax Ay
rectangle output.cgns cells_i cells_j zones_i Lx Ly periodic_x
clustered-rectangle output.cgns cells_i cells_j zones_i Lx Ly cluster_x cluster_y strength periodic_x
periodic-channel output.cgns cells_i cells_j cells_k zones_i zones_k Lx Ly Lz wall_cluster_strength [origin_y]
```

- `periodic-square`/`warped-periodic-square` 生成 2×2 原生多块、x/y 双周期网格。
- 扭曲方形采用连续解析映射；幅值必须满足生成器的正 Jacobian 条件。
- `clustered-rectangle` 向给定 x/y 内点光滑加密；`strength` 越大目标点附近越密。
- `periodic-channel` 生成 x/z 周期、y 向有上下物理壁的三维通道；壁面加密强度 0 为均匀网格。可选 `origin_y` 是 y 区间起点，默认 0；例如 `Ly=2,origin_y=-1` 产生 `[-1,1]`。
- `invalid-one-sided` 专用于失败测试，故意生成单向连接，严禁物理解算。

生成器覆盖测试网格，不是通用网格转换器。外部网格工具产生的 CGNS 仍须满足本节契约，并首先通过 `wcns_run --dry-run`。

## 7. 配置文件语法与必填清单

### 7.1 语法

```text
# 这是独占一行的注释
key = value
```

规则如下：

- 每个有效行必须恰好有一个 `=`；
- 键和值两侧空白会去掉；值中没有引号、变量替换、算式或单位换算；
- 注释只在去除行首空白后第一个字符是 `#` 时生效；**不支持行尾注释**；
- 布尔值只能小写 `true` 或 `false`；枚举也区分大小写；
- 列表用英文逗号分隔，不允许尾逗号、连续逗号或空项；
- 路径不加引号。当前解析器按整行文本处理，因此路径本身不能包含 `=`。

错误示例：

```text
run.cfl = 0.2 # CFL
output.field.quantities = rho,u,
source.enabled = TRUE
```

### 7.2 无条件必填键

以下键即使某功能关闭也必须存在：

```text
schema_version
case.name
mesh.path
algorithm.profile
algorithm.flux_difference（可省略；默认 `profile`）
algorithm.reconstruction
algorithm.reconstruction_variables
algorithm.riemann
gas.gamma
reference.velocity
reference.density
reference.temperature
reference.length
reference.viscosity
partition.mode
partition.allow_idle_ranks
partition.max_load_ratio
partition.min_cells_per_active_direction
initial.type
boundary.default
source.enabled
run.mode
run.viscous
run.cfl
run.max_steps
output.directory
output.allow_existing
output.dimensional
output.field.enabled
output.history.enabled
output.statistics.enabled
output.checkpoint.enabled
```

此外 `gas.molar_mass` 和 `gas.specific_gas_constant` 必须恰好给一个。非定常必须有正的`run.t_end`；启用流场必须有非空 `output.field.quantities`；启用统计必须有非空 `output.statistics.quantities`；启用源项必须有非空 `source.models`。

## 8. 逐组填写配置

### 8.1 算例和路径

| 键 | 要求与行为 |
|---|---|
| `schema_version` | `1`，或 v2 开发分支的 `2` |
| `case.name` | 非空；输出文件名中非字母数字、`-`、`_` 字符会替换为 `_` |
| `mesh.path` | 非空；相对路径以配置文件所在目录为基准 |
| `restart.path` | 可选；相对路径同样以配置文件目录为基准 |
| `output.directory` | 相对路径以启动程序时的当前工作目录为基准，不是配置目录 |

为了避免路径基准混淆，推荐把网格与配置放在一个算例目录中，对 `mesh.path` 写短相对路径，
对 `output.directory` 写从固定启动目录可定位的绝对路径或明确的仓库根相对路径。

schema 2 的层流兼容配置为：

```text
turbulence.model = none
time.integrator = ssprk3
preconditioner.type = none
```

schema 2 可将模型段改为：

```text
turbulence.model = sa_neg
turbulence.prandtl = 0.9
turbulence.wall_treatment = resolved
turbulence.sa.farfield_nu_tilde_ratio = 3
turbulence.sa.source_treatment = local_implicit
```

SA-neg 必须与 `run.viscous = true`和至少一个 no-slip wall 同时使用；可把源项处理改为
`explicit`。当前 `local_implicit` 只是模型源项的局部对角更新，不是 LU-SGS。模型活动后，
定常停止判定会同时要求五个平均流残差和 `nu_tilde` 残差通过，checkpoint 会保存
`NuTilde` 及其参考残差。可输出 `nu_tilde,mu_t_over_mu,sa_production,sa_destruction,`
`wall_distance,sa_negative_branch`。SA-neg、SST-2003m 与实验级标准 k-epsilon 已实现；短
smoke 不能当作定量物理验证。阶段 Y 的真实 RANS
验收固定使用 TMR 平板、发展槽道和 no-plenum hump；资产散列、三级网格、阈值和命令见
[`cases/validation/sst_tmr`](../cases/validation/sst_tmr/README.md)，阶段状态以
[`stage-y-acceptance.md`](v2.0.0/stage-y-acceptance.md) 为准。可用
`wcns_extract_rans_profile <field.cgns> <target-x> <output.txt>` 从真实 CGNS 顶点壁面位置抽取
二维 RANS 剖面。

阶段 Z 候选支持三维非定常 `smagorinsky|scale_similarity|mixed_smagorinsky_similarity|`
`dynamic_smagorinsky|wale`。LES 强制 `run.viscous=true`、`time.integrator=lu_sgs`、BDF2
双时间和 `preconditioner.type=none`；二维、steady、SSPRK3 与无黏组合在启动前拒绝。完整且
可复制的模型专属键见 [`config-schema-2-draft.md`](v2.0.0/config-schema-2-draft.md)。其中测试
滤波比固定为 2，van-Driest 名称虽可识别，但在生产 wall-$y^+$ 场接入残差前明确拒绝。
schema 1 继续使用原配置，不应添加上述键。

### 8.2 算法选择

| 键 | 可选值 | 建议 |
|---|---|---|
| `algorithm.profile` | `phenglei_wcns`, `scmm6_wcns` | 同一网格可分别运行两次比较，不可混用组件 |
| `algorithm.flux_difference` | `profile`, `conservative_two_point` | 默认 `profile`；强激波正性困难时才显式选用守恒两点差分 |
| `algorithm.reconstruction` | `zero_order`, `linear5`, `weno_js`, `weno_z`, `mdcd_linear`, `mdcd_hybrid` | 间断优先从 `weno_z`/`mdcd_hybrid` 开始；`zero_order` 主要用于调试和高耗散基线 |
| `algorithm.reconstruction_variables` | `conservative`, `primitive`, `characteristic` | 强间断通常用 `characteristic` |
| `algorithm.riemann` | `rusanov`, `hllc`, `roe` | Rusanov 更耗散；HLLC/Roe 分辨率更高 |

六种重构都保持六点标量模板和三层 cell-centered ghost。`zero_order` 不缩小模板：同一面六点为 `(q[j-2],q[j-1],q[j],q[j+1],q[j+2],q[j+3])` 时，左值严格取第 3 点 `q[j]`，右值严格取第 4 点 `q[j+1]`。非法重构状态会按确定性策略回退，HLLC/Roe 的非法中间状态也会回退；历史文件记录累计回退数。回退不是静默成功，数量异常增大时应检查网格、CFL、初边值和正性。

`algorithm.flux_difference=profile` 严格使用 profile 冻结的高阶通量差分：
`phenglei_wcns` 使用其 D4/D2 路径，`scmm6_wcns` 使用其 D6/D4 路径。
`conservative_two_point` 则对相邻的左右真实面数值通量作差；它仍使用所选 profile 的
坐标度量、同一界面重构和同一 Riemann 求解器，并保持有限体积意义下的共享面守恒，
但通量散度不再是 profile 的高阶 D4/D6 算子。该模式是强激波出现负密度或负内能时的
显式稳定性选项，不能把其结果表述成“纯 profile 高阶通量差分”结果。该键会进入摘要和
restart signature，修改后不能直接续算旧检查点。

`linear5` 主要供线性回退、光滑基线和算法测试使用；普通有激波计算不把它或 `zero_order` 当作高分辨率首选方案。

尖后缘曲线网格可显式调整高阶 Jacobian 与有限体积参考体积的一致性门槛：

```text
geometry.metric.fallback = strict
geometry.metric.maximum_reference_relative_difference = 0.35
```

默认值分别为 `strict` 和 `0.20`。程序逐单元比较高阶 Jacobian 与正的有限体积参考体积；
任一值非有限、非正，或相对差超过门槛时，`strict` 在推进前失败。仅
`phenglei_wcns` 可选 `phenglei_finite_volume` 回退，且日志会报告替换单元数；该回退改变数值
离散身份，不应作为修补坏网格的通用手段。Family II NACA0012 的尖尾缘三层网格经只读检查
得到最大相对差 `0.349915/0.264342/0.165909`，因此该算例族冻结为 `strict + 0.35`，没有发生
有限体积替换。非默认设置进入 restart signature，续算时必须一致。

### 8.2.1 SSPRK 事后稳健化（v1.1）

```text
robustness.enabled = true
robustness.max_local_recomputations = 3
robustness.max_step_retries = 4
robustness.time_step_reduction = 0.5
robustness.minimum_time_step = 1e-12
```

该功能默认关闭；关闭时继续使用 v1.0 的直接 SSPRK3 提交路径。开启后，每个 RK 子步先在独立
缓冲区形成候选态，检查五个守恒分量有限，并检查 `rho/p/T` 不低于统一 `NumericalFloors`、
比内能为正。候选不合法时不会夹断状态，也不会覆盖当前正式状态，而是把失败单元标记为
troubled cell。

程序先找出真实残差中能影响该单元的直接面通量，再增加一层转置离散支持保护：凡残差会
读取这些直接面的真实单元，其完整面支持也并入本轮请求。对守恒两点差分，直接支持是相邻的
`i`、`i+1` 两面，保护层自然形成一个单元宽的局部带；对 profile 差分则直接读取与残差计算
相同的 `LineOperators` 行。连接附近还包含 PH 的
`i-1...i+2` 或 SCMM6 的 `i-2...i+3` 连接 halo 支持。多块/MPI 情况下 receiver 把请求反向
传播给共享面的唯一 owner，owner 对所有请求取最大级别并计算一份权威通量，再沿正常方向
发布。完整原理、系数公式、伪代码与终止性说明见 [`../算法补充.md`](../算法补充.md) 11.2.3。

有效降阶序列为“用户方案 → 同重构 primitive → `linear5/primitive` →
`zero_order/conservative + Rusanov`”，相邻重复组合自动删除；每轮每面最多升一级，同一 RK
阶段只升不降。超过局部重算上限后拒绝整步，从未修改的 `U^n` 恢复，并按
`dt <- time_step_reduction * dt` 重试；只有接受的时间步才增加 step/time 和产生正常输出。
五个配置键都会进入 restart signature，修改后不能直接续算旧检查点。

MDCD 的两个线性谱控制参数已经暴露：

```text
algorithm.mdcd.disp = 0.0463783
algorithm.mdcd.diss = 0.01
```

`disp` 对应 `WcnsParameters::mdcd_dispersion`，`diss` 对应 `mdcd_dissipation`；二者同时作用于 `mdcd_linear` 与 `mdcd_hybrid` 的线性部分。两值必须有限，并满足 `disp>0`、`0<=diss<disp`、`3*disp+9*diss<1`。省略时使用上面的默认值。即使当前选择不是 MDCD，出现的键仍被解析、验证、写入摘要和 restart signature；因此不要在不同重启段任意改值。

当前生产配置仍**没有**低 Mach 预处理键，也没有暴露 WENO epsilon、MDCD 传感器常数、Roe 熵修正、数值 floor 或无粘边界强约束开关。底层 `strong_boundary_face_state` 默认且实际为 `true`；标准 `wcns_run` 不能通过 `.wcns` 关闭它。

### 8.3 气体与无量纲参考量

`gas.gamma>1`。用摩尔质量时：

```text
gas.gamma = 1.4
gas.molar_mass = 0.029
```

用比气体常数时：

```text
gas.gamma = 1.4
gas.specific_gas_constant = 287.0
```

两者不能同时存在。程序使用 `R0=8.314 J/(mol*K)`，并由 `R=R0/M` 或反向关系得到另一项。

五个参考量必须是有量纲有限正数：

```text
reference.velocity = U_ref
reference.density = rho_ref
reference.temperature = T_ref
reference.length = L_ref
reference.viscosity = mu_ref
```

程序导出：

```text
Re = rho_ref * U_ref * L_ref / mu_ref
Ma = U_ref / sqrt(gamma * R * T_ref)
p_ref = rho_ref * U_ref^2
t_ref = L_ref / U_ref
```

禁止配置 `Re`、`Ma`、`reference.reynolds` 或 `reference.mach`。内部状态使用 `rho/rho_ref`、`u/U_ref`、`T/T_ref`、`p/(rho_ref*U_ref^2)`。因此理想气体无量纲关系为 `p=rho*T/(gamma*Ma^2)`；`rho=1,T=1` 并不意味着 `p=1`。

#### 层流输运（v1.1）

旧配置不写任何 `transport.*` 键时严格迁移为：

```text
transport.model = constant
transport.prandtl = 0.72
transport.constant.viscosity_ratio = 1.0
```

常黏度比是 `mu_const*/mu_ref`。Sutherland 模式写为：

```text
transport.model = sutherland
transport.prandtl = 0.72
transport.sutherland.reference_viscosity_ratio = 1.0
transport.sutherland.temperature = 110.4
```

最后一行单位为 K，也可唯一地替换成
`transport.sutherland.temperature_ratio = S*/T_ref`。两种表示若物理等价，内部配置与重启签名
相同；两者同时出现会失败。Sutherland 模式必须显式给参考黏度比和一个温度常数，常黏度与
Sutherland 专属键不能混用。所有值必须有限且为正。

内部公式为 `mu(T)=mu_Tref_ratio*T^(3/2)*(1+S)/(T+S)` 和
`chi=mu/[(gamma-1)*Ma^2*Pr]`；`1/Re` 仍只由粘性散度和时间步限制施加。启动摘要会在初场或
重启状态上报告全局 `T_range` 与真实 `mu_range`，`viscosity`/`mu_w` 输出也使用这个有效模型。
修改任一输运参数后旧检查点会被拒绝；缺省 v1.0 输运检查点允许按固定默认语义恢复。

### 8.4 MPI 分区

| 键 | 约束 | 含义 |
|---|---|---|
| `partition.mode` | 三种枚举 | `zones_only` 只用原 zone；`auto_split`/`force_split` 可二分 zone |
| `partition.allow_idle_ranks` | `true\|false` | 是否允许某些 rank 无叶块 |
| `partition.max_load_ratio` | `>=1` | 最大 rank 单元负载与平均负载的目标比 |
| `partition.min_cells_per_active_direction` | 整数 | 每个叶块每个活动方向最少单元数 |

当前 `auto_split` 与 `force_split` 在实现中走同一可切分流程；`force_split` 是为后续更明确策略保留的配置名，不能理解为“即使没有并行需要也必然增加叶块”。

`zones_only` 在 zone 少于 rank 且不允许 idle 时失败。另两种模式先把较大叶块确定性二分到足够覆盖 rank，再在负载比超过目标时继续切分，最多探索到约 `4*ranks` 个叶块。负载只按单元数估计，不计通信面、边界或异构硬件。

profile 的硬下限为：`phenglei_wcns>=4`，`scmm6_wcns>=5`。若 rank 数超过按该下限可形成的最大叶块数，且 `allow_idle_ranks=false`，启动明确失败。通常取 8 或更大更稳妥。

### 8.5 初场

所有初场参数都是无量纲量。只保留所选初场相关的键，避免以后切换类型时误用旧参数。

#### `uniform`

```text
initial.type = uniform
initial.rho = 1.0
initial.u = 0.2
initial.v = 0.0
initial.w = 0.0
initial.temperature = 1.0
```

默认 `rho=1,u=v=w=0,T=1`。可用 `initial.pressure` 代替温度；如果两者都给，当前实现优先压力，但手册要求只给一个以保持语义清楚。

#### `sod_x`

```text
initial.type = sod_x
initial.x0 = 0.5
initial.left_rho = 1.0
initial.left_u = 0.0
initial.left_v = 0.0
initial.left_p = 1.0
initial.right_rho = 0.125
initial.right_u = 0.0
initial.right_v = 0.0
initial.right_p = 0.1
```

`x<x0` 用左态，否则用右态；w 恒为零。以上值也是内建默认值。

#### `quadrant_riemann`

使用 `x0,y0` 和 `ne/nw/sw/se` 四组 `rho,u,v,p`。判断规则是东侧 `x>=x0`、北侧`y>=y0`。内建默认状态为：

| 象限 | rho | u | v | p |
|---|---:|---:|---:|---:|
| NE | 1.5 | 0 | 0 | 1.5 |
| NW | 0.5323 | 1.206 | 0 | 0.3 |
| SW | 0.138 | 1.206 | 1.206 | 0.029 |
| SE | 0.5323 | 0 | 1.206 | 0.3 |

完整键名可直接按下列模板填写：

```text
initial.type = quadrant_riemann
initial.x0 = 0.5
initial.y0 = 0.5
initial.ne_rho = 1.5
initial.ne_u = 0.0
initial.ne_v = 0.0
initial.ne_p = 1.5
initial.nw_rho = 0.5323
initial.nw_u = 1.206
initial.nw_v = 0.0
initial.nw_p = 0.3
initial.sw_rho = 0.138
initial.sw_u = 1.206
initial.sw_v = 1.206
initial.sw_p = 0.029
initial.se_rho = 0.5323
initial.se_u = 0.0
initial.se_v = 1.206
initial.se_p = 0.3
```

人工验收 case01 展示了完整配置与结果。

#### `isentropic_vortex`

```text
initial.type = isentropic_vortex
initial.x0 = 5.0
initial.y0 = 5.0
initial.beta = 5.0
initial.background_u = 1.0
initial.background_v = 1.0
initial.period_x = 10.0
initial.period_y = 10.0
```

周期域必须显式给正的 `period_x/period_y`，使初场用最短周期距离；0 表示普通非周期距离。程序根据 `gamma` 和导出的 Ma 计算温度、密度并检查正性。case03 给出扭曲周期网格的一周期示例。

#### `couette`

用 `eta=(y-y0)/(y1-y0)`：

```text
u = lower_velocity + (upper_velocity-lower_velocity)*eta
    + velocity_curvature*eta*(1-eta)
T = lower_temperature + (upper_temperature-lower_temperature)*eta
    + temperature_curvature*eta*(1-eta)
```

对应配置键为：

```text
initial.type = couette
initial.y0 = 0.0
initial.y1 = 1.0
initial.lower_velocity = 0.0
initial.upper_velocity = 1.0
initial.velocity_curvature = 0.0
initial.lower_temperature = 1.0
initial.upper_temperature = 1.0
initial.temperature_curvature = 0.0
initial.pressure = 1.0
```

压力恒定，密度由状态方程得到。默认 `y0=0,y1=1,lower_velocity=0,upper_velocity=1`，两个曲率为 0，端温默认 1。压力省略时用 `1/(gamma*Ma^2)`。

#### `poiseuille`

```text
initial.type = poiseuille
initial.y0 = 0.0
initial.y1 = 1.0
initial.centerline_velocity = 1.0
initial.temperature = 1.0
initial.temperature_curvature = 0.0
initial.pressure = 1.0
```

速度 `u=4*Ucenter*eta*(1-eta)`，横向速度为零。温度使用手册模板对应的四次多项式修正；压力恒定、密度由状态方程闭合。维持周期通道流还必须配置 `pressure_gradient` 源项。case02 展示了三维均匀/壁面加密网格，但加密算例在 3600 步被人工强制停止，不能作为已收敛基线。

#### `turbulent_channel`

该类型专用于三维双壁周期槽道的非定常初场：

```text
initial.type = turbulent_channel
initial.y0 = -1.0
initial.y1 = 1.0
initial.x0 = 0.0
initial.z0 = 0.0
initial.period_x = 6.283185307179586
initial.period_z = 3.141592653589793
initial.re_tau = 180.0
initial.bulk_velocity = 1.0
initial.bulk_velocity_plus = 15.481978793165828
initial.perturbation_amplitude = 0.05
initial.rho = 1.0
initial.temperature = 1.0
```

平均速度用对称 Reichardt 类复合壁律，并用一个高次小修正使中心线导数为零；
三个速度叠加可解析的低阶 x/z 周期扰动，扰动在两面壁上为零且 x-z 面平均为零。
`period_x/period_z`、`re_tau`、`bulk_velocity`、`bulk_velocity_plus` 必须为正，
`perturbation_amplitude` 只允许 `[0,0.5]`。修改壁律或参考尺度时必须重新计算
`bulk_velocity_plus`、体系 Re 和驱动体积力。完整公式、稀疏网格限制和可执行示例见
[`case05`](../cases/manual/case05_3d_turbulent_channel/README.md)。
case05 以初始体积平均速度 `U_ref=U_b,0` 缩放，所以无量纲
`bulk_velocity=1`，而 `bulk_velocity_plus=U_b,0/u_tau` 仍为 15.4819787932。
因此参考 Reynolds 数是 `Re_b^(h)=U_b^+*Re_tau`，体积力是
`a_x*h/U_b,0^2=1/(U_b^+)^2`。若使用其他速度尺度，必须联动修改这些值，不能照搬。

#### `linear_conduction`

速度为零，`T=lower_temperature+(upper_temperature-lower_temperature)*eta`，压力恒定。默认下/上温为 1/2。

#### `manufactured_periodic`

`beta` 默认 0.01，背景速度默认 `(0.2,0.1,0)`；在 x/y（以及三维 z）用固定 `2*pi` 三角函数构造光滑非均匀场。它用于程序回归，不会自动匹配任意域周期或任意制造源。

#### `double_mach_reflection`

该类型只实现经典 Woodward--Colella 二维无粘问题，并固定采用 `gamma=1.4`、入射 Mach 10、激波与 x 轴夹角 60 度、激波足点默认 `x0=1/6`。标准区域为 `[0,4]x[0,1]`：

```text
initial.type = double_mach_reflection
initial.x0 = 0.16666666666666667
```

激波前状态为 `(rho,u,v,p)=(1.4,0,0,1)`，激波后状态为`(8,8.25*cos(30deg),-8.25*sin(30deg),116.5)`。初始激波线为`x=x0+y/tan(60deg)`；激波后状态位于其左侧。该初场必须配合下文专用边界，且 `run.viscous=false`、`source.enabled=false`。非二维网格、非 1.4 的 gamma 或只设置初场而没有专用边界都会在推进前失败。

### 8.6 物理边界

可选类型及行为：

| 类型 | ghost/真实面行为 | 需要的数据 |
|---|---|---|
| `farfield` | 目标态与内部态按特征方向组合 | 目标态；省略时取初场在原点的状态 |
| `inflow` | 同一特征边界算法，超声速入流完全取目标态 | 目标态；省略时取初场在原点的状态 |
| `outflow` | 无目标时外推；有目标时使用特征组合 | 目标态可选 |
| `slip_wall` | 反射相对壁面法向速度 | 壁速可选，默认零 |
| `symmetry` | 与静止滑移壁相同的法向反射 | 通常无需数据 |
| `no_slip_adiabatic_wall` | ghost 反射三速度，真实粘性面强制无滑移和零法向温度梯度 | 壁速可选 |
| `no_slip_isothermal_wall` | ghost/真实粘性面施加无滑移与给定壁温 | 壁温必需；程序有回退值但建议显式给 |
| `double_mach_reflection` | 经典算例的左固定激波后态、顶面移动激波和底面分段入流/滑移壁 | 只与同名初场配套；无需普通目标态字段 |

目标态写法：

```text
boundary.inlet.type = inflow
boundary.inlet.rho = 1.0
boundary.inlet.u = 0.2
boundary.inlet.v = 0.0
boundary.inlet.w = 0.0
boundary.inlet.temperature = 1.0
```

也可把最后一行换为 `boundary.inlet.pressure=...`，但温度和压力不能同时给。一旦出现任意目标态字段，`rho` 与恰好一个 `temperature|pressure` 必须存在；缺省速度分量按零处理。

壁面写法：

```text
boundary.default = outflow
boundary.bottom.type = no_slip_isothermal_wall
boundary.bottom.wall_velocity_x = 0.0
boundary.bottom.wall_velocity_y = 0.0
boundary.bottom.wall_velocity_z = 0.0
boundary.bottom.wall_temperature = 1.0
```

patch 名必须与 CGNS 名完全一致。自动分区会保留原 patch 名。只填面 ghost 的三层法向数据；边和角 ghost 不填，且物理边界 ghost 上只有由边界条件得到的原始量、压力及守恒量可用，坐标、度量、梯度和其他二级量不可读取。

当前无粘真实边界面在重构后总执行强约束；粘性壁面速度/温度或热流约束也总是强制。配置 schema 尚未暴露 `strong_boundary_face_state`。

阶段 X 的二维 NACA numerical-analysis 对照可在 farfield 上叠加固定点涡：

```text
boundary.point_vortex.enabled = true
boundary.point_vortex.lift_coefficient = 1.09125
boundary.point_vortex.center_x = 0.25
boundary.point_vortex.center_y = 0.0
boundary.point_vortex.chord = 1.0
```

它按给定的 `lift_coefficient` 一次性构造远场目标速度，不读取或反馈当前计算载荷；正升力采用
顺时针环量约定。只允许二维算例，且只改变 `farfield` 目标态的速度，密度与温度仍取原远场值。
弦长必须为正，点涡中心不得落在远场面中心；非法组合在推进前或首次边界构造时明确失败。
默认关闭，启用后的参数全部进入 restart signature。公式和符号约定见
[`../算法补充.md`](../算法补充.md) 12.11。

经典双马赫反射的边界配置必须按 CGNS patch 名显式写成：

```text
boundary.default = outflow
boundary.left.type = double_mach_reflection
boundary.bottom.type = double_mach_reflection
boundary.top.type = double_mach_reflection
```

左边界恒取激波后态；顶边界使用`x_s(y,t)=x0+y/tan(60deg)+10*t/sin(60deg)`判断移动激波两侧；底边界在`x<x0`取激波后态，在`x>=x0`按静止滑移壁反射法向速度；右边界由 default 保持出流。专用算子只用真实边界顶点求面中心，不构造或读取 ghost 坐标，三层 ghost 的物理状态在每个 SSPRK 子步按该子步时间刷新。

### 8.7 源项

关闭源项：

```text
source.enabled = false
```

此时必须省略 `source.models`，所有已出现的 `source.*` 数值也必须为 0。开启时：

```text
source.enabled = true
source.models = body_force,pressure_gradient
source.body.ax = 0.0
source.body.ay = -1.0
source.body.az = 0.0
source.pressure_gradient.x = 0.01
source.pressure_gradient.y = 0.0
source.pressure_gradient.z = 0.0
```

| 模型名 | 参数 | 守恒源含义 |
|---|---|---|
| `uniform_conservative` | 五个 `source.uniform.*` | 直接加到五分量守恒方程 |
| `body_force` | `ax,ay,az` | 动量源 `rho*a`，能量源 `(rho*u)·a` |
| `pressure_gradient` | `x,y,z` | 保存 `G=-grad(p)`，动量源 `G`，能量源 `u·G` |
| `manufactured` | 五个 `source.manufactured.*` | 幅值乘 `1+x+y(+z)+t` |

多个模型同址相加并在每个 SSPRK 子步按该子步时间计算。模型名不可重复，二维 z 动量总源必须为零。所有源项值是当前无量纲控制方程中的值，程序不解析带单位表达式。

### 8.8 运行模式和停止

公共配置：

```text
run.mode = unsteady
run.viscous = false
run.cfl = 0.2
run.max_steps = 1000000
run.t_end = 0.3
run.max_wall_time = 0
```

- `run.viscous=false` 解热完全 Euler；`true` 解层流 Navier–Stokes。
- CFL 必须为正。稳定上限依网格、状态、重构和粘性尺度而变；默认值不是稳定性保证。
- `run.t_end` 是无量纲时间。量纲时间为 `t*t_ref`。
- `run.max_wall_time` 单位秒；正值要求 `output.checkpoint.enabled=true`。

定常停止配置：

```text
steady.min_steps = 100
steady.check_interval_steps = 10
steady.consecutive_checks = 3
steady.reference_floor = 1e-30
steady.l2_absolute = 1e-12
steady.l2_relative = 1e-8
steady.linf_enabled = true
steady.linf_absolute = 1e-11
steady.linf_relative = 1e-8
```

在 `step % check_interval_steps == 0` 时检查。第一次检查冻结五分量参考 `L2/Linf`。每个分量必须满足 `L2绝对阈值 OR L2相对阈值`，并在启用时满足 `Linf绝对 OR Linf相对`；全部分量连续通过指定次数且达到 `min_steps` 才停止。

停止优先级是：数值失败 → steady 收敛或 t_end → 用户信号 → 墙钟 → 最大步数。SIGINT/SIGTERM 和墙钟只在完整步结束时生效。不要用任务管理器“结束进程”代替一次 Ctrl+C，因为强制杀进程无法保证原子提交最终检查点。

### 8.9 输出总设置与调度

```text
output.directory = output/my_case
output.allow_existing = false
output.dimensional = false
```

`allow_existing=false` 时只要目录存在即失败，保护历史结果。`true` 允许目录存在并允许同名最终文件被替换，但历史/统计文件是本次运行重写而不是追加。重启推荐写到新目录；如必须复用目录，应先备份并理解覆盖行为。

每类输出都有：

```text
every_steps = 0
every_time = 0
explicit_times = 0.1,0.25,0.5
write_initial = false
write_final = true
```

在键前加对应前缀 `output.field|history|statistics|boundary|checkpoint`。五种事件取并集；同一`(step,time)` 去重。`explicit_times` 必须非负、严格递增。`every_time` 和显式时刻会让时间步裁剪到事件时刻；定常时这里的 time 是伪时间。

## 9. 输出文件逐项说明

### 9.1 流场快照

```text
output.field.enabled = true
output.field.format = cgns
output.field.every_steps = 100
output.field.write_initial = true
output.field.write_final = true
output.field.quantities = rho,u,v,w,p,T
```

格式可为 `cgns|tecplot|both`。文件名类似：

```text
my_case.field.step00000100.time1p250000000eM02.cgns
my_case.field.step00000100.time1p250000000eM02.dat
```

内建字段：

| 名称 | 含义 | 量纲恢复尺度 |
|---|---|---|
| `rho` | 密度 | `rho_ref` |
| `u,v,w` | 速度 | `U_ref` |
| `p` | 压力 | `rho_ref*U_ref^2` |
| `T` | 温度 | `T_ref` |
| `rho_u,rho_v,rho_w` | 动量密度 | `rho_ref*U_ref` |
| `rho_E` | 总能量密度 | `rho_ref*U_ref^2` |
| `sound_speed` | 声速 | `U_ref` |
| `mach` | 局部速度模/声速 | 1 |
| `total_enthalpy` | 单位质量总焓 | `U_ref^2` |
| `entropy_proxy` | `p/rho^gamma` 代理量 | 当前保持无量纲 |
| `viscosity` | 当前输运模型黏度 | `mu_ref` |
| `jacobian` | 计算到物理空间 Jacobian | `L_ref^dimension` |

CGNS 输出在 rank 0 把运行时叶块重新拼回输入原 zone，写 cell-centered `FlowSolution` 和时间descriptor。Tecplot 是 ASCII ordered zone，坐标取 cell center。当前场输出只复制坐标和场，不承诺完整保留输入 BC/connectivity，因此它是后处理文件，不是正式检查点，也不应替代原网格。

### 9.2 残差历史

```text
output.history.enabled = true
output.history.format = txt
output.history.every_steps = 1
output.history.write_initial = true
output.history.write_final = true
```

禁止设置 `output.history.quantities`；历史 schema 固定。文件名为`<case>.history.r<ranks>.txt|dat`。TXT 首行以 `#` 开头，依次包含：step、time、dt、CFL、
wall time、总 L2、五分量 L2、参考 L2、归一化 L2、五分量 Linf、参考 Linf、归一化 Linf、连续通过次数、重构回退数、Riemann 回退数、稳健化 level 0--3 owner 面数、troubled-cell 数、局部重算数、整步 retry 数、proposed/accepted dt、候选最小 `rho/p/T/e`、本步是否检查残差、停止原因。初始行 `dt=0`；重启后的首行也不要当成本步时间步长。稳健化关闭时级别计数和重试计数为零，尚无候选最小值时写 `nan`。

### 9.3 全场统计

```text
output.statistics.enabled = true
output.statistics.format = txt
output.statistics.every_steps = 10
output.statistics.quantities = total_mass,total_momentum_x,total_energy
```

可选 `total_mass,total_momentum_x,total_momentum_y,total_momentum_z,total_energy`。这些值使用原 zone 守恒积分权重、正 Jacobian 和 MPI 全局归约，避免人工切分接口重复计数。`output.dimensional=true` 时还乘相应场尺度和 `L_ref^dimension`；二维结果表示单位出平面厚度下的积分，应按二维模型解释。

schema 2 的非定常计算可对上述已选量启用接受步时间累计：

```text
statistics.time.enabled = true
statistics.time.start = 1.0
statistics.time.end = 5.0
statistics.time.every_steps = 1
statistics.time.weight = accepted_dt
```

统计区间端点按物理步与区间的实际重叠时间裁剪；LU-SGS 内迭代、失败尝试和同一
`(step,time)` 重放不采样。最终 `.time_statistics.r<ranks>.txt` 保存每个量的加权 mean/RMS、
Favre mean/RMS 以及所有唯一量对的 Reynolds/Favre covariance。累计器 v2 状态随 checkpoint
保存，允许改变合法 rank 数后连续累计。

三维 Poiseuille/槽道流还可打开多个 x-z 层监测：

```text
output.statistics.enabled = true
output.statistics.xz_planes.enabled = true
output.statistics.xz_planes.cell_j_indices = 0,11,23,35,47
```

索引是每个**原始 CGNS zone** 的零基 cell-j 索引，不是顶点索引，也不是自动分区后的叶块局部索引。开启后程序按给定次序为每个 `j` 自动追加两列：`xz_mean_u_j<j>` 为面积加权`<u>_xz=sum(A*u)/sum(A)`；`xz_mass_flow_x_j<j>` 为`sum(A*rho*u)`。后者是 x 方向质量流动在 x-z 壁平行面的面积积分指标，不是穿过该面的法向通量；真正的通道截面流量应使用下述法向为 x 的 y-z 面接口。

所有原 zone 必须为三维、每个索引在各 zone 中都有效，并且同一索引对应的所有单元中心 y 坐标必须在容差内共面。面积使用该 cell-j 层上下两个 J 面面积的平均；人工 MPI 切分时通过 `PartitionLeaf` 映射回原 zone，再由 MPI 求和，因此不会重复累计 ghost 或接口。无量纲输出时两列分别为速度和`rho*u*L^2`；量纲输出分别乘 `U_ref` 和`rho_ref*U_ref*L_ref^2`。关闭功能可以保留索引列表但不会生成列；开启时无需、也不要把自动列名手工重复写进`output.statistics.quantities`。

真实 y-z 截面平均速度和 x 向质量流量使用目标 x 坐标配置：

```text
output.statistics.yz_planes.enabled = true
output.statistics.yz_planes.target_x_coordinates = 0.0,3.141592653589793
```

对每个目标，程序选择与其重合或正 x 侧最近的常 x 单元中心面，并按配置顺序自动追加
`yz_mean_u_plane<n>` 和 `yz_mass_flow_x_plane<n>`。前者为面积加权平均速度，后者为
`integral(rho*u dA_yz)`，即穿过截面的真实 x 向质量流量。程序只支持三维、几何上常 x
的平面；不会在任意曲面上做隐式插值。这两类自动列也不应重复写入
`output.statistics.quantities`。

三维等温槽道还可开启两面壁摩擦统计：

```text
output.statistics.channel_walls.enabled = true
output.statistics.channel_walls.lower_patch = bottom
output.statistics.channel_walls.upper_patch = top
output.statistics.channel_walls.half_height = 1.0
```

开启后自动追加 `channel_wall_shear_lower`、`channel_wall_shear_upper`、
`channel_wall_shear_mean`、`channel_friction_velocity`、`channel_re_tau`。两个 patch 必须分别是
J-lower/J-upper 的 `no_slip_isothermal_wall`，网格必须为三维平面 x-z 壁，且
`run.viscous=true`。壁切应力是面积加权的流向切应力绝对值；摩擦速度用面平均壁密度，
`channel_re_tau` 还使用壁温对应黏性和所配半高。这些列是瞬时空间统计，不是累积时间平均。

### 9.4 边界面与载荷

```text
output.boundary.enabled = true
output.boundary.format = tecplot
output.boundary.every_steps = 100
output.boundary.write_initial = true
output.boundary.write_final = true
output.boundary.patches = cylinder
output.boundary.quantities = p_w,T_w,mu_w,Cp,Cf,q_wall,traction_x,traction_y
output.boundary.reference_pressure = 17.8571428759968
output.boundary.reference_density = 1.0
output.boundary.reference_velocity_x = 1.0
output.boundary.reference_velocity_y = 0.0
output.boundary.reference_velocity_z = 0.0
output.boundary.reference_area = 1.0
output.boundary.reference_length = 1.0
output.boundary.moment_center_x = 0.0
output.boundary.moment_center_y = 0.0
output.boundary.moment_center_z = 0.0
output.boundary.drag_direction_x = 1.0
output.boundary.drag_direction_y = 0.0
output.boundary.drag_direction_z = 0.0
output.boundary.lift_direction_x = 0.0
output.boundary.lift_direction_y = 1.0
output.boundary.lift_direction_z = 0.0
output.boundary.tangent_direction_x = 1.0
output.boundary.tangent_direction_y = 0.0
output.boundary.tangent_direction_z = 0.0
```

三维算例可再给出与 drag/lift 构成右手系的单位展向方向和严格递增的面心分箱边界：

```text
output.boundary.span_direction_x = 0
output.boundary.span_direction_y = 0
output.boundary.span_direction_z = 1
output.boundary.span_bin_edges = 0,0.5,1
```

`patches` 和 `quantities` 均为禁止重复的逗号列表。逐面内建量还包括
`pressure_traction_x/y/z`、`viscous_traction_x/y/z` 和 `traction_x/y/z`。RANS 无滑移壁还可输出
`wall_distance,friction_velocity,wall_y_plus,wall_y_plus_class`；最后一项依次用 0/1/2/3 表示
$y^+\le5$、$5<y^+<30$、$30\le y^+\le300$、$y^+>300$。这些量只作诊断，不改变
`resolved`/`wall_function` 的显式选择；无 RANS 模型或非无滑移 patch 请求它们会失败。`Cp/Cf` 使用配置中的
`q_inf=0.5*rho_inf*|u_inf|^2`；参考动压必须大于统一阈值，方向必须为单位向量且 drag/lift
正交。无粘运行请求 `Cf,q_wall` 或黏性牵引会在写文件前失败。

无粘边界的压力（以及显式请求时的温度）默认采用当前 profile 的高阶内部迹；若该迹在激波附近
出现非有限或不大于对应数值 floor 的输出侧超调，则仅对写出值退回该面的最近内部真实单元值。
两者都无效时运行失败。这个保护不修改求解状态、残差、通量、时间步或 restart signature；未请求
`T_w,mu_w` 的压力/载荷输出也不会额外计算温度迹与输运系数。

面面积元使用全局守恒权重，不按面数平均；运行时切分通过原 zone 索引恢复后，在 root 进行
确定性排序和查重。逐面文件名为 `<case>.boundary.r<ranks>.step....txt|dat`，几何列固定包含
原 patch/zone 索引、全局面索引、面心、面积和流体域外法向。固定载荷历史
`<case>.loads.r<ranks>.txt` 保存压力、黏性、总力/力矩以及 `Cd/Cl/Cm`。二维力按单位展向长度
解释并在文件头写 `force_per_unit_span=true`；启用展向分箱后，每次事件另写
`<case>.spanwise_loads.r<ranks>.step....txt`。面按面心
投影完整归入一个半开区间（末箱含右端点），分箱压力/黏性/总载荷之和应回收到整体载荷。
`output.dimensional=true` 时坐标、面积、牵引、
力和力矩分别按 `L_ref`、`L_ref^(d-1)`、`rho_ref U_ref^2`、
`rho_ref U_ref^2 L_ref^(d-1)` 和 `rho_ref U_ref^2 L_ref^d` 恢复量纲；`wall_distance` 与
`friction_velocity` 分别按 $L_{ref}$、$U_{ref}$ 恢复，$y^+$ 和分类保持无量纲，其余系数不变。

### 9.5 检查点

检查点固定保存五个无量纲守恒量，不使用 `quantities` 键。每次事件生成：

```text
my_case.checkpoint.step00001000.time....cgns
my_case.checkpoint.latest.cgns
```

带 step/time 的文件是事件快照；`latest` 是滚动副本。还保存格式版本、step/time/dt、网格签名、数值签名、定常参考残差与连续计数。

### 9.6 manifest 与临时文件

正常进入最终化后，rank 0 写 `<case>.manifest.r<ranks>.txt`，其中包含程序版本、Git 提交、编译器、Release/Debug、MPI 数、配置/分区摘要、网格/重启签名、最终状态、停止原因、末个接受步的稳健化级别/重试/最小状态诊断和成功提交的文件列表。审查结果时先看 manifest，再看 history。

输出先写 `.tmp`，关闭成功后改名。启动阶段异常可能没有 manifest；I/O 中断可能留下 `.tmp`，它只是诊断残留，不能当成有效结果。

## 10. 串行与 MPI 运行

串行：

```powershell
build-user-serial\wcns_run.exe --config path\case.wcns
```

MPI：

```powershell
mpiexec -n 4 build-user-mpi\wcns_run.exe --config path\case.wcns
```

配置只由 rank 0 读取并广播，随后所有 rank 校验摘要 digest；网格元数据和分区也要求一致。求解 halo 使用 MPI，但 CGNS 读写和原 zone 重组当前不是并行 I/O。rank 增多不保证更快，尤其小网格和频繁全场输出可能由通信/根进程 I/O 主导。

更换 rank 数前先执行 dry-run：

```powershell
mpiexec -n 8 build-user-mpi\wcns_run.exe --config path\case.wcns --dry-run
```

检查摘要中的每个 `leaf=...,range=...,owner=...`，确认没有不可行分区。正式可复现性检查应在同一网格/配置下运行 1 rank 和目标 rank，并用验证器比较最终场。

## 11. 检查点重启操作

### 步骤 1：源计算启用检查点

```text
output.checkpoint.enabled = true
output.checkpoint.every_steps = 1000
output.checkpoint.write_final = true
```

### 步骤 2：确认检查点已完整提交

查看源运行 manifest 中是否列出 `.checkpoint.latest.cgns`，并确认源停止原因。不要从 `.tmp` 或数值失败状态恢复。

### 步骤 3：复制配置并设置新输出目录

```text
restart.path = ../run-a/output/my_case.checkpoint.latest.cgns
output.directory = output/run-b
```

可以改变 rank 数、合法的叶块分区、输出设置、`run.max_steps` 和非定常 `run.t_end`。也可改变case 名。初场配置仍是 schema 必填，但恢复时不会用于覆盖检查点状态。

### 步骤 4：先 dry-run，再续算

```powershell
mpiexec -n 4 build-user-mpi\wcns_run.exe --config run-b\restart.wcns --dry-run
mpiexec -n 4 build-user-mpi\wcns_run.exe --config run-b\restart.wcns
```

程序要求 profile、重构、Riemann、气体、参考量、边界数据、源项、黏性开关和网格签名兼容。当前网格签名覆盖 base/zone 名称、维数、尺寸和坐标；不要依赖它发现所有 BC/connectivity 语义变化，实际重启应保持原网格文件不变，只改变运行时分区。

历史/统计不会把源文件自动拼接到新文件。分析连续轨迹时按 checkpoint 的 step/time 合并两次运行的序列，并去掉重复的重启初始行。

## 12. 常用算例操作配方

### 12.1 自由流和新网格最低验收

1. 生成或导入网格；设 `initial.type=uniform`。
2. 所有物理边界设合适远场/滑移壁，周期方向使用 CGNS 连接。
3. 分别选择两套 profile，短推 10–100 步。
4. 输出初末场、全部统计；检查回退数、场差和守恒量漂移。
5. 至少比较串行与目标 MPI rank。

### 12.2 二维 Riemann 问题

使用 `quadrant_riemann`、`outflow`、`run.mode=unsteady`、特征重构和 HLLC/Roe，输出初末场。完整 256² 均匀/局部加密操作见 [`case01`](../cases/manual/case01_2d_riemann/README.md)。

### 12.3 三维 Poiseuille

使用 `periodic-channel` 网格、上下无滑移等温壁、`poiseuille` 初场、`source.models=pressure_gradient`、`run.viscous=true` 和 `run.mode=steady`。完整设置及一次未完成的壁面加密运行分析见 [`case02`](../cases/manual/case02_3d_poiseuille/README.md)。

### 12.4 扭曲网格等熵涡

使用 `warped-periodic-square`、周期初场距离、两套独立 profile 各运行一周期；用`wcns_compare_metric_profiles` 和 `field-error` 比较。完整实测见
[`case03`](../cases/manual/case03_2d_vortex/README.md)。

### 12.5 经典双马赫反射

用 `wcns_generate_release_cgns rectangle ... 960 240 ... 4.0 1.0 false` 生成`[0,4]x[0,1]`结构网格；配置 `double_mach_reflection` 初场和 left/bottom/top 三个同名专用边界，以 WENO-Z 特征重构和 HLLC 从 `t=0` 推进到 `t=0.2`。完整命令、配置解释和输出判读见 [`case04`](../cases/manual/case04_2d_double_mach_reflection/README.md)。

### 12.6 \(Re_\tau=180\) 非定常湍流槽道

使用 x/z 双周期 `periodic-channel` 网格、y 向两面等温无滑移壁、
`turbulent_channel` 复合壁律加低模态扰动初场、`body_force` 定常体积力、
`run.mode=unsteady`，并开启 \(x=0,\pi\) 附近的 y-z 截面流量/平均速度和两壁摩擦统计。当前 36×48×36 网格只通过
4-rank、5 步工程可行性卡口，未作湍流统计/DNS 验收。完整公式、配置、命令、实测结果和
Linux 迁移前检查见 [`case05`](../cases/manual/case05_3d_turbulent_channel/README.md)。

### 12.7 二维圆柱低速与高超声速绕流

用 `wcns_generate_release_cgns cylinder-o` 生成带首尾周期连接的多块 O 网格。低速圆柱使用
可压缩层流 Navier--Stokes、绝热无滑移壁和远场边界，可通过 Re=20/40 与 Re=100/200
分别观察稳定对称尾迹和非定常涡脱落；Mach 5 钝体功能检查使用 Euler、滑移壁和远场边界。
当前粗网格结果只作定性验收。v1.1 已用求解器权威边界面迹直接输出压力/黏性牵引、热流，
并对全局边界权重积分得到升阻力、力矩和系数；Case07 保留旧后处理结果仅供历史对照。完整参数、
命令、实际结果、图像和限制见
[`case07`](../cases/manual/case07_2d_cylinder/README.md)。

## 13. 独立验证工具

最常用命令：

```powershell
wcns_validate_release_case finite final.cgns
wcns_validate_release_case compare serial.cgns mpi.cgns 1e-12
wcns_validate_release_case field-error initial.cgns final.cgns
wcns_validate_release_case tecplot-consistency final.cgns final.dat 1e-12
wcns_compare_metric_profiles mesh.cgns
wcns_inspect_structured_mesh mesh.cgns Wall
```

`wcns_inspect_structured_mesh` 的可选第二参数是需要详细检查的壁面 patch 名；省略时检查所有
CGNS wall 类型 patch。首层高度是壁面线与第一层网格线的中点间距，仅是几何量，不能在缺少
Re、物性和壁面剪切速度时直接解释为 $y^+$。

解析算例还支持 `uniform`、`vortex`、`sod`、`diagonal-symmetry`、`viscous-profile`、`poiseuille-profile`、`uniform-source`、`derived` 和 `nonzero`。完整参数表和 Python 矩阵驱动见 [`release-validation.md`](release-validation.md)。验证器独立通过 CGNS API 重读文件，不直接信任求解器内存结果。

## 14. 运行前、运行中和运行后检查表

### 14.1 运行前

- Release 构建及 CTest 零失败；MPI 运行库来自同一安装。
- 网格是 ADF-CGNS、Structured、2D/3D，边界范围和双向 1-to-1 连接正确。
- 网格和所有流动输入已无量纲化；参考量为有限正数。
- 配置无重复/未知键，gas 二选一，Re/Ma 未直接输入。
- 2D 的 w 和所有 z 源为零；初始密度、温度/压力为正。
- patch 名与 CGNS 完全一致；周期用 connectivity。
- 输出目录不存在或已明确允许覆盖；墙钟限制有检查点。
- 目标 rank 的 dry-run 成功，分区摘要合理。

### 14.2 运行中

- 观察 step/time/dt 是否推进，残差是否有限。
- 观察重构/Riemann 回退是否突然持续增长。
- 定常算例查看归一化残差和 `consecutive`，而非只看 total L2。
- 非定常算例确认 time 将精确命中事件和 t_end。
- 长算例确认检查点按计划更新；磁盘空间足够容纳原 zone 全场输出。

### 14.3 运行后

- 记录进程退出码；打开 manifest 确认 `stop_reason`。
- 检查最终 step/time、Git commit、MPI 数和文件清单。
- 用 `finite` 检查最终场；按算例做解析误差/对称性/守恒检查。
- 检查 statistics 漂移、history 回退计数和数值失败痕迹。
- MPI 计算与串行参考比较；保存配置、原网格、manifest、日志和验证输出。

## 15. 常见故障与处理

| 报错或现象 | 原因 | 处理 |
|---|---|---|
| `unknown configuration key` | 拼写、大小写或当前 schema 未支持 | 对照模板；不要随意添加规划中的键 |
| `must contain exactly one '='` | 行尾又有 `=` 或格式错误 | 注释独占一行；每行只保留一次赋值 |
| gas 要求 exactly one | 两个气体键都给或都没给 | 只保留摩尔质量或比气体常数之一 |
| Re/Ma cannot be configured | 直接输入派生量 | 删除该键，通过五个参考量控制 |
| `coordinate ... outside` | CGNS 范围越界 | 修复网格生成器的 1-based PointRange |
| `edge or corner, not a face` | BC PointRange 没覆盖完整面 | 在 CGNS 中写面范围 |
| connectivity/transform 错误 | 供体名、范围、变换或反向连接错误 | 修复双向共形 1-to-1 定义 |
| partition infeasible/idle | rank 太多或叶块太窄 | 减 rank、允许 idle、降低安全下限但不得低于硬下限、加密网格 |
| `output directory already exists` | 保护策略生效 | 使用新目录；确认后才设 `allow_existing=true` |
| unknown field/statistic | 名称未注册或重复 | 检查内建列表；自定义量需重新编译并注入 registry |
| checkpoint signature differs | 数值配置或网格不兼容 | 使用原配置/网格；不要绕过签名 |
| 退出码 2 | 最大步、墙钟或信号停止 | 看 manifest；增加资源/步数或从检查点续算 |
| 退出码 3 | 非有限状态、正性或算子失败 | 检查初边值、网格、CFL、源项和回退计数 |
| MPI 挂起或 rank 不一致 | 环境混用、网格连接问题、某 rank I/O 失败 | 用 1 rank 和 2 rank 最小复现，检查首个失败 rank 日志 |
| 结果“没有收敛但程序结束” | 达到 max_steps | `maximum_steps` 不是定常收敛，必须增加步数或修正设置 |

## 16. 当前功能边界

当前开发分支可运行单组分热完全理想气体、层流常比热、常黏度/Sutherland 输运、显式
SSPRK3、LU-SGS、BDF2 双时间、Weiss--Smith 低 Mach 预处理、结构共形网格和内建源项。
SA-neg、SST-2003m、实验级标准 k-epsilon 及五种三维 LES 已实现；阶段 Z 的本机证据仅是公式、
解析/微型场、重启和并行离散卡口，完整 HIT/槽道/能谱和长时间统计仍待服务器验证。化学反应、
通用表达式源项、动态插件和涡量/Q 等体派生输出仍未实现。
默认输运为 `Pr=0.72` 和 `mu/mu_ref=1` 的常黏度。

这些限制不能通过写一个未知配置键绕过。需要扩展时按开发手册同时修改数据结构、严格 parser、验证、摘要/重启签名、生产装配、测试、模板和文档。

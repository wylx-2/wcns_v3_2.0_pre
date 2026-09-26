# WCNS 用户自定义开发指南

本文面向需要修改或扩展 WCNS 的开发者，对应 `v2.0_pre` 的 schema 1/2 统一源码。目标不是只告诉读者“改哪个文件”，而是说明一次扩展必须穿过哪些数据、验证、并行、重启、输出和测试路径，避免新增代码在串行小算例中可运行、到多块/MPI/重启时失效。

精简的 `WCNS_v2.0_pre` 目录按发布要求不携带开发仓库中的 `tests/`、`cases/`、`examples/`、
人工算例结果和阶段记录。本文中涉及这些目录的回归方法仍用于说明扩展应达到的验证层级；
需要复现完整矩阵时，应使用 `wcns` 开发仓库。精简目录只能用用户自备网格和配置做端到端检查。

先按 [`user-manual.md`](user-manual.md) 完成串行与 MPI 构建，并阅读[`算法补充.md`](../算法补充.md) 和 [`known-limitations.md`](known-limitations.md)。本项目目前没有稳定的对外 ABI；“接口”是源码扩展点，不是无需重编译的动态插件接口。

## 1. 扩展前的不可破坏契约

任何算法或模型扩展都必须保持以下约束，除非该任务明确包含一次经过设计和验收的契约升级。

1. 内部单元索引是 0-based 半开区间；CGNS 是 1-based，转换后必须检查上下界。
2. 二维仍使用三维索引容器和五分量守恒状态，但 k extent 为 1、w/z 动量为零。
3. WCNS 重构当前冻结为六点标量模板：首偏移 -2、点数 6、halo 宽度 3。
4. 块 ID、rank ID 和本地容器下标不是同一概念，不得互相替代。
5. 原生连接和运行时切分连接都要支持同 rank 与跨 rank；周期连接还要正确旋转矢量/通量。
6. 一个共享连接面的数值通量只由确定性 owner 计算，再通信给另一侧；接收侧不可独立重算。
7. `phenglei_wcns` 和 `scmm6_wcns` 是不可交叉的完整 profile bundle。
8. 物理边界 ghost 上只有边界条件构造的物理量、压力和守恒量有效；边/角 ghost、ghost坐标、度量、梯度及任意派生量不可读。
9. 密度、压力、温度、Jacobian、面面积和重构尺度必须经过统一 `NumericalFloors`/有限性检查；不允许用无记录的截断或零值掩盖错误。
10. 所有 rank 必须以相同顺序进入集合通信；某个 rank 的异常必须转成一致失败，不能让其他rank 停在 MPI 调用中。
11. 会改变数值轨迹或状态解释的配置必须进入摘要和重启签名；所有输入键进入配置 digest。
12. 新功能必须同时覆盖 2D/3D 适用性、物理边界、原生多块、运行时切分、串行/MPI 和重启影响。

## 2. 代码架构和一次时间步的数据流

生产入口在 `src/app/wcns_run.cpp`。启动和推进关系如下：

```text
.wcns 配置 --严格解析/广播--> CaseConfig
CGNS 元数据 -----------------> StructuredPartitionPlan
CGNS zone + 分区叶块 --------> StructuredMesh + LocalBlockSet
连接/所属 rank --------------> DistributedTopology / halo plan
原 zone 坐标 + profile ------> MetricField，再按叶块切片分发
初场或 checkpoint -----------> 每个真实单元的守恒场
边界配置 --------------------> BlockBoundaryDataMap
                                      |
                                      v
每个 SSPRK3 子步：更新原始量 -> 守恒 halo -> 物理 ghost
                -> 左右面重构 -> 物理面强约束 -> Riemann 通量
                -> 共享面通量通信 -> profile 通量散度
                -> 粘性项（可选）+ 源项 -> residual -> RK 更新
                                      |
                                      v
SimulationDriver -> 停止判据 -> OutputSchedule -> 场/历史/统计/checkpoint/manifest
```

主要数据对象：

| 对象 | 位置 | 职责 |
|---|---|---|
| `StructuredBlock` | `include/wcns/mesh/structured_block.hpp` | 坐标、度量、流场、边界和连接的本地块容器 |
| `StructuredMesh` | `include/wcns/mesh/structured_mesh.hpp` | 全局块描述及拓扑 |
| `LocalBlockSet` | `include/wcns/parallel/block_distribution.hpp` | 当前 rank 拥有的块 |
| `MetricField` | `include/wcns/mesh/high_order_metrics.hpp` | profile 绑定的单元坐标、J、三向面面积矢量 |
| `CaseConfig` | `include/wcns/runtime/case_config.hpp` | 严格配置的唯一结构化表示 |
| `InviscidWcnsConfig` | `include/wcns/solver/inviscid_wcns_solver.hpp` | 重构、Riemann、无粘边界、源项配置 |
| `SimulationDriver` | `include/wcns/runtime/simulation_driver.hpp` | SSPRK 外层步进、时间裁剪、停止和 observer |
| `RuntimeOutputManager` | `include/wcns/runtime/output_manager.hpp` | 输出事件、历史/统计、manifest 和原子提交 |

## 3. 推荐开发流程

### 3.1 建立可重复基线

```powershell
git status --short --branch
cmake -S . -B build-dev-serial -G "MinGW Makefiles" `
  -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_MPI=OFF `
  -DWCNS_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-dev-serial --parallel 4
ctest --test-dir build-dev-serial --output-on-failure
```

再建立独立 MPI 目录并运行 CTest。保存基线的 Git commit、测试结果和一个代表性算例 manifest。
当前仓库不启用外部 CI，因此开发者必须在合并前完成与改动范围相称的本机串行/MPI 自动卡口。
已有未提交修改属于用户工作，新增功能不能覆盖或顺带格式化无关文件。

### 3.2 格式与静态清理

C/C++ 采用仓库根目录 `.clang-format`，Python 采用 `pyproject.toml` 中的 Black 配置；行尾和
基础编辑约束见 `.gitattributes`、`.editorconfig`。本机清理命令为：

```powershell
$cpp = git ls-files '*.cpp' '*.hpp' '*.h' '*.c' '*.cc' '*.cxx'
clang-format -i --style=file $cpp
$python = git ls-files '*.py'
black $python
clang-format --dry-run --Werror --style=file $cpp
black --check $python
```

机械格式改动应单独提交。删除“冗余代码”前必须证明它不是公开扩展点、条件编译路径或故意
保留的独立验证入口；清理后至少执行受影响测试，正式版本前重新执行完整串行/MPI 回归。

### 3.3 把扩展拆成小卡口

建议顺序：

1. 写数学/接口约定和适用范围；
2. 增加数据结构与严格验证，暂不接生产路径；
3. 写最小单元测试和失败测试；
4. 接入串行生产路径；
5. 接入多块/MPI、周期变换和共享面 owner 路径；
6. 更新配置摘要、restart signature、模板和文档；
7. 运行串行/MPI CTest；
8. 用正式 `wcns_run` 做具体算例，保存配置、日志、manifest 和独立验证结果；
9. 每个被认可的小卡口单独提交，提交中不要混入算例生成物或用户未审核修改。

### 3.4 完成定义

“代码能编译”不是完成。一个扩展至少应有：合法输入测试、非法输入测试、解析或制造结果、自由流/守恒检查、2D/3D 适用性、物理边界、原生连接、运行时切分、1/多 rank 等价、输出选择、重启兼容/拒绝、README/模板更新和 Release 构建实测。

## 4. 增加一个配置键

配置 parser 不是自动反射，新增键必须完整经过以下步骤。

1. 在合适的 config 结构中增加带安全默认值的成员，例如`include/wcns/runtime/case_config.hpp`。
2. 把完整键名加入 `src/runtime/case_config.cpp::fixed_keys()`。按 patch 动态命名的边界键需更新 `dynamic_boundary_key()` 属性集合。
3. 在 `CaseConfig::from_text()` 中用 `require` 或 `optional_*` 解析。枚举写专用解析函数；布尔只接受 `true|false`；实数必须有限。
4. 在对应 `validate()` 中检查数值范围、互斥/依赖条件、维数约束和功能开关状态。
5. 在 `summary()` 中输出用户可审查的最终值。
6. 如果该键影响离散算子、物理状态、边界、源项或恢复后轨迹，把它加入`CaseConfig::restart_signature()` 或相应子配置的 `restart_signature()`。
7. 把它传到真正消费该值的生产对象；只解析但不使用属于错误。
8. 在 `tests/test_case_config.cpp` 增加默认、合法、非法、未知/重复键和摘要/签名测试。
9. 更新 `examples/full_case_template.wcns`、用户手册和相关算例。

配置 digest 对规范化后的全部输入键值做 FNV-1a，用于各 rank 一致性，不等同于重启兼容签名。输出调度、目录等可不进入重启签名；数值算法参数必须进入。

## 5. 增加或修改界面重构算法

### 5.1 公共接口

接口位于 `include/wcns/solver/wcns_reconstruction.hpp`：

```cpp
class IReconstructionScheme {
public:
    virtual ~IReconstructionScheme() = default;
    virtual std::string_view name() const noexcept = 0;
    virtual StencilRequirement stencil_requirement() const noexcept = 0;
    virtual Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const = 0;
};
```

当前 registry 强制 `first_offset=-2,point_count=6,halo_width=3`。如果新方法需要其他模板宽度，不能只改这个检查；还必须同时改块 ghost 分配、CGNS 叶块读取、物理 ghost、状态 halo、分区最小宽度、面模板提取、测试和所有 profile 边界闭合。这应作为独立架构阶段处理。

### 5.2 六点模板算法接入步骤

1. 在 `src/solver/wcns_reconstruction.cpp` 或新的 solver 源文件实现`IReconstructionScheme`。算法名只能用小写字母、数字和下划线。
2. `TraceSide::Right` 必须与左侧镜像方向一致；可复用现有 `orient_stencil` 思路，不能简单用同一权重直接作用于未反转模板。
3. 所有输入先检查有限性；`context.scale` 和 `context.parameters` 用于无量纲 smoothness，不能绕过统一 epsilon/scale floor。
4. 在 `ReconstructionRegistry::with_builtins()` 注册工厂：

```cpp
result.register_scheme("my_scheme", [] {
    return std::make_unique<MyScheme>();
});
```

5. 当前 `reconstruct_thermodynamic_face()` 在内部创建 `with_builtins()`，所以只在应用层另建registry 不会生效。要让标准求解器识别新算法，必须把工厂加入该 built-in 注册路径，或先重构求解器使 registry 可注入。
6. `algorithm.reconstruction` 本身保存字符串，不需增加枚举即可选择注册名；若还维护`ReconstructionKind/reconstruction_name()`，则同步增加枚举分支。
7. 算法常数若要从配置输入，按第 4 章扩展 `WcnsParameters`、parser、summary 和签名。

### 5.3 正性与回退契约

重构结果必须通过守恒/原始态转换和 floor。特征重构失败时先回退同算法 primitive，再回退 `linear5`，最后一阶；不要在新算法内部悄悄夹断状态而不记录。若要改变回退链，应同步修改`ReconstructionDiagnostics`、历史输出、重启签名和专项测试。

### 5.4 必需测试

- 常数保持、左右镜像、线性/高阶多项式精度；
- 光滑正弦网格收敛阶和间断非振荡性；
- conservative/primitive/characteristic 三种变量空间；
- 非法特征基和非物理解的回退事件、位置、计数；
- 扭曲网格自由流、物理边界面、原生多块和人工切分面；
- 1/2/4 rank 最终场与回退统计一致。

## 6. 增加 Riemann 求解器

### 6.1 实现接口

`include/wcns/solver/riemann_solver.hpp` 的 `IRiemannSolver::solve()` 接收左右压力原始态、单位法向、气体和 floors，返回**单位面积**守恒通量与谱半径。面积由上层乘入。

```cpp
class MyRiemann final : public wcns::IRiemannSolver {
public:
    std::string_view name() const noexcept override { return "my_riemann"; }

    wcns::RiemannResult solve(
        const wcns::PressurePrimitiveState& left,
        const wcns::PressurePrimitiveState& right,
        wcns::Normal3 unit_normal,
        const wcns::GasModel& gas,
        const wcns::NumericalFloors& floors) const override;
};
```

返回结果必须满足：通量全部有限、谱半径有限且非负、`requested_solver` 为注册名、无回退时`used_solver` 相同且 path 为空。有回退时必须给非 `None` 原因和连续的
`from_solver -> to_solver` 路径，最终到达 `used_solver`。

### 6.2 注册与配置

在 `RiemannSolverRegistry::with_builtins()` 注册。当前无粘和粘性 WCNS solver 构造函数都在内部创建 built-in registry，所以标准入口必须修改这一注册路径，或重构构造函数接受注入的registry。`algorithm.riemann` 是字符串，注册后即可由 schema 选择。

若有新参数：扩展 `RiemannSolverParameters`、验证、summary/restart signature 和配置 parser。低 Mach 预处理不是简单开一个布尔值：需要定义预处理矩阵、适用 Mach、稳态/非定常一致性、时间尺度和回退行为，并建立低 Mach 解析/收敛验收后才可暴露配置。

### 6.3 测试

至少覆盖相同状态物理通量、接触/激波、法向旋转不变性、二维 z 分量、熵修正、非物理中间态回退、强壁面无穿透、周期旋转连接通量变换、自由流和串并行共享面唯一计算。

## 7. 增加算法 profile 或修改度量

仅增加重构/Riemann 不需要新 profile。profile 用于绑定一整套几何和线性算子；新 profile 必须作为不可拆 bundle 设计。

需要审查的入口至少包括：

- `include/wcns/mesh/algorithm_profile.hpp`：profile 与各组件 enum；
- `src/mesh/algorithm_profile.cpp`：期望组件、字符串工厂和禁止混用验证；
- `src/mesh/linear_operators.cpp`：插值、通量导数和单边闭合；
- `src/mesh/high_order_metrics.cpp`：单元坐标、Jacobian、面面积矢量和几何诊断；
- `src/mesh/conservation_weights.cpp`：残差/统计积分权重；
- `src/runtime/structured_partition.cpp`：该模板所需的叶块最小活动方向宽度；
- `src/app/wcns_run.cpp::initialize_partitioned_metrics()`：在原 zone 上计算后切片/打包/分发；
- 面通量 halo 层数、共享面 owner、metric pack/unpack 和 restart signature。

度量约定固定为计算步长 1，`J=partial(x,y,z)/partial(xi,eta,zeta)`，二维有面积、三维有体积意义，必须为正。人工 MPI 切面不能变成几何单边界；高阶度量应先在完整原 CGNS zone 上计算，再切片给运行时叶块。物理边界外没有可用坐标，不得通过伪造坐标把中心模板伸入 ghost。

验收包括笛卡尔精确值、解析扭曲网格收敛、离散几何守恒、正 J、2D 所有 z 分量初始化、3D 三方向、原生连接两侧一致、切分前后逐值一致、pack/unpack 无损和两套旧 profile 回归不变。

## 8. 增加源项模型

当前源项不是多态插件，而是 `SourceModelKind + SourceTermConfig + evaluate switch`。接入步骤：

1. 在 `include/wcns/physics/source_terms.hpp` 增加 enum 和参数存储；复杂参数可建独立结构。
2. 在 `source_model_name()`、`validate_kind()`、`SourceTermConfig::validate/summary/restart_signature()` 增加分支。
3. 在 `src/runtime/case_config.cpp` 增加允许键、字符串解析和数值读取。
4. 在 `SourceTermRegistry::evaluate()` 实现每单位无量纲物理体积的五分量守恒源。
5. 明确源是状态依赖、坐标依赖还是时间依赖；当前函数已提供 `U,x,t,dimension`。
6. 二维必须保证 z 动量源严格为零；能量源要与动量功一致。
7. 新参数和模型顺序进入摘要/重启签名，关闭开关时非零参数必须拒绝。
8. 若模型需要气体、参考量或其他场，扩展 registry 上下文，而不是读取全局变量。

测试源项自身的解析值、多个模型叠加、每个 SSPRK stage time、均匀场解析演化、守恒积分、2D/3D、MPI 等价和 restart 连续性。若是制造源，必须把“初场/解析场”和“源”分别推导，不能因为二者名字相同就假定自动匹配。

## 9. 增加或修改边界条件

边界扩展至少跨越以下层次：

1. `include/wcns/mesh/topology.hpp::BoundaryType` 增加类型；
2. `src/runtime/case_config.cpp` 增加字符串解析、名称输出、动态物理数据键和组合验证；
3. 如需读取 CGNS 原类型，更新 `src/io/cgns_reader.cpp` 的类型映射；注意生产入口最终仍以配置的 default/patch override 为准；
4. 扩展 `BoundaryPhysicalDataConfig` 和运行时 `BoundaryData`；
5. 在 `src/app/wcns_run.cpp::make_boundary_data()` 完成无量纲数据闭合和默认策略；
6. 在 `PhysicalGhostStateOperator` 的 `make_ghost()` 构造三层面 ghost 的`rho,u,v,w,T`，再由同一气体模型计算 p 和守恒量；
7. 在 `apply_inviscid_boundary_face_state()` 定义重构后的无粘真实面强/弱约束；
8. 在 `src/solver/viscous_boundary.cpp` 定义粘性真实面状态与梯度约束；
9. 把新增物理数据写入 summary 和 restart signature；
10. 更新配置模板、CGNS patch 名示例和边界测试。

不要给物理 boundary ghost 生成虚假坐标或度量。需要边界法向时使用真实面 `FaceMetric`；需要高阶单边导数时只使用真实域可用数据和 profile 的单边闭合。边/角 ghost 当前不填，新算法若读取它们必须先提出并实现清晰的多边界优先级，而不是依赖未初始化内存。

周期不是普通物理边界类型扩展。周期必须通过 `GridConnectivity1to1` 进入拓扑和 halo 路径，处理坐标平移、旋转下的速度/动量/梯度/通量变换及反向连接一致性。

边界测试至少包含静止/移动滑移壁、无滑移绝热/等温壁、入流/远场亚声速与超声速特征选择、出口有/无目标态、两侧法向符号、扭曲面法向、三层 ghost、无粘壁通量、粘性剪切/热流、切分 patch 继承和 MPI。

## 10. 增加气体、输运或其他物理模型

### 10.1 当前气体模型

`GasModel` 只表示常 gamma 的热完全理想气体。所有温度/压力/守恒转换集中在`thermodynamics.hpp/.cpp`。增加真实气体、多组分或变比热会影响状态变量、声速、焓、Riemann、特征矩阵、时间步、边界、输出和 checkpoint，不适合只在一个函数中加 switch。

推荐先把 EOS/热力学改成明确的模型接口或封闭 variant，并保持转换函数为唯一入口。新增模型必须有状态域验证、正性 floor、无量纲关系、声速/焓一致性和序列化签名。

### 10.2 当前输运模型

`TransportConfig` 的常黏度和 Sutherland、Prandtl 已完整接入 schema 1、生产装配、配置摘要、
manifest 与重启签名。修改或增加输运模型时：

1. 在 `CaseConfig` 增加 transport 配置与键；
2. 解析 `constant|sutherland` 及参数并验证正值；
3. 把**同一份** config 同时传给 `ViscousWcnsConfig::transport` 和输出
   `QuantityContext::transport`，否则求解黏度与输出黏度不一致；
4. 把 transport summary 纳入真正被 checkpoint 使用的 `CaseConfig::restart_signature()`；
5. 更新 Re/Pr/Sutherland 解析解、时间步稳定性、黏度/壁面输出和配置迁移测试。

若增加湍流或多方程模型，还需扩展流场分量、halo payload、checkpoint 和输出，不应假装为一个等效黏度标量完成。

## 11. 初始化新的算例流场

### 11.1 接入步骤

1. 在 `InitialConditionConfig::validate()` 的 `valid_types` 加入新名字。
2. 如有新参数，把 `initial.<name>` 加入 `fixed_keys()`；否则严格 parser 会在启动时拒绝。
3. 在验证中检查必需参数、正性、维数和参数组合，而不是全部依赖内建默认值。
4. 在 `src/runtime/flow_initializer.cpp` 写返回 `TemperaturePrimitiveState` 的纯函数。
5. 在 `FlowInitializer::evaluate()` 分派新类型。
6. 使用 `pressure_primitive()`、`temperature_primitive()` 和`thermodynamic_conservative()` 做闭合，不复制另一套状态方程。
7. `initialize_block()` 只遍历真实单元；不要初始化 ghost。halo 与物理边界在空间算子前填充。
8. 对周期解析场使用明确的周期距离，不能根据扭曲网格坐标包围盒猜周期。
9. 更新模板、用户手册和一个正式算例。

### 11.2 函数骨架

```cpp
TemperaturePrimitiveState my_initial_state(
    const InitialConditionConfig& config,
    Real x, Real y, Real z,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension)
{
    const Real rho = config.parameter("rho", 1.0);
    const Real u = /* analytic u(x,y,z) */;
    const Real v = /* analytic v(x,y,z) */;
    const Real w = dimension == 3 ? /* analytic w */ : 0.0;
    const Real temperature = /* positive nondimensional T */;
    TemperaturePrimitiveState state {rho, u, v, w, temperature};
    static_cast<void>(pressure_primitive(
        state, gas, reference, floors, dimension));
    return state;
}
```

测试坐标分区边界两侧、恰好位于间断/中心的点、默认参数、非法温度/压力、2D z 分量、多个叶块初始化一致和解析输出。重启会跳过初场写入；因此修改初场不会改变已有 checkpoint 的状态。

## 12. 增加新的流场输出量

### 12.1 接口和描述符

实现 `IFieldQuantity`：

```cpp
class MyFieldQuantity final : public wcns::IFieldQuantity {
public:
    MyFieldQuantity()
    {
        descriptor_.name = "my_quantity";
        descriptor_.location = wcns::TopologyLocation::Cell;
        descriptor_.dimensional_unit = "1";
        descriptor_.scale = wcns::QuantityScale::Dimensionless;
    }

    const wcns::QuantityDescriptor& descriptor() const override
    {
        return descriptor_;
    }

    wcns::Real evaluate_cell(
        const wcns::StructuredBlock& block,
        const wcns::MetricField& metric,
        wcns::Index3 index,
        const wcns::QuantityContext& context) const override
    {
        // 只返回该真实单元上的有限无量纲值。
        return 0.0;
    }

private:
    wcns::QuantityDescriptor descriptor_;
};
```

当前只接受 `Cell` 位置。可选尺度为 Dimensionless、Density、Velocity、Pressure、Temperature、Momentum、Energy、SpecificEnergy、Viscosity、LengthPower。量纲输出时 registry 自动乘尺度。

`descriptor.dependencies` 当前用于检查未知依赖、重复依赖和环，并不会把已计算依赖值传给`evaluate_cell()`；实现仍需从 block/metric/context 自行计算，或以后扩展缓存接口。

### 12.2 注入标准入口

`ProductionFieldWriter` 支持传入 registry，但当前 `wcns_run` 使用默认参数。修改启动装配：

```cpp
auto field_registry = wcns::FieldQuantityRegistry::create_builtin();
field_registry.register_quantity(
    std::make_shared<MyFieldQuantity>());

wcns::ProductionFieldWriter field_writer(
    mpi, config, plan, local_blocks, metrics,
    quantity_context, mesh_name, std::move(field_registry));
```

然后配置 `output.field.quantities` 加入 `my_quantity`。还要把新 `.cpp` 加入 CMake 对应库。

### 12.3 梯度型二级量的特别要求

涡量、散度、Q、壁面剪切等不能在 `evaluate_cell()` 中无条件读取三方向 ghost：物理边界的ghost 坐标/度量/梯度无效，输出时也没有独立的“刷新所有派生 halo”契约。可选设计是：

- 只用真实单元和 profile 单边闭合计算边界附近值；
- 在输出前增加明确的同步/precompute 阶段，把二级量写入有版本的 `TopologyField`；
- 对壁面量建立 face/boundary quantity 接口，而不是冒充 cell-centered 量。

任何方案都要测试原生接口和人工 MPI 切面的值一致，且不可把分区切面误当物理边界。

### 12.4 输出测试

测试唯一名称、重复注册、未知/循环依赖、有限性、无量纲/量纲尺度、CGNS 字段名、Tecplot 列、多 zone 重组、串并行逐值一致和未知配置选择明确失败。

## 13. 增加新的全场统计量

实现 `IStatisticQuantity`：

```cpp
class MyStatistic final : public wcns::IStatisticQuantity {
public:
    MyStatistic()
    {
        descriptor_.name = "my_statistic";
        descriptor_.location = wcns::TopologyLocation::Cell;
        descriptor_.scale = wcns::QuantityScale::Dimensionless;
    }

    const wcns::QuantityDescriptor& descriptor() const override
    {
        return descriptor_;
    }

    wcns::Real evaluate(const wcns::StatisticContext& context) const override
    {
        wcns::Real local = 0.0;
        for (const auto& block : context.local_blocks.blocks()) {
            // 只累计本 rank 拥有的真实单元；需要积分时使用全局守恒权重和 J。
        }
        return context.mpi.sum(local);
    }

private:
    wcns::QuantityDescriptor descriptor_;
};
```

在应用层注入：

```cpp
auto statistics = wcns::StatisticRegistry::create_builtin();
statistics.register_quantity(std::make_shared<MyStatistic>());

wcns::RuntimeOutputManager output(
    mpi, config, plan, checkpoint.mesh_signature(),
    &statistic_context, event_writer, std::move(statistics));
```

每个统计事件所有 rank 都会调用 `evaluate()`，所以实现中集合通信的次数和顺序必须与名称选择无关并在所有 rank 一致。只遍历 owned real cells，不能累计 ghost。体积积分使用`GlobalConservationWeights` 与对应 `MetricField::jacobian()`，避免原 zone 物理边界闭合权重或人工切分重复计数。

当前 registry 在 `output.dimensional=true` 时还会自动乘 descriptor 的场尺度和`L_ref^dimension`。若统计量已经在实现中自行恢复量纲，descriptor 应保持 Dimensionless，否则会二次缩放；更推荐实现始终返回无量纲值，由 descriptor 统一缩放。

测试解析积分、网格加密、原生多 zone/单 zone 切分等价、1/多 rank 归约、量纲尺度、重复名、未知选择、输出列顺序和 restart 前后时间序列拼接。

## 14. 增加新的输出格式或观察量类别

字段格式位于 `ProductionFieldWriter::write()`，历史/统计和 manifest 位于 `RuntimeOutputManager`。新增格式需：

1. 扩展配置枚举、parser、名称、summary 和模板；
2. 明确 cell/face 位置、zone 顺序、坐标、单位和时间元数据；
3. 先写同目录 `.tmp`，完整关闭后再按 `allow_existing` 原子提交；
4. 只由确定 rank 写文件，但所有 rank 仍以一致顺序参加 gather/reduction；
5. 返回最终文件路径，让 manifest 记录；
6. 写独立重读器测试，不只检查“文件存在”；
7. 保持原 zone 重组，不把 MPI 叶块布局泄漏成用户物理网格，除非格式明确声明为分区输出。

若新增随时间统计类别，优先复用 `OutputSchedule` 和 observer；不要把 I/O 直接塞进求解器子步。

## 15. 修改停止机制或时间推进

停止判据在 `StopController`，时间步裁剪和 observer 通知在 `SimulationDriver`，SSPRK3 状态更新在 solver/time integrator。新增停止条件必须定义：

- steady/unsteady 哪些模式生效；
- 与数值失败、正常完成、信号、墙钟、最大步数的优先级；
- 是否要求安全 checkpoint；
- manifest 名称和进程退出码；
- MPI 如何归约为全局一致决定；
- restart 需要保存哪些内部计数。

修改时间推进需重新核对源项 stage time、输出精确时间事件、t_end 不越界、checkpoint dt、粘性稳定限制和残差定义。隐式、本地时间步或双时间步不能复用“非定常物理时间”语义而不区分内外迭代。

## 16. 修改 CGNS 读取、连接或自动分区

### 16.1 CGNS 读取

新增 CGNS 特性时，元数据和实际块读取必须使用同一验证。所有 `cgsize_t` 转换到 `int` 前检查范围；所有 Vertex PointRange 转换到 0-based 后检查 `vertex_extent`，再推导 cell/face range 并检查 `cell_extent`。不要信任供体名、transform 或边界定位。

### 16.2 连接

连接需验证接收/供体计数、轴置换、方向、互逆性和唯一覆盖。周期刚体变换必须是有限的正确旋转/平移，并在状态 halo、梯度和面通量上用相符变换。消息 tag 必须确定且无冲突。

### 16.3 分区

自动分区当前是按单元数的确定性递归二分。改进负载模型时必须保持：同一元数据在所有 rank产生相同 digest、叶块不重叠且完整覆盖原 zone、物理 patch 正确切片、跨原连接正确切片、兄弟连接成对、每个活动方向满足 profile 下限、rank owner 稳定。

分区变化不应改变原 zone 上的度量、守恒积分和最终重组字段。用同一逻辑网格验证：单 zone 运行时切分、原生多 zone、不同 rank 数最终结果一致。

## 17. 检查点格式与兼容性开发

检查点版本在 `src/runtime/checkpoint.cpp`。新增持久状态时：

1. 明确它是否可由恢复后的守恒量重建；不可重建则写入文件；
2. 增加版本或向后兼容读取规则，绝不把缺失字段默认为看似合理的零；
3. 数值模型变化加入 restart signature；
4. 原 zone 存储，不绑定当前叶块/rank；
5. 保持无量纲守恒状态为恢复权威；
6. 测试 1→2、2→1、2→4 rank 和不同合法分区；
7. 测试网格、profile、气体、边界、源项、模型参数不匹配时在推进前失败；
8. 连续 N 步与 K 步 checkpoint + N-K 步恢复逐场比较。

当前 mesh signature 主要覆盖 base/zone 名、维数、extent 和坐标，不应把它当成所有 CGNS 拓扑语义的密码学证明。若扩展连接/边界语义，考虑把规范化拓扑也纳入签名。

## 18. 测试如何加入工程

### 18.1 普通单元测试

1. 在 `tests/test_<feature>.cpp` 写一个清晰命名的入口函数；
2. 在测试函数附近用简短注释说明它验证的目的、输入和主要失败条件；
3. 把源文件加入 `tests/CMakeLists.txt` 的 `wcns_unit_tests`；
4. 在 `tests/test_main.cpp` 声明并调用入口；
5. 保证测试没有外部顺序依赖、不会写源码目录、失败时抛出包含上下文的异常。

### 18.2 独立/MPI 测试

需要独立进程、CGNS fixture 或 MPI 的测试单独建 executable/add_test。MPI 测试至少运行 2 rank，对可能死锁的失败路径设置合理 TIMEOUT。临时目录必须位于构建树；清理脚本只允许删除带自身标记的目录。

### 18.3 生产入口算例测试

最终验收必须调用 `wcns_run --config`，不能用简化测试 solver 代替。驱动器应保存每条命令、退出码、stdout/stderr、配置、网格、manifest 和独立验证结果。比较时排除 wall time 等天然非确定列，但不能排除物理场、step/time、残差和回退计数。

### 18.4 常用命令

```powershell
cmake --build build-dev-serial --parallel 4
ctest --test-dir build-dev-serial -R wcns.unit --output-on-failure
ctest --test-dir build-dev-serial --output-on-failure

cmake --build build-dev-mpi --parallel 4
ctest --test-dir build-dev-mpi --output-on-failure
```

提交前再做全新 Release 构建，避免旧对象或 CMake cache 隐藏依赖遗漏。

## 19. 自定义扩展审查清单

### 配置和文档

- [ ] 名称、默认值、单位、范围和互斥条件已定义；未知/重复键仍失败。
- [ ] 生产代码实际消费该值；summary 和日志可看见最终选择。
- [ ] 影响数值轨迹的值进入 restart signature。
- [ ] 完整模板、用户手册、算法文档和已知限制同步。

### 数值和物理

- [ ] 推导与无量纲方程一致；2D/3D 分量和单位一致。
- [ ] 常数保持、解析值、精度/收敛和正性已测试。
- [ ] 回退、floor 或限制器全部可诊断，没有静默夹断。
- [ ] 物理边界只读取可用数据，边/角 ghost 未被意外访问。

### 多块和 MPI

- [ ] 原生连接、周期连接、人工切分面和同 rank/跨 rank 均覆盖。
- [ ] 共享面只计算一次，旋转/方向变换正确。
- [ ] 所有集合通信顺序一致，失败路径不会死锁。
- [ ] 1/多 rank 结果及诊断在规定容差内等价。

### 输出和重启

- [ ] 新字段/统计量可选择、量纲缩放正确、未知名失败。
- [ ] 文件经独立读回验证，原子提交和 manifest 记录完整。
- [ ] restart 连续性通过，不兼容配置在推进前拒绝。
- [ ] 不把 MPI 叶块布局错误地当成物理原 zone。

### 工程

- [ ] 串行/MPI、Debug/Release 中相关目标编译无警告升级。
- [ ] 单元、失败、端到端和人工算例均有记录。
- [ ] Git 提交只包含本小目标；没有覆盖用户未提交文件。

## 20. 常见开发错误

- 只把名字加入 parser，忘记 built-in registry，导致配置能读但 solver 报 unknown。
- 只在一个 rank 上注册统计量，其他 rank 进入不同集合通信顺序。
- 新重构声称需要更宽模板，却仍用三层 ghost 和旧分区下限。
- 在人工 MPI 切面两侧各算一次 Riemann 通量，破坏逐位守恒和串并行等价。
- 为物理 ghost 推造坐标/度量，让高阶中心差分越过真实边界。
- 给粘性 solver 改了 TransportConfig，却让 `viscosity` 输出仍使用默认模型。
- 新数值参数没有进入 restart signature，旧 checkpoint 被错误接受。
- 在输出 quantity 中读取“上次残留”的 halo 或派生场，没有版本和预输出同步。
- 用 `maximum_steps` 退出当作定常收敛，或只看 total L2 不看五分量判据。
- 测试只调用底层函数，没有通过正式配置、生产装配、输出和重启路径。

按照本指南完成扩展后，仍应以具体算例报告作为最终证据。报告必须公开说明对程序源码的所有修改，不得把为算例添加的算法、边界、初场或验证逻辑隐藏在结果目录中。

## 21. 实战一：把 MDCD 的 `disp` 和 `diss` 暴露给配置文件

本节记录已经落地的第一次用户自定义开发。对应提交为 `cafdfc3`。目标不是简单增加两个
parser 键，而是让它们沿“配置文本 → 类型化对象 → 数值核 → 摘要/重启签名 → 测试 → 用户
文档”完整闭环。

### 21.1 先固定外部契约

最终配置写法为：

```text
algorithm.mdcd.disp = 0.0463783
algorithm.mdcd.diss = 0.01
```

二者是可选实数键。省略时保持以上默认值；`disp` 映射 `WcnsParameters::mdcd_dispersion`，`diss` 映射 `WcnsParameters::mdcd_dissipation`。它们共同影响 `mdcd_linear` 和 `mdcd_hybrid` 中使用的MDCD线性候选。选择 WENO、`linear5` 或 `zero_order` 时不会消费其结果，但任何显式写出的非法值仍会在启动阶段失败，不能让“休眠配置”携带无效数据。

当前合法域为

\[
\mathrm{disp}>0,\qquad
0\le\mathrm{diss}<\mathrm{disp},\qquad
3\,\mathrm{disp}+9\,\mathrm{diss}<1,
\]

且二者都必须有限。最后一个不等式保证实现中相应系数组合留在已定义的稳定参数区域。不要在 parser 中静默截断、取绝对值或回退默认值；用户输入错误必须携带参数上下文明确失败。

### 21.2 数据结构为什么放在 `WcnsParameters`

数值参数已经由 `ReconstructionContext::parameters` 统一送入每个`IReconstructionScheme::reconstruct_scalar`。因此扩展已有`include/wcns/solver/wcns_reconstruction.hpp` 中的 `WcnsParameters`，可以避免给 MDCD 另开一条旁路，也保证单元测试可直接构造不同参数。配置层的 `ReconstructionConfig::nonlinear` 复用同一结构；正式装配无需复制两份默认值。

设计检查顺序如下：

1. 参数有唯一内部所有者，不在 parser、MDCD-LINEAR 和 MDCD-HYBRID 三处各存一份；
2. 默认值在 C++ 结构中初始化，配置只覆盖，不依赖配置模板才能工作；
3. `WcnsParameters::validate()` 是数值层最终防线，直接调用底层 API 也无法绕过验证；
4. 配置摘要和重启签名从最终类型化参数生成，而不是回显未经解析的文本。

### 21.3 parser 的实际修改路线

在 `src/runtime/case_config.cpp` 中进行了以下动作：

1. 把两个精确键加入固定 schema 集合；否则严格配置会先报 unknown key；
2. 在解析完重构名称和变量空间后调用 `optional_real`；
3. 以 `result.reconstruction.nonlinear` 当前值作为缺省值；
4. 继续走 `CaseConfig::validate()` → `InviscidWcnsConfig::validate()` →  `WcnsParameters::validate()`，不在 parser 重复数学规则。

核心数据流等价于：

```cpp
parameters.mdcd_dispersion = optional_real(
    entries, "algorithm.mdcd.disp", parameters.mdcd_dispersion);
parameters.mdcd_dissipation = optional_real(
    entries, "algorithm.mdcd.diss", parameters.mdcd_dissipation);
```

如果以后增加另一项可选重构实数，照搬代码前必须回答三个问题：它是否应对所有重构验证；默认值是否属于算法定义；改变它是否会改变 checkpoint 后续数值轨迹。第三个问题为“是”时必须进入重启签名。

### 21.4 数值核、摘要和重启兼容

`src/solver/wcns_reconstruction.cpp` 中 MDCD 的两条实现都从传入参数读取最终值，不再依赖局部硬编码。`WcnsParameters::validate()` 同时检查有限性和上述三个不等式。

`ReconstructionConfig::summary()` 输出最终 `mdcd_dispersion`、`mdcd_dissipation`；`restart_signature()` 也包含它们。结果是：

- 日志能复核实际采用的值；
- 输入省略时记录默认值，而不是记录“空”；
- 修改参数后，旧 checkpoint 即使网格和其他算法相同也会被拒绝；
- MPI 各 rank 对完整配置 digest 做一致性检查，不能出现各 rank 系数不同。

### 21.5 如何验证这项开发

`tests/test_case_config.cpp` 覆盖以下行为：

- `disp=0.04,diss=0.005` 能准确读入；
- summary/restart signature 能看到最终值；
- 令 `diss==disp` 会确定性失败；
- 原有不写新键的生产配置继续使用默认值。

`tests/test_inviscid_reconstruction.cpp` 直接改变 `WcnsParameters`，确认 MDCD 数值核实际消费参数，而不是仅由 parser 接受但计算忽略。复现最小测试：

```powershell
cmake --build build-dev-serial --parallel 4
ctest --test-dir build-dev-serial -R wcns.unit --output-on-failure
```

人工配置测试还应做三次独立运行：默认参数、一个合法自定义组合、一个非法组合。前两次比较启动摘要及重构结果差异；第三次必须在读初场/推进前失败。不要只看“配置能读”。

### 21.6 同步的文档位置

- `examples/full_case_template.wcns`：给出可取消注释的两个键和默认值；
- `docs/user-manual.md`：给出含义、合法域和重启影响；
- `docs/stage-l-design.md`、`算法补充.md`：同步算法符号及约束；
- 本指南：保留真实开发路线，供新增数值参数时复用。

## 22. 实战二：加入保持六点契约的零阶重构

这项开发也在提交 `cafdfc3` 中完成。需求中特别强调“兼容六点模板而不缩小模板”，因此它
不是改变整个残差算子的 halo 宽度，而是新增一个接收相同六点输入、只使用中心相邻两点的
策略。

### 22.1 固定模板编号，避免左右值错一位

对面 \(j+1/2\)，通用调用者提供：

```text
数组下标       0       1       2       3       4       5
物理点       q[j-2]  q[j-1]   q[j]   q[j+1] q[j+2] q[j+3]
```

零阶结果严格定义为

\[
q^L_{j+1/2}=q_j=\mathrm{stencil}[2],\qquad
q^R_{j+1/2}=q_{j+1}=\mathrm{stencil}[3].
\]

这里“第三/第四点”是自然语言的一基编号，对应 C++ 下标 2/3。它仍报告
`StencilRequirement{first_offset=-2, point_count=6, halo_width=3}`，因此：

- 不改变块最小宽度和自动分区约束；
- 不改变物理边界三层 ghost 约定；
- 不为一个低阶策略产生第二套残差/通信调用链；
- 可以与其他六点策略在同一测试和诊断框架内替换。

### 22.2 实际增加一个策略需要改哪些位置

1. 在 `ReconstructionKind` 增加 `ZeroOrder`；
2. 在 `src/solver/wcns_reconstruction.cpp` 实现 `ZeroOrderScheme`；
3. `name()` 返回唯一规范名 `zero_order`；
4. `stencil_requirement()` 返回通用六点要求；
5. `reconstruct_scalar()` 先调用 `checked_stencil`，确保长度、有限性等共同前置契约成立；
6. 在 `ReconstructionRegistry::with_builtins()` 注册工厂；
7. 在 `reconstruction_name()` 补齐枚举到字符串映射；
8. 在配置模板、用户手册、算法文档和顶层能力说明加入第六种重构。

实现核心保持刻意简单：

```cpp
const auto values = checked_stencil(stencil);
return side == TraceSide::Left ? values[2] : values[3];
```

没有另加 limiter，也不把它伪装成一阶外推。它是分片常数的左右单元迹值；Riemann 求解器仍负责由这两个状态产生面通量。

### 22.3 注册表和生产装配如何连通

配置中的

```text
algorithm.reconstruction = zero_order
```

先作为规范字符串进入 `ReconstructionConfig`，生产装配再从 built-in registry 创建策略。未知拼写必须由 registry 报错，不能回退到 WENO 或 linear。通用面循环只依赖`IReconstructionScheme`，因此没有在 Euler 残差中增加 `if (zero_order)` 分支。

降阶诊断把 `zero_order` 视作用户主动选择的线性/低阶基线，而不是“高阶失败后又回退到linear5”。否则每个面都会产生虚假的降阶计数。真正的非法热力学状态和 Riemann 回退仍按公共诊断记录。

### 22.4 测试为什么必须使用非对称数据

`tests/test_inviscid_reconstruction.cpp` 用左右不对称的六个值创建 `zero_order`，验证：

- registry 能按 `zero_order` 创建；
- `point_count==6`；
- 左值精确等于下标 2；
- 右值精确等于下标 3；
- 枚举名称映射正确；
- 错误模板长度仍失败。

若用常数模板，取错任何点都可能通过；若用对称模板，也可能掩盖左右反转。因此新增模板算法时，应专门构造“每个位置可辨识”的数据，例如 `{11,22,33,44,55,66}`。

端到端 x-z 统计 smoke 配置也选择 `zero_order`，这额外证明正式配置、生产 registry、完整残差装配和输出路径都接受它。这个 smoke 的目的不是证明零阶有高精度，而是防止“单元测试注册了、生产程序没注册”的断链。

### 22.5 用户何时可以使用它

`zero_order` 适合：

- 检查 Riemann、边界和通信而尽量去掉高阶插值影响；
- 做强耗散基线；
- 构造端到端装配测试。

它不适合用来宣称 WCNS 的光滑高阶精度，也不应作为双马赫反射等正式结果的默认算法。切换算法后必须换输出目录并保留配置/manifest，不能覆盖高阶结果。

## 23. 实战三：可开关、多截面的 x-z 平面监测

主体实现提交为 `eb56dcd`，正式入口端到端测试提交为 `93bb56a`。需求面向泊肃叶流和槽道湍流：按若干 y 层持续监测流向平均速度以及 x 向质量通量密度在 x-z 平面上的积分。

### 23.1 先澄清“截面”和“流量”定义

本功能选择原始 CGNS zone 的一个 cell-centered J 层。每个选中层输出：

\[
\overline{u}_{xz}(j)=
\frac{\sum_{c\in j}A_c u_c}{\sum_{c\in j}A_c},
\qquad
M_x(j)=\sum_{c\in j}A_c\rho_cu_c.
\]

面积近似取该单元上下两个 J 面面积的算术平均：

\[
A_c=\frac12\left(|\boldsymbol S_{j-1/2}|+|\boldsymbol S_{j+1/2}|\right).
\]

在平直笛卡尔通道中这就是单元对应的 x-z 面积。`M_x` 的 SI 量纲是 kg/s，但 x-z 面的法向是 y，所以它不是“穿过 x-z 面的法向质量流率”；它是壁平行层上对流向质量通量密度`rho*u` 的面积积分。这个量适合比较不同 y 层、不同时间或不同网格的流向输运强度。若要真正测量入口 y-z 截面的穿面流量，应另建 y-z 统计并使用面法向通量，不能改标签不改公式。

### 23.2 外部配置契约

在普通统计配置之后增加：

```text
output.statistics.enabled = true
output.statistics.quantities = total_mass,total_energy
output.statistics.xz_planes.enabled = true
output.statistics.xz_planes.cell_j_indices = 0,11,23,35,47
```

索引是原始 CGNS zone 的零基单元 J 索引，不是顶点编号、CGNS 一基编号、MPI rank 局部编号或自动切分后的叶块编号。打开子开关时：

- 父级 `output.statistics.enabled` 必须为 true；
- 列表不能为空；
- 每项必须非负且唯一；
- 每个原 zone 必须是三维；
- 每个原 zone 都必须包含所有请求的 J 索引；
- 同一索引在全域的单元中心 y 坐标必须共面到实现容差。

每个 N 自动追加两个规范列名：

```text
xz_mean_u_jN
xz_mass_flow_x_jN
```

用户不需要、也不应在 `output.statistics.quantities` 重复写这些名字；自动追加后重复选择会由registry 校验拒绝。关闭子开关时不注册、不计算、不产生集合通信，也不改变原统计文件。

### 23.3 配置数据结构和启动前验证

`include/wcns/runtime/case_config.hpp` 增加 `XzPlaneStatisticsConfig`，并作为`OutputConfig::xz_planes` 的成员。它保存 `enabled` 和 `cell_j_indices`，负责父开关、空列表、非负和重复值等不依赖网格的验证。

`src/runtime/case_config.cpp` 完成四件事：

1. 把两个键加入严格 schema；
2. 用 `optional_bool` 和 `optional_integer_list` 解析；
3. 开启时调用 `xz_plane_statistic_names` 自动补列；
4. 在 output summary 中打印开关及索引列表。

网格范围不能在纯文本解析时验证，因为这时尚未读取 CGNS。`src/app/wcns_run.cpp` 在建立`StructuredPartitionPlan` 后、分配场数组及推进前调用 `validate_xz_plane_statistics`。这保证越界和二维误用尽早失败，同时验证依据是原 zone 元数据而不是某个 rank 的局部叶块。

### 23.4 从原 zone 索引映射到 MPI 叶块

自动分区可能把一个原 zone 切成多个 `PartitionLeaf`。统计不能假设“一 zone 等于一 block”。
对请求的 `target_j`，每个 rank 对本地块执行：

1. 由 `block.id()` 查到对应 leaf；
2. 判断 `target_j` 是否在半开区间`[leaf.cells.begin.j, leaf.cells.end.j)`；
3. 若不在，跳过；若在，令`local_j = target_j - leaf.cells.begin.j`；
4. 遍历该本地层的 i/k 单元，累加面积、面积乘速度或面积乘 `rho*u`；
5. 所有 rank 按完全相同顺序做 `sum/min/max` 集合通信。

这样一个原 zone J 层即使被 i、j、k 多次切分也只覆盖一次。不能直接用本地 `j=N`：当叶块从原 zone 的 J=24 开始时，本地 j=0 对应的全局层是 24。

### 23.5 几何和物理合法性检查

实现从真实单元的守恒量读取 `rho`、`rho*u`，从 cell coordinates 与两个真实 J 面度量读取 y 和面积。它检查：

- 面积有限且严格为正；
- 密度有限且高于公共 density floor；
- 速度或 `rho*u` 有限；
- 请求层在 MPI 汇总后具有非零全局面积；
- 全局 `y_max-y_min` 不超过`1e-10*(1+max(abs(y_min),abs(y_max)))`。

最后一项意味着该功能当前只支持真正平面的 x-z 层。一般扭曲网格的同一 J 层若 y 随 x/z 变化，会被明确拒绝，而不是把曲面结果误标成平面。若未来要支持曲面平均，应新增不同名称、用一致曲面面积定义并记录几何语义，不应放宽容差偷换概念。

### 23.6 统计 registry 和量纲恢复

`src/runtime/quantity_registry.cpp` 提供三个公共入口：

- `xz_plane_statistic_names`：确定性产生列名；
- `validate_xz_plane_statistics`：做原 zone/索引验证；
- `register_xz_plane_statistics`：按索引注册两个闭包。

`src/app/wcns_run.cpp` 先创建 built-in `StatisticRegistry`，再在开关开启时注册平面量，最后交给`RuntimeOutputManager`。所有 rank 都依据广播后的同一配置注册，保证集合通信调用序列一致。

这次扩展还在 `QuantityDescriptor` 增加 `integration_length_power`。原因是不能仅凭`QuantityScale::Momentum` 推断一个积分过多少个空间维度：

- `xz_mean_u_jN` 恢复量纲时只乘 `U_ref`，长度幂为 0；
- `xz_mass_flow_x_jN` 乘 `rho_ref*U_ref*L_ref^2`，长度幂为 2。

如果遗漏长度平方，无量纲输出可能看似正确，但 `output.dimensional=true` 会系统性错误。新增面积、体积、线或面统计时，都要明确基础物理尺度和积分长度幂。

### 23.7 测试分层

`tests/test_case_config.cpp` 验证开关、列表、自动列名、重复索引和父开关关系。
`tests/test_quantity_registry.cpp` 构造解析均匀场，验证名称、公式、合法范围和二维拒绝。

串行端到端测试由三步 CTest fixture 组成：

1. `wcns.generate_xz_plane_statistics_mesh` 生成 8×8×8 CGNS；
2. `wcns.run.xz_plane_statistics.smoke.serial` 用正式 `wcns_run` 和 `zero_order` 推进一步；
3. `wcns.check.xz_plane_statistics.serial` 用
   `tests/check_xz_statistics.py` 读取真实 statistics 文件，检查初/终记录的列和 0.2 解析值。

MPI 构建另以独立输出目录执行 `wcns.run.xz_plane_statistics.smoke.2` 和 `wcns.check.xz_plane_statistics.2`。两个 rank 会切分同一原 zone，检查器仍要求面积平均和质量积分为 0.2，从而直接覆盖原 zone J 索引到叶块局部索引的映射及 `MPI_Allreduce` 路径。

运行：

```powershell
ctest --test-dir build-dev-serial -R "xz_plane_statistics|xz_statistics" --output-on-failure
```

端到端测试很重要，因为只有它能同时发现“parser 已加但生产 registry 未注册”“列有了但调度没写初/终值”“无量纲值对但输出量纲错”等跨模块问题。

### 23.8 在泊肃叶算例中的实际落点

`cases/manual/case02_3d_poiseuille/uniform_36x48x36.wcns` 和 `wall_clustered_36x48x36.wcns` 已开启 J=`0,11,23,35,47`。这五层用于观察近壁至中心的速度和流向输运分布。

必须注意：case02 目录已有的长时结果生成于本功能之前，旧 statistics 文件没有这些列。配置升级不等于旧结果自动升级；只有重新运行并保存新 manifest/统计文件后，才能验收截面监测。该事实已明确写进 case02 报告，避免文档和历史证据冲突。

### 23.9 扩展成其他截面的正确步骤

若要增加 y-z 入口流率，不应复制名称后仍用 J 面面积。应：

1. 定义原 zone I 面或 cell-I 层的索引语义；
2. 根据物理问题决定用单元层近似还是直接用真实 I 面通量；
3. 真正穿面流率使用 `rho*(u dot n)*A` 并固定法向符号；
4. 处理轴交换、反向和多 zone 连接；
5. 给它独立的配置前缀与列名；
6. 明确平面/曲面的支持范围；
7. 做串行、自动切分、多 rank 和量纲恢复测试。

## 24. 实战四：经典双马赫反射的初场与特殊边界

实现提交为 `0f4c0fa`。这项扩展横跨 physics、初始化、配置、边界、时间推进、CGNS 算例和测试；只加入一个初场函数而使用普通 outflow 边界并不能构成经典双马赫反射问题。

### 24.1 把经典模型做成单一真源

新增 `include/wcns/physics/double_mach_reflection.hpp` 和 `src/physics/double_mach_reflection.cpp`。`DoubleMachReflection` 统一提供：

- `shock_foot()`；
- `shock_x(y,time)`；
- `upstream_state()`；
- `post_shock_state()`；
- `exact_state(x,y,time)`，仅作为初场/外边界给定；
- `restart_signature()`；
- `validate(gamma,dimension)`。

它冻结经典域中使用的 Mach 10、60 度和两侧状态：

\[
U_0=(1.4,0,0,0,1),\quad
U_1=(8,8.25\sqrt3/2,-4.125,0,116.5),
\]

并使用

\[
x_s=x_0+y/\sqrt3+20t/\sqrt3.
\]

把公式放在单一 physics 类中，避免初场、顶边界和测试各复制一套常数后发生漂移。模型要求二维和 `gamma=1.4`；未来若要任意 Mach/角度，必须从激波关系推导相容状态，不能简单开放六个独立配置数值。

### 24.2 将新初场接入正式初始化链

`src/runtime/flow_initializer.cpp` 在 `initial.type == "double_mach_reflection"` 时：

1. 从 `initial.x0` 取激波足，省略时为 1/6；
2. 建立并验证 `DoubleMachReflection`；
3. 使用每个真实单元已有的 cell coordinates；
4. 在 `t=0` 调用 `exact_state` 选择压力原始量；
5. 通过统一气体模型转换温度原始量和守恒量；
6. 不填物理 ghost，ghost 留给边界算子统一处理。

配置层允许 `initial.*` 动态实数参数，因此不需要为 `x0` 另造固定 parser 键，但必须在模型创建处给出相同默认值，并让它进入 summary/restart signature。

### 24.3 新边界类型和配置一致性

`BoundaryType` 增加 `DoubleMachReflection`，外部规范字符串为`double_mach_reflection`。`CaseConfig::validate()` 建立双向约束：

- 专用初场至少要有一个同名专用边界；
- 出现专用边界时，初场也必须是同名模型；
- 算例必须无粘且关闭源项；
- gamma 必须为 1.4；
- 边界不能携带普通 inflow 目标态、壁温或非零壁速度。

正式算例进一步把 `left/bottom/top` 三面全部配置为专用边界。代码保留“至少一个”的通用解析检查，是因为纯配置阶段还没读到真实 CGNS patch 集合；读取网格后，patch 名和方向检查会阻止错误布置推进。人工 case 验收仍必须逐一确认三面和右 outflow，不能把最小 parser 约束当作完整物理证明。

### 24.4 边界数据如何装配

`BoundaryData` 增加可选 `DoubleMachReflection` 对象。`src/app/wcns_run.cpp` 在遍历本地物理patch 时，对专用类型用配置中同一个 `x0` 建模；`BoundaryData::validate()` 保证只有专用类型携带该对象，普通边界不能意外读取它。

三种边界操作为：

1. i-lower：固定激波后态；
2. j-upper：按真实面中心 `(x,y)` 和当前时间，从移动激波线选择前/后态；
3. j-lower 且 `x<x0`：固定激波后态；
4. j-lower 且 `x>=x0`：调用公共 stationary slip-wall 镜像；
5. 其他轴/侧：明确抛出配置错误。

右边界不使用专用模型，继承 `outflow`。周期边界仍必须由 connectivity 表达，不能伪装成物理 patch。

### 24.5 为什么必须把 SSPRK 子步时间送到边界

顶部激波位置随时间改变。如果只在完整步开始更新边界，SSPRK3 的中间 stage 会使用过期边界，时间精度和冲击位置都不一致。因此本次修改把 `Real time` 沿以下调用链传递：

```text
SSPRK stage time
  -> PhysicalGhostStateOperator::fill(..., time)
  -> compute_inviscid_face_fluxes(..., time)
  -> apply_inviscid_boundary_face_state(..., face_coordinates, time)
```

每个 stage 的 ghost 和重构后的强面状态看到同一个物理时间。时间必须有限且非负；这个检查在物理 ghost 填充入口完成。以后增加任何移动壁、脉动入口或时变远场，都必须复用 stage时间路径，不能从全局变量猜测当前时间。

### 24.6 严守“物理 ghost 只有物理量有效”

`boundary_face_coordinates(block,patch,face)` 只平均真实边界面顶点：二维平均两个切向顶点，三维平均四个。它不外推 ghost 坐标。计算法向仍使用真实边界面度量。

对每个真实边界面和三层 ghost：

1. 找到对应镜像内部真实单元；
2. 用专用状态或公共滑移壁逻辑产生 `rho,u,v,w,T`；
3. 由相同气体模型得到 `p` 和五分量守恒量；
4. 只写这三组物理状态数组；
5. 不写、也不读取 ghost 坐标、Jacobian、面度量、梯度或其他二级量；
6. 不填边/角 ghost，因为当前方向分裂六点面模板不会访问它们。

这与项目已冻结的物理边界契约一致。若未来算法需要斜向二维模板而访问角 ghost，应先重新
设计并验证边角边界合成规则，不能让当前算例偷偷读取未定义数据。

### 24.7 重构后强边界处理

三层 ghost 先参与与内点相同的六点重构。随后当`InviscidBoundaryOptions::strong_boundary_face_state=true`（当前生产默认且配置不可关闭）时，真实边界面外侧迹值再次按专用规则约束：左/顶/底部入流段给定状态，底部壁段做反射迹值。

这样 ghost 影响高阶插值，而最终 Riemann 面状态仍严格满足经典边界。如果未来公开“是否强约束”开关，必须将开关加入 schema、summary、restart signature，并分别验证开/关模式；本次开发没有暗中增加未完成的用户开关。

### 24.8 重启和配置签名

`double_mach_reflection_v1` 及 `x0` 进入 `CaseConfig::restart_signature()`；边界类型本来也进入签名。改变激波足或把边界换成普通类型会拒绝旧 checkpoint。非定常目标时间和输出调度不一定属于数值状态兼容签名，但恢复后的边界在每个 stage 使用恢复物理时间，因此移动激波位置可以连续。

### 24.9 测试如何覆盖数学和调用链

`tests/test_double_mach_reflection.cpp` 覆盖：

- 两侧经典状态常数；
- 激波线在 y/t 方向的位置；
- 初场分类；
- 左、上、下分段边界；
- 三层 ghost；
- 改变时间后顶部 ghost 随激波移动；
- 非法 gamma/维数/边界方向失败。

`tests/data/double_mach_reflection_case.wcns.in` 与 CTest `wcns.run.double_mach_reflection.smoke.serial` 使用正式 `wcns_run`、生成的 CGNS、多块装配和一步推进，防止只测试 physics 小类而漏接生产入口。

手工复现局部验收：

```powershell
ctest --test-dir build-dev-serial -R "double_mach_reflection|wcns.unit" --output-on-failure
```

### 24.10 正式人工算例

`cases/manual/case04_2d_double_mach_reflection/` 包含：

- `double_mach_960x240.wcns`：4×1 域、960×240、WENO-Z characteristic/HLLC、`t=0.2`；
- `run_case04.py`：生成网格、正式 dry-run 和 MPI 运行；
- `README.md`：公式、逐项配置、命令、输出和验收判据。

该目录没有预制完整结果。开发完成与数值算例通过是两级证据：单元/端到端 smoke 证明程序路径连通，960×240 完整周期运行和可视化/网格/rank 对比才是后续人工物理验收。

### 24.11 将来增加另一种“初场 + 特殊边界”模型的模板

按以下顺序工作可避免遗漏：

1. 在 physics 模块写单一模型类，统一常数、公式、合法域和签名；
2. 在 `FlowInitializer` 只初始化真实单元；
3. 增加边界枚举、字符串解析和名称输出；
4. 扩展 `BoundaryData`，保证类型和数据一一对应；
5. 明确每个允许 face axis/side 的 ghost 规则；
6. 需要坐标时只用真实面坐标；需要时间时传 SSPRK stage 时间；
7. 明确重构后强约束规则；
8. 让所有改变数值轨迹的模型参数进入 restart signature；
9. 加模型单测、边界单测、配置失败测试和生产入口 CGNS smoke；
10. 建独立人工 case，公开所有源码修改、配置、命令、输出与未完成验收项；
11. 同步用户手册、算法文档、模板、已知限制和顶层 README。

## 25. 四项实战开发的追踪矩阵与复现卡口

| 需求 | 主要生产文件 | 配置入口 | 核心自动证据 | 实际提交 |
|---|---|---|---|---|
| MDCD `disp/diss` | `wcns_reconstruction.*`、`case_config.*` | `algorithm.mdcd.*` | case-config + reconstruction unit | `cafdfc3` |
| 六点零阶重构 | `wcns_reconstruction.*` | `algorithm.reconstruction=zero_order` | registry/stencil unit + x-z smoke | `cafdfc3` |
| 多 x-z 截面监测 | `case_config.*`、`quantity_registry.*`、`wcns_run.cpp` | `output.statistics.xz_planes.*` | unit + 正式输出三段 CTest | `eb56dcd`、`93bb56a` |
| 经典双马赫反射 | `double_mach_reflection.*`、initializer、boundary、driver | `initial.type` 与同名 boundary | physics/boundary unit + 正式 CGNS smoke | `0f4c0fa` |

建议最终一次性执行：

```powershell
cmake --build build-dev-serial --parallel 4
ctest --test-dir build-dev-serial --output-on-failure

cmake --build build-dev-mpi --parallel 4
ctest --test-dir build-dev-mpi --output-on-failure
```

验收不能只看最后一行“100% passed”。还要保留 CMake 配置、编译器、MPI 实现、Git 提交、CTest 完整输出，并逐项检查：

- 改 MDCD 参数会改变摘要和重启签名，非法组合启动即失败；
- `zero_order` 保持六点/三层 halo 契约，左右确取下标 2/3；
- x-z 开关关闭零开销，开启时多层列名、解析值、MPI 归约和量纲恢复正确；
- 双马赫顶部边界使用 stage 时间且不读 ghost 几何，正式 case 能到 `t=0.2`；
- 串行与 MPI 全量旧测试没有回归；
- 文档没有宣称未真正运行的人工算例已经通过。

若其中任何一项失败，只提交已经独立成立的小目标；修复后重新验收。历史 case 结果与新配置
不一致时必须标明生成版本并重算，禁止通过手工补列或修改结果文本伪造证据。

## 26. v1.1 热路径工作区与发布维护

### 26.1 结构缓存键和版本

v1.1 的求解器不再在每次残差内重新装配不变对象。缓存身份由

```text
(mesh_signature, partition_digest, profile, dimension, local_layout)
```

共同决定：每个本地块/有效轴的 `LineOperators`、face/operand/gradient/viscous halo plan、
block registry、面通量/梯度字段、SSPRK 初始状态和通信消息缓冲都必须与该键一致。不得只按
网格尺寸复用，因为相同尺寸可能具有不同拓扑、profile 或 owner 布局。

字段/描述符的单调 `version` 是时序一致性，不属于结构键。一次残差开始时统一推进版本，
plan 的 `set_version` 只能改变消息头版本，不能改变 pairs、owner、tag 或 payload 长度。工作
字段改写前用 NaN 重置；消费者逐项验证 profile、extent 和版本。若改动网格、分区、profile、
维数或本地布局，必须显式重建完整 workspace，不能在旧缓冲上局部修补。

### 26.2 API 与性能不变量

- 生产残差优先使用预构造 `LineOperators` 和写入既有字段的 `*_into` 入口；返回值 API 只作
  兼容薄包装，不能重新成为热路径默认。
- SSPRK 每次整步覆盖复用初始状态，Q 阶段稳健模式的未修改 `U^n` 回滚契约仍优先于复用。
- MPI exchanger 在初始化时按 descriptor 精确分配缓冲；交换只做 pack/post/wait/version
  check/unpack。没有经过内部区/连接邻域证明，不得把同步交换改名为“通信计算重叠”。
- root 输出每次只持有当前 zone/quantity payload；CGNS 检查点和场量写完即释放。修改 I/O
  时必须用独立重读同时检查 schema 和数值，不能只看文件是否存在。

涉及这些路径的提交至少运行 allocation probe、串行/MPI 数值等价、旧版本/错误长度/NaN
注入以及对应性能组。详细公式、分配减少率和扩展效率定义见《算法补充》11.6；冻结阈值见
开发仓库的 `docs/v1.1.0/stage-t-design.md`。

### 26.3 发布步骤

1. 在独立 stage 分支提交设计、实现和机器可读证据；禁止把运行输出或用户未跟踪目录加入。
2. 从空串行/MPI 目录配置 Release、编译、全量 CTest，并执行算法规格校验。
3. v2.0_pre 用 `tools/package_v2_pre.py` 从 `HEAD` 已跟踪 payload 生成无算例源码目录；核验
   `WCNS_SOURCE_REVISION` 和包内 `PACKAGE_CONTENTS.sha256`。
4. 从没有 `.git` 的解包目录重新构建、安装、运行和重启，核对 manifest 的版本/来源提交。
5. 运行发布矩阵、错误路径和固定性能协议，写阶段验收报告；候选标签只指向报告提交。
6. 候选后的人工卡口通过前，不得合并 `main` 或创建正式版本标签。

v2.0_pre 的完整发布卡口见 [`v2.0-pre-capability-matrix.md`](v2.0-pre-capability-matrix.md) 和
[`v2.0-pre-validation.md`](v2.0-pre-validation.md)；历史 v1.1.0 卡口保存在开发仓库。

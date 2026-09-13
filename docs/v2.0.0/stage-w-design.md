# v2.0.0 阶段 W 设计：通用 RANS/LES 基础设施

状态：**实现完成，等待阶段 W 自动验收结果回填。** 本阶段不实现任何具体 RANS/LES 闭合；
生产 registry 只注册 `none`，因此 `sa_neg`、SST、k--epsilon 和各 LES 名称即使能被识别，
也会在启动前明确拒绝。这样可以先冻结公共数据和离散契约，又不会误报模型可用。

## 1. 范围与不变量

阶段 W 建立模型场、模型接口、标量通量、边界/容许性、MPI、壁距、黏性耦合、源 Jacobian、
配置、输出与重启的共同底座。以下不变量必须保持：

1. 平均流始终是五个 Euler 守恒分量；模型变量存于独立 `TurbulenceFieldSet`，不得以第六个
   Euler 分量进入特征分解或 Riemann 求解器。
2. `model=none` 不分配模型场、不增加热路径临时分配，原五分量通量调用的是显式零闭合重载。
3. 模型标量和平均流共享块划分、连接索引变换及非阻塞 halo 机制，但标量按自身分量数打包。
4. 任一非有限或违反 descriptor 下界的候选整体拒绝；不进行无日志 clipping。
5. schema 1 的摘要、restart signature 和运行语义保持不变。schema 2 在 W 只允许
   `none + ssprk3 + none`，其余已规划组合在实现阶段到来前 fail closed。

## 2. 数据所有权与生命周期

`TurbulenceFieldDescriptor` 固定字段名、角色、量纲、严格正性和下界；descriptor 顺序同时是
内存分量顺序、MPI 消息顺序、输出注册顺序及重启 payload 顺序。签名使用全部 descriptor 元数据，
恢复时先比较签名和长度，再一次性写入并复验容许性。

`TurbulenceFieldSet` 使用一个 cell-major `Field<Real>` 保存任意少量模型量，与块的真实区及 ghost
范围一致。空集合持有零分量且不分配数据；`reset` 先在临时对象中完整校验名称唯一性和布局，
再提交新存储，避免半初始化状态。interior 的 checkpoint 打包采用
`component -> k -> j -> i` 的确定顺序，ghost 永不持久化。

模型生命周期冻结为：配置解析与拒绝未知组合 → registry 构造模型 → 读取 descriptor 并分配场 →
初始化真实单元 → 交换连接 halo/施加物理边界 → 形成面通量/源/闭合 → 检查候选 → 接受后输出或
checkpoint。阶段 W 通过 `none` 和 dummy 数据逐层验证该契约；具体模型在 X/Y/Z 接入 driver。

## 3. 标量输运、边界和事务

离散公式的唯一规格见《算法补充》第 12.2 节。面接口接收平均流权威质量通量、左右比值、非负
扩散系数和法向梯度，返回可分别诊断的 `advective`、`diffusive` 与
`net = advective - diffusive`。上风方向只由质量通量符号决定；左右模型量不会改变 Euler 波速。

物理边界公共规则包括零阶外推与 Dirichlet 镜像。模型可在后续阶段根据壁面/远场语义产生边界值，
但不得另造一套连接 halo。候选更新在临时向量中计算并按 descriptor 逐项检查，全部通过后才返回；
异常时原状态不变。逐面降阶和整步拒绝继续沿用第 11.2 节事务，本阶段仅冻结模型量必须参与同一
接受决定，具体模型的重构回退在 X/Y 接入。

W 的制造验证在周期均匀网格上使用面中心二阶插值与中心扩散。对

$$
\phi_m=2+\sin((m+1)x),\quad
S_m=a(m+1)\cos((m+1)x)+\mathcal D(m+1)^2\sin((m+1)x)
$$

验证 1/2/3 分量的稳态残差二阶收敛；线性耗散 `dphi/dt=-0.7 phi` 验证 SSPRK3 三阶。该测试
验证公共面通量的符号、分量无关性和源组合，不代替后续 SA/SST 的高阶 WCNS MMS。

## 4. MPI 与连接

现有 `BlockFieldRegistry/HaloExchanger` 本来就是任意固定分量的标量容器。W 将模型场直接注册到
该路径，并在原生多块轴置换网格上逐项检查 1/2/3 分量、1/2/4 rank 的 donor interior 与 receiver
ghost 完全相等。标量没有方向分量，因此连接 transform 只映射索引，不旋转数值。

Case06 同时需要两个基础兼容修正：

- 对 `CellDimension=2, PhysicalDimension=3`，只有全部有限且位于同一常 Z 平面的网格进入现有
  XY 二维度量；非平面嵌入立即拒绝。
- reciprocal 同 zone 自连接按连接几何规范键选一对，仍生成两个有向 exchange；排序键加入连接
  名称以保持同块双方确定顺序。

这两项有独立单元测试，并让 Case06 dry-run 越过阶段 V 的两个失败签名。随后暴露的“一个块面
同时含物理壁和尾迹连接”非张量边界权重缺口不属于 W 的通用湍流底座，已登记为 AB 正式算例
前置回归；不得用忽略部分连接或修改网格绕过。

## 5. 壁距算法

壁距公式、规范物理键和 tie-break 见《算法补充》第 12.3 节。二维 wall patch 的相邻顶点形成
segment；三维每个结构面四边形沿固定 `q00-q11` 对角线拆为两个 triangle。退化或非有限几何元
立即拒绝。只收集配置解释后的无滑移绝热/等温壁，滑移、对称、远场、连接和运行时人工面排除。

各 rank 将几何元序列化为固定十个 `Real`（类型加三个三维点），通过 `MPI_Allgatherv` 收集，
再规范化、排序和去重。因此每个 rank 构造完全相同的 BVH。BVH 以质心最大展宽轴稳定二分，叶
节点最多八个元；查询先用 AABB 平方距离剪枝，再调用点到 segment/triangle 的精确最近点算法。
等距时选择规范全局表的较小序号。测试覆盖平板/槽道、三角面内外最近点、跨块去重、等距选择、
非有限输入和 1/2/4 rank 一致性；圆/球的系统网格误差验证随具体模型算例在 X/Y 扩展。

## 6. 黏性闭合与源 Jacobian

`TurbulenceViscousContribution` 以对称六分量应力、三个能量热流分量及涡黏度组成。唯一黏性
Cartesian 通量入口先形成原层流通量，再按同一正扩散号约定加入模型应力；能量增量为

$$
\Delta F^E_j=\widetilde u_i\tau^t_{ij}+q^t_j.
$$

二维禁止 `xz/yz/z-energy` 通量，全部值必须有限且 `mu_t>=0`。显式零贡献与旧五参数重载逐位
相等，防止 `none` 基线漂移。

模型源接口同时返回 $S_t$ 和行主序 $J_t$，隐式对角严格构造
$A_{t,diag}=\alpha I-J_t$；尺寸、有限性与正时间对角均检查。AA 阶段 LU-SGS 直接消费该块，
无需猜测模型源符号。

## 7. 配置、输出和重启骨架

schema 2 在生产 parser 中启用三个阶段 W 固定键：`turbulence.model`、`time.integrator`、
`preconditioner.type`；RANS 可选 `turbulence.prandtl/wall_treatment` 也完成类型解析，但非 `none`
模型仍因尚未注册而拒绝。schema 1 明确拒绝 v2 键，避免静默忽略。

`FieldQuantityRegistry::register_turbulence_field` 依据 descriptor 映射无量纲、速度平方、运动/动力
黏度、逆时间、耗散和长度尺度，并从块的独立模型场读取。没有活动模型时不注册这些量。
`capture/restore_turbulence_restart` 冻结 descriptor 签名和真实区 payload；实际 CGNS checkpoint
在 `none` 下继续使用 v1 格式，X 首次出现生产模型场时再增加 CGNS 数组和版本迁移，旧 checkpoint
仅能恢复到 `none`。schema 2 的模型/时间/预处理选择已经进入 summary 与 restart signature。

## 8. 自动验收矩阵与通过规则

阶段 W 必须全部满足后才允许创建候选标签：

1. 独立串行/MPI Release 配置、构建、安装和完整 CTest 全过；算法规格检查全过。
2. `none` 零闭合逐位相等，schema 1 release matrix、输出、checkpoint/restart 和分配探针无回归。
3. 1/2/3 标量制造残差观测阶大于 1.95，SSPRK3 线性源观测阶大于 2.9。
4. 1/2/4 rank 模型 scalar halo 与 wall primitive 全局结果一致。
5. 壁距、黏性应力/能量、源 Jacobian、非法值和 schema 2 fail-closed 单元测试通过。
6. Case06 原文件哈希不变、`cgnscheck` 无 error；dry-run 必须精确到达已登记的部分面权重缺口，
   若出现其他错误则阶段失败。
7. `git diff --check` 通过，工作树仅保留两处已登记用户目录；不提交构建物或长算结果。

依据项目负责人 2026-09-13 对 V→W 连续执行且中间不设人工审批的授权，V 通过后已直接进入 W。
W 自动卡口通过后可创建 `v2.0.0-w-candidate.1` 并非快进合入 `release/v2.0.0`；阶段 X 未获连续
授权，必须在 W 合入后停在人工判断点。

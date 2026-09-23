# v2.0.0 阶段 Z 设计冻结：LES 与可重启统计

状态：**Z0 设计冻结；实现与自动验收进行中。** 基线为 `release/v2.0.0` 提交
`b9da7a4`，阶段分支为 `stage/v2.0.0-z`。本文件只授权轻量解析场、公式单元和微型三维
算例；HIT、正式槽道、完整能谱及长时间统计继续列为服务器 S 级 `pending`。

## 1. 范围与不可越过的限制

阶段 Z 完成 Z0--Z11：确定性三维滤波、五种 SGS 闭合、SGS 动量/能量通量、能量诊断、壁面/
截面量、按接受物理步加权的统计状态、重启以及 BDF2 双时间 LU-SGS 接口。LES 仅允许三维、
黏性、非定常、`lu_sgs+bdf2`；非法组合必须在 LES 工作区分配前拒绝。

本机卡口不把短 smoke 称作物理湍流验证，不运行长期槽道或 HIT。用户的 Case05/Case06 手工目录、
已有结果、根目录 `output/` 和 `tmp/` 均不读取、不修改、不暂存。

## 2. Z0 数值定义

### 2.1 网格滤波和边界

1. 网格/测试滤波均使用逐逻辑方向可分离三点核
   $w_{-1},w_0,w_1=(1/4,1/2,1/4)$；三维张量积权重和严格为 1。
2. 已交换的块连接与 MPI ghost 属于完整 stencil，不能按本地块边界截断。物理边界使用最近合法
   单元常值延拓后仍应用完整固定核；因此常量精确保持，边界不因 rank 或遍历顺序改变。该闭合
   仅是一阶边界闭合，必须在 manifest 标识为 `nearest_constant_extension`。
3. 守恒量先过滤；速度和二阶速度矩使用 Favre 形式。测试滤波比固定为
   $\widehat\Delta/\Delta=2$。
4. $\Delta=C_\Delta V_c^{1/3}$，默认 $C_\Delta=1$。滤波核仍沿三个逻辑方向作用；首版不以
   最短边静默缩小 $\Delta$。同时输出网格各向异性诊断，超过 8 的配置只具受限支持等级。
5. 测试滤波和动态平均所需 halo 固定为一层；初始化时检查现有 halo 足够。两级操作之间先交换
   中间标量，所以总离散支持不因运行时切分而改变。

### 2.2 模型常数和保护

| 项 | 冻结值/范围 |
|---|---|
| `les.smagorinsky.cs` | 默认 0.17，合法 `[0,0.3]` |
| `les.smagorinsky.wall_damping` | `none`（默认）或显式 `van_driest` |
| van-Driest | $f_d=[1-\exp(-y^+/26)]^2$，$\mu_{sgs}$ 乘 $f_d$ |
| `les.similarity.cb` | 默认 1，任意有限值 |
| `les.wale.cw` | 默认 0.325，合法 `(0,1]` |
| `les.sgs_prandtl` | 默认 0.9，正有限值 |
| dynamic average | 固定 `local_box3_tensor` |
| dynamic denominator floor | $\epsilon_M=10^{-20}$，正有限值 |
| dynamic clipping | $C_d\in[-0.05,0.09]$；负值保留并单独诊断 |

van-Driest 只在已定义壁面距离和同一壁面 registry 给出的 $y^+$ 可用时生效；否则启动失败，
不能用单元编号、应变率或经验距离替代。动态系数先形成逐点 $L:M$ 和 $M:M$，再以同一
三点张量核分别平均后相除、保护并裁剪。纯尺度相似模型不附加隐藏涡黏性。

### 2.3 应力、热通量和符号

模型返回的是《算法补充》12.12 定义的 $\tau_{ij}^{sgs,d}$。求解器黏性正通量使用
$-\tau_{ij}^{sgs,d}$；涡黏性分支因此贡献 $+2\mu_{sgs}S_{ij}^d$，相似分支贡献
$-C_BL_{ij}^d$。SGS 热通量只由涡黏性部分给出，系数为
$\mu_{sgs}/[(\gamma-1)M_\infty^2Pr_{sgs}]$，与现有无量纲温度梯度约定一致。

$\Pi_{sgs}=-\tau_{ij}^{sgs,d}S_{ij}^d$；分别累计 $\max(\Pi,0)$ 和
$\min(\Pi,0)$。动态负 $C_d$ 与相似应力允许局部 backscatter，但所有应力和热通量必须有限；
裁剪事件及上下限命中次数属于诊断，不改变模型身份。

## 3. 统计和几何定义

1. 普通统计以接受物理步的实际 $\Delta t$ 为权；Favre 统计以 $\rho\Delta t$ 为权。禁止以
   LU-SGS 内迭代、被拒绝尝试、输出重读或同一 `(step,time)` 的重复事件更新。
2. 冻结在线状态为 `sample_count, W, mean, M2/covariance, first_step/time,
   last_step/time`。RMS 为总体量 $\sqrt{M_2/W}$；协方差也除以 $W$。
3. 合并两个分区/重启段时采用 Chan--Golub--LeVeque 的确定性加权合并式；MPI 以全局块 ID、
   分箱 ID 固定顺序合并，不依赖归约树。
4. checkpoint 保存累加器版本、配置身份、全部状态和最后接受事件；相同事件重放必须无操作，
   不兼容的统计身份必须拒绝。
5. 三维方向采用显式右手系 $(e_D,e_L,e_Y)$；分箱由冻结的物理翼展坐标边界定义。一个物理面
   只由全局面身份累计一次，压力与黏性载荷分开保存后再求和。
6. 壁面 $u_\tau,y^+,\mu_{model}/\mu$ 使用《算法补充》12.9 的同一牵引、壁面密度/黏度和
   真实几何距离；未定义热学分母时 `St/Nu` 明确为 `not_applicable`。

## 4. 实现分层

- `les_filter`：固定核、物理边界闭合、Favre 矩和中间 halo 契约；
- `les_models`：逐点张量、Smagorinsky、Bardina/混合、Germano--Lilly、WALE 和能量传递；
- `turbulence_model`：配置、registry、restart signature 和面通量适配；
- `weighted_statistics`：去重的接受步状态、加权均值/RMS/协方差/Favre 与可序列化状态；
- 求解器/runtime：三维启动前验证、SGS 工作区、输出/边界 registry、checkpoint 与 LU-SGS
  谱半径。

生产路径不得用测试专用数据绕过上述层次。`C_s=0` 必须走精确零 SGS 通量，保证与 `none`
残差逐位一致（除诊断字段外）。

## 5. 阶段 Z 自动验收卡口

1. 规格/格式：`git diff --check`、算法规格检查、配置文档与 restart signature 一致；
2. 公式：五种模型逐点参考、旋转/纯剪切/近壁 WALE、`C_s=0`、热通量和正负 $\Pi$；
3. 滤波：常量、线性、单 Fourier 模式、物理边界以及原生块/运行时切分的独立卷积对照；
4. 配置：二维、steady、无黏、SSPRK3 LES 在分配 LES 时间/统计工作区前拒绝；所有数值键进入
   summary/restart signature；
5. 统计：解析加权均值/RMS/协方差/Favre、重复/拒绝/伪迭代不计数、连续与序列化重启一致；
6. 求解器：均匀流、解析三维梯度和受限周期盒；1/2/4 rank 残差、滤波支持、动态系数和统计
   在冻结容差内一致；
7. 隐式：每个 LES 模型完成微型 BDF2 双时间运行，内残差达到冻结容差，统计样本数只随接受
   物理层增加；
8. 资源：本机单项默认不超过 60 秒、整套 L 级不超过 20 分钟；HIT/正式槽道/长 TGV 明确
   `pending_server`，不能用短 smoke 顶替。

全部自动卡口通过后创建不可移动候选标签 `v2.0.0-z-candidate.1` 并停在人工判断。人工重点检查
滤波边界、相似/动态 backscatter、SGS 与数值耗散区分、统计重启接缝和受限支持声明。

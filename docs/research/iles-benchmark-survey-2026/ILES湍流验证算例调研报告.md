# 粗网格 ILES 湍流验证算例调研报告

研究背景 参考配置 可信数据与验证方法

调研日期：2026 年 10 月 3 日　适用对象：WCNS 多块结构网格可压缩流程序

## 1 研究结论与范围

建议在已有槽道湍流、周期山、SD7003 翼型和超声速压缩拐角的基础上，首先增加 Taylor–Green 涡、衰减均匀各向同性湍流和低雷诺数后台阶，随后增加方管、零压梯度边界层及激波与各向同性湍流相互作用。这些算例分别隔离能量耗散、近壁输运、入口湍流、分离再附和激波捕捉问题，能够解释复杂算例中的误差来自哪里。圆管、逆压梯度边界层、混合层、射流、壁面隆起和圆柱尾流则用于检验跨流动类型的适用性。

本报告详细综述 19 类算例，其中 4 类已在程序中实现，15 类为新增候选。Gaussian bump 等近期研究作为经典算例的后续发展讨论。每类均说明物理背景、可复现配置或实验设计、参考结果及数据入口、应比较的统计量、对湍流计算的指导意义和研究进展。这里的 ILES 指保留分子黏性、没有显式 SGS 模型和壁面模型，由空间离散、通量、滤波或限幅产生隐式小尺度耗散的三维非定常计算。

用户所说的“精确结果”，在湍流验证中应理解为有明确来源和精度说明的实验或充分分辨 DNS 统计。本报告区分原始数值表、论文表值、读图近似和高分辨率 LES；没有取得原始数组的项目明确标出。DNS 仍有离散、计算域和采样误差，实验仍受探针分辨率、壁面条件和入口状态影响。不存在一组可跨不同工况直接使用的通用“真值”。

本次只开展文献调研、已有资料核对和报告制作，没有修改求解器、生成新增生产算例或启动流场计算。有关结构网格和验证顺序的建议属于本报告的工程判断，不是文献已验证的 WCNS 配置。

## 2 算例体系与选择依据

表 1 按主要物理机制归类。DNS 与实验表示主参考类型；“补充 LES”不能替代独立验证证据。表内工况用于选取参考数据，具体边界条件和数据限制见各节。

| 编号 | 算例 | 主参考工况 | 主测试目标 | 参考类型 |
|---|---|---|---|---|
| C01 | 平面槽道 | Reτ 180 及 550–5200 | 壁面应力与湍流维持 | DNS |
| C02 | 周期山 | Reₕ 1400／2800／10595 | 周期驱动与分离再附 | 分工况 DNS／LES |
| C03 | SD7003 | Re꜀ 60000，α 4° | 分离泡与自然转捩 | 实验，补充 ILES |
| C04 | 24° 压缩拐角 | Ma 2.25，Reδ₀ 15800 | 激波分离与再附加热 | 高分辨率计算 |
| C05 | Taylor–Green 涡 | Re 1600，Ma 0.1 | 转捩及动能耗散 | 谱方法 DNS |
| C06 | 衰减各向同性湍流 | CBC 5.08 cm 格栅 | 能谱与衰减率 | 实验 |
| C07 | 强迫各向同性湍流 | JHTDB 1024³ | 稳态能量通量 | DNS |
| C08 | 光滑圆管 | Reτ 180–5200 等 | 曲率与摩擦系数 | DNS、实验 |
| C09 | 光滑方管 | Reτ 约 150–1055 | 应力各向异性和二次流 | DNS |
| C10 | 零压梯度边界层 | Reθ 约 670–4060 | 入口发展与壁面输运 | DNS |
| C11 | 逆压梯度边界层 | β 1 与约 39 | 非平衡与压力历史 | DNS，补充 LES |
| C12 | 可压缩平板边界层 | Ma 2–10，冷热壁 | 密度加权与热通量 | DNS |
| C13 | 平面混合层 | 时间发展／空间发展 | 卷吸与剪切层增长 | DNS、实验 |
| C14 | 圆射流 | Ma 0.9，Reᴅ 10⁶ | 入口状态、混合与噪声 | 实验，补充 LES |
| C15 | 后台阶 | Reₕ 5100／约 36000 | 固定分离与恢复 | DNS／实验 |
| C16 | NASA 壁面隆起 | Re꜀ 约 9.3×10⁵ | 光滑表面分离 | 实验 |
| C17 | 圆柱尾流 | Reᴅ 3900 | 脱涡、转捩与尾迹 | 实验、DNS 文献 |
| C18 | 平面激波与各向同性湍流 | 上游 Reλ 40 等 | 激波后湍动能放大 | DNS |
| C19 | 斜激波与平板边界层 | Ma 2.28，转角 8° | 激波低频与入口影响 | DNS、实验 |

TGV 是确定性初值触发的转捩基准，并非已经充分发展的统计稳态湍流。C03 和 C17 也同时检验转捩机制。将它们与成熟湍流槽道、强迫 HIT 和湍流边界层组合，才能避免仅在某一种失稳机制上得到好结果。

对无壁面模型的粗网格 ILES，高雷诺数近壁流存在固有成本限制。缩减网格而保持无滑移边界，并不等于自动获得合理的壁面模型。应从有 DNS 的低至中等雷诺数开始；高雷诺数实验算例用来识别适用边界，不能仅凭平均压力吻合宣称壁面应力和热流可靠。

## 3 程序已实现的四类算例

### C01 平面槽道湍流

**背景与设置。** 压力驱动的双平板槽道把入口与几何复杂性移除，是研究近壁条带、准流向涡、摩擦阻力和 Reynolds 应力输运的标准问题。流向、展向周期，上下壁无滑移；必须区分恒压梯度和恒流量驱动。现有基线目标 Reτ=180，初始体积平均 Ma=0.1，计算域为 2πh×2h×πh，36×48×36 单元，壁面法向加密。这一较小周期盒与大盒 DNS 的能容尺度并不完全一致。[本地配置与依据](../../../cases/manual/case05_3d_turbulent_channel/README.md)

**参考结果。** 本地采用 ERCOFTAC Case 032 保存的 Kim–Moin–Moser 数据，标称 Reτ=180。由本地 simul1.dat 离散剖面积分得到 Uᵦ⁺≈15.5172，由此计算 Cf≈0.008306；这是离散表的派生值，不是原论文给出的无限精度常数。另行下载的 Lee–Moser 2015 文件 LM_Channel_0180_mean_prof.dat 实际 Reτ=182.088，ν=0.00035，Uᵦ=1，uτ=0.0637309；据文件头计算 Cf≈0.00812326。两组不能因文件名含“180”就混为同一参考。[R01](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case032) [R02](https://turbulence.oden.utexas.edu/channel2015/content/README_2015.html)

**测试量与意义。** 最少比较 U⁺、三分量脉动均方根、负 Reynolds 剪应力、Cf、实际 Reτ、上下壁对称性、流向及展向能谱。总剪应力应符合充分发展动量平衡，体力做功、黏性耗散及数值耗散应闭合。仅拟合对数律可能掩盖法向脉动不足和非物理大尺度能量堆积。对低马赫 WCNS，这是检查通量耗散是否压制湍流的直接算例。

**已有证据与进展。** 本地 final-campaign-report 记录了不同算法分支合计约 227 万推进步及 DNS 对比，并非只有短步测试；其版本号为早期 1.0.0／1.0.1。该记录不自动验证后续 v2.5 的所有路径。末段判稳表仍显示分支间摩擦和壁面对称性差异，适合继续用统一统计窗重评，而不宜只报单个末时刻。现代参考已经扩展至 Reτ≈5200，并提供统计标准误及应力预算；应优先使用这些误差信息，而非把所有偏差归为 ILES 模型误差。[R02](https://turbulence.oden.utexas.edu/channel2015/content/README_2015.html)

### C02 周期山流动

**背景与设置。** 平滑山丘诱导的分离剪切层、回流泡和再附使本算例比槽道更接近工程分离流，同时保持流向和展向周期，不需要人为合成入口湍流。现有 v2.4 配置依据用户论文第 5.3.4 节，Reₕ=1400、Maᵦ=0.1，域长高宽 9h×3.035h×4.5h，88×48×24=101376 单元，通过山顶截面流量反馈调整均匀驱动力。该三向网格分解、等温壁和 PI 参数是程序补充选择。[本地实施依据](../../periodic-hill-v2.4.md)

**参考结果必须按雷诺数分开。** Fröhlich 等的 Reₕ=10595 常用数据是高分辨率 LES，公开平均剖面、应力和预算；再附位置约 xᵣ/h=4.6–4.7。NASA/TMR 的另一数据集是 Balakumar 的可压缩 DNS，Reₕ=2800、Ma=0.2，801×351×513 网格，分离约 xₛ/h=0.233、再附约 xᵣ/h=5.50。二者都不能直接用于当前 Reₕ=1400 的定量验收。Breuer 等 2009 给出跨雷诺数的实验及 DNS／LES 比较，可作为寻找同工况数据的主线；本次尚未获取与现有 Reₕ=1400 完全匹配的原始 DNS 数组。[R03](https://tmbwg.github.io/turbmodels/Other_LES_Data/2dhill_periodic.html) [R04](https://tmbwg.github.io/turbmodels/Other_DNS_Data/2dhill_periodic_compress.html) [R05](https://www.sciencedirect.com/science/article/pii/S0045793008001126)

**测试量与意义。** 比较下壁有符号 Cf、Cp、分离和再附、回流强度、平均二维速度及全部六个应力分量；保留现有 x/h=0.05、2、4、8 的剖面，并按参考库补齐其他截面。流量控制误差和压力梯度历史也是必需输出。山后剪切层中应力增长过慢、再附过迟，可能源于入口处已周期循环的湍流过弱、网格分辨率不足或耗散过强，不能只通过加大驱动力修正。

**研究进展与结构网格建议。** 后续研究已经从单一山形扩展到几何参数变化和数据驱动闭合评估。对本程序，先在相同物理工况下做三向加密与展宽敏感性，再改变山形；曲率和块接口处的度量一致性会直接影响局部应力。现有算例可作为“无入口不确定性的分离流”主基准，但同工况原始数据仍是待补项。

### C03 SD7003 低雷诺数翼型

**背景与设置。** SD7003 上表面出现层流分离、剪切层失稳、三维转捩和再附，适合检验 ILES 能否自行产生合理转捩。现有配置为 Re꜀=60000、Ma=0.2、α=4°、展宽 0.2c、远场约 30c、绝热无滑移壁，8 块结构网格共 1720320 单元。参考 UIUC 坐标，但尾缘圆钝化和节点分布为工程重建。高阶工作坊的相关工况 Ma=0.1，应与本程序 Ma=0.2 区分。[本地资料](../../chapter5-v2.5-research.md) [R06](https://cfd.ku.edu/hiocfd/case_c3.3.html)

**可核验参考结果。** Galbraith–Visbal 2010 是 ILES 文献，不是 DNS。其零来流湍流度表值为 xₛ/c=0.23、xₜ/c=0.55、xᵣ/c=0.65，最大泡高约 0.030c。用户文档列出的 IAR 实验为 0.33／0.57／0.63，TU-BS 为 0.30／0.53／0.62，AFRL 为 0.18／0.47／0.58。不同风洞的湍流度和设施影响明显，不能将这些点拼成同一条“精确”曲线。TU-BS 在原论文中 Tu=0.08%，用户文档约写为 0.1%。[R07](https://doi.org/10.2514/6.2010-4737)

**测试量与意义。** 比较上、下表面 Cp，上表面 Cf，分离／转捩／再附位置，泡高，升阻力及俯仰力矩均值和 RMS，x/c=0.1、0.5、0.95 的速度与应力，剪切层速度频谱和展向相关长度。转捩位置需要固定识别标准，不能把 Cf 零点当成转捩点。过量数值耗散会延迟三维失稳，过大的入口噪声会提前转捩，两者有可能在升力上相互抵消。

**数据限制与后续研究。** 已有参考包含论文表值及带像素误差的曲线数字化，不含完整原始实验时序。应以实验检验可信范围，以文献 ILES 检查同类方法差异。后续重点是来流湍流度、展宽、尾缘形状及扰动可重复性，不宜仅增加总网格数。v2.5 短算测得尾缘小尺度使显式时间步极小，生产效率与物理统计都尚需目标服务器评估；本报告未新做长算。

### C04 超声速压缩拐角

**背景与设置。** 24° 拐角迫使来流压缩，形成激波诱导分离、压力平台、再附及局部热流峰。现有配置沿用用户文档：Ma=2.25、Reδ₀=15800、δ₀=0.0006096 m、T∞=170 K、Tw=323 K，拐点 x/δ₀=100、展宽 6δ₀。层流相似入口经体力扰动触发转捩；12 块结构网格共 8847360 单元。体力幅值和波形为程序补充设置，必须先验证拐角前湍流状态。[本地资料](../../chapter5-v2.5-research.md)

**参考结果及边界。** Porter–Poggie 2019 及作者算例页区分约 30 亿单元高分辨率计算与约 4000 万单元长时间计算，后者适合研究低频，不能把二者视作同一数据序列。前期 2017 文献的展宽为 10δ₀，与现有 6δ₀ 不同。已有论文图数字化大致给出 xₛ/δ₀≈88、xᵣ/δ₀≈107，仅适合趋势核对。精确热流、压力 RMS、全部 Favre 应力和长时间谱数组本次仍未取得。[R08](https://doi.org/10.1063/1.5078938) [R09](https://engineering.purdue.edu/~jpoggie/ramp/index.html)

**测试量与意义。** 入口首先检查 x/δ₀=70 附近的实际壁面单位及 x/δ₀=80 的变换速度、厚度、形状因子和应力。交互区比较有符号壁面切向 Cf、流向分量 Cfx、Cp、热流／Stanton 数、压力和温度 RMS、分离面积、激波足位置 PDF 与功率谱，以及壁面压力的两点相关和相干性。倾斜壁面上的 Cf 与 Cfx 不等价。

**研究进展与应用判断。** 研究重心已从平均分离长度扩展到低频运动的上游／下游联系、间歇性和热边界影响。对 ILES，这同时检验低耗散湍流推进和局部激波稳定化；只看激波没有振荡远远不够。应记录传感器、限幅及降阶触发区域，检查其是否大范围侵入湍流边界层。该例属于后期综合验证，宜先通过 C12、C18 和 C19 的单独检查。

## 4 自由湍流与数值耗散基准

### C05 Taylor–Green 涡

**背景与配置。** 三维大尺度涡经拉伸产生细尺度，随后衰减；初值确定、三方向周期，适合区分离散格式的耗散与色散。高阶工作坊 C3.5 指定 Re=1600、立方域 [−πL,πL]³；可压缩分支为 Ma=0.1、γ=1.4、Pr=0.71、常物性。初速度采用标准正弦余弦 TGV 场，压力采用对应空间扰动。特别要注意原配置要求初始温度均匀、密度随压力变化，不能同时另设密度均匀而仍声称完全同一初场。时间范围为 tU₀/L=0–20。[R10](https://cfd.ku.edu/hiocfd/case_c3.5.pdf)

**实际参考数组。** 本次取得工作坊 data.tgz 中 spectral_Re1600_512.gdiag，文件明确标为 512³ 谱方法结果。以下为原数组的四舍五入显示，K 以 U₀² 归一化，耗散 ε 以 U₀³/L 归一化，涡量平方积分量 Ω 为半个体积平均涡量平方。原数据每 0.01 输出，最后记录在 19.99；不将它改写成 t=20。[R11](https://cfd.ku.edu/hiocfd/data.tgz)

| 无量纲时间 | K | ε | Ω |
|---|---:|---:|---:|
| 0 | 0.12500000 | 0.000468743 | 0.375000 |
| 5 | 0.11843575 | 0.004127097 | 3.301699 |
| 8 | 0.09828988 | 0.010372895 | 8.298340 |
| 8.97 | 0.08678956 | 0.012857528 | 10.286127 |
| 10 | 0.07447917 | 0.011272054 | 9.017462 |
| 19.99 | 0.02156736 | 0.001873496 | 1.498795 |

8.97 是此离散数组中 ε 最大值所在时刻，不是对连续真解的无限精度极值。采用可压缩分支时，应把压力膨胀和密度扰动的影响纳入动能收支，再与不可压缩谱参考解释差异。

**测试量与 ILES 意义。** 同时比较 K、负时间导数、分子黏性耗散、Ω、能谱和 t=8 的固定截面涡量结构。只看总动能曲线可能出现“数值耗散补偿分子耗散不足”的假吻合。对周期盒，可以用完整动能预算残差估计净数值耗散，但它不是普适、局部、恒正的等效 SGS 黏度。

**实施与进展。** 建议先用 32³、64³、128³ 三档做方法比较，并对同一网格改变块划分，检查接口造成的耗散变化；这些是本报告建议，不能当作已通过验证的分辨率。TGV 广泛用于高阶方法评估，但晚期仍可能存在非均匀性和各向异性，因此还必须加入 C06／C07；无黏 TGV 的结论也不能直接替代 Re=1600 的黏性基准。

### C06 Comte–Bellot 与 Corrsin 衰减各向同性湍流

**实验背景与设置。** 格栅后的近似各向同性湍流是验证小尺度耗散、能谱输运和衰减规律的经典实验。1971 年研究的主要格栅间距 M=5.08 cm，来流 U₀=10 m/s，格栅实度约 0.34，并通过下游弱收缩改善各向同性。常用测量时标 U₀t/M 为 42、98、171；论文还含 2.54 cm 格栅，不能混用。原论文表 2、3 提供一维及三维能谱，表 4 给出整体量。[R12](https://doi.org/10.1017/S0022112071001599)

| U₀t/M | u′ RMS cm/s | ε cm²/s³ | Reλ | Kolmogorov 尺度 cm |
|---|---:|---:|---:|---:|
| 42 | 22.2 | 4740 | 71.6 | 0.029 |
| 98 | 12.8 | 633 | 65.3 | 0.048 |
| 171 | 8.95 | 174 | 60.7 | 0.066 |

上表为 5.08 cm 格栅的原论文表值，保留原有效位数；例如 ε=4740 cm²/s³ 对应 0.474 m²/s³。Reλ 最后一点为 60.7，不采用部分二手资料中的 60.3。[原始全文与表 4](https://courses.washington.edu/mengr544/handouts/comtebellot-corrsin-jfm-71.pdf)

**计算对应。** 通常不直接解析实体格栅，而在三周期盒中生成与第一测站能谱一致的无散随机速度场，再计算到后两测站。必须保存随机种子、盒长与波数单位、低波数补偿方式和初始离散投影；首次滤波之后的实际 K、谱和无散误差都要检查。风洞的空间测量转为时间衰减需要声明 Taylor 冻结假设及使用的对流速度。

**测试量与意义。** 比较三时刻 E(k)、各分量能量、衰减率、积分尺度、Reλ 和各向异性。该例比 TGV 更接近已存在湍流的演化，能发现某格式对随机宽带场过度耗散，或把能量堵塞在截止波数附近。无需近壁网格，适合当前程序先行增加。

**研究进展及限制。** 后续工作广泛将实验谱用于 LES 初值与合成湍流验证，但实验一维谱与三维球壳谱有不同定义。应做多个随机实现，比较集合统计及低波数敏感性，避免只选一组种子使衰减曲线吻合。历史实验的有限测量范围和噪声修正也意味着不能用尾部谱值制定没有误差容限的硬阈值。

### C07 强迫均匀各向同性湍流

**背景与配置。** 在周期盒低波数持续注入能量，可形成统计稳态的级联，避免衰减问题中时刻对齐的困难。JHTDB 1024³ 数据采用 2π 立方域、ν=0.000185，保持波数模长不大于 2 的模态总能量，形成稳定大尺度强迫。模拟时间步 0.0002，存储间隔 0.002。不能将恒能低模态强迫替换为任意随机体力后仍宣称精确复现。[R13](https://turbulence.pha.jhu.edu/docs/README-isotropic.pdf)

**参考统计及时间窗。** 原始短窗 t=0–2.048 的统计为 K=0.695、ε=0.0928、u′=0.681、Reλ=433。扩展数据说明对 t=0–10.056 的平均为 K=0.705、ε=0.103、u′=0.686、Reλ=418，积分尺度约 1.364，约覆盖五个大涡周转时间。这一差异说明“数据库 Reλ=433”并非任意时间窗都严格成立。本次保存了原始 README，尚未下载体量巨大的三维场。[R13](https://turbulence.pha.jhu.edu/docs/README-isotropic.pdf)

**测试量与意义。** 除平均 K、输入功率和耗散平衡外，应比较 E(k)、积分尺度、二阶和三阶结构函数、速度梯度偏度／平坦度、压力统计及能量通量。具有 −5/3 区间并不单独证明耗散正确；谱幅值、通量及大尺度统计必须同时合理。可压缩程序应控制湍流马赫数，区分声学模态与旋度模态。

**结构网格与研究进展。** 笛卡尔多块周期网格很容易构造，难点在强迫与频域统计。JHTDB 已扩展到更高分辨率并更新访问工具；部分超大数据集只包含少数快照，不能拿它们计算长时间频谱。使用接口插值或有限差分梯度时，还应承认这与原谱导数不同，数据库场的“非零数值散度”不必然意味着 DNS 本身不无散。[R14](https://turbulence.idies.jhu.edu/database)

## 5 壁面湍流的扩展验证

### C08 光滑圆管湍流

**背景与配置。** 圆管提供封闭横截面、轴向周期和恒定压力梯度或恒流量驱动，可检验槽道结果能否迁移到曲率与轴向大尺度不同的壁湍流。定义 Reτ 时以管半径 R 为长度，体积雷诺数 Reᵦ 则通常以直径 2R 为长度。Yao 等 2023 的 DNS 使用轴向长度 10πR，覆盖至 Reτ≈5200，采用谱方法并量化统计不确定性。Pirozzoli 等 2021 的另一系列覆盖至约 6000，典型长度 15R。[R15](https://doi.org/10.1017/jfm.2022.1013) [R16](https://doi.org/10.1017/jfm.2021.727)

**参考结果与数据入口。** 2023 DNS 的均值、应力、压力和谱数据通过 Texas Data Repository 的 turbpipe 数据空间公开。建议选择其中低 Reτ 数据集建立首例，记录实际 Reᵦ、Reτ、长度、驱动和采样时间后再导入数组。本次确认论文与数据入口，未逐个下载仓库文件。摩擦应明确是 Cf 还是 Darcy 摩阻系数；后者为前者的 4 倍。不同 DNS 系列在摩阻和近壁脉动峰上有可辨差异，不能把某条经验摩阻律当成比 DNS 更精确的真值。[数据仓库](https://dataverse.tdl.org/dataverse/turbpipe)

**测试量与 ILES 意义。** 检查压降、体积流量、摩阻、径向速度剖面、三方向应力、近壁脉动峰、压力 RMS、轴向及周向谱。轴线附近应保持速度和几何正则性。若槽道表现良好而圆管摩阻偏差大，应先排查曲线网格度量、中心区域拓扑、周向分辨率和轴向盒长，再讨论模型普适性。

**实施与进展。** 适合中心方块加外围 O 型块的结构拓扑，避免退化轴线单元；统计应在物理径向上归并，不按块索引简单平均。近年高 Re DNS 继续讨论对数律和外层峰值，其拟合常数仍存在数据集与尺度选择差异。本程序宜先使用低 Reτ、充分近壁分辨的管流，暂不以极高雷诺数拟合常数作为粗网格 ILES 的合格条件。

### C09 光滑方管与第二类二次流

**背景与配置。** 四个无滑移壁围成方形截面，流向周期并由压力驱动。角区应力各向异性产生 Prandtl 第二类二次流，典型平均横截面有八个涡胞。Pirozzoli 等的 DNS 覆盖以半边长及周向平均壁剪定义的 Reτ 约 150–1055。平均横向速度仅为体积平均速度的约 1%–2%，因此需要比普通平均轴向速度更长的统计时间。[R17](https://doi.org/10.1017/jfm.2018.66)

**参考结果与数据状态。** 应从论文同一 Re 数据获取完整截面 U、V、W、应力及周向壁剪。Modesti 等的动量／阻力分解表明，在所研究条件下二次流贡献约占总摩阻 6%；这是指定分解下的结果，不是所有方管的固定比例。本次核实论文和表述，未取得可以立即用于误差计算的完整截面原始数组，实施前应联系作者或确认作者公开数据版本。[R18](https://doi.org/10.1017/jfm.2018.391)

**测试量与指导意义。** 除总体 Cf，还要比较二次流矢量、流向涡量、八涡拓扑、沿周向壁剪以及角平分线和壁面中线上的应力。标量涡黏性 RANS 难以完整描述驱动这种二次流的应力差异；ILES 则应通过解析非定常结构产生它。若横向平均速度接近采样噪声，不能凭图像像“八个涡”就判定通过，应报告时间分块置信区间。

**实施与进展。** 方形截面比圆管容易构造多块结构网格，能有效区分“总压降正确”和“各向异性正确”。2026 年 9 月 Verolini、Xiao、Pirozzoli 将研究扩展到极高 Re 的 WMLES，关注二次流强度的尺度变化；该研究使用壁面模型，不能直接证明无壁面模型粗网格 ILES 同样适用。[R19](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/turbulent-flow-in-a-square-duct-at-extreme-reynolds-number/C3DAB1AC62E975591CFD93D57C53E31E)

### C10 零压梯度平板湍流边界层

**背景与设置。** 与周期槽道不同，边界层沿下游增厚，入口湍流生成和对流发展成为独立误差来源。标准配置为平板无滑移壁、恒定外缘速度、展向周期、足够远的顶部和非反射出口；入口可由上游前驱场或经验证的回收再标度方法提供。应保持边界层外缘近零压梯度，而不是仅给一个恒定出口压力。

**参考结果。** KTH 的 Schlatter–Örlü 2010 DNS 数据库公开多站剖面与预算，名义 Reθ 从约 670 到 4060。实际下载的 vel_1410_dns.prof 文件头为 Reθ=1420.960、Reδ*=2030.877、Reτ=492.2115、H₁₂=1.429229、Cf=0.003884512。应以文件头的实际值匹配，不将文件名“1410”当精确雷诺数。该文件还给出均速、三分量 RMS、剪应力、压力及涡量统计。[R20](https://www.mech.kth.se/~pschlatt/DATA/) [原始剖面](https://www.mech.kth.se/~pschlatt/DATA/vel_1410_dns.prof)

**测试量与意义。** 沿 x 比较 Cf、δ₉₉、位移厚度、动量厚度、形状因子和实际 Reθ；在匹配 Reθ 的测站比较 U⁺、应力和谱，不能只在相同几何 x 对比。该例能够检查入口扰动是否在下游异常衰减，统计是否由入口强迫而非自主湍流维持，也是后台阶、隆起、压缩拐角前驱流的必要基础。

**实施与进展。** 多块网格可沿流向分段、壁法向拉伸，但接口处应避免尺度跳变。KTH 页面同时列有更高 Re 的 LES 和早期数据，不能将所有下载项统称为 DNS。首轮建议选中等 Reθ 剖面做空间发展验证，并同时保留入口及若干上游剖面；调整入口直至下游 Cf 匹配而不检查应力谱，可能掩盖合成入口的缺陷。

### C11 逆压梯度湍流边界层

**背景与配置。** 外缘速度沿流向降低后，湍流边界层变厚、应力外峰增强，并可接近分离。它比几何分离问题更便于控制压力历史。应给定完整的外缘压力或速度分布，采用可信湍流入口，并监测 Clauser 参数 β；β 是位移厚度乘外缘压力梯度再除以壁面剪应力。只指定测站的 Reθ 和 β，不能完整定义该流动。

**DNS 基准。** Kitsios 等 2017 的自相似强 APG DNS 总体 Reθ 范围约 570–13800，自相似区域为 10000–12300，平均 β≈39，并与 β≈1 和零压梯度结果比较。强 APG 下应力、产生和耗散出现明显外峰，适合测试从近壁湍流向自由剪切层行为的转变。本次获得论文与配置依据，未获得其完整原始剖面数据库；不另造某个峰值或通用再附位置。[R21](https://doi.org/10.1017/jfm.2017.549)

**测试量与 ILES 意义。** 需比较 Cf、H、外缘速度和 β 的流向历史，内外尺度下的均值／应力、压力 RMS、产生项及外峰位置。强 APG 接近分离时 uτ 变小，单用壁面单位会放大或扭曲误差，应同时采用外缘速度和厚度归一化。该例可检验在槽道调好的离散耗散是否在非平衡区域仍能保留大尺度结构。

**研究进展与实施顺序。** Bobke 等 2017 的高分辨率 LES 强调压力梯度历史：同一局部 β 和 Re 也可能有不同统计。其证据类型为 LES，不能写成另一组 DNS。对无壁面模型粗网格计算，宜先选择较温和、低 Re 的 APG，再以强 APG 为难度上限；复制 β≈39 的完整高 Re 近壁解析计算不应是首个新增目标。[R22](https://doi.org/10.1017/jfm.2017.236)

### C12 可压缩平板湍流边界层与冷热壁

**背景与配置。** 高速壁面湍流即使没有激波，也存在显著密度、黏度及温度变化。该例可将热流和可压缩效应从压缩拐角中独立出来。Larsson 课题组公开了 8 组 Ma=2–10 的 DNS 平均数据，覆盖冷壁、近绝热状态和加热壁，Reδ₂ 约 2000，采用 γ=1.4 的量热完全气体与幂律黏度。实际输入文件一并提供，定义和热壁比值应逐例读取。[R23](https://larsson.umd.edu/data/)

**参考数据及限制。** 数据含沿程壁面量、边界层厚度和下游法向剖面，比只有变换速度图更适合验证。页面推荐引用 Kumar–Larsson 2025。此次核实到数据入口，未完成该压缩包下载，故没有列出未经读取的 Cf、St 或壁温比数值。更早的 Volpiani 等 DNS 数据可作补充，应保持黏度律与归一化方式一致。[R24](https://doi.org/10.1103/2w9d-ys9k)

**测试量与意义。** 比较 Reynolds 与 Favre 均值的差别、全部 Favre 应力、温度／密度／压力 RMS、壁面剪应力和热通量、湍流热通量、总焓剖面，以及速度变换后的偏差。Van Driest 变换不是所有强冷壁条件下的精确规律；只比较变换后的 U⁺，可能掩盖实际温度和热流错误。Stanton 数必须写明参考温差、焓差及热流正号，尤其不能在接近绝热温差为零时机械归一化。

**进展与实施判断。** 2025 年相关研究表明，壁面模型的动量精度与热流精度可能明显不同；引用它的 DNS 作为参考，不意味着本程序要加入该壁面模型。建议先复现 Ma≈2 的一组冷／暖壁配对工况，再进入激波算例。高马赫量热完全气体只是研究隔离条件，不能直接推广到有振动激发、化学反应的真实高焓流。

## 6 自由剪切流

### C13 平面混合层

**背景与设置。** 两股不同速度流体之间的剪切层经历 Kelvin–Helmholtz 失稳、涡配对、三维破碎和自相似增长，是研究卷吸、混合与可压缩增长抑制的经典问题。应先区分时间发展与空间发展：前者通常流向和展向周期、层厚随时间增长；后者有分隔板／入口、层厚沿下游增长，不能直接互换增长率。

**DNS 与实验参考。** Rogers–Moser 1994 的三个不可压缩 DNS 使用来自 Reθ=300 湍流边界层的初始速度场，并改变叠加的二维扰动；部分工况达到可视厚度雷诺数约 20000 的自相似阶段。强初始扰动改变了涡配对和统计演化。这一研究本身说明任意双曲正切均速加白噪声不等于复现其 DNS。本次仅核实论文配置与结果概要，未取得逐点统计数组。[R25](https://ntrs.nasa.gov/citations/19970022168)

**可压缩扩展。** Papamoschou–Roshko 1988 的实验以不同气体和两侧马赫数研究增长率降低，涵盖两股来流 Ma 约 0.2–4。选取具体实验组时须同时给出速度比、密度比、温度、对流马赫数及厚度定义，不使用单一经验增长常数替代实验。对同种气体的简化工况，可用速度差除以两侧声速之和表征对流马赫数，但异密度、异气体的转换应服从原文定义。[R26](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/compressible-turbulent-shear-layer-an-experimental-study/8603B210FE5A2550E5C5784F5D8D48A6)

**测试量与意义。** 需比较动量厚度、涡量厚度、增长率、平均剖面、应力峰及预算、卷吸率、谱与跨层相关。若现有方程没有被动标量，先完成速度验证；不能将可视涡卷吸当作已经验证分子混合。它对剪切层人工耗散、入口扰动及远场反射尤其敏感，直接关联 SD7003 分离剪切层和射流计算。

**实施与进展。** 推荐先建立可重复的低马赫时间发展算例，使用平直结构网格和固定种子，再做空间发展及可压缩扩展。研究由可视大涡逐步转向自相似统计、初始条件记忆和可压缩预算；长时间或宽域是否达到相似状态，必须由增长率及无量纲剖面共同判断。

### C14 亚声速圆射流

**背景与配置。** 圆喷口下游有环形剪切层、势核、卷吸、三维涡破碎及自相似远场，并产生宽带噪声。Brès 等 2018 研究等温 Ma=0.9、以出口直径定义 Reᴅ=10⁶ 的收缩加直管喷嘴，结合实验验证 LES，强调喷口边界层状态。其数值方法在喷嘴内使用壁面模型；不能把该网格直接用于无壁面模型计算而期待相同代价和精度。[R27](https://doi.org/10.1017/jfm.2018.476)

**参考结果及获取方式。** 参考包含 PIV／速度测量和声学结果，应优先提取同一喷口状态下的中心线速度衰减、径向均值及应力、势核长度和声压谱。原论文全文已保存，但本次未获得这些量的完整原始数组，因此不指定未经核实的通用势核长度、衰减常数或声压级。该高 Re 工况的 LES 是辅助解释，实验才是独立物理验证依据。[作者公开全文](https://atowne.com/wp-content/uploads/2018/07/BresEtAl_2018_JFM_jetLES.pdf)

**测试量与意义。** 除平均速度和扩张率，还需入口边界层厚度、形状因子、湍流强度、环向模态能量、两点相关、近场压力谱；若验证噪声，还要统一测角、距离、频率归一化和声学传播处理。层流喷口与湍流喷口可能得到截然不同的初始卷起、混合及噪声，仅匹配 Ma 和 Re 不够。

**网格与进展。** 多块结构网格可用中心方块加环形块，沿剪切层与势核加密，并配置足够宽远场与缓冲区。小尺度湍流及声波对数值耗散的要求不同，应分别验证。研究已由单点谱扩展到波包、相干结构及谱 POD；这些分析应建立在足够长且均匀采样的压力／速度时序上。对当前无壁面模型程序，圆射流宜列为中后期目标，先解决可信入口和出流声反射。

## 7 分离再附与钝体尾迹

### C15 后台阶流动

**背景与配置。** 后台阶把分离位置固定在几何突变处，减少光滑表面分离点的不确定性，可集中检验剪切层增长、再附、角区回流及恢复。应分成两个独立基准：Le–Moin–Kim 1997 DNS 为 Reₕ=5100、扩张比 1.20；Driver–Seegmiller 实验经 TMR 修正后的 Reₕ 约 36000，上游 Reθ 约 5000，Ma 约 0.128、入口边界层厚度约 1.5h。不要把后者继续标成早期部分资料中的 50000。[R28](https://doi.org/10.1017/S0022112096003941) [R29](https://tmbwg.github.io/turbmodels/backstep_val.html)

**可核验结果。** 低 Re DNS 的平均再附长度为 6.28h，并报告约 0.06 的剪切层特征 Strouhal 数。高 Re 实验的再附位置为 xᵣ/h=6.26±0.10，来自壁面摩擦测量解释。两值接近不代表应力、近壁回流和入口条件相同。TMR 公开 Cp、Cf 及多个速度／应力剖面；本次已经下载 Cp 和带不确定度的 Cf 表。稀疏 Cf 点的简单线性插值约得 6.279，不应取代作者给出的 6.26±0.10。[R28](https://doi.org/10.1017/S0022112096003941) [R29](https://tmbwg.github.io/turbmodels/backstep_val.html)

**测试量与意义。** 比较再附长度、负 Cf 峰、压力恢复、回流区速度和应力、剪切层厚度及 x/h=1、4、6、10 等剖面。低 Re DNS 表明再附后相当远处仍未完全恢复普通平板壁律；不应把远下游偏离对数律一概判为错误。对 ILES，后台阶比周期山更敏感于入口湍流是否具有正确长度尺度。

**结构网格与进展。** 可用台阶上下游 H 型块实现，是适合当前程序的新增分离基准。首轮推荐 Reₕ=5100 的 DNS 路线，待前驱边界层成熟后再挑战实验高 Re 工况。顶壁角度、无滑移／滑移条件和扩张比必须匹配具体实验；TMR 的二维 RANS 设置不能直接充当三维 ILES 的入口与展宽配置。

### C16 NASA 壁面隆起与光滑表面分离

**背景与实验设计。** NASA wall-mounted hump 采用 Glauert–Goldschmied 型隆起，先加速后形成强逆压梯度，后缘附近发生分离。基线无流动控制工况弦长约 420 mm，U∞≈34.6 m/s、Ma≈0.1；TMR 使用 Re꜀=936000，Greenblatt 等论文报告约 929000。入口 x/c=−2.14 附近边界层厚度约 0.08c。端板和风洞堵塞效应会改变外缘压力，应采用对应几何或经说明的等效顶部轮廓。[R30](https://tmbwg.github.io/turbmodels/nasahump_val.html)

**实验数据。** TMR 提供无控制的 Cp、Cf、入口剖面、下游速度／应力及 PIV。已下载 Cf 表的两个原始记录为：x/c=1.09，Cf=−0.00038、不确定度 0.00004；x/c=1.15，Cf=0.00049、不确定度 0.00004。两点线性插值给出再附约 1.116c，仅为当前表的派生估计。分离区前缘附近测点较稀，不能用跨越很长距离的两点线性插值伪造精确分离点。[原始 Cf 数据](https://tmbwg.github.io/turbmodels/Nasahump_validation/noflow_cf.exp.dat)

**测试量与意义。** 比较 Cp、带误差条的 Cf、再附、剪切层中心位置、回流强度和 x/c=0.65、0.8、0.9、1.0、1.1、1.2、1.3 的剖面。与后台阶不同，本例分离由压力梯度与湍流共同决定；它更直接考查应力对分离泡长度的反馈。压力场吻合而再附错误，通常说明仅捕捉了外流几何效应。

**实施与新进展。** 曲面挤出的多块结构网格可实现该例，但三维展宽、风洞顶部及入口发展不能省略。它的 Re 较高，当前无壁面模型路线成本较大。2025–2026 年 Gaussian bump 研究进一步考查光滑分离的三维性、低频运动和侧壁约束；作者预印本指出有限展宽结构可能解释部分实验／计算差异。因此，小展宽周期计算不一定代表有限风洞，即使中心线均值接近。[R31](https://arxiv.org/abs/2512.13582) [R32](https://arxiv.org/abs/2603.05765)

### C17 雷诺数 3900 圆柱绕流

**背景与实验设计。** Reᴅ=3900 位于亚临界圆柱绕流区，壁面边界层与分离剪切层的转捩、周期脱涡和三维尾迹相互作用。Parnaudeau 等 2008 同时开展热丝、PIV 和 LES。其实验具有有限阻塞和长径比；二维计算无法代表这里的三维尾迹湍流，LES 部分也不能误称 DNS。[R33](https://doi.org/10.1063/1.2957018)

**原论文表值。** 该论文表 II 的热丝测量 St=0.208±0.002；PIV 回流区长度 Lᵣ/D=1.51，中心线最小平均速度 Umin/U꜀=−0.34。论文高分辨率 LES 分别给出 St=0.208±0.001、Lᵣ/D=1.56、Umin/U꜀=−0.26，应单独列为数值对照。Lᵣ 从圆柱后表面起算，而部分图的 x 从圆心起算，比较时需处理 0.5D 的坐标差。[作者公开论文](https://www.researchgate.net/profile/Dominique-Heitz/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900/links/0912f50a3828de19d1000000/Experimental-and-numerical-studies-of-the-flow-over-a-circular-cylinder-at-Reynolds-number-3900.pdf)

**测试量与意义。** 需要升阻力均值／RMS、St、表面 Cp、分离角、尾流回流长度、中心线速度、横向剖面与应力、展向相关和速度谱。St 正确但应力和回流长度错误仍可能发生，不能只用脱涡频率一个数判定通过。对于曲面 WCNS，该例还暴露网格非正交、尾迹块接口与数值扰动引起的转捩差异。

**网格与进展。** 可采用贴体 O 型近场与尾迹块；展宽及边界条件应跟选定文献保持一致并做加宽检查。已有文献对紧邻圆柱的统计存在显著差别，来源包括入口扰动、PIV 空间分辨率、阻塞和统计长度。因此本报告推荐以一个设施的成套数据为主，不混合不同实验中最容易吻合的指标。流动控制等后续研究也应先通过无控制基线。

## 8 激波与湍流相互作用

### C18 平面激波与均匀各向同性湍流

**背景与配置。** 让来流湍流穿过平均平面激波，在没有壁面分离、入口边界层和复杂几何的条件下检验激波后湍动能放大、各向异性、声学及熵扰动。Larsson–Lele 2009 的 DNS 研究包含上游 Reλ≈40、平均激波马赫数约 1.3–6、湍流马赫数约 0.16–0.38 的系列，后续研究进一步考查强度与雷诺数影响。每个系列工况应独立指定，不从这些范围自由拼接参数。[R34](https://doi.org/10.1063/1.3275856) [R35](https://doi.org/10.1017/jfm.2012.573)

**数据与统计定义。** 作者公开 ASCII 平均数据压缩包，含普通均值、Favre 均值及应力。页面明确应力单位为速度平方，没有再乘密度，且由有限三维快照平均得到，仍有可见采样噪声。本次核实下载入口，但压缩包下载未成功，故不虚列湍动能放大倍数。正式复现前应锁定一个数据文件并读取头部，确认激波位置、上游谱和归一化。[R23](https://larsson.umd.edu/data/)

**测试量与 ILES 意义。** 先用无扰动激波检查 Rankine–Hugoniot 跳跃及守恒，再检查湍流相互作用。湍流平均跳跃还包含脉动输运的贡献，不应强制与无湍流关系完全相同。主要比较三个法向应力、湍动能增益、各向异性、密度／压力 RMS、激波起伏与下游恢复谱。同步输出激波离散厚度和传感器，检查激波捕捉耗散是否污染上下游涡结构。

**实施与进展。** 计算域可采用激波法向开边界、横向周期的笛卡尔分块网格，避免曲面和壁面干扰。难点在可信且可持续的上游湍流供给、非反射边界和平均激波位置控制。它非常适合在压缩拐角之前使用：如果此例已经严重压低激波后的旋度与应力，复杂拐角中的误差就不宜仅归因于近壁网格。

### C19 斜激波与湍流平板边界层

**背景与配置。** 外部斜激波照射平板边界层，可在几何平直的壁面上研究分离、反射激波和低频运动。Bernardini 等 2023 DNS 采用 Ma=2.28、激波发生器转角 8°、域大小 70δ₀×12δ₀×6.5δ₀。无滑移壁的温度设为上游恢复温度，入口使用湍流生成，上边界局部施加激波跳跃。这里的恒定恢复温度壁不能在代码中不加核对地换成零热流壁。[R36](https://doi.org/10.1017/jfm.2022.1038)

**可核验结果。** 原文表 1 的 DNS 上游 Reθ=6882、以壁面黏度定义的 Reδ₂≈3900、Reτ≈1100、Cf=1.98×10⁻³；对照实验 Reθ=5100，并非完全相同雷诺数。表 2 给出 DNS 相互作用长度 Lint/δ₀=3.30、分离尺度 Lsep/δ₀=2.16；实验相应为 4.18 和约 3.40。后者部分来自长度换算，不能视作与 DNS 同精度的直接测量。文中低频峰约 fLint/U∞=0.04。[R36](https://doi.org/10.1017/jfm.2022.1038)

**定义与测试意义。** 原文对分离尺度采用包含间歇分离区的特定定义，不完全等于“平均 Cf 负值区间长度”；必须先统一判据。比较入口均值／应力、Cp、Cf、压力 RMS、平均及瞬时分离面积、激波足 PDF、频谱、跨展向相干和条件平均。已发表较长时间 DNS 显示低频活动具有间歇性，用很短时窗的单个 FFT 峰值评价误差并不可靠。

**实施与进展。** 该例以平直网格隔离激波入射条件，适合作为压缩拐角的交叉检查。后续微型斜坡控制研究增加了三维涡和侧向调制，应先完成无控制基线。当前无壁面模型程序宜选择较低 Re、数据完整的同系列条件启动，再向上述 Reτ≈1100 工况推进；本报告没有把该高分辨率 DNS 网格数量作为可在本机长算的建议。

## 9 怎样用这些算例验证粗网格 ILES

### 9 1 先定义比较对象

ILES 的有效耗散由重构、Riemann 通量、低马赫修正、传感器、时间推进、网格及边界共同决定。切换其中任一项，都可能改变隐式模型。因此一次验证记录应包括程序版本／源码哈希、完整数值选项、网格哈希、编译器、并行分区和精度；不能只写“turbulence.model=none”。不同算例可有合理的激波开关，但需要在事前确定使用原则，避免逐个调节耗散直到匹配参考。

DNS 的全部小尺度方差不等于粗网格上的解析方差。验证二阶统计时，建议同时提供两种比较：与未滤波 DNS 比较总量偏差；在获得三维场时，对 DNS 施加明确的空间滤波或单元体积平均，再比较解析尺度统计。对于有限体积量，参考取样应考虑体积平均；对曲线网格，应说明物理空间的滤波与插值。ILES 没有唯一已知的显式滤波核，不应把任意一次滤波后的吻合当成严格证明。

实验的 PIV 窗口、热丝长度和传感器带宽也起滤波作用。应保留测量分辨率，尽量让计算统计与实验观测算子一致。仅能取得论文曲线时，数字化误差和原实验误差分别保留；图像坐标的多位小数不能变成物理精度。

### 9 2 统一统计定义

| 统计项 | 应保存的内容 | 容易混淆的地方 |
|---|---|---|
| 时间平均 | 开始时刻、累计物理时间、样本数 | 变时间步应按实际 dt 加权 |
| 普通应力 | 一阶矩和全部六个速度二阶矩 | 先局部乘积后平均 |
| Favre 应力 | 密度、一阶质量通量、密度加权二阶矩 | 与普通 Reynolds 应力分开 |
| 壁面量 | 压力、切向剪应力、法向热流和壁温 | 物理切向与坐标分量不同 |
| 归一化 | 参考速度、密度、黏度及长度 | 半高／全高、半径／直径不同 |
| 谱 | 采样间隔、窗函数、窗长、重叠、PSD 单位 | 空间谱与时间谱不能直接混用 |
| 分离与再附 | 判据、零点插值、坐标原点、误差 | 瞬时反流与平均 Cf 零点不同 |
| 守恒与耗散 | 质量、动量、能量、源项功及边界通量 | 开域流动不能照搬周期盒预算 |

普通协方差等于速度乘积的平均减去两个平均速度的乘积。Favre 均值是密度加权速度平均除以平均密度，Favre 应力须用相同密度权重形成。若先做展向速度平均，再对其计算方差，会丢失展向湍流能量，尤其损害周期山、翼型和压缩拐角的比较。

Reτ 使用摩擦速度和明确的特征高度／半径，实际值由统计壁剪获得。Cf 为有符号切向剪应力除以选定参考动压，Cp 为压力相对基准值除以该动压。斜壁热流和剪应力须沿物理法向／切向计算。将壁面求导后处理与求解器黏性面通量作一致性核对，比只检查输出公式更有价值。

谱分析必须有足够均匀的物理时间采样。变步长历史先作有说明的重采样和抗混叠处理，不能直接按样本编号 FFT。既要给主频，也要给频率分辨率、可用低频区间与置信区间。升阻力、激波位置、压力及流量可以连续高频保存，全三维场则按储存预算较疏输出。

### 9 3 分清统计误差与分辨率误差

先剔除启动、转捩与大尺度调整阶段，再用若干不重叠时间块检查均值、二阶量和关键壁面指标的漂移。充分发展判据应至少覆盖流量／驱动力、摩擦、动能及一项脉动统计。大量相邻时间步并非大量独立样本；应依据自相关时间或分块统计估计有效样本数。对于低频激波运动，需要覆盖许多低频周期，而不是只覆盖许多声学 CFL 步。[R37](https://doi.org/10.1063/1.4866813)

建议每个主指标给出时间分块置信区间和加长统计窗的稳定性；实验误差、参考 DNS 统计误差与本次采样误差应分别列出。初步工程筛查可将主要整体量的几百分点偏差列为关注区间，但本报告不把统一“误差小于 5%”设为所有算例的合格线：弱二次流、零点、热流峰和频谱需要各自的绝对／相对误差尺度。

网格研究至少改变三个方向的分辨率并检查周期盒尺寸；必须记录第一层中心 y⁺、流向／展向间距、局部网格拉伸及关键剪切层厚度内的单元数。对高阶方法，还要区分单元数、自由度和实际谱分辨率。同样 100 万单元的不同方法、不同块拓扑或不同近壁分配，不具有自动可比的有效分辨率。

ILES 加密时隐式耗散也改变，统计误差未必单调。可以报告趋势和有限分辨率范围内的可信度；在未确认渐近收敛、统计误差足够小之前，不机械套用 Richardson 外推或给出“网格无关”的结论。最有说服力的证据是：同一数值设置在多个独立物理机制下，同时改善均值、应力、谱和预算。

### 9 4 建议的实施和验证顺序

| 阶段 | 优先对象 | 新增工作重点 | 进入下一阶段的依据 |
|---|---|---|---|
| A | TGV 与衰减 HIT | 可重复初值、能谱和动能预算 | 多分辨率耗散与谱趋势可信 |
| B | 已有槽道与方管 | 壁剪、应力、二次流及统计误差 | 湍流持续、近壁与各向异性合理 |
| C | 零压梯度边界层和低 Re 后台阶 | 前驱入口、开边界与再附统计 | 入流状态及剖面同时符合参考 |
| D | 周期山、圆柱、SD7003 | 曲线多块接口、展宽和转捩 | 分离／尾流与二阶量一致 |
| E | 可压缩边界层和激波 HIT | Favre 统计、热流、局部激波耗散 | 无激波热流与无壁面激波响应可信 |
| F | 斜激波及压缩拐角 | 长时间低频统计、再附加热 | 入口、均值、脉动及热流分别验证 |

圆管可与 B 并行作为曲率测试；强迫 HIT 适合在 A 后增加强迫与稳态统计；混合层可与 C 的入口开发配合。高 Re NASA 隆起、强 APG 和圆射流宜在基础能力可靠且服务器资源评估完成后安排。这里的顺序是研发建议，不是本次已启动的计算计划。

## 10 参考数据清单与可追溯性

本次随报告保存若干体量较小的公开原始数据和来源清单。sources/manifest.json 记录获取日期、来源 URL、文件大小、SHA-256，以及未成功获取的项目。原始文件不改写列名、数值或不确定度；报告中的派生值用独立脚本重算。大体量三维数据库和原作者未公开的时序没有下载。

| 数据 | 本次状态 | 可立即核验的内容 |
|---|---|---|
| 工作坊 TGV 规范与 data.tgz | 已保存并读取 | K、ε、Ω 时间历程及离散峰值 |
| CBC 1971 原论文 | 已保存并核对表 4 页面 | 三测站 RMS、耗散、尺度与 Reλ |
| Lee–Moser 180 均值及应力表 | 已保存并读取文件头 | 实际 Reτ、U、应力剖面 |
| KTH vel_1410_dns.prof | 已保存并读取 | Reθ、H、Cf、速度及应力 |
| NASA 后台阶 Cp／Cf | 已保存 | 曲线、Cf 不确定度 |
| NASA 隆起 Cp／Cf | 已保存 | 曲线、Cf 不确定度和再附估计 |
| JHTDB 1024³ README | 已保存并读取 | 原始／扩展时间窗及统计 |
| 圆管 turbpipe 仓库 | 已确认入口 | 具体数据文件尚未逐项下载 |
| 方管、APG、混合层 | 已核实论文 | 完整逐点数组仍待获取 |
| 圆射流 2018 与斜激波 2023 论文 | 已保存全文 | 配置及公开表值；没有完整原始时序 |
| Larsson 激波 HIT 数据包 | 已确认入口但下载失败 | 尚不能声称已持有原始数组 |
| Larsson 可压缩边界层数据包 | 已确认入口 | 尚未下载并逐文件核对 |
| 原有 SD7003 与压缩拐角数据 | 引用既有目录 | 表值／图像数字化，精度有限 |

TMR 已于 2026 年 2 月 24 日迁移到 TMBWG 的 GitHub 站点，报告采用新地址。该站同时包含 RANS 代码校核、实验、DNS 和 LES，使用时必须识别具体栏目；两个 RANS 程序互相吻合属于实现校核，不是湍流物理真值。[R38](https://tmbwg.github.io/turbmodels/index.html)

进一步扩展可考虑粗糙壁槽道、旋转槽道、曲率边界层、三维交汇流，以及 Rayleigh–Taylor／Richtmyer–Meshkov 混合。它们分别检验粗糙度、旋转、曲率、角区流动和变密度混合，但会引入新的几何／体力／多组分或初始界面问题。本轮优先保留单组分黏性流中与当前程序最直接相关的 19 类，不将这些扩展仅凭名称加入核心验收集。

## 11 参考文献与数据入口

以下按正文编号排列。研究进展检索截至 2026 年 10 月 3 日；明确标为预印本的内容不冒充已经同行评审的结论。页面获取失败不代表资料不存在，具体下载状态以第 10 节及本地清单为准。

R01　Kim, J., Moin, P., Moser, R. 1987. Turbulence statistics in fully developed channel flow at low Reynolds number. Journal of Fluid Mechanics 177, 133–166. 数据整理见 [ERCOFTAC Case 032](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case032)。本地 simul1.dat 为当前槽道主对照。

R02　Lee, M., Moser, R. D. 2015. Direct numerical simulation of turbulent channel flow up to Reτ≈5200. Journal of Fluid Mechanics 774, 395–415. [DOI 10.1017/jfm.2015.268](https://doi.org/10.1017/jfm.2015.268)；[数据及说明](https://turbulence.oden.utexas.edu/channel2015/content/README_2015.html)。

R03　Fröhlich, J., Mellen, C. P., Rodi, W., Temmerman, L., Leschziner, M. A. 2005. Highly resolved large-eddy simulation of separated flow in a channel with streamwise periodic constrictions. Journal of Fluid Mechanics 526, 19–66. [DOI 10.1017/S0022112004002812](https://doi.org/10.1017/S0022112004002812)；[TMR LES 数据](https://tmbwg.github.io/turbmodels/Other_LES_Data/2dhill_periodic.html)。

R04　Balakumar, P. 2015. DNS/LES Simulations of Separated Flows at High Reynolds Numbers. AIAA 2015-2783. [DOI 10.2514/6.2015-2783](https://doi.org/10.2514/6.2015-2783)；[Reₕ=2800 可压缩 DNS 数据页](https://tmbwg.github.io/turbmodels/Other_DNS_Data/2dhill_periodic_compress.html)。

R05　Breuer, M., Peller, N., Rapp, C., Manhart, M. 2009. Flow over periodic hills—Numerical and experimental study in a wide range of Reynolds numbers. Computers & Fluids 38, 433–457. [出版者页面](https://www.sciencedirect.com/science/article/pii/S0045793008001126)。

R06　International Workshop on High-Order CFD Methods. C3.3 SD7003 airfoil benchmark. [算例规范](https://cfd.ku.edu/hiocfd/case_c3.3.html)；[UIUC 原始翼型坐标](https://m-selig.ae.illinois.edu/ads/coord/sd7003.dat)。

R07　Galbraith, M. C., Visbal, M. R. 2010. Implicit Large Eddy Simulation of Low-Reynolds-Number Transitional Flow Past the SD7003 Airfoil. AIAA 2010-4737. [DOI 10.2514/6.2010-4737](https://doi.org/10.2514/6.2010-4737)。

R08　Porter, K. M., Poggie, J. 2019. Selective upstream influence on the unsteadiness of a separated turbulent compression ramp flow. Physics of Fluids 31, 016104. [DOI 10.1063/1.5078938](https://doi.org/10.1063/1.5078938)。

R09　Poggie 课题组. Compression ramp benchmark and publications. [作者算例页](https://engineering.purdue.edu/~jpoggie/ramp/index.html)；[2017 年前期论文](https://engineering.purdue.edu/~jpoggie/papers/AIAA-2017-0533.pdf)。

R10　International Workshop on High-Order CFD Methods. C3.5 Direct Numerical Simulation of the Taylor–Green Vortex at Re=1600. [原始规范](https://cfd.ku.edu/hiocfd/case_c3.5.pdf)。

R11　同一工作坊. Taylor–Green vortex reference data. [谱方法参考包](https://cfd.ku.edu/hiocfd/data.tgz)。本次读取文件 spectral_Re1600_512.gdiag。

R12　Comte-Bellot, G., Corrsin, S. 1971. Simple Eulerian time correlation of full- and narrow-band velocity signals in grid-generated, isotropic turbulence. Journal of Fluid Mechanics 48, 273–337. [DOI 10.1017/S0022112071001599](https://doi.org/10.1017/S0022112071001599)；[论文全文](https://courses.washington.edu/mengr544/handouts/comtebellot-corrsin-jfm-71.pdf)。

R13　Johns Hopkins Turbulence Databases. Forced isotropic turbulence dataset extended. [参数及两个统计窗的原始 README](https://turbulence.pha.jhu.edu/docs/README-isotropic.pdf)。

R14　Johns Hopkins Turbulence Databases. [数据与访问工具入口](https://turbulence.idies.jhu.edu/database)；[原站数据目录](https://turbulence.pha.jhu.edu/datasets.aspx)。

R15　Yao, J., Rezaeiravesh, S., Schlatter, P., Hussain, F. 2023. Direct numerical simulations of turbulent pipe flow up to Reτ≈5200. Journal of Fluid Mechanics 956, A18. [DOI 10.1017/jfm.2022.1013](https://doi.org/10.1017/jfm.2022.1013)；[turbpipe 数据仓库](https://dataverse.tdl.org/dataverse/turbpipe)。

R16　Pirozzoli, S., Romero, J., Fatica, M., Verzicco, R., Orlandi, P. 2021. One-point statistics for turbulent pipe flow up to Reτ≈6000. Journal of Fluid Mechanics 926, A28. [DOI 10.1017/jfm.2021.727](https://doi.org/10.1017/jfm.2021.727)。

R17　Pirozzoli, S., Modesti, D., Orlandi, P., Grasso, F. 2018. Turbulence and secondary motions in square duct flow. Journal of Fluid Mechanics 840, 631–655. [DOI 10.1017/jfm.2018.66](https://doi.org/10.1017/jfm.2018.66)。

R18　Modesti, D., Pirozzoli, S., Orlandi, P., Grasso, F. 2018. On the role of secondary motions in turbulent square duct flow. [DOI 10.1017/jfm.2018.391](https://doi.org/10.1017/jfm.2018.391)。

R19　Verolini, F., Xiao, M., Pirozzoli, S. 2026. Turbulent flow in a square duct at extreme Reynolds number. Journal of Fluid Mechanics，2026 年 9 月 10 日在线发表. [出版者全文](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/turbulent-flow-in-a-square-duct-at-extreme-reynolds-number/C3DAB1AC62E975591CFD93D57C53E31E)。该研究使用 WMLES。

R20　Schlatter, P., Örlü, R. 2010. Assessment of direct numerical simulation data of turbulent boundary layers. Journal of Fluid Mechanics 659, 116–126. [作者维护的数据与文献目录](https://www.mech.kth.se/~pschlatt/DATA/)。

R21　Kitsios, V., Sekimoto, A., Atkinson, C., et al. 2017. Direct numerical simulation of a self-similar adverse pressure gradient turbulent boundary layer at the verge of separation. Journal of Fluid Mechanics 829, 392–419. [DOI 10.1017/jfm.2017.549](https://doi.org/10.1017/jfm.2017.549)；[作者预印本](https://arxiv.org/abs/1708.03039)。

R22　Bobke, A., Vinuesa, R., Örlü, R., Schlatter, P. 2017. History effects and near equilibrium in adverse-pressure-gradient turbulent boundary layers. Journal of Fluid Mechanics 820, 667–692. [DOI 10.1017/jfm.2017.236](https://doi.org/10.1017/jfm.2017.236)。

R23　Larsson Computational Turbulence Laboratory. [激波与湍流及可压缩边界层数据页](https://larsson.umd.edu/data/)。含数据头定义与下载链接。

R24　Kumar, V., Larsson, J. 2025. Improved heat flux modeling for high-speed wall-modeled large eddy simulation. Physical Review Fluids 10, 124605，2025 年 12 月 24 日发表. [DOI 10.1103/2w9d-ys9k](https://doi.org/10.1103/2w9d-ys9k)。本文采用的是作者页公开 DNS 参考，不是要求移植该壁面模型。

R25　Rogers, M. M., Moser, R. D. 1994. Direct simulation of a self-similar turbulent mixing layer. Physics of Fluids 6, 903–923. [NASA 原始记录](https://ntrs.nasa.gov/citations/19970022168)；[DOI 10.1063/1.868325](https://doi.org/10.1063/1.868325)。

R26　Papamoschou, D., Roshko, A. 1988. The compressible turbulent shear layer an experimental study. Journal of Fluid Mechanics 197, 453–477. [DOI 10.1017/S0022112088003325](https://doi.org/10.1017/S0022112088003325)。

R27　Brès, G. A., Jordan, P., Jaunet, V., et al. 2018. Importance of the nozzle-exit boundary-layer state in subsonic turbulent jets. Journal of Fluid Mechanics 851, 83–124. [DOI 10.1017/jfm.2018.476](https://doi.org/10.1017/jfm.2018.476)。

R28　Le, H., Moin, P., Kim, J. 1997. Direct numerical simulation of turbulent flow over a backward-facing step. Journal of Fluid Mechanics 330, 349–374. [DOI 10.1017/S0022112096003941](https://doi.org/10.1017/S0022112096003941)。

R29　Driver, D. M., Seegmiller, H. L. 1985. Features of a reattaching turbulent shear layer in divergent channel flow. AIAA Journal 23, 163–171. [DOI 10.2514/3.8890](https://doi.org/10.2514/3.8890)；[TMR 修正工况与实验数据](https://tmbwg.github.io/turbmodels/backstep_val.html)。

R30　Greenblatt, D., Paschal, K. B., Yao, C.-S., Harris, J., Schaeffler, N. W., Washburn, A. E. 2006. Experimental investigation of separation control Part 1 Baseline and steady suction. AIAA Journal 44, 2820–2830. [DOI 10.2514/1.13817](https://doi.org/10.2514/1.13817)；壁剪研究 Naughton, Viken, Greenblatt 2006，[DOI 10.2514/1.14192](https://doi.org/10.2514/1.14192)；[TMR 实验数据页](https://tmbwg.github.io/turbmodels/nasahump_val.html)。

R31　Klopsch, R., Fuchs, L. M., Rigas, G., Oberleithner, K., von Saldern, J. G. R. 2025. Spectral analysis of attached and separated turbulent flows over a Gaussian-shaped bump. [作者预印本 arXiv 2512.13582](https://arxiv.org/abs/2512.13582)。此处按所核验预印本版本引用。

R32　Manohar, K. H., Annamalai, H., Williams, O., Morton, C., Martinuzzi, R. J. 2026. Unsteadiness in turbulent separated flow over a three-dimensional Gaussian bump. [作者预印本 arXiv 2603.05765](https://arxiv.org/abs/2603.05765)。正文仅用于说明低频与三维机制的发展，不作为新的已完成验收数据集。

R33　Parnaudeau, P., Carlier, J., Heitz, D., Lamballais, E. 2008. Experimental and numerical studies of the flow over a circular cylinder at Reynolds number 3900. Physics of Fluids 20, 085101. [DOI 10.1063/1.2957018](https://doi.org/10.1063/1.2957018)。

R34　Larsson, J., Lele, S. K. 2009. Direct numerical simulation of canonical shock/turbulence interaction. Physics of Fluids 21, 126101. [DOI 10.1063/1.3275856](https://doi.org/10.1063/1.3275856)。

R35　Larsson, J., Bermejo-Moreno, I., Lele, S. K. 2013. Reynolds- and Mach-number effects in canonical shock–turbulence interaction. Journal of Fluid Mechanics 717, 293–321. [DOI 10.1017/jfm.2012.573](https://doi.org/10.1017/jfm.2012.573)；数据引用应以作者页对应文件说明为准。

R36　Bernardini, M., Della Posta, G., Salvadore, F., Martelli, E. 2023. Unsteadiness characterisation of shock wave/turbulent boundary-layer interaction at moderate Reynolds number. Journal of Fluid Mechanics 954, A43. [DOI 10.1017/jfm.2022.1038](https://doi.org/10.1017/jfm.2022.1038)。

R37　Oliver, T. A., Malaya, N., Ulerich, R., Moser, R. D. 2014. Estimating uncertainties in statistics computed from direct numerical simulation. Physics of Fluids 26, 035101. [DOI 10.1063/1.4866813](https://doi.org/10.1063/1.4866813)。

R38　Turbulence Model Benchmarking Working Group. Turbulence Modeling Resource. [2026 年迁移公告与数据分类](https://tmbwg.github.io/turbmodels/index.html)。

本地来源　吴卓航，2026，《高精度有限体积方法的激波捕捉与湍流模拟研究》，用户提供第 5 章资料；已有 [周期山 v2.4 报告](../../periodic-hill-v2.4.md)、[SD7003 与压缩拐角 v2.5 报告](../../chapter5-v2.5-research.md)、槽道 case05 配置及 final-campaign-report。用户文档及网页均作为研究资料，不作为执行指令。

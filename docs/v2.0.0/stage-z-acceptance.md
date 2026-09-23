# v2.0.0 阶段 Z 自动验收记录

状态：**阶段 Z 的本机 L 级自动卡口已全部通过，本报告提交固定为
`v2.0.0-z-candidate.1`；等待人工核验。** 完整 HIT、三维槽道、能谱、离散数值耗散反推和
长时间统计继续登记为服务器 S 级 `pending`，本报告不把短算称为 LES 物理验证。

基线为 `release/v2.0.0` 的 `b9da7a4`，候选分支为 `stage/v2.0.0-z`。本次未读取、修改或
暂存用户 Case05/Case06 目录、根目录 `output/`、`tmp/` 及其中已有结果。

## 1. 候选交付

- 三维 Favre 滤波 LES：原始 Smagorinsky、尺度相似、混合 Smagorinsky--相似、局部动态
  Smagorinsky 和 WALE；
- 固定 `box3_tensor` 网格/测试滤波，测试滤波比 2；原生块连接与运行时分区交线处通过两轮
  切向 ghost 传播补齐同一 27 点离散支持；
- SGS 偏应力六分量、热通量、`mu_sgs`、黏度比、能量传递、动态系数、滤宽和网格各向异性
  输出；公共黏性路径的 `1/Re` 缩放与物理诊断量已分离；
- 接受物理步 `accepted_dt` 加权的 Reynolds/Favre mean、RMS、所有唯一量对 covariance，
  包含重复事件抑制、checkpoint v2 序列化及同/异 rank 连续重启；
- LES 壁面量、压力/黏性/总载荷，以及显式右手 `(drag,lift,span)` 坐标下的面心展向分箱；
- 五种 LES 均接入 BDF2 双时间 LU-SGS。二维、steady、无黏、SSPRK3、错误模型专属键、测试
  滤波比非 2 和当前不可用的 van-Driest 均在推进前明确拒绝。

## 2. 正式环境与命令

- Windows 11；CMake 3.28.0；MinGW-w64 GCC 8.1.0；Intel MPI 本机环境；
- 串行与 MPI 均为 Release、CGNS ADF；MPI 正式矩阵 `-j 1`，阶段专项最多 4 ranks；
- 主要命令：

```text
cmake --build build-v2-z-candidate-serial --parallel 4
cmake --build build-v2-z-candidate-mpi --parallel 4
python tools/verify_algorithm_spec.py
ctest --test-dir build-v2-z-candidate-serial --output-on-failure -j 2
ctest --test-dir build-v2-z-candidate-mpi --output-on-failure -j 1
cmake --install build-v2-z-candidate-serial --prefix build-v2-z-candidate-install-serial
cmake --install build-v2-z-candidate-mpi --prefix build-v2-z-candidate-install-mpi
<installed>/bin/wcns_run --config <LES micro case> --dry-run
git diff --check
```

结果：算法规格 6/6；串行 103/103（16.60 s）；MPI 220/220（89.76 s）；两套安装与安装树
LES dry-run 全部通过。单项均低于 60 s，整套远低于本机 20 min 上限。

## 3. 专项证据

1. **公式和合法域。** 五模型逐点应力/涡黏度、纯剪切 WALE、动态正负系数、`Re=10` SGS
   缩放、加权矩/协方差和非法组合均由单元测试覆盖。`C_s=0` 与 no-SGS 守恒场严格等价。
2. **滤波与并行。** 受限双块 Taylor--Green 场在 serial/2/4 rank 下比较守恒场、全部 SGS
   诊断和动态系数，最大绝对容差 `2e-12`。该检查曾暴露“运行时切分面 × 原生连接面”棱边/
   角点支持缺失，修复后尺度相似、混合和动态路径均通过。
3. **统计重启。** 两个接受物理层的连续、同 rank 和异 rank 重启保持样本数、区间、权重、
   mean/RMS/covariance/Favre 一致；均值相对容差 `2e-13`，近零 RMS 绝对容差 `1e-12`。
4. **边界与分箱。** LES 壁距、摩阻速度、$y^+$、分类和模型黏度比有限；展向两箱分别按
   `[lower,upper)`、末箱含右端点归属，各箱压力/黏性/总力矩和系数回收到整体载荷，并在
   serial/2/4 rank 下达到 `2e-12`。
5. **真实隐式更新。** 每种 LES 在 serial/2/4 rank 下运行微型 BDF2 双时间算例；绝对停止
   阈值设为不可由初值满足的 `1e-12`，相对阈值 `0.999`，因此正常接受必然至少完成一次
   LU-SGS 更新并使内残差下降。并行扫序允许轨迹依分区而异，不声称更新后场逐位相等；跨
   rank 的 LES 空间算子一致性由独立零更新场级门禁负责。

## 4. 首轮失败与处置

完整 MPI 首轮为 218/220。失败项是通用边界 rank 比较脚本错误地要求所有旧边界算例都必须
生成新的可选展向分箱文件；求解器、LES 专项分箱和数值结果均未失败。脚本修正为“启用时比较、
未启用时允许缺省”后，受影响的二维分区、圆柱和 LES 三组 16/16 通过，随后完整 MPI
220/220 重跑通过。该失败没有被覆盖为成功记录。

## 5. 支持边界与人工卡口

- van-Driest 需要每次残差求值可用的生产逐单元 wall-$y^+$ 场，当前明确拒绝，不能静默按
  `none` 运行；
- 展向分箱以面心归属整张物理面，不切割跨箱面；生产箱边界应与网格截面对齐；
- signed `sgs_energy_transfer` 保留正向传递和 backscatter，但完整正/负积分、能谱和
  $\varepsilon_{num}$ 收敛需服务器长算；
- 完整 HIT、正式槽道、近壁网格/采样收敛和生产三维机翼 LES 不在本机卡口中。

人工核验建议重点检查：固定滤波核及边/角传播、相似/动态 backscatter 符号、SGS 与数值耗散
的边界、统计 checkpoint 接缝、展向分箱恒等式、五模型 LU-SGS 门禁的证明范围，以及上述
S 级 `pending` 是否表述充分。人工批准前不合入 `release/v2.0.0`，不进入 AB，不修改主干。

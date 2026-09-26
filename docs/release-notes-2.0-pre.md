# WCNS v2.0_pre 版本说明

## 定位

`v2.0_pre` 是内部前瞻源码版本：原定 v2 通用求解、湍流模型、隐式推进、低马赫和统计功能已
进入统一程序，但验证证据限定为公式、单元测试、制造解及本机微型二维/三维算例。它不是
Case06、长期 LES 或目标三维翼型的生产物理验证版本，也不对这些问题给出精度承诺。

## 主要功能

- 二维/三维结构多块 CGNS、共形一对一连接、MPI 自动切分、检查点及合法异 rank 重启；
- Roe、Rusanov、HLLC，WCNS-JS/WCNS-Z/MDCD 等重构与两套几何 profile；
- Euler、层流 Navier--Stokes、SA-neg、modified SST-2003m 和实验性标准 k--epsilon；
- 原始 Smagorinsky、尺度相似、混合、动态 Smagorinsky 和 WALE；
- SSPRK3、定常 LU-SGS、非定常 BDF2/BDF1 双时间 LU-SGS；
- Roe+LU-SGS 的 Weiss--Smith 低马赫预处理；
- 壁面 `Cp/Cf/q_wall/u_tau/y+`，压力/黏性/总载荷、三轴力矩和展向载荷分箱；
- RANS/LES 模型字段、SGS 应力/热流/能量传递、动态系数、滤宽及网格各向异性诊断；
- 空间积分/截面统计，以及只按接受物理步 `accepted_dt` 加权并可重启的 Reynolds/Favre
  mean、RMS 和 covariance。

## 统计支持边界

逐项能力、字段、壁面量、空间统计和接受步统计清单见
[`v2.0-pre-capability-matrix.md`](v2.0-pre-capability-matrix.md)。

时间统计只接受已经注册、在当前物理/模型组合中可计算的标量统计量。伪时间迭代、被拒绝步和
重复 `(step,time)` 不采样；checkpoint 保存统计身份与累加器并支持合法异 rank 续算。未知量、
重复量、二维误用三维截面或模型不匹配均在启动/首次计算前报错，不输出占位零值。

内建频谱、自相关置信区间、自动平稳性判定和任意曲面插值截面尚未提供。用户如需频谱，应从
明确采样间隔的历史/探针序列独立后处理并报告窗函数、分段、频率分辨率和统计时长。

## 验证边界

- 本机执行串行及 2/4 rank 的完整 Release 自动矩阵和小型有限性/重启测试；
- Case06 897×257 只读取和 dry-run，执行时间步数为 0；
- 不运行大型 TMR 网格、HIT、长期槽道/hump 或目标三维翼型；
- SA/SST/k--epsilon/LES 的工程精度、网格无关性和长期统计充分性均保持“未由本版本验证”；
- SST 的 NACA0012 公开参考是 SSTm，而程序实现是 SST-2003m，不能当作严格同模型真值。

## 发布目录

`WCNS_v2.0_pre` 是精简源码目录，不含 `cases/`、`examples/`、`tests/`、网格、运行结果、构建
缓存或 Git 元数据。目录内 `PACKAGE_CONTENTS.sha256` 可核验每个文件，
`WCNS_SOURCE_REVISION` 记录来源提交。

完整的本机命令、串行/MPI 计数和证据边界见
[`v2.0-pre-validation.md`](v2.0-pre-validation.md)。

本版本继续只用于本机独立开发和内部评估；没有恢复外部 CI，也没有解决 WCNS 自有代码许可。

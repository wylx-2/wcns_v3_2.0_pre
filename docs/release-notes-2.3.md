# WCNS v2.3 版本说明

## 主要变化

- 标准 Roe、Li--Gu `roe_all_speed` 和 Rieper `roe_all_speed_rieper` 改用代数等价的
  Roe--Pike 闭式波强度与耗散恢复，不再逐面显式构造通用特征矩阵；Weiss--Smith 预处理
  Roe 保留专用矩阵实现。
- Riemann 求解器把每侧守恒状态、通量、声速、法向速度和总焓集中派生一次，供 HLL、
  HLLC、Rusanov、Roe 族及回退路径复用。
- 六点重构模板改为一次检查后同时计算左右状态，MDCD-HYBRID 复用归一化模板；公开接口
  的输入验证和运行时非有限值回退保持不变。
- 清理 Euler 基础通量、Rusanov、总焓和 Roe 特征基中的嵌套重复转换与校验。
- 新增矩阵实现对闭式 Roe--Pike 的公式级等价测试，并完成串行、MPI 和短时性能回归。
- 版本号、完整配置模板、Linux 指南、发布契约和确定性源码打包流程升级为 v2.3。

## 兼容性

v2.3 不新增或删除数值配置项，restart signature 与 v2.2 保持一致；v2.2 checkpoint 在算法
参数不变时可用 `restart.mode=strict` 续算。若在续算时切换 Riemann、重构、时间推进或其他
核心算法，仍须使用 `restart.mode=algorithm_change`。

## 性能结果与边界

同机 Release 小规模基准相对修改前改善约 6.4%–22.0%，详见
[`performance-optimization-v2.3.md`](performance-optimization-v2.3.md)。本次没有执行长期
Re_tau=180 槽道、Case06 大网格或目标三维翼型；v2.3 仍为内部、物理验证受限版本。

自动卡口与构建记录见 [`v2.3-validation.md`](v2.3-validation.md)。

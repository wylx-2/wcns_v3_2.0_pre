# WCNS v2.2 版本说明

## 主要变化

- 审查确认 v2.1 `roe_all_speed` 与指定报告第五章 Li--Gu 公式不一致；v2.2 将同名算法改为
  Li--Gu 特征值修正、压力跳跃修正和可配置 $M_{ref},c_1,c_2$。
- v2.1 的 Rieper 型实现改名为 `roe_all_speed_rieper`，用于复现旧结果。
- 新增 `hll` 和 `roe_rotated`，从而覆盖 `riemann_solver_tensor_report.pdf` 明确列出的 Roe、
  HLL、HLLC、旋转 Roe 四种求解器。
- `roe_rotated` 支持二维和三维确定性正交基；速度跳跃退化时安全回到面法向 Roe 耗散。
- Riemann 配置摘要和 restart signature 升级，Li--Gu 三参数写入数值身份。
- 更新算法、配置、用户、开发、Linux 和已知限制文档，并补充求解器逐项覆盖报告。

## 兼容性

v2.2 的 `roe_all_speed` 与 v2.1 同名算法数值含义不同。已有 v2.1 checkpoint 若要切换到
Li--Gu 算法，第一次必须使用 `restart.mode=algorithm_change`，新分支后续参数不变时再使用
`strict`。若要继续旧算法，应把配置改为：

```text
algorithm.riemann = roe_all_speed_rieper
```

## 验证边界

本版本只执行公式级单元测试、完整本机 CTest、小网格 SSPRK3、重启和无 Git 发布包重建。
没有运行长期 Re_tau=180 槽道、Case06 大网格或目标三维翼型，也没有把旋转 Roe 与激波传感器
组成混合通量。发布状态仍是内部、物理验证受限版本。

详细算法见 [`riemann-solvers-v2.2.md`](riemann-solvers-v2.2.md) 和
[`all-speed-roe.md`](all-speed-roe.md)，自动卡口结果见
[`v2.2-validation.md`](v2.2-validation.md)。

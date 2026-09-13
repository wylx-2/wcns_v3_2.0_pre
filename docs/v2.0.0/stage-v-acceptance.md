# v2.0.0 阶段 V 自动验收

结论：**通过（2026-09-13）**。本结论只表示规格、输入身份和基线已冻结，不表示 RANS、LES、
LU-SGS、低 Mach 或 Case06 物理计算已经实现。

## 1. 范围与 Git

- 基线/测试源码提交：`f72df3952f80687cd05e7dba5074dee8680ad5e2`；阶段分支：
  `stage/v2.0.0-v`。
- V 未修改生产 C/C++；交付为设计、schema 2 草案、Case06 只读 smoke 配置、三维 intake 模板
  和本报告。
- 用户目录 `cases/manual/case06_2d_naca0012/` 与
  `cases/manual/case05_3d_turbulent_channel/results/lowmach-longrun-segment01/` 均未暂存。

## 2. 环境与构建

| 项目 | 结果 |
|---|---|
| OS | Microsoft Windows NT 10.0.26200.0 |
| CPU | Intel64 Family 6 Model 151 Stepping 2，12 logical CPU |
| CMake / generator | 3.28.0 / MinGW Makefiles |
| C/C++ | MinGW-w64 GCC/G++ 8.1.0 |
| MPI | Intel MPI，`C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe` |
| 串行 Release | 独立配置、构建、安装；CTest 60/60，19.00 s |
| MPI Release | 独立配置、构建、安装；CTest 108/108，32.72 s |
| 公式规格 | `python tools/verify_algorithm_spec.py`，6/6 |

两套配置均显式启用 CGNS 和测试；MPI 构建覆盖现有 1/2/4/8-rank 回归。CTest 中的 release
matrix、输出/重启、故障注入、Case07 边界、黏性与守恒测试全部通过。由于 V 不改生产源码，
固定机器性能参考继续使用已入库的 `docs/v1.1.0/stage-p-performance-baseline.json` 与阶段 U
复核值，W 完成后另做前后对照。

## 3. Case06 只读 intake

| 项目 | 结果 |
|---|---|
| 文件 | `naca0012_897x257_str.cgns`，5,545,984 byte |
| SHA-256（前/后） | `ADBEC980FB1C1CC2956DCFDFB918ED25F7159F2287C958DB5056C926D725351E` / 相同 |
| `cgnscheck -v` | exit 0；0 error，6 warning，与 intake 登记一致 |
| 串行 dry-run 探针 | 精确退出：`self-connectivity is not supported by stage D` |
| 2-rank dry-run 探针 | 各 rank 精确退出：`2D metrics currently require two-dimensional physical space` |

上述两个退出是阶段 V 从真实输入识别并冻结的既有兼容缺口，不是“成功运行”的替代说法。
它们进入 W 的平面二维/自连接基础回归，并在 AB 以前保持为 Case06 必修项；其他错误会使 V
卡口失败。原 CGNS 未被转换、写回或纳入 Git。

## 4. 自动卡口逐项结论

1. 文档 whitespace 和《算法补充》规格检查通过；公式来源与实现名已唯一登记。
2. 串行/MPI 干净 Release 的 configure/build/install/CTest 全过。
3. schema 2 草案覆盖兼容映射、合法/非法组合、数值摘要与 restart signature 规则；生产 parser
   仍保持 schema 1，符合 V 阶段边界。
4. `model=none` 基线由未修改的 v1.1.0 源码和本轮完整输出/重启测试锁定；不存在新模型场。
5. 工作树例外仅为上述两个已登记用户目录；验收文档不保存大场文件。

依据项目负责人 2026-09-13 对 V→W 的连续授权，V 候选可非快进合入 `release/v2.0.0`，随即进入
W；中间不设置人工等待。

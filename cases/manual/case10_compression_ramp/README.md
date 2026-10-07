# 24° 压缩拐角：WCNS v2.5 多块结构网格 ILES

Ma=2.25、Re_δ₀=15800、δ₀=0.0006096 m、T∞=170 K、Tw=323 K，展宽 6δ₀；层流相似解入口，局部体力促转捩，关闭湍流模型。生产 CGNS 含 12 块、8847360 单元（相对文档 −0.10%）。[网格图](grids/mesh_xy.png)、[调研报告](report/chapter5-v2.5-research.md)、[使用和统计定义](report/chapter5-v2.5-implementation.md)、[验证记录](report/v2.5-validation.md)。

| 文件 | 用途 |
|---|---|
| `production.wcns` / `grids/ramp_coarse.cgns` | 服务器 ILES，统计窗 t=1000…4000，以 δ₀/U∞ 计 |
| `smoke.wcns` / `grids/ramp_smoke.cgns` | 18432 单元功能检查，最多 3 步，120 秒步间墙钟限制 |
| `production_restart.wcns` | 严格恢复流场和累计统计，新输出目录 |
| `reference/` | 入口表及生成参数、参考 Cp/Cfx/Van Driest 数字化曲线和原图 |
| `validation/` | 几何、短算、MPI 和改变进程分配的续算证据 |

在此目录启动（源码包与本目录相邻）：

```bash
sha256sum -c PACKAGE_CONTENTS.sha256
mpiexec -n 4 ../WCNS_v2.5/build-mpi/wcns_run --config smoke.wcns
# maximum_steps 的预期返回码为 2；重复测试须选择新的 output.directory。
```

生产网格采用默认严格度量阈值且全部通过。极小 smoke 网格专用阈值为 0.9，仅作程序功能检查；不能把 smoke 的分离、激波或摩阻作为物理解。

正式计算只在服务器进行。先核查体力条带后是否形成正确的入射湍流边界层，比较 x/δ₀=70 的 Δx+/Δy+/Δz+ 以及 x/δ₀=80 的 Van Driest 剖面；幅值 A=5 和波形是记录明确的补充选择，未经长算验证。初始场不是平衡激波/边界层，t=1000 也不自动保证启动过程消失。低频分析需增加足够长的时间窗。

```bash
python tools/analyze_chapter5.py --case ramp --prefix output/ramp-production/ramp-iles --output analysis --start 1000 --span-cells 72 --reference reference
```

输出 Cp/Cfx、热流、Reynolds/Favre 应力、边界层参数、分离/再附、固定探针、PSD 和压力相干；详细符号与单位见使用说明。新的统计 CSV 始终无量纲。现有 wall-pressure-gradient 位置是激波运动代理量，不能当作原文精确激波脚定义。参考曲线为带读图误差的数据，完整原始压力/热流时间序列未取得。

网格直接可用；可通过 `python tools/prepare_chapter5_meshes.py ramp --case-dir .` 重建二维映射，再用源码包内 `wcns_generate_chapter5_cgns grids/ramp_coarse.bin NEW.cgns` 写新 CGNS；`wcns_audit_chapter5_metrics NEW.cgns` 严格逐块检查。重新生成入口表时用 `python tools/generate_ramp_inlet.py --output-directory reference-new --table reference-new/ramp_inlet_table.inc`，不要无意修改正在运行程序所用的表。详见 [参考说明](reference/README.md)。

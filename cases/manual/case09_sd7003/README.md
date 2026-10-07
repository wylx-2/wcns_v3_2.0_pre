# SD7003：WCNS v2.5 多块结构网格 ILES

Ma=0.2、Re_c=60000、攻角 4°、展宽 0.2c、绝热无滑移壁，关闭湍流模型。生产 CGNS 含 8 块、1720320 单元（相对文档主网格 +0.53%）。[网格图](grids/mesh_xy.png)、[完整调研报告](report/chapter5-v2.5-research.md)、[使用和统计定义](report/chapter5-v2.5-implementation.md)、[验证记录](report/v2.5-validation.md)。

| 文件 | 用途 |
|---|---|
| `production.wcns` / `grids/sd7003_coarse.cgns` | 服务器显式 ILES 基线，统计窗 t=20…60 |
| `smoke.wcns` / `grids/sd7003_smoke.cgns` | 保留二维网格，展向 8 层，最多 3 步与 120 秒步间墙钟限制 |
| `production_restart.wcns` | 严格恢复流场和累计统计，新输出目录 |
| `reference/` | UIUC 几何、原文表值、带误差的参考曲线 CSV、来源和原图 |
| `validation/` | 最终短算、MPI、续算和几何证据，不是收敛物理结果 |

在此目录启动（下面假设源码包与本目录相邻；开发仓库请改可执行路径）：

```bash
sha256sum -c PACKAGE_CONTENTS.sha256
mpiexec -n 4 ../WCNS_v2.5/build-mpi/wcns_run --config smoke.wcns
# maximum_steps 的预期返回码为 2；重复测试须选择新的 output.directory。
```

**请先评估计算成本。** smoke 测得 CFL=0.1 时 Δt≈5.86×10⁻⁸ c/U∞，完整 t=60 可能需要约十亿步。生产配置可运行，但长算经济性尚未验证；服务器先做资源/步长 pilot，再决定是否继续。SD7003 的自然转捩对离散耗散和微扰敏感，本次没有计算出收敛气动力。

服务器得到统计后：

```bash
python tools/analyze_chapter5.py --case sd7003 --prefix output/sd7003-production/sd7003-iles --output analysis --start 20 --span-cells 40 --reference reference
```

默认探针只是预置位置，应根据 pilot 的 δ99 将三个剖面处探针移到半边界层厚度，另开统计段。严格续算不允许悄悄更换探针或统计窗；`algorithm_change` 模式用于显式重置统计。后处理依赖 NumPy≥2.0、SciPy、Matplotlib。新统计 CSV 始终无量纲。

网格可直接使用；重新生成二维映射可执行 `python tools/prepare_chapter5_meshes.py sd7003 --case-dir .`。然后使用源码包构建的 `wcns_generate_chapter5_cgns grids/sd7003_coarse.bin NEW.cgns`，生成器拒绝覆盖现有 CGNS。`wcns_audit_chapter5_metrics NEW.cgns` 逐块严格审计。参考读图数据的有效精度见 [参考说明](reference/README.md)。

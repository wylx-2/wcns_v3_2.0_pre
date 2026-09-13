# SA-neg 独立验证输入

本目录保存阶段 X 的小型配置、散列清单和验收脚本；NASA TMR 原始网格/参考数据体积较大且
可能更新，不提交 Git，由 `tools/fetch_sa_tmr_assets.ps1` 下载到指定缓存后用
`tools/verify_sa_tmr_assets.py` 校验。未经 SHA-256 校验的文件不得用于验收。

验收分为 SA MMS、TMR finite flat plate、TMR NACA0012 validation family 和 TMR NACA0012
numerical Family II。最粗网格只做读取/有限性 smoke。用户 Case06 不属于这里，也不能替代这些
独立网格。详细工况、阈值、网格层级与人工检查见 `docs/v2.0.0/stage-x-design.md`。

目录约定：

```text
<cache>/
  references/       # Cf/CD/CL/Cp 等 TMR 文本数据
  archives/         # 原始 zip 和生成器
  grids/            # 经明确转换并记录散列的结构网格
  results/          # solver 输出和机器可读自动报告
```

任何 PLOT3D/CGNS 转换必须保留原始坐标、C-grid wake 连接、边界索引和 double precision，并在
报告中同时记录原文件和转换文件散列。TMR 网格内携带的 BC 仅供参考；实际 boundary map 必须按
对应 `.nmf` 和 stage X 冻结工况复核。

## 规范处理流程

从仓库根目录执行，缓存目录必须位于 Git 忽略路径中：

```powershell
powershell -ExecutionPolicy Bypass -File tools\fetch_sa_tmr_assets.ps1 -CacheDirectory build-v2-x-assets
python tools\verify_sa_tmr_assets.py --manifest cases\validation\sa_tmr\manifest.json --cache build-v2-x-assets
cmake --build build-v2-x-mpi-mingw --parallel 4
```

NACA0012 validation 原始 formatted PLOT3D 转换参数固定为：

| 网格 | wall first/last vertex index |
|---|---:|
| 225x65 | 49 / 177 |
| 449x129 | 97 / 353 |
| 897x257 | 193 / 705 |

示例：

```powershell
build-v2-x-mpi-mingw\wcns_convert_tmr_p2d_to_cgns.exe input.p2dfmt output.cgns 49 177
python tools\prepare_sa_tmr_smoke.py --template cases\validation\sa_tmr\naca0012_validation_smoke.wcns.in --grid output.cgns --output-directory build-v2-x-assets\results\naca-225 --config build-v2-x-assets\naca-225.wcns
build-v2-x-mpi-mingw\wcns_run.exe --config build-v2-x-assets\naca-225.wcns --dry-run
build-v2-x-mpi-mingw\wcns_run.exe --config build-v2-x-assets\naca-225.wcns
```

平板网格由散列通过的 TMR Fortran 生成器产生；交互输入网格 level 6、5、4，
生成的 unformatted `grid.p3d` 立即重命名并转换：

```powershell
build-v2-x-mpi-mingw\wcns_convert_tmr_flatplate_p3d_to_cgns.exe grid-level6.p3d flatplate-level6.cgns
python tools\prepare_sa_tmr_smoke.py --template cases\validation\sa_tmr\finite_flat_plate_smoke.wcns.in --grid flatplate-level6.cgns --output-directory build-v2-x-assets\results\flatplate-level6 --config build-v2-x-assets\flatplate-level6.wcns
build-v2-x-mpi-mingw\wcns_run.exe --config build-v2-x-assets\flatplate-level6.wcns --dry-run
build-v2-x-mpi-mingw\wcns_run.exe --config build-v2-x-assets\flatplate-level6.wcns
```

模板中的 `run.max_steps=1` 只是读取、拓扑、边界和首步有限性 smoke，绝对不是物理验收。
定量验收必须复制配置、使用新输出目录、运行到残差与载荷联合收敛，然后按
`manifest.json` 阈值进行多网格比对。当前全局 SSPRK3 步长为 $10^{-9}$ 量级，在项目负责人
批准阶段依赖调整前不开始伪收敛长算。

## 当前自动结果

- 7/7 原始资产散列通过；
- flat-plate level 6/5/4 均通过转换、dry-run 和首步有限性；
- NACA0012 225/449/897 均通过转换、dry-run 和首步有限性；
- NACA0012 225 的 1/2/4-rank 场最大绝对差为 $4.27\times10^{-14}$；
- level 8 平板仅有 7 个法向顶点，被高阶 Jacobian 严格检查拒绝；它仅是官方
  网格族中的最粗 smoke 层，不在 level 6/5/4 定量验收层级内。

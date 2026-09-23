# SA-neg 独立验证输入

本目录保存阶段 X 的小型配置、散列清单和验收脚本；NASA TMR 原始网格/参考数据体积较大且
可能更新，不提交 Git，由 `tools/fetch_sa_tmr_assets.ps1` 下载到指定缓存后用
`tools/verify_sa_tmr_assets.py` 校验。未经 SHA-256 校验的文件不得用于验收。

2026-09-14 起本机只运行不超过 32,768 cells 的 L 级 smoke/小算例，单例最多 2,000 accepted
steps/1,800 s；449/897、三级充分收敛和严格表面曲线比较是服务器 S 级。897 仅有一次累计不超过
2 小时的 LU-SGS 稳定参数诊断例外，不得转为本机物理长算。完整限制见
[`docs/v2.0.0/local-validation-policy.md`](../../../docs/v2.0.0/local-validation-policy.md)。

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
服务器 S 级定量验收必须复制配置、使用新输出目录、运行到残差与载荷联合收敛，然后按
`manifest.json` 阈值进行多网格比对。当前全局 SSPRK3 步长为 $10^{-9}$ 量级，在项目负责人
批准的顺序调整下，必须等前置 AA 的层流/SA LU-SGS 自动卡口通过并回合到 X 后，
才开始 X-B 验收。前置计划已获人工批准；本机只执行 L 级，任何 smoke 仍不得当作物理验收。

AA 已通过后，定量配置由同一模板以 `acceptance` 模式生成；该模式固定 LU-SGS、严格残差
判停、每 1000 步 checkpoint 和每 20 步边界/载荷输出，避免手工复制时遗漏恢复能力：

```powershell
python tools\prepare_sa_tmr_smoke.py --template cases\validation\sa_tmr\naca0012_validation_smoke.wcns.in --grid build-v2-x-assets\grids\n0012_225-65_2d.cgns --output-directory build-v2-x-assets\results\naca225-a10-local --config build-v2-x-assets\naca225-a10-local.wcns --mode local --case-name tmr-naca225-a10-local --angle-deg 10 --cfl 1 --max-steps 500 --max-wall-time 1800
python tools\check_validation_budget.py --config build-v2-x-assets\naca225-a10-local.wcns --kind probe --dimension 2 --cells 14336
```

服务器 `acceptance` 模式必须另行显式给出已批准的 `--max-wall-time`，且生成后的配置不通过本机
预算检查器；不得在本机启动。

Family II 的压缩包、no-PV 与带 PV 参考数据也由 manifest 固定散列。只解出
`n0012familyII.6/.5/.4.p2dfmt.gz`，按 NMF 的 `[49,177]`、`[97,353]`、`[193,705]`
壁面顶点范围转换为 225/449/897 三层 CGNS 后，执行：

TMR 将该二维网格解释为 $x$--$z$ 平面并报告 `CMy`；WCNS 将 PLOT2D 第二坐标映射到
$+y$ 并输出右手系 `Cm_z`，因此比较时固定使用 `CMy=-Cm_z`，报告同时保留原始 `Cm_z`。

```powershell
python tools\verify_sa_tmr_assets.py --manifest cases\validation\sa_tmr\manifest.json --cache build-v2-x-assets --include-derived
```

no-PV 是主分支。带 PV 对照在生成命令末尾增加 `--point-vortex-cl 1.09125`；公式、正号约定
和拒绝条件见《算法补充》12.11。

## 当前自动结果

- 13/13 原始资产散列通过；Family II 三层派生 CGNS 也已登记并通过散列；
- flat-plate level 6/5/4 均通过转换、dry-run 和首步有限性；
- NACA0012 225/449/897 均通过转换、dry-run 和首步有限性；
- NACA0012 225 的 1/2/4-rank 场最大绝对差为 $4.27\times10^{-14}$；
- level 8 平板仅有 7 个法向顶点，被高阶 Jacobian 严格检查拒绝；它仅是官方
  网格族中的最粗 smoke 层，不在 level 6/5/4 定量验收层级内。

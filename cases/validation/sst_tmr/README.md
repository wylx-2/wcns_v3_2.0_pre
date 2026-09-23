# SST-2003m 独立验收输入

本目录定义阶段 Y 的真实 RANS 验收，不把公式/MMS、一步 smoke 或用户 Case06 当作物理验证。

自 2026-09-14 起分为本机 L 级和服务器 S 级：本机主例仅为 35x25 平板，补充例从 41x129
槽道/103x28 hump 二选一，均受 2,000 accepted steps、1,800 s 单例上限约束；其余网格和充分
定常收敛只在服务器执行。L 级判断软件通路和定性行为，不替代 S 级物理精度。详见
[`docs/v2.0.0/local-validation-policy.md`](../../../docs/v2.0.0/local-validation-policy.md)。
模型身份固定为本项目《算法补充》12.1、12.2、12.5 所述 modified SST-2003m；TMR 页面中
名称为 SSTm 的平板、槽道和 hump 数据只作为具有独立来源的工程物理基线，报告必须明确这一
模型标识差异。

三组 S 级卡口分别是：零压梯度湍流平板的三级网格趋势、高 Reynolds 数发展槽道的壁面与对数层
统计、NASA no-plenum wall-mounted hump 的逆压梯度分离/再附。所有远程输入和派生 CGNS
散列均冻结在 `manifest.json`；缓存位于 Git 忽略的 `build-v2-y-assets`，不得提交大网格和原始
输出。

## 资产与转换

```powershell
powershell -ExecutionPolicy Bypass -File tools\fetch_sst_tmr_assets.ps1 -CacheDirectory build-v2-y-assets
python tools\verify_sa_tmr_assets.py --manifest cases\validation\sst_tmr\manifest.json --cache build-v2-y-assets
python tools\verify_sa_tmr_assets.py --manifest cases\validation\sst_tmr\manifest.json --cache build-v2-y-assets --include-derived
```

`wcns_convert_tmr_p2d_to_cgns` 的整数参数模式用于平板；`channel` 模式把两个 J 面写为解析无滑移
壁；`hump` 模式把下壁写为无滑移壁、上壁写为滑移面。转换前必须按官方 NMF 复核索引，转换后
必须先 dry-run 并比对 manifest 中的 SHA-256。

剖面比较的列契约不可互换：平板 `analyze_sst_flatplate.py --uplus-reference` 必须传
`references/flatplate/sst-upyp_cfl3d.dat`，该文件两列为 `log10(y+)`、`u+`；
`u+y+.dat` 的两列是未取对数的 `y+`、`u+`，只作为独立一般壁律参考。槽道
`sst_u+y+KM_cfl3d.dat` 的三列依次为 `u+`、`log10(y+)`、`KM`，分析器按第 2 列为横轴、
第 1 列为纵轴读取。任何列交换都必须使验收失败，不能通过调宽误差阈值补偿。

这些 SST 验收均使用 `resolved` 壁面处理。场文件的单元量 `wall_y_plus` 只服务壁面函数路径，
在 resolved 路径中不作为验收数据；边界面文件则必须输出 `mu_w,wall_distance,`
`friction_velocity,wall_y_plus,wall_y_plus_class`。其中边界 `wall_y_plus` 由权威壁面黏性牵引、
第一层真实几何距离、最近单元密度和冻结 Reynolds 数计算。平板和槽道分析器还必须从最终
$C_f$ 和剖面独立重算 $y^+$，不能把场文件的零诊断列误报为首层 $y^+=0$。

`wall_y_plus_class` 是便于筛选的数值枚举：0 表示 $y^+\le5$，1 表示 $5<y^+<30$，2 表示
$30\le y^+\le300$，3 表示 $y^+>300$；它不参与推进，也不得触发壁面处理自动切换。

## 本机 L 级生成与预算检查

```powershell
python tools\prepare_sst_tmr_case.py --template cases\validation\sst_tmr\flat_plate_smoke.wcns.in --grid build-v2-y-assets\tmr-flatplate-35x25.cgns --output-directory build-v2-y-assets\local-flatplate-output --config build-v2-y-assets\local-flatplate.wcns --mode local --case-name y-local-sst-flatplate-35x25 --cfl 0.05 --relaxation 0.5 --max-steps 500 --max-wall-time 1800
python tools\check_validation_budget.py --config build-v2-y-assets\local-flatplate.wcns --kind probe --dimension 2 --cells 816
```

实际派生网格名以散列清单为准。服务器 `acceptance` 模式必须显式给出已批准的
`--max-wall-time`；不得在本机运行该模式。

## 服务器 S 级自动判定

每个 S 级定量算例必须：残差至少下降 3 个数量级、末 10 个边界输出的载荷/摩阻相对跨度不超过
0.2%、末 500 个接受步无 `k/omega` 下界投影、全程投影率不超过 0.5%，且 `omega` 投影恒为
零。程序停止原因为 `maximum_steps` 只表示数值推进结束，不等于验收通过。

- 平板：35x25、69x49、137x97 使用相同工况和停止标准；细网格 `x=0.97` 的 `Cf` 相对 TMR
  理论/公开结果不超过 5%，相邻两级变化不超过 3%，并报告 `u+`、`y+`、`mu_t/mu`。
- 槽道：以官方网格在 `x=500` 比较 `Cf`、摩阻速度、`u+` 和 `mu_t/mu`，`Cf` 误差不超过
  5%，对数层拟合 `kappa` 必须位于 `[0.38,0.44]`。
- hump：103x28、205x55、409x109 同一停止标准；最细层分离/再附分别位于 manifest 冻结区间，
  `Cp/Cf` 相对 CFL3D 曲线的 L2 误差分别不超过 12%/25%，三级结果需呈可解释趋势。

真实运行的命令、版本、rank、输入散列、墙钟、停止原因、残差、投影计数和所有曲线误差统一
写入阶段 Y 自动报告。Y 的 L 级全部自动卡口通过后可创建明确标注 S 级状态的候选标签并停在
人工验收点；S 级未通过时不得升级物理支持声明。

解析壁面逐面恒等式使用：

```text
python tools/verify_resolved_wall_units.py --boundary <latest-boundary> --reynolds <Re>
```

脚本从输出的切向黏性牵引与 $y^+$ 分别反推壁面密度，二者相对差必须 `<=2e-12`，同时严格
核验单位法向和四档分类。

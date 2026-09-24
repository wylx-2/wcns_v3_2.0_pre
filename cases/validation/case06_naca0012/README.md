# Case06 NACA0012 冻结输入与分层验收

本目录保存阶段 AB 的可跟踪配置模板、NASA TMR 数据清单和运行生成器。用户原始
`cases/manual/case06_2d_naca0012/naca0012_897x257_str.cgns` 仍是未跟踪输入，不复制、不改写、
不纳入源码包。公开参考数据下载到 Git 忽略的缓存，并在使用前逐文件校验 SHA-256。

## 1. 冻结工况

主工况取 NASA TMR 数值分析定义：二维、全湍流、$M_\infty=0.15$、$Re_c=6\times10^6$、
$T_\infty=540\,R=300\,K$、$Pr=0.72$、$Pr_t=0.9$、绝热无滑移壁，弦长/单位翼展参考面积均为
1，力矩中心为 $(0.25,0,0)$。首要验收点是 $\alpha=10^\circ$；$0^\circ$ 用于对称性与零升力
诊断，$15^\circ$ 只在首要工况通过后作为高升力扩展，不作为本机长算。

`sa_neg` 使用 $\widetilde\nu_\infty/\nu_\infty=3$。`sst_2003m` 使用 TMR NACA 页面给出的
$I=0.052\%$、$\mu_t/\mu=0.009$，换算为
$k^*=4.056\times10^{-7}$、$\omega^*=270.4$、$L_t/c=4.300130725961134\times10^{-6}$。
公式和反算恒等式见根目录《算法补充》12.17，并由 `wcns.case06_preparation` 自动测试。

TMR SA 曲线可作为 SA-neg 非负稳态分支的直接工程对照；负分支占比仍必须单独报告。项目实现
是 SST-2003m，而公开 NACA0012 页面曲线是 SSTm，二者不是同一模型。SSTm 数据只用于发现
量级、方向或曲线形状错误，报告不得把容差通过写成 SST-2003m 的严格同模型验证。

## 2. 参考数据层次

- `n0012*sa.dat` / `n0012*sst.dat`：897×257 上的 CFL3D 参考 CFD，用于代码到代码的积分量和
  `Cp/Cf` 比较，不是真值；这些文件的升力工况启用了固定点涡远场修正。
- `CLCD_Ladson_expdata.dat`：$Re=6\times10^6$ 的 tripped 力数据，是全湍流力的主要实验对照。
- `CP_Gregory_expdata.dat`：上表面压力分辨率较好，但 $Re\approx2.88\times10^6$ 且来自数字化
  图线，只作压力形状对照，不能与 Case06 的 Reynolds 数混同。
- TMR 没有本算例的实验 `Cf`；`Cf` 只能对参考 CFD、网格趋势和壁面单位恒等式评估。

无点涡分支是自洽主解。只有服务器上的同源 TMR 对照才使用 `--reference-point-vortex`，其
$C_{L,pv}$ 由 manifest 按模型和攻角读取；程序不会从迭代载荷反馈更新点涡。两个分支必须分开
命名、分开报告，禁止把固定参考升力当成求解预测。

## 3. 本机 L 级

本机只使用公开 113×33（3,584 cells）网格，执行 dry-run 或最多 20 个接受步、300 s 的短探针。
它只证明配置、拓扑、有限性、载荷方向、模型量合法域和短趋势，不证明物理收敛。生成例：

```powershell
python tools\prepare_case06_naca0012.py --model sa_neg --grid build-v2-x-assets\n0012_113-33_2d.cgns --output-directory build-v2-ab-assets\sa-local-output --config build-v2-ab-assets\sa-local.wcns --case-name case06-sa-local --mode local --angle-deg 10 --max-steps 20 --max-wall-time 300
python tools\check_validation_budget.py --config build-v2-ab-assets\sa-local.wcns --kind smoke --dimension 2 --cells 3584
build-v2-z-candidate-serial\wcns_run.exe --config build-v2-ab-assets\sa-local.wcns --dry-run
```

本机不得对 897×257 Case06 执行时间步；该文件只允许哈希、网格检查和 dry-run。

## 4. 完整功能后的服务器 S 级

服务器任务推迟到 AC/AD 通用三维功能完成之后、AE 发布候选之前。先运行
$\alpha=10^\circ$：SA-neg 使用 225/449/897 三级趋势，SST-2003m 先用 225 排错再运行 Case06
897 主解。两模型都先无点涡，只有主解通过后才运行固定参考 $C_L$ 的点涡对照。$0^\circ$ 和
$15^\circ$ 不得抢占首要工况预算。

每个服务器配置必须用 `--mode server` 显式给出 `--max-steps` 和 `--max-wall-time`，并在提交前
记录 commit、二进制、网格/配置散列、rank、内存、磁盘、检查点和提前失败条件。例：

```powershell
python tools\prepare_case06_naca0012.py --model sa_neg --grid <verified-897.cgns> --output-directory <server-output> --config <server-config.wcns> --case-name case06-sa-a10 --mode server --angle-deg 10 --max-steps 50000 --max-wall-time <approved-seconds>
```

必须同时满足残差、载荷窗口、质量守恒、`y+`、投影/降阶、MPI/restart 和参考曲线卡口；
`maximum_steps` 或 `maximum_wall_time` 结束不等于收敛。manifest 内的数值容差是在运行前冻结的
故障探测阈值，实验差异仍需人工解释。

## 5. 资产获取

```powershell
powershell -ExecutionPolicy Bypass -File tools\fetch_case06_tmr_assets.ps1 -CacheDirectory build-v2-ab-assets
python tools\verify_sa_tmr_assets.py --manifest cases\validation\case06_naca0012\manifest.json --cache build-v2-ab-assets
```

缓存、网格和求解输出均不得提交 Git；只提交 manifest、模板、生成器和小型验收摘要。

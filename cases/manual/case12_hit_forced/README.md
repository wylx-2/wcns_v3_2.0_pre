# WCNS v2.6 JHTDB 强迫均匀各向同性湍流

三维、常黏度、无显式湍流模型的周期盒。默认 `production.wcns` 为64³八块网格，`coarse32.wcns` 与 `refine128.wcns` 为分辨率对照，`smoke.wcns` 为16³最多3步。

计算域[0,2π]³，ν=0.000185，Re=5405.405405，Ma_ref=0.1。初始场匹配截断的公开 JHTDB 能谱。`jhtdb_shells` 在每个接受步后将两个整数波数壳层能量设为0.30与0.13，分别对应0.5≤|m|≤1.5、1.5<|m|≤2.5。均匀内能增长率由 thermostat 扣除。

## 运行

```sh
mpiexec -n 4 /path/to/wcns_run --config smoke.wcns
mpiexec -n 4 /path/to/wcns_run --config production.wcns --dry-run
mpiexec -n 4 /path/to/wcns_run --config production.wcns
python tools/analyze_hit.py --case forced --reference reference --output output/hit-forced-production --result comparison
```

正式配置计划终止时间80、平均窗口20至80，必须在服务器上检查平稳性再决定是否延长。`constant_power.wcns` 的功率为0.103，用于强迫方法敏感性研究，不是与主配置完全相同的DNS强迫。返回码2可表示正常达到 max_steps。输出目录不得与旧目录重用。

`production_restart.wcns` 含需要替换的检查点占位符；`tools/prepare_hit_restart.py` 可生成严格续算配置。跨MPI进程数续算已经短测。统计平均随检查点保存；分块误差估计需要通过 --history 补入旧段历史文件。

## 参考与统计

`reference/spectrum.txt` 与 `ener_Re_time.txt` 是公开原始数值文件。后者0至10.04的去重梯形平均动能为0.705309882769；它与资料中不同时间窗口的K≈0.695或0.705不可混淆。两版README、文件来源和SHA256均保留。

程序输出能谱、强迫能量与功率、动能、应力、梯度耗散、尺度、偏斜度、峰度、一维谱和二阶空间相关。比较应以统计分布和平均为主；不同相位不能逐点重现DNS时间序列。粗网格未解析耗散尺度，不能把已解析 ε_viscous 或 Reλ 直接当DNS全场量。

完整背景、强迫推导、统计定义、并行结构及使用限制见 [详细报告](report/HIT算例与v2.6实现报告.md)。`validation` 是解析与短步一致性证据，不是长时间湍流验证结果。

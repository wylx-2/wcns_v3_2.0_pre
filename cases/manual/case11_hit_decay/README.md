# WCNS v2.6 CBC 自由衰减均匀各向同性湍流

2026-10-07 更新：已新增“预演化 → Kang 逐模态重匹配 → 正式时间归零”的路线。使用 [prepared_production.wcns](prepared_production.wcns)（64³）、`prepared_coarse32.wcns`、`prepared_refine128.wcns`；本机仅测试 [prepared_smoke.wcns](prepared_smoke.wcns)。原配置保留为直接随机初场。详见 [初始化更新说明](../../../docs/hit-v2.6/HIT初始化与预演化说明.md)。

三维、常黏度、无显式湍流模型，使用 Comte-Bellot 与 Corrsin 1971 的42、98、171测站参考谱。默认 `production.wcns` 为64³八块网格，`coarse32.wcns` 与 `refine128.wcns` 用于网格敏感性；`smoke.wcns` 为16³最多3步。

物理盒长0.5588 m，速度尺度0.2718933614489327 m/s，ν=1.5e-5 m²/s，计算域[0,1]³，Ma_ref=0.1。精确采样时间为0、0.13841843855582028、0.31885676024465742。完整实验参考 K 为1、0.332440548656765、0.162532464897330；粗网格初始谱截断后不重新放大至K=1。

## 运行

```sh
mpiexec -n 4 /path/to/wcns_run --config smoke.wcns
mpiexec -n 4 /path/to/wcns_run --config production.wcns --dry-run
mpiexec -n 4 /path/to/wcns_run --config production.wcns
python tools/analyze_hit.py --case decay --reference reference --output output/hit-decay-production --result comparison
```

正式运行仅建议在计算服务器执行。返回码2可表示正常达到冒烟测试的 max_steps，需查看日志。每次使用新输出目录。`production_restart.wcns` 含必须替换的检查点占位符；推荐使用 `tools/prepare_hit_restart.py`，具体命令见报告。

## 初始化与统计

旧配置按原始谱作单位换算、对数插值和球形N/3截断，以固定种子生成共轭对称无散初始场。初始密度和温度为1，平均速度0，随风洞流动的坐标系中无需加入10 m/s平均风。旧配置未经预演化，早期能量传递与声学调整可能影响实验匹配。新 `prepared_*` 配置在此基础上新增准备阶段；正式配置的准备时长为 0.10381382891686522，准备数据独立写入 `preparation/`，不计入正式统计。

输出动能、Reynolds与Favre应力、谱、梯度耗散、涡量、尺度、梯度高阶矩、沿x的一维谱和二阶相关。所有量的定义与有效比较范围见 [详细报告](report/HIT算例与v2.6实现报告.md)。CBC应比较三个指定时刻，不能用全程平均代替测站数据。

`reference` 内含原始 MAT、表3的CSV、无量纲输入谱、表4 JSON 与下载来源清单。`tools/prepare_hit_reference.py` 可从原始CSV重新生成谱；续算期间不要改写输入文件，因为 strict 检查包括文件内容散列。`validation` 保留本次短测记录，未完成三个测站的长计算。

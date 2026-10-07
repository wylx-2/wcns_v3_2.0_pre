# WCNS v2.6：CBC 实验谱九算例服务器计算包

生成日期：2026-10-07。本包包含当前工作区源码快照（含 MPI 傅里叶初始化、预演化及谱重匹配更新），不是早先发布的旧 v2.6 ZIP。程序版本保持 2.6；准确源码身份见 `source/WCNS_SOURCE_REVISION`、`campaign.json` 和 `PACKAGE_CONTENTS.sha256`。本地仅做短时验证，正式九组结果需在服务器计算。

## 1. 计算矩阵和固定设置

运行顺序如下，同一时刻只运行一个 MPI 算例。

| 顺序 | 算例 ID | 总网格单元数 | Riemann 求解器 |
|---|---|---:|---|
| 1 | cbc_n032_roe | 32³ = 32,768 | Roe |
| 2 | cbc_n032_rusanov | 32³ | Rusanov |
| 3 | cbc_n032_roe_all_speed | 32³ | Roe all-speed |
| 4 | cbc_n064_roe | 64³ = 262,144 | Roe |
| 5 | cbc_n064_rusanov | 64³ | Rusanov |
| 6 | cbc_n064_roe_all_speed | 64³ | Roe all-speed |
| 7 | cbc_n128_roe | 128³ = 2,097,152 | Roe |
| 8 | cbc_n128_rusanov | 128³ | Rusanov |
| 9 | cbc_n128_roe_all_speed | 128³ | Roe all-speed |

所有正式网格为 `[0,1]³` 周期均匀立方体，沿三个方向各分为两块，共 **8 个结构块**。网格已经随包提供，不需要服务器另行生成。`inputs/grids/hit16.cgns` 仅用于短测。

| 配置项 | 数值 |
|---|---|
| `algorithm.profile` | `scmm6_wcns` |
| `algorithm.reconstruction` | **`mdcd_linear`** |
| `algorithm.reconstruction_variables` | `characteristic` |
| `algorithm.mdcd.disp` / `.diss` | 0.0463783 / **0.001** |
| `time.integrator` / `run.cfl` | `ssprk3` / **0.3** |
| `turbulence.model` / `preconditioner.type` | `none` / `none` |
| `run.viscous` / `transport.model` | `true` / `constant` |
| γ / Pr | 1.4 / 0.72 |
| `robustness.enabled` | `false` |
| 初始条件 | CBC 谱、种子 20261003、保留整数波矢满足 `0 < |m| ≤ N/3` 的球形区域 |
| 时间上限 | 预演化 0.10381382891686522；正式衰减 0.31885676024465742 |
| 步数保护上限 | 每阶段 10,000,000 步 |

Roe all-speed 使用当前程序的 Li–Gu 实现，`reference_mach=0.1`、`dissipation_scale=0.02`、`pressure_coefficient=0.05`；普通 Roe、Rusanov 配置不包含这些专属参数。三种通量的稳定时间步不同，以同一物理时刻比较，不能以相同步数比较。脚本不自动降低 CFL 或更换格式；程序已有的通量可容许性回退仍可能发生，须保留日志中的相关计数用于分析。

`hit.type=decay`，没有外加强迫、温控或显式 SGS 模型。保留分子黏性，数值耗散承担 ILES 的未解析尺度耗散。

## 2. CBC 准备流程和对比时刻

本包使用 `hit.initialization=shell_spectrum`，按如下流程计算每一个组合：

1. 按 CBC 第 42 测站谱生成随机无散速度场。
2. 用该组合指定的 WCNS/通量方案预演化至 `t_pre*=0.10381382891686522`。
3. 对预演化速度作横向投影并逐傅里叶模态重匹配 CBC 初始谱，保留横向相位；密度、内能保留。
4. 正式时间、步数及 HIT 统计清零；新 `t*=0` 对应 CBC 第 42 测站。
5. 继续衰减至第 98、171 测站，保存瞬时统计与全流场。

同一网格、同一随机种子的三个组合具有相同准备初场和目标谱；由于各自使用不同通量预演化，**重匹配后的瞬时相位不相同**。本轮比较的是完整准备协议下的计算结果。不能把某一通量的准备检查点直接作为另一通量的 strict 重启文件，也不能把三个结果的差异全部解释成同一瞬时初场下的通量差异。跨网格比较也有截断尺度差异；目前仅一个随机种子，不是多样本统计集合。

本流程参考 [Bae 的 CBC 初场说明](https://web.stanford.edu/~hjbae/CBC) 的预演化与 Kang 重匹配思路；原实现采用 QR 模型，本包按要求采用无显式湍流模型的 ILES。这里是算法适配，不声称逐点复现原 MATLAB 初场。

取实验网格宽度 `M=0.0508 m`、来流 `U0=10 m/s`，参考量 `Lref=0.5588 m`、`Uref=0.27189336144893272 m/s`、`ν=1.5e-5 m²/s`，得到 `Re=10128.9340`。参考密度 1 kg/m³，理想气体常数取 1，参考温度 5.2804285714285708，使参考马赫数为 0.1。温度是所选可压缩数值模型的缩放量，不能理解为风洞实际温度。

\[
t^*=(t'-42)\frac{M U_{ref}}{U_0 L_{ref}},\qquad t'=tU_0/M.
\]

| CBC 测站 t′ | 正式 t* | 自第 42 测站起的物理时间 / s | 实验总动能 K/Uref² |
|---:|---:|---:|---:|
| 42 | 0 | 0 | 1.0000000000 |
| 98 | 0.13841843855582028 | 0.28448 | 0.3324405487 |
| 171 | 0.31885676024465742 | 0.65532 | 0.1625324649 |

输入实验谱、三个测站数据和数据出处均在 `inputs/reference/`。`CBC_exp.mat` 为公开转录数据，原始来源为 Comte-Bellot & Corrsin (1971), *Simple Eulerian time correlation of full- and narrow-band velocity signals in grid-generated, ‘isotropic’ turbulence*, JFM 48, 273–337，表 3、4。下载地址与 SHA256 见 `download_manifest.json`。

初始化不把截断之外的能量强行补回。粗网格解析动能通常小于实验全谱动能；重匹配后的离散壳谱还受到壳内模态计数影响。后续比较必须同时给出全谱实验值与同带宽参考值。`preparation/*_rematch.csv` 含重匹配前/后谱及离散目标，供检查初始化误差。

## 3. 本地上传步骤（Windows PowerShell）

上传文件为 `WCNS_v2.6_CBC_9cases_20261007.tar.gz` 及同名 `.sha256` 文件。以下 `your_user@your_server`、`/your/work/path` 均需替换为自己的服务器账号和工作目录。

```powershell
cd D:\program\WCNS_v3
Get-FileHash .\WCNS_v2.6_CBC_9cases_20261007.tar.gz -Algorithm SHA256
Get-Content .\WCNS_v2.6_CBC_9cases_20261007.tar.gz.sha256
ssh your_user@your_server "mkdir -p /your/work/path"
scp .\WCNS_v2.6_CBC_9cases_20261007.tar.gz .\WCNS_v2.6_CBC_9cases_20261007.tar.gz.sha256 your_user@your_server:/your/work/path/
ssh your_user@your_server
```

也可用 WinSCP/MobaXterm 按二进制模式上传这两个文件。不要用编辑器转换包内源码、配置或脚本的换行符。

## 4. 服务器解压和编译（Linux Bash）

```bash
cd /your/work/path
sha256sum -c WCNS_v2.6_CBC_9cases_20261007.tar.gz.sha256
tar -xzf WCNS_v2.6_CBC_9cases_20261007.tar.gz
cd WCNS_v2.6_CBC_9cases_20261007
python3 scripts/verify_package.py
```

需要 CMake ≥3.20、支持 C++20 的编译器（建议 GCC 11 或更新）、MPI 开发环境、Python ≥3.9、Bash 和 make/Ninja。FFT 的本地变换采用随包 FFTW，分布式转置采用程序的 MPI slab 实现；不需要额外安装 FFTW-MPI。CGNS 和 FFTW 源码压缩包都已提供，CMake 不需要从网络下载它们。

在集群上先按站点说明加载兼容的编译器/MPI 模块；以下只检查已经加载的环境：

```bash
cmake --version
c++ --version
mpicxx --version
mpiexec --version
python3 --version
JOBS=2 bash build_linux.sh
```

若需要指定编译器：`JOBS=4 bash build_linux.sh -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++`。不要关闭 MPI/FFTW/CGNS。编译器或 MPI 实现改变时，使用新的编译目录并重新验证，不能混用其他 MPI 实现的 launcher。主程序是 `build/wcns_run`；编译日志保存在 `build-logs/`。

本包保留所有当前求解器源码及 CMake 所需工具源码，但未打包开发工作区的完整测试/其他算例数据。`WCNS_BUILD_TESTS=OFF` 是有意设置。

## 5. 先做服务器短测

下面命令运行九种配置标签的 **16³ 缩短流程**，每个组合的预演化与正式阶段均只到 0.001，最多各 30 步，用于验证编译、MPI、三种通量、重匹配和输出流程。结果进入独立的 `results_smoke/`，不会被当成正式结果。

```bash
bash run_all.sh --smoke --np 8
cat results_smoke/summary.csv
```

正式网格的初始化/分区检查可单独执行：

```bash
bash run_all.sh --dry-run --case cbc_n032_roe --np 8
```

`--dry-run` 不推进时间，但仍需要加载网格、构建几何指标并初始化，128³ 应在有足够内存的计算节点执行。要检查所有真实分辨率，可去掉 `--case`。短测成功仅说明流程工作，不能证明正式长时间稳定性或实验精度。

## 6. 顺序运行全部九算例

**在已分配的计算节点/资源内运行**。`--np` 为每个算例的 MPI 进程数，不是九算例总进程数。默认每组 8 个 MPI 进程，线程数默认为 1；应按节点内存与资源分配调整。网格块可自动切分，不要求进程数正好等于块数；过多进程导致每方向子块不足 6 单元时会拒绝分区。

```bash
bash run_all.sh --np 8
```

可对不同网格分配不同进程数（下列数字仅为调用示例，必须与已申请资源匹配）：

```bash
bash run_all.sh --np32 8 --np64 16 --np128 32
```

脚本会按表中顺序执行 `mpiexec -n P build/wcns_run --config ...`。默认不设额外墙钟上限，直到九组完成或暂停。每组保留自己的输出和日志。结果位置：

```text
results/
  campaign_state.json           # 续算状态；请保留
  summary.csv                  # 九组状态总表
  cbc_n032_roe/
    segment-0001/
      config.wcns              # 本段实际执行配置，含绝对路径
      run.log                  # 完整 MPI/求解器输出与命令
      output/
        preparation/           # 准备阶段历史、检查点、谱重匹配记录
        cbc_n032_roe_hit_*.csv  # 正式阶段统计
        ...manifest.r8.txt
        ...checkpoint.*.cgns
        ...field.*.cgns
    segment-0002/              # 若续算，使用新的输出目录
```

独占服务器无调度器、需要断开 SSH 时，可以在 `tmux` 中运行；或：

```bash
nohup bash run_all.sh --np 8 > campaign-console.log 2>&1 < /dev/null &
echo $! > campaign.pid
tail -f campaign-console.log
```

不要在集群登录节点用 nohup 启动长计算。已有调度系统时使用作业队列。

## 7. Slurm 与墙钟额度

包内 `job.slurm` 是单节点 8 进程、24 小时的模板。根据集群要求编辑 partition、account、节点/进程数、模块加载和时间，再提交：

```bash
sbatch job.slurm
squeue -u "$USER"
```

模板给 runner 设置 `--job-seconds 85800`，预留 300 秒供最后检查点和退出，作业请求为 86400 秒。更改 Slurm 时间必须同步修改该预算，并按文件系统速度增大检查点余量。预算覆盖整轮本次调用，不是每组各自获得一个完整时限。下一次 `sbatch job.slurm` 自动跳过已完成组并继续未完成组。

默认在 Slurm 分配内使用 `mpiexec`。若站点规定使用 srun 且已配置匹配的 MPI 插件，可将执行行改成 `bash run_all.sh --launcher srun --np "$SLURM_NTASKS" --job-seconds 85800`；站点额外参数可写成 `--launcher-arg=--mpi=pmix` 等，按站点实际要求选择。PBS/其他队列也可在已分配资源内调用同一个 `run_all.sh`，设置小于作业时限的 `--job-seconds`。

## 8. 中断、续算与失败处理

- 正式终点达到且三个实验时刻的历史记录齐全才标记 `complete`。仅 `--dry-run` 成功标记 `validated`，不会混入正式结果。
- 退出码 2（墙钟、用户信号、最大步数）标记 `pending`，停止调度后续组合。再次运行**同一命令/同一结果目录**自动续算。准备检查点继续准备，正式检查点继续正式阶段，不重复已完成的重匹配。
- 数值失败或启动错误标记 `failed`，保留日志并继续尝试其他组合；整批最终返回非零。检查日志后可用 `--retry-failed` 从随机初始化重新尝试失败组合。它不会自动修改 CFL、通量或其他物理参数。
- 若达到最大步数，上限是阶段绝对步数，重复原命令不会产生有效推进。需检查原因并有意识地调整配置/管理重启；脚本不会无止境自动重试。当前保护值为 1000 万步。
- 每次调用获得独立的 `segment-NNNN`，已有数据不覆盖；一个结果目录同一时间只允许一个 runner。完整保留所有段，以免丢失早期实验测站。
- 输入、可执行文件 SHA256 和包目录身份均记录。计算开始后不要移动目录、改源码重编译或改九份配置继续同一结果树；脚本会拒绝身份不符。MPI 进程数可变，配置允许的 strict 重分区重启可用。
- 墙钟预算是首选的正常暂停方式。Ctrl+C/SIGTERM 会尝试传递信号，但不同 MPI launcher 对信号的处理不同；强制结束/调度器 SIGKILL 可能只能退回先前周期检查点。遇到未完整写出的检查点，先保留现场并检查日志，不要删除其他检查点。

常用命令：

```bash
cat results/summary.csv
tail -n 50 results/cbc_n032_roe/segment-0001/run.log
bash run_all.sh --np 8 --job-seconds 3600
# 上一轮暂停后，再运行同一条命令即可继续。
```

## 9. 输出统计与后续比较

统计采用空间体积平均及 HIT 专用时间积分；预演化和正式统计分开。正式阶段每 10 步采样，并在 42/98/171 实验测站强制输出；每 100 次采样刷新统计文件，阶段结束也写出。全场在三个实验时刻输出。检查点每 0.05 个代码时间输出，并在初始、最终及实验时刻保存。

- `*_hit_history.csv`：K、K_mass、速度均值/雷诺应力/Favre 应力、密度/压力/温度涨落、黏性耗散、涡量平方、散度、压力膨胀、积分尺度、Taylor 尺度、Reλ、Kolmogorov 尺度、Mt、螺旋度、梯度偏斜度/平坦度、总质量/总能量及 Parseval 误差等。
- `*_hit_spectrum.csv`：三维壳谱、无散/膨胀分量和分速度能谱。
- `*_hit_spectrum_1d.csv`：一维纵向/横向谱 F11、F22、F33；`*_hit_correlation.csv`：相关函数及二阶结构函数。
- `*_hit_means.csv`、`*_hit_spectrum_mean.csv`：积分窗口信息及平均谱。衰减湍流对比实验以测站瞬时空间统计为主，不能拿整个时间窗均值代替测站值。
- `*_hit_metadata.txt`、manifest、日志：阶段标记、FFT 实现、源版本、分区、配置身份和停止原因。

后续将拼接各段**正式阶段**数据、按 `(step,time)` 去除重启边界重复记录，并核验三站输出。依次比较 E(k)、F11/F22、解析 K 衰减率、Reλ 与尺度；检查各向同性、散度和守恒误差。代码输出的壳谱满足 `K = Σ E_shell Δk`，不是把离散 E 列直接求和。实验全谱值、按相同解析带宽截断的值以及离散初始化目标需分别列出。

`epsilon_balance_residual` 同时包含空间数值耗散、时间离散和统计差分误差，不能直接标成“精确 SGS 耗散”；耗散相关指标在粗网格上尤其需要结合能量预算和网格收敛解释。不同时间段的累计均值也不能直接再做无权平均。

## 10. 打包结果并传回本地

在本轮作业结束/暂停且 runner 不再运行时执行：

```bash
bash collect_results.sh
```

会生成 `CBC9_results_时间戳_statistics.tar.gz` 及 `.sha256`：包含**所有九组、所有续算段**的统计 CSV、阶段记录、实际配置、运行日志、manifest、状态表、参考数据和编译日志。默认不包含 CGNS 大文件，足以开展第一轮谱与统计分析。无论是否全部完成，状态表都会如实保留。

若还需要检查/可视化流场或在另一台服务器恢复计算：

```bash
bash collect_results.sh --full
```

此包额外包含所有 CGNS 全场及检查点，可能达到数 GB 或更多；传输前查看 `du -sh results` 和生成包大小。统计包不具备恢复求解所需的全场检查点。

Windows PowerShell 下载示例（将文件名替换成实际生成的名称）：

```powershell
New-Item -ItemType Directory -Force D:\program\WCNS_v3\CBC9_return
scp your_user@your_server:/your/work/path/WCNS_v2.6_CBC_9cases_20261007/CBC9_results_时间戳_statistics.tar.gz D:\program\WCNS_v3\CBC9_return\
scp your_user@your_server:/your/work/path/WCNS_v2.6_CBC_9cases_20261007/CBC9_results_时间戳_statistics.tar.gz.sha256 D:\program\WCNS_v3\CBC9_return\
Get-FileHash D:\program\WCNS_v3\CBC9_return\CBC9_results_时间戳_statistics.tar.gz -Algorithm SHA256
Get-Content D:\program\WCNS_v3\CBC9_return\CBC9_results_时间戳_statistics.tar.gz.sha256
```

核对 SHA256 后，保留原包并告知本地文件位置，即可继续做九组统一对比。请回传整个统计包，不要只挑最后一个续算段或图片。

本次验证范围和实测结果见同目录 `validation.txt`。

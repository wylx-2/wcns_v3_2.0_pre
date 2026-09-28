# WCNS v2.1 Linux 服务器操作指南

本指南面向不含网格/结果的 `WCNS_v2.1` 精简源码目录。它说明上传、完整性核验、串行/MPI
编译、小规模预检、安全停止与重启。大型 Case06、长期湍流和目标三维翼型仍需另行制定资源
预算与验收方案；不要在登录节点直接运行计算。

## 1. 上传前准备

在开发机上从 `WCNS_v2.1` 的父目录打包，保留 UTF-8 文件名：

```text
tar -czf WCNS_v2.1.tar.gz WCNS_v2.1
scp WCNS_v2.1.tar.gz user@server:/work/user/wcns/
```

程序包不含 `cases/`、网格或运行结果。配置和 CGNS 网格应另行上传到独立算例目录，避免把
输入、输出和源码混在一起。WCNS 自有代码尚无对外许可证；服务器副本仍限内部使用。

## 2. 解包与完整性核验

```bash
cd /work/user/wcns
tar -xzf WCNS_v2.1.tar.gz
cd WCNS_v2.1
sha256sum -c PACKAGE_CONTENTS.sha256
cat WCNS_SOURCE_REVISION
```

散列必须全部显示 `OK`。若文件系统或传输工具改写了文件、缺少文件或出现额外来源不明的
文件，应重新上传，不要在损坏目录上继续编译。`WCNS_SOURCE_REVISION` 用于把服务器结果追溯
到开发仓库提交。

## 3. 环境要求

- CMake 3.20 或更高；
- 支持 C++20 的 GCC 或 Clang；建议 GCC 10 以上或集群当前受支持版本；
- MPI 构建需要同一工具链下的 MPI C++ 包装器和运行器，例如 OpenMPI 或 MPICH；
- 默认 CGNS 4.4.0 ADF 源码已位于 `third_party/cgns/`，配置与编译不需要联网或 HDF5；
- 精简包没有开发测试和 Python 验收脚本，`WCNS_BUILD_TESTS` 保持关闭。

使用模块系统的示例应按服务器实际模块名调整：

```bash
module purge
module load gcc/12 cmake/3.26 openmpi/4.1
cmake --version
c++ --version
mpicxx --version
mpirun --version
```

不要混用编译器和 MPI ABI，例如用 GCC 编译程序却在运行时加载由另一套编译器构建的 MPI。

## 4. 独立构建与安装

串行 Release：

```bash
cmake -S . -B build/serial \
  -DCMAKE_BUILD_TYPE=Release \
  -DWCNS_ENABLE_MPI=OFF \
  -DWCNS_BUILD_TESTS=OFF \
  -DWCNS_INSTALL_EXAMPLES=OFF
cmake --build build/serial --parallel 4
cmake --install build/serial --prefix install/serial
```

MPI Release：

```bash
cmake -S . -B build/mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DWCNS_ENABLE_MPI=ON \
  -DWCNS_BUILD_TESTS=OFF \
  -DWCNS_INSTALL_EXAMPLES=OFF \
  -DMPI_CXX_COMPILER="$(command -v mpicxx)"
cmake --build build/mpi --parallel 4
cmake --install build/mpi --prefix install/mpi
```

配置日志必须确认找到了预期 MPI。安装后检查动态库来源：

```bash
ldd install/mpi/bin/wcns_run
install/serial/bin/wcns_run --help || test "$?" -eq 1
```

`--help` 按程序约定打印帮助后返回 1；这不是构建失败。构建并行度应服从登录节点政策和作业
配额，不应直接取整机核心数。

## 5. 算例目录与小规模预检

建议目录彼此隔离：

```text
/work/user/wcns/WCNS_v2.1/           # 只读源码和构建
/work/user/wcns/cases/my_case/       # 配置与网格
/scratch/user/wcns/my_case/run-001/  # 本次输出
```

把配置中的 `mesh.file` 和 `output.directory` 改为服务器可访问路径。先做串行和少 rank 预检：

```bash
install/serial/bin/wcns_inspect_structured_mesh /work/user/wcns/cases/my_case/mesh.cgns
install/serial/bin/wcns_run --config /work/user/wcns/cases/my_case/case.wcns --dry-run
mpirun -np 2 install/mpi/bin/wcns_run \
  --config /work/user/wcns/cases/my_case/case.wcns --dry-run
```

人工核对程序版本、来源提交、网格维数/块数、总单元数、边界类型、参考 Mach/Reynolds 数、
湍流模型、Riemann 求解器、时间推进、rank 分区、输出目录和预计内存。随后只推进极少步，检查
残差、容许性修复计数、模型量、壁面量、检查点和输出文件均有限且符合预期。正式大计算前，
至少比较 1/2/4 rank 小规模终态或关键积分量。

## 6. Slurm 作业模板

以下模板只演示资源绑定；账户、分区、模块和路径必须按服务器修改：

```bash
#!/usr/bin/env bash
#SBATCH --job-name=wcns-preflight
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.out

set -euo pipefail
module purge
module load gcc/12 cmake/3.26 openmpi/4.1

source_root=/work/user/wcns/WCNS_v2.1
case_file=/work/user/wcns/cases/my_case/case.wcns
export OMP_NUM_THREADS=1

srun --ntasks="${SLURM_NTASKS}" --cpu-bind=cores \
  "${source_root}/install/mpi/bin/wcns_run" --config "${case_file}"
```

WCNS 当前以 MPI 为主，`OMP_NUM_THREADS=1` 可防止数学库或运行库意外超额使用线程。先在单节点
完成预检，再依据网格块数、每 rank 内存和 I/O 规模决定是否多节点运行。

## 7. 安全停止、检查点与续算

长作业必须启用检查点，并把 `run.max_wall_time` 设为比调度器时限短 5--10 分钟的秒数。
程序在完整步边界响应该上限或 `SIGINT/SIGTERM`，写入安全检查点后以退出码 2 结束。不要使用
`SIGKILL`；强制终止不能保证检查点原子完成。

停止后核对 manifest 的 `stop_reason=wall_time_checkpoint` 或
`stop_reason=user_signal_checkpoint`，并确认 `*.checkpoint.latest.cgns` 存在且没有 `.tmp`
后缀。续算配置使用：

```text
restart.path = ../run-001/output/my_case.checkpoint.latest.cgns
restart.mode = strict
output.directory = ../run-002/output
output.allow_existing = false
```

先对续算配置执行 `--dry-run`。网格和数值签名必须兼容；rank 数可以在合法分区范围内改变。

若服务器任务要从现有场切换重构、Riemann、湍流闭合、LU-SGS 或低 Mach 参数，应另建分支目录
并设置 `restart.mode=algorithm_change`。该模式仍要求同一网格，但重置 BDF 历史、定常收敛状态
和时间统计；RANS 字段布局不同时按新模型初始化。先进行 dry-run 和少量步数验收，再申请完整
资源，不能把这种分支解释为原计算的无缝续算。

## 8. 从 case05 t=250 分叉到显式 all-speed Roe

发布包提供 `examples/channel_retau180_from_t250_all_speed_roe.wcns`。把它复制到
case05 算例目录，保留 SCMM6、线性 MDCD、`diss=0.001` 和 SSPRK3，只替换真实网格/
checkpoint 路径、唯一输出目录和墙钟：

```bash
cp examples/channel_retau180_from_t250_all_speed_roe.wcns \
  /work/user/wcns/cases/case05/case05-all-speed-from-t250.wcns
```

必须确认下列配置：

```text
algorithm.profile = scmm6_wcns
algorithm.reconstruction = mdcd_linear
algorithm.riemann = roe_all_speed
algorithm.mdcd.diss = 0.001
run.cfl = 0.3
restart.path = /absolute/path/to/t250.checkpoint.latest.cgns
restart.mode = algorithm_change
```

`run.t_end` 是绝对终止时间；从 t=250 再计算 50 应写 `300`。首次切换算法必须用
`algorithm_change` 并采用全新 `case.name/output.directory`。先分别执行串行和 4-rank dry-run，
然后复制为 smoke 配置，把 `run.max_steps` 改为 2--5 并再次使用独立输出目录。启动摘要必须同时
出现 `riemann_solver=roe_all_speed` 和 `time(integrator=ssprk3`；首步时间应接续约 250，且
守恒量有限、Riemann 回退未异常增长。正式运行与后续 `strict` 续段见
[`all-speed-roe.md`](all-speed-roe.md)。

必须同时上传生成 t=250 checkpoint 时的原始网格，不能重新生成同名网格。已核验种子的 mesh
signature 是 `18099232003167909757`；若 dry-run 报 `checkpoint mesh signature differs`，应
更正 `mesh.path`，不得绕过签名检查。

## 9. 正式计算前卡口

只有以下项目全部通过后才提交大型作业：

1. 包散列、来源提交、编译器/MPI/动态库记录完整；
2. 网格检查和串行/2-rank dry-run 通过，边界及模型配置经人工复核；
3. 小步运行没有非有限值、异常 floor repair、异常壁面量或输出覆盖；
4. 检查点安全停止和异 rank 续算至少演练一次；
5. 作业的核时、内存、磁盘、检查点频率和最大墙钟有明确上限；
6. Case06、长期 RANS/LES 和三维翼型分别有独立物理验收指标，不能把本机微型门禁当作工程
   精度结论。

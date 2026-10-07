# 最小构建包运行说明

入口是 `wcns_run --config <文件> [--dry-run]`。
配置中的 `mesh.path` 和谱文件路径相对于配置文件所在目录解析；
`output.directory` 相对于启动进程的工作目录解析，正式计算建议使用绝对路径。

MPI rank 数由 `mpiexec -n P` 指定。`partition.mode=auto_split` 可以继续拆分结构块，
因此进程数不必等于原始块数；每方向的最小单元数约束仍必须满足。

HIT 配置使用 `initial.type=hit`；`hit.type=decay` 或 `forced`；
`hit.initialization=shell_spectrum` 或 `analytic_random_phase`。
CBC 准备阶段需要 `hit.spectrum_file` 和正的 `hit.preparation.time`，
先预演化、重匹配，再清零正式时间和统计。HIT 的所有方向必须周期。

续算使用 `restart.path`、`restart.mode=strict`，并指定新的输出目录。
准备检查点继续准备；正式检查点不会重复已完成的准备阶段。
已有检查点可以在符合分区约束的不同 rank 数下恢复。

返回码 0 表示达到正常停止条件；2 表示步数、墙钟或用户停止，需要检查 manifest；
3 表示数值失败。不能把预演化结束或返回码 2 当成正式算例已经完成。

通用配置模板和许可文件随包提供。具体物理参数、网格及统计对比设置请使用主仓库的
案例说明；不要用模板默认值代替物理验证算例的参考配置。

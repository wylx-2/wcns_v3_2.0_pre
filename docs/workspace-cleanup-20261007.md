# 工作区清理与 v2.6 交付整理

按用户确认的范围，保留主仓库源码、最新算例与计算结果，删除旧版本副本及可重建产物。
远程同步目标为 `v2-pre`，即 `wylx-2/wcns_v3_2.0_pre` 的 `main` 分支。

## 清理结果

- 清理 39 个明确列出的旧目录/归档及构建、验证产物，逻辑文件体积合计 2,039,352,637 字节。
  此统计包括本次恢复网格验证时新建的临时文件，不等同于清理前原有工作区的净减少量。
- 旧 v2.2—v2.6 源码发布副本和独立案例副本已删除；案例真值保留在主仓库 `cases/manual/`。
- 清理旧 CMake 构建树、临时 Python 依赖缓存及重复依赖解包目录。
- 14 个结果/当前 CBC 包目录的文件清单、大小和修改时间摘要在清理前后完全一致。
- 为待删除目录中的数据建立了 746 条保留校验记录；其中 616 个文件复制到本地
  `maintenance/20261007/retained-data/`，共 51,826,306 字节，其余在主目录已有相同 SHA256 文件。
- 旧 `wcns-v1.0-state-only` 工作树已移除，提交仍由主仓库 `release/v1.0.1` 分支保存。
- 旧独立 `wcns_v3_release` 仓库清理前已创建并验证完整 Git bundle，保存在本地维护记录目录。

原有 `WCNS_v2.6_CBC_9cases_20261007` 目录、上传压缩包及校验码原样保留，方便核对正在服务器上运行的版本。
`PHengLEI`、`mesh_test` 和原有计算结果不属于本轮旧副本清理范围。
主仓库内的历史源码、文档、旧发布脚本和示例仍保留；不会因已经被 Git 跟踪就直接删除。

## 最新交付方式

最小包由 [package_minimal.py](../tools/package_minimal.py) 从明确的 Git 提交生成：

```bash
python tools/package_minimal.py --ref HEAD --output ../dist/WCNS_v2.6_minimal_20261007
```

求解器和所有物理模型代码保留，独立包关闭开发工具与测试，包含 CGNS、FFTW 离线源码依赖。
构建及安装方法见 [最小包说明](minimal-build-readme.md)，验证范围见 [验证记录](minimal-build-validation.md)。

新增案例的大网格以无损 XZ 压缩保存于 `cases/grid-archives/`，避免把 210 MB 单文件直接提交到 GitHub。
17 个路径共用 13 个内容寻址归档；新克隆后运行：

```bash
python tools/restore_case_grids.py
```

恢复程序验证压缩文件及解压网格的 SHA256，拒绝覆盖内容不同的已有文件。
原有本地网格保持原位置。计算结果继续保留本地，新增结果和构建产物由 `.gitignore` 排除，
此前已跟踪的验证资料仍保留在 Git 历史及当前工作树中。

详细本地清单、数据校验记录和旧仓库 bundle 位于工作区根目录的 `maintenance/20261007/`。

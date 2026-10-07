# WCNS v2.6 最小构建包

本包包含最新求解器源码、CGNS/FFTW 离线依赖、许可证及必要配置模板。
所有物理模型和已实现算例入口保留；网格、计算结果、开发测试和网格处理工具在主仓库中维护。
精确版本见 `WCNS_SOURCE_REVISION`，文件完整性见 `PACKAGE_CONTENTS.sha256`。

## Linux 构建

需要 CMake ≥3.20、支持 C++20 的编译器、MPI 开发环境和 make/Ninja。
先按服务器要求加载编译器与 MPI 模块，然后在包根目录执行：

```bash
sha256sum -c PACKAGE_CONTENTS.sha256
bash build_linux.sh
```

默认编译 MPI、FFTW、CGNS，Release 模式，2 个构建任务。
可用 `JOBS=4 bash build_linux.sh` 调整编译并行度。
等价的 CMake 命令为：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DWCNS_ENABLE_MPI=ON -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_FFTW=ON \
  -DWCNS_BUILD_TESTS=OFF -DWCNS_BUILD_TOOLS=OFF
cmake --build build --target wcns_run --parallel 2
cmake --install build --prefix install --component WCNS
```

CGNS 和 FFTW 均从 `third_party/` 的本地源码包构建，无须下载。
FFT 是本地 FFTW 变换加 MPI slab 转置；不要求另装 FFTW-MPI。
**必须保留 `WCNS_BUILD_TOOLS=OFF`**，因为本包有意省去辅助工具源码。
安装使用 `--component WCNS`，仅安装程序及配套文档，不要求额外编译上游 CGNS 命令行工具。

## 运行

```bash
mpiexec -n 8 ./build/wcns_run --config /path/to/case.wcns --dry-run
mpiexec -n 8 ./build/wcns_run --config /path/to/case.wcns
```

网格和配置可直接使用已有 CBC 九算例包或主仓库中的案例。
本包不包含预编译的服务器二进制文件；编译与执行应使用同一种 MPI 环境。
`--dry-run` 检查配置、分区、几何指标与初始化，不推进时间。

## 版本与计算范围

程序版本仍为 v2.6，含周期山、SD7003、压缩拐角、自由衰减/强迫 HIT、
CBC 预演化重匹配及解析谱初始化。CBC、解析谱初始化可通过参数切换。
本次整理只改变构建和交付方式，没有改变求解器数值公式。
通用模板保留原文件名 `examples/full_case_template_v2.3.wcns`；HIT 等完整专用配置在主仓库中。

完整源码、测试、报告和案例：<https://github.com/wylx-2/wcns_v3_2.0_pre>

压缩包是按白名单从指定 Git 提交生成的源码快照；不会包含 `.git`、临时构建目录、
计算检查点或用户本地环境。仅构建验证与短时测试在本机执行，正式长计算应在服务器完成。

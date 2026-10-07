# Third-party notices

## CGNS 4.4.0

WCNS 随仓库保存 `third_party/cgns/CGNS-4.4.0.zip`，其 SHA-256 为
`72AD499936089102C1D77E3CDD78A623D61D1C6AD94C4DF3B23F36DB9DDFAD6A`。当前构建仅使用
静态 ADF 后端，并关闭 HDF5、Fortran、共享库和 CGNS 自带测试。归档来自 PHengLEI 随附的
未修改 CGNS 4.4.0 发布包；CGNS 的许可证通知如下：

> This software is provided 'as-is', without any express or implied warranty.
> In no event will the authors be held liable for any damages arising from the
> use of this software.
>
> Permission is granted to anyone to use this software for any purpose,
> including commercial applications, and to alter it and redistribute it
> freely, subject to the following restrictions:
>
> 1. The origin of this software must not be misrepresented; you must not claim
>    that you wrote the original software. If you use this software in a
>    product, an acknowledgment in the product documentation would be
>    appreciated but is not required.
> 2. Altered source versions must be plainly marked as such, and must not be
>    misrepresented as being the original software.
> 3. This notice may not be removed or altered from any source distribution.

CGNS 原文还注明该许可证借自 zlib/libpng License，并取代 CGNS 软件此前使用的 LGPL；
参见归档内未修改的 `license.txt`。

该通知是归档根目录 `license.txt` 的完整实质条款。归档自身同时保存原始 `license.txt`，应与
任何包含 CGNS 源码的再分发一同保留。

## PHengLEI 参考关系

WCNS 的算法和数据结构设计参考 PHengLEI 的公开实现思路，文档中标明了具体对照位置；当前
仓库清单未包含 PHengLEI 自有源文件或其预编译库。上述 CGNS 原始归档是唯一从 PHengLEI
随附第三方目录复制的文件，仍按 CGNS 自身许可证分发。

## FFTW 3.3.11

The v2.6 spectral backend optionally links FFTW. The original source archive, SHA256, COPYRIGHT and COPYING are in `third_party/fftw/`. Upstream is licensed under GPL v2 or later. WCNS implements its own MPI transpose layer; native FFTW MPI APIs are not used. `WCNS_ENABLE_FFTW=OFF` builds the internal radix-2 alternative.

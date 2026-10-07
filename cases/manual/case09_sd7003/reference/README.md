# SD7003 参考数据

`separation_table.csv` 为原文/用户文档的有限精度表值。Galbraith & Visbal：AIAA 2010-4737，https://doi.org/10.2514/6.2010-4737。`sd7003-uiuc.dat` 来自 https://m-selig.ae.illinois.edu/ads/coord/sd7003.dat，网格构造另含明确的尾缘圆钝化。

`*_digitized.csv` 是用户提供文档图 5.24、5.25 的黑色参考方点数字化，**不是原作者原始数组**。`source_*` 保存取点原图，`*_digitization_qa.png` 标记点位置，`*_provenance.json` 记录轴标定和 SHA-256。±4 原图像素是读图误差估计，不是源数据误差界限。部分遮挡点未保留，没有用红色“本文结果”或插值补齐。

速度图横轴是切向速度 / U∞，纵轴是壁法向距离 / c；文件名给出 x/c 站位。Cp 文件含上下表面点，未强行给前缘重叠点分配表面身份。Cf 为上表面。Tu 数值有设施和四舍五入差异；不可把实验设施数据混为一条精确真解。

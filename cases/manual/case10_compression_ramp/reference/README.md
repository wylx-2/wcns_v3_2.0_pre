# 压缩拐角参考数据

Porter & Poggie (2019), Physics of Fluids 31, 016104，https://doi.org/10.1063/1.5078938。作者算例页 https://engineering.purdue.edu/~jpoggie/ramp/index.html 区分细网格和长时间粗网格数据。2017 前期论文 https://engineering.purdue.edu/~jpoggie/papers/AIAA-2017-0533.pdf 的展宽为 10δ₀，本任务采用用户文档的 6δ₀。

`Cp_digitized.csv`、`Cfx_digitized.csv`、`van_driest_x80_digitized.csv` 是用户文档图 5.34、5.33 中黑色参考方点数字化，**不是原作者精确数组**。原图、标点 QA、像素标定和 SHA-256 一同保存。坐标范围列依据 ±4 原图像素；只包含清晰可识别的参考点，不插值补造遮挡数据。原始压力/热流/激波时间序列和完整应力数组本次未取得。

`inlet_similarity.csv/json` 是本程序生成的压缩性层流入口及其 BVP 参数/残差，**不是湍流 DNS 验证结果**。该剖面加体力转捩后，必须通过 x/δ₀=70、80 的发展检查。体力幅值/波形不在这些参考资料中完整给出，v2.5 的适配不能冒充原激励。

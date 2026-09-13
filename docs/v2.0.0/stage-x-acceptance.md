# v2.0.0 阶段 X 自动验收

结论：**阶段实现和非长算卡口已完成，定量物理卡口未关闭，禁止创建 X 候选标签。**

阶段 W 已获人工验收，X 以 NASA TMR SA-neg 作为首个湍流模型。验收始终分开
“模型公式/软件通路”、“独立公开规范算例”和“项目 Case06”；用户 Case06 未被读取、
改写或用作 SA 验收替代物。

## 1. 已关闭卡口

| 卡口 | 结果 | 证据边界 |
|---|---|---|
| 模型身份 | 通过 | `sa_neg` 固定为 TMR SA-neg modified-vorticity 1(c)；正分支保留 $f_{t2}$，trip 源和压缩修正关闭 |
| 公式级单元测试 | 通过 | 常数、$f_{v1}/f_{v2}/f_{t2}/f_n/f_w$、正/负分支、产生/破坏/cross diffusion、源 Jacobian 逐点比较 |
| SA 制造解算子 | 通过 | 32/64/128 周期标量对流—扩散—源离散的 $L_1/L_2/L_\infty$ 相邻网格观测阶均 $>1.9$ |
| 生产求解通路 | 通过 | 守恒量 $\rho\widetilde\nu$、初场/边界/halo、壁距、对流/扩散/源、平均流黏性和热流闭合、PH/SCMM 已接线 |
| 配置与输出 | 通过 | schema 2 严格校验；注册 `nu_tilde`、`mu_t_over_mu`、产生/破坏、壁距和负分支；模型残差进入 history 和停止判定 |
| checkpoint/restart | 通过 | 活动 SA 使用 checkpoint v2，保存 `NuTilde`、descriptor 和残差参考；连续运行/中途续算的平均流与 `nu_tilde` 一致，串行容差 $10^{-13}$、2-rank 容差 $10^{-12}$ |
| TMR 资产身份 | 通过 | 7 个公开原始资产 URL 和 SHA-256 固定在 `manifest.json`；未校验文件不能入场 |
| 规范网格转换 | 通过 | 平板 level 6/5/4 与 NACA0012 225/449/897 均以 double precision 转换，保留壁面索引和 C-grid 尾迹自连接 |
| 启动与首步有限性 | 通过 | 上述 6 个验收网格均通过 dry-run 和实际首步；平板 level 8 因仅 7 个法向顶点被高阶 Jacobian 严格拒绝，不转化为验收网格 |
| C-grid 与 MPI | 通过 | 新增“同块面物理壁面+尾迹自连接”权重/通量 halo 回归；NACA 225 的 1/2/4-rank 场最大差 $4.27\times10^{-14}$ |
| 本机完整回归 | 通过 | 串行 Release CTest 64/64；MPI Release CTest 116/116；算法规格 6/6；串行/MPI 安装及安装树 dry-run 通过 |

## 2. 尚未关闭的强制卡口

- [ ] TMR finite-flat-plate level 6/5/4 充分定常收敛后，$C_f$、阻力、速度剖面和首层 $y^+$ 进入冻结容差；
- [ ] TMR NACA0012 validation 225/449/897 在 $\alpha=0,10,15^\circ$ 收敛后呈系统网格趋势，最细层 `CL/CD/Cp/Cf` 进入冻结容差；
- [ ] TMR NACA0012 numerical Family II 的 `CL/CD/CM`、无 point-vortex 主分支和带 point-vortex 对照完成；
- [ ] 记录定常残差/载荷联合停止、负分支占比、无隐藏 floor/裁剪和可追溯 manifest。

## 3. 卡口发现的阶段依赖

规范近壁网格的全局显式 SSPRK3 稳定步长实测为：

| 网格 | 首步 `dt` |
|---|---:|
| flat plate level 6 / 5 / 4 | $8.67\times10^{-9}$ / $3.44\times10^{-9}$ / $1.36\times10^{-9}$ |
| NACA0012 225 / 449 / 897 | $4.56\times10^{-9}$ / $1.82\times10^{-9}$ / $6.78\times10^{-10}$ |

这些数值不是数值失败，而是 resolved-wall 网格在“全局显式步长”下的预期刚性。
按当前阶段顺序，本地伪时间和 LU-SGS 要到 AA 才实现；在此之前用本机全局 SSPRK3
完成六组多网格定常收敛不具备现实可行性。放宽残差、使用最粗 smoke、减少网格层级或
以 Case06 替代都会破坏已冻结的 X 验收规范，本报告不采用这些做法。

## 4. Git 决策

- 已验证的实现、测试、转换工具和文档可在 `stage/v2.0.0-x` 上形成可追溯提交；
- 不创建 `v2.0.0-x-candidate.1`，不合入 `release/v2.0.0`，不开始阶段 Y；
- 需由人工批准阶段依赖调整后，才能继续长算并关闭 X 候选卡口。

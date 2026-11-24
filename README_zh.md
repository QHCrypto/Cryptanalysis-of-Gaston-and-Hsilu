# Gaston 与 Hsilu 论文配套代码

本文件夹对应论文 **Cryptanalysis of Gaston and Hsilu**。通讯作者为 Meicheng Liu，邮箱为 liumeicheng@iie.ac.cn。

目标期刊：**Frontiers of Computer Science**。

通讯作者单位：中国科学院信息工程研究所，网络空间安全防御全国重点实验室，北京，中国；中国科学院大学网络空间安全学院，北京，中国。

公开仓库：https://github.com/QHCrypto/Cryptanalysis-of-Gaston-and-Hsilu

代码可用性：论文直接引用本公开仓库，不使用出版社分配的补充材料编号。内容为 Gaston 与 Hsilu 的差分、差分线性、原像和 cube 分析源代码、核验脚本、必要的小型输入及部分计算记录。论文元数据保存在 `ARTICLE_METADATA.json`。

这是轻量代码包，保留源代码、必要的小型输入、结果摘要和统一核验入口。按作者要求，不包含完整 DRAT 证明证书、大型匹配证书和可重新生成的多项式中间文件。原研究目录中的文件保持不变。

## 运行方法

需要 Python 3.10 以上、NumPy、支持 C++17 的 GCC/G++、make 和 GMP 开发库。以下核验不需要 Gurobi：

```bash
python3 run_checks.py --profile quick
python3 run_checks.py --profile datf
```

请在本 README 所在文件夹执行。程序在临时副本中编译、运行，将日志和汇总写到 `validation/日期时间_配置/`，不覆盖归档输入。成功时各组显示 `PASS`。

`quick` 包含参数和模式核验、差分路径核验、Gaston 支撑空间核验、差分线性基础测试和抽样、三轮 Hsilu 条件及具体分支枚举、低轮 cube 实验、Huang 口径的复杂度计算。`datf` 单独重新计算两个四轮端点的 DATF Algorithm 5 模型值。

## 文件对应关系

以下路径均位于 `research/` 中；完整入口表见英文 [README](README.md)。

| 内容 | 主要目录 |
|---|---|
| 原始置换和模式接口 | `hsilu_ref/`、`gaston_ref/`、`nist_compatibility_20261109/` |
| Gaston 差分界 | `gaston_active_pattern/`、`Gaston_差分路径重量下界/` |
| Hsilu 重量 80 的三轮路径 | `Hsilu_trail/search_below82_20261117/` |
| 差分线性分析 | `差分线性/gaston_dl/`、`差分线性/hsilu_dl/`、`差分线性/路径分段汇总_20261027/` |
| Gaston 原像分析 | `gaston分析/` |
| Hsilu 二轮原像分析 | `代数攻击/hsilu_2r/`；当前图对应 `hsilu_no_start.svg` |
| Hsilu 三轮原像分析 | `代数攻击/Hsilu3r_new/`，具体核验在 `validation_current_prefix/` |
| Gaston cube | `GastonCube/`、`Cube/gaston_rank_certificate/` |
| Hsilu cube | `Cube/hsilu_quartic_screen/`、`Cube/src/hsilu_7r_batch_audit.py` |

包根目录的 `scripts/preimage_accounting.py` 是当前复杂度口径的计算入口，三个结果依次为 `2^118.06`、`2^106`、`2^126`。按正文采用的 Huang 口径计入预计算和高斯消元。这些数值是模型内的复杂度估计。

## 核验范围

- Hsilu 重量 80 是一条已验证的具体路径，不表示已证明全局最优。
- Gaston 支撑空间核验会复核已有枚举摘要和零空间计算；轻量包不重新执行全部历史枚举，也不包含完整证明证书。
- 三轮 Hsilu 的 `2^23` 分支枚举只针对保存的具体条件系统，不是执行完整 128 位原像求解。
- Gaston 四轮 cube 求和会重新执行；五轮只用已有完整求和记录重建矩阵并核对秩 125，不重跑 `2^31` 点求和。六轮的秩假设仍按正文处理。
- Hsilu 七轮的完整批量计算只做计数分析与缩小实例核验，不实际执行完整规模。
- 差分线性的模型估计、有限样本观察与严格恒等式分别记录。DATF 对 Hsilu 四轮端点返回零，不表示实际相关为零。

MILP 搜索另需 Gurobi 及有效许可；可通过 `GUROBI_HOME` 指定本机安装位置。可选 SAT 搜索还需 `python-sat` 及相应求解器。包内没有求解器安装文件、许可证或虚拟环境。

`SOURCE_MANIFEST.json` 保存原路径和原文件哈希；`PORTABILITY_CHANGES.json` 记录整理副本的适配；`SHA256SUMS` 用于检查分发文件完整性。运行 `sha256sum -c SHA256SUMS` 可核对。内部旧 README 保留了历史探索信息，当前稿件对应关系以本说明和英文 README 为准。

核验所需的参数表、轮常数、模式初始值和差分路径数据位于 `research/Gaston_Hsilu_CTP/verification_inputs/paper_parameters.tex`，代码可独立运行。

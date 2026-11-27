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

`quick` 包含参数和模式核验、差分路径核验、Gaston 支撑空间核验、差分线性基础测试和抽样、三轮 Hsilu 条件及具体分支枚举、低轮 cube 实验、Huang 口径的复杂度计算。本轮增加对 Hsilu 重量80所选相消核心的 C++ 完整扩展，枚举全部16384个兼容中间输出。

`datf` 已同步为当前论文的两轮差分前缀加权：DATF Algorithm 3 和 HNW26 在 G4、H4 上均分别得到 `2^-8`、`2^-9`。HNW26 采用作者公开 PRESENT 代码（`getBias_Opt2`/`getBias_Opt4`）已有的联合末端处理，将各活跃列的掩码合并后经 `L^T` 回传；该处理归于原作者，不作为本文新方法。基本 Algorithm 1 逐列乘积读出得到的两个 `0` 单独保留为对照。首轮四个分支全部计入，第二轮的完整差分分布通过仿射像消元求和，保留差分坐标之间的依赖。计算设置和独立核验见 [当前四轮计算说明](research/差分线性/joint_dl_20261125/README.md)。

旧完整四轮 Algorithm 5 可通过 `python3 run_checks.py --profile datf-legacy` 复现；`3/2048` 和 `0` 仅作为历史模型结果保留。

## 文件对应关系

以下路径均位于 `research/` 中；完整入口表见英文 [README](README.md)。

| 内容 | 主要目录 |
|---|---|
| 原始置换和模式接口 | `hsilu_ref/`、`gaston_ref/`、`nist_compatibility_20261109/` |
| Gaston 差分界 | `gaston_active_pattern/`、`Gaston_差分路径重量下界/` |
| Hsilu 重量 80 的三轮路径 | `Hsilu_trail/search_below82_20261117/` |
| 当前四轮差分线性计算 | `差分线性/joint_dl_20261125/`；包根目录 `scripts/current_dl_results.py` |
| 差分线性搜索及历史记录 | `差分线性/gaston_dl/`、`差分线性/hsilu_dl/`、`差分线性/路径分段汇总_20261027/` |
| Gaston 原像分析 | `gaston分析/` |
| Hsilu 二轮原像分析 | `代数攻击/hsilu_2r/`；当前图对应 `hsilu_no_start.svg` |
| Hsilu 三轮原像分析 | `代数攻击/Hsilu3r_new/`，具体核验在 `validation_current_prefix/` |
| Gaston cube | `GastonCube/`、`Cube/gaston_rank_certificate/` |
| Hsilu cube | `Cube/hsilu_quartic_screen/`、`Cube/src/hsilu_7r_batch_audit.py` |

包根目录的 `scripts/preimage_accounting.py` 是当前复杂度口径的计算入口，三个时间结果依次为 `2^118.54`、`2^106.97`、`2^126.23`，机器可读结果见 `CURRENT_PREIMAGE_RESULTS.json`。在正文的预计算、高斯消元与候选计费模型内保留 `C_s*(1+T_s)+2^64`；候选数量和求解成本公式不变。这些数值是模型内的时间估计，未新增系数构造或候选处理的计费项。

## 本次同步的结果

Hsilu 三轮差分路径搜索算法和状态无需修改：置零节点为4、9、10、12，轮重量为 `(16,14,50)`，活跃 S 盒数为 `(8,7,22)`。新核验入口重跑所选核心的扩展，与独立逐列 DDT 验证相互核对；未重跑全部历史搜索。

当前差分线性七组结果见包根目录 `CURRENT_DL_RESULTS.json`。旧编号 G3a、G4a、H4a 分别对应当前 G3、G4、H4。

| ID | 轮数 | 实验相关度 | DATF Algorithm 3 | HNW26 逐轮近似 | 求值方式 |
|---|---:|---:|---:|---:|---|
| G2 | 2 | 1 | 1 | 1 | 整体 |
| G3 | 3 | -2^-2.42 | -3/16 | -3/16 | 整体 |
| G4 | 4 | 2^-7.96 | 2^-8 | 2^-8 | 整体，两轮差分前缀 |
| G5 | 5 | 未测量 | 2^-59 | 2^-59 | 拼接 (1,2,2) |
| H3 | 3 | 1 | 1 | 1 | 整体 |
| H4 | 4 | 2^-9.01 | 2^-9 | 2^-9 | 整体，两轮差分前缀 |
| H5 | 5 | 未测量 | 2^-63 | 2^-63 | 拼接 (0,2,3) |

原始实验计数保持不变：G2、G3、H3 各 `2^22` 对，G4 为 `2^26` 对，H4 为 `2^28` 对。使用未舍入的实验相关度计算 G3、G4、H4 的数据复杂度指数，分别为4.83、15.93、18.03；G5、H5 的拼接估计对应118、126。模型内精确有理数求值不等于真实固定置换相关度的严格证明。

当前附录顺序为：A 模式接口证书、B 立方计算、C Hsilu 重量 80 路径、D Gaston 差分界、E 差分线性端点、F 原像配置。内部历史记录可能沿用旧编号。

## 核验范围

- Hsilu 重量 80 是一条已验证的具体路径，不表示已证明全局最优。
- Gaston 支撑空间核验会复核已有枚举摘要和零空间计算；轻量包不重新执行全部历史枚举，也不包含完整证明证书。
- 三轮 Hsilu 的 `2^23` 分支枚举只针对保存的具体条件系统，不是执行完整 128 位原像求解。
- Gaston 四轮 cube 求和会重新执行；五轮只用已有完整求和记录重建矩阵并核对秩 125，不重跑 `2^31` 点求和。六轮的秩假设仍按正文处理。
- Hsilu 七轮的完整批量计算只做计数分析与缩小实例核验，不实际执行完整规模。
- 差分线性的模型估计、有限样本观察与严格恒等式分别记录。当前 Hsilu 四轮 DATF 和采用作者联合末端处理的 HNW26 前缀加权估计均为 `2^-9`；基本逐列读出和旧 Algorithm 5 的零值分别作为对照与历史记录保留。

MILP 搜索另需 Gurobi 及有效许可；可通过 `GUROBI_HOME` 指定本机安装位置。可选 SAT 搜索还需 `python-sat` 及相应求解器。包内没有求解器安装文件、许可证或虚拟环境。

`SOURCE_MANIFEST.json` 保存原路径和原文件哈希；`PORTABILITY_CHANGES.json` 记录整理副本的适配；`SHA256SUMS` 用于检查分发文件完整性。运行 `sha256sum -c SHA256SUMS` 可核对。内部旧 README 保留了历史探索信息，当前稿件对应关系以本说明和英文 README 为准。

核验所需的参数表、轮常数、模式初始值和差分路径数据位于 `research/Gaston_Hsilu_CTP/verification_inputs/paper_parameters.tex`，代码可独立运行。

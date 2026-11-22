# Hsilu 路径搜索：PySAT + Kissat

PySAT 负责生成 CNF 和基数约束，外部 **Kissat 4.0.4** 负责求解。默认使用论文《Multiple Rows Mixers and Hsilu》的 Algorithm 3、Table 9 参数，以及从标准 χ₅ 定义直接枚举得到的精确 DDT／LAT 支集。

当前差分粘贴区按要求排除了输入差分 `11111`：实际启用 **31 条子句**，保留 301 个合法差分对；原完整 DDT 有 317 个，其中被删除的 16 个恰好都是输入为 `11111` 的差分对。其余输入的关系保持精确，没有增加非法差分对。代码中先放完整的 46 条版本并全部注释，随后单独放当前启用的 31 条版本；线性 66 条不变。这个限制只作用于实际建模的 χ 层，末端未建模 χ 的输入 x2 不额外限制。

这是完整差分关系的受限子集。路径模型直接使用当前启用的子句，不将手动 CNF 与完整或受限差分表进行对应检查。

## S 盒 CNF 与路径建模已分开

独立生成工具在 `sbox_cnf/`，详细操作见 [sbox_cnf/README_zh.md](sbox_cnf/README_zh.md)。常用流程：

1. 运行 `python3 sbox_cnf/generate.py`，生成局部差分／线性 CNF。
2. 将 `sbox_cnf/generated/paste_into_model.py` 的全部内容手动复制到 `hsilu_sat.py` 顶部 `BEGIN MANUAL SBOX CNF` 至 `END MANUAL SBOX CNF` 的区域。
3. 照常用 `./run.sh` 搜索路径，不需要额外的排除输入配置或 CNF 对应检查。

默认 χ₅ 子句已预填。生成器不依赖 PySAT 或 Kissat，也不会改写建模文件；路径模型不导入生成器，只实例化粘贴区里的子句。建模文件中的 DDT／LAT 表计算保留用于精确重量计算，不再用于生成 CNF。局部文字 1…5 是输入行 0…4，6…10 是输出行 0…4，Python 元组中不包含 DIMACS 行末的 0。

要复现当前的 31 条差分子句，使用 `python3 sbox_cnf/generate.py --exclude-diff-input 0b11111 --output sbox_cnf/generated_no11111`，复制该目录下的代码块。省略排除参数生成的是完整的 46 条差分关系，仍输出到 `generated/`，不会自动改变当前粘贴区。

手动切换为 46 条时，将 46 条区块取消注释、31 条区块全部注释；切回 31 条时反向操作。两组子句只启用一组，不需要修改其他配置。

## 本次确定的三层计数模型

```text
x0 ── χ0 ──> y0 ── L0 ──> x1 ── χ1 ──> y1 ── L1 ──> x2
             ↑ 第 0 层计数  ↑ 第 1 层计数              ↑ 第 2 层计数
```

`--rounds 3` 默认建立 **两个 χ 层、两个线性层**。活跃盒数统计 **y0、x1、x2** 的非零列：第一轮统计 χ0 输出（第一轮线性层输入），其余轮次仍统计 χ 输入；x2 是第三个 χ 层的输入，不为第三个 χ 层创建 SAT 变量或约束。

此位置约定用于差分和线性模式的活跃盒等值、上下界、指定活跃位置和总活跃数，也用于结果中的 `active`。`active_count_states` 明确记录相应状态名。活跃比特约束 `--bits`、`--min-bits`、`--max-bits` 及 `active_bits` 保持原来的 x0、x1、x2 位置（`reference` 拓扑首项为 y0），由 `bit_count_states` 标明。精确差分／线性重量仍根据每个 χ 的输入输出计算。

标准 χ5 是置换，因此合法差分传播中每列的输入、输出同时为零或同时非零；计数移到 y0 后，合法路径集合不变，但活跃变量直接连接到 y0，CNF 的传播行为和求解耗时可能变化。修改前的实验结果及模型快照保留原计数方式。

状态是五个 64 位字，行顺序为 x[0]…x[4]，列 z 是每个字的第 z 个最低有效位编号（0…63）。文本图从高位到低位输出，因此 z=0 位于最右边。χ 的局部编号采用行 0 为最低位：

```text
χ(a)[i] = a[i] XOR ((NOT a[(i+1) mod 5]) AND a[(i+2) mod 5])
```

前向线性层依论文使用左旋。差分关系是 `x[r+1] = L(y[r])`；线性掩码关系是 `y[r] = transpose(L)(x[r+1])`，转置电路使用右旋。常数加法不改变差分和线性支集，也不改变平方相关重量；程序不报告包含轮常数的相关符号。

## 直接运行

```bash
cd /home/yxb/分析/Hsilu_trail/pysat_kissat

# 三层活跃盒数分别固定为 6、6、27；没有默认时间上限
./run.sh --mode diff --active 6 6 27 --prefix results/diff_6_6_27

# 线性路径示例：固定活跃盒数；这只是查询，不保证存在解
./run.sh --mode linear --active 16 2 6 --prefix results/linear_16_2_6

# 每层下界、上界，和总活跃盒数的等值约束
./run.sh --mode linear --min-active 1 1 1 --max-active 16 2 6 \
  --total-active 24 --prefix results/linear_total24

# 只导出 CNF，不启动求解
./run.sh --mode diff --active 6 6 27 --export-only --prefix results/diff_6_6_27

# 需要限制运行时间时，显式给出秒数
./run.sh --mode diff --active 6 6 27 --timeout 600 --prefix results/timed
```

每次查询产生 `前缀.cnf`、`前缀.log` 和 `前缀.json`。求解器运行日志写入 `.log`，可以在另一终端运行 `tail -f results/diff_6_6_27.log`。SAT 的 JSON 包括各状态、活跃盒数、活跃比特数和精确重量；UNSAT／timeout 的 JSON 记录状态，不会保留上次同名前缀的旧路径。

第一轮默认非零，约束第一轮计数状态 y0 至少有一个活跃列。已有正的活跃计数时不额外添加冗余非零子句。需要研究全零情况时可显式指定 `--allow-zero`，Python 对应 `first_nonzero=False`。保留 `--nonzero` 兼容原命令。

其他约束示例：

```bash
# 第 1 层只有列 0、7、20 活跃，其他列均不活跃
./run.sh --active-positions 1:0,7,20 --prefix results/positions

# 固定 χ 输入，五个字按行排列；整数使用 Python 风格十六进制前缀
./run.sh --fix x0=0x1,0,0,0,0 --prefix results/fixed_input

# 固定某层线性输入（也就是 χ 输出）
./run.sh --fix y1=0,0,0x81,0,0 --prefix results/fixed_middle

# 活跃比特数的等值／下界／上界分别由 --bits、--min-bits、--max-bits 给出
./run.sh --bits 1 11 83 --prefix results/bits
```

计数列表长度必须等于 `--rounds`；多个 `--fix`、`--active-positions` 可以叠加。`--full` 建立完整的 n 个 χ 层和 n 个线性层，仍统计 n 个 χ 输入。

## 精确重量输出

SAT 搜索条件与原程序一样，主要是活跃数和固定状态条件，**没有把重量最小化作为优化目标**。

- 差分重量为 `−log2(DDT[a,b]/32)`。
- 线性重量按论文 §2 的定义为 `−log2((LAT[a,b]/32)^2)`，这里 LAT 保存 Walsh 和。它是平方相关重量。
- `modeled_sbox_weights` 给出实际建模的各 χ 层精确重量。
- `terminal_sbox_completion` 在得到路径后，根据最后的输入状态逐列选择绝对 DDT／LAT 值最大的合法输出，附上这个未建模 χ 层的输出与重量。
- `completed_trail_weight` 是加入该末层补全后的总重量。它是当前这条路径补全的重量，不是全局最优值。对 χ₅ 差分，每个固定输入的非零 DDT 值相等；线性模式选择固定输入下重量最小的末层输出。

末层补全只是求解后的表查询，不增加第三个 χ 的 CNF，也不影响活跃数搜索。原参考拓扑额外省略了首层 χ，故该拓扑不报告完整路径总重量。

## Python 接口

使用已经安装的解释器：`../tools/sat/venv/bin/python`。

```python
from ascon import AsconDifferential, AsconLinear

A = AsconDifferential(nrounds=3)
A.set_nactive_sboxes_at_round(0, 6)
A.set_nactive_sboxes_at_round(1, 6)
A.set_nactive_sboxes_at_round(2, 27)
A.output_prefix = 'results/python_diff'
answer = A.has_trail()  # 可传 timeout=秒数
if answer is True:
    print(A.trail())
    print(A.result())
elif answer is False:
    print('UNSAT')
else:
    print(answer)  # timeout / interrupted；不要写 if answer: 来判断 SAT
```

保留常用原接口名称，同时支持直接导入 `HsiluDifferential`、`HsiluLinear`：

| 功能 | 接口 |
| --- | --- |
| 每轮活跃盒等值／上下界 | `set_nactive_sboxes_at_round`、`set_min_nactive_sboxes_at_round`、`set_max_nactive_sboxes_at_round` |
| 每轮活跃比特等值／上下界 | `set_nactive_bits_at_round`、`set_min_nactive_bits_at_round`、`set_max_nactive_bits_at_round` |
| 精确指定所有活跃位置 | `set_active_sboxes_at_round(r, positions)` |
| 固定 χ 输入／输出 | `set_input_nonlinear_state_at_round`、`set_output_nonlinear_state_at_round` |
| 固定 L 输入／输出 | `set_input_linear_state_at_round`、`set_output_linear_state_at_round` |
| 固定轮输入／输出 | `set_input_round_state`、`set_output_round_state` |
| 固定总活跃盒数 | `has_trail(nactive_sboxes=k)`；这是等值，不是上界 |
| 读取各层值与计数 | `input_nonlinear_values`、`output_nonlinear_values`、`input_linear_values`、`output_linear_values`、`active_sboxes_at_round`、`nactive_sboxes_at_round`、`nactive_bits_at_round`、`nactive_sboxes` |
| 输出结果／时间 | `trail()`、`trail_at_round(r)`、`result()`、`solving_time()` |
| 导出 DIMACS | `export_cnf(path, nactive_sboxes=None)` |

设置状态的参数是 `(r, word0, word1, word2, word3, word4)`。默认截断模型中 `input_values()` 返回 x0，`output_values()` 返回最后的 x[n−1]，`trail()` 只显示实际连接的状态。没有 χ 的末轮不能调用 `trail_at_round(n−1)` 取得 χ 输出，应读取最后的输入或 JSON 末层补全。

STP 专有表达式接口不兼容：`solver` 现在是 PySAT `CNF`，`vars('x', r)` 等返回 5×64 个正整数文字，`nactive_sboxes_var(r)` 返回 64 个活跃文字，`additional_constraints` 接收 CNF 子句列表。原有直接操作 STP BitVec 的自定义约束需要改写；原内部 z／zz 辅助位向量不再提供。每次 `has_trail` 启动独立 Kissat 进程，不是增量 Kissat；临时总数约束不会累积到下一次查询。

## 差分线性层：合并只使用一次的 XOR 中间量

默认 `--xor-encoding binary`：每个字级 XOR 创建 64 个输出变量，每一位使用 4 条子句；字旋转仅重排变量编号。默认启用等价变量合并，线性层最终结果直接作为下一状态的变量；使用 `--no-merge-equalities` 可恢复原来的相等子句连接。

可选 `--xor-encoding inline` 会把不是层输出、且只有一个后续使用位置的中间字代入后续 XOR。多处使用的中间字保留，避免重复展开。合并时同时处理旋转量相加取模 64，以及重复 XOR 项的抵消。

论文 Hsilu 参数中，z0…z4 表示线性层输入，z5…z17 表示内部字。默认合并 z7、z8、z10 后，相关表达式为：

```text
z11 = z1 XOR ROL(z2,54) XOR z6
z12 = z11 XOR ROL(z0,14) XOR ROL(z1,44)
z15 = z13 XOR ROL(z5,10) XOR ROL(z3,1)
```

每层从 13 个二输入 XOR 变成 7 个二输入 XOR 和 3 个三输入 XOR，省去 3×64 个中间变量。两个差分线性层共省去 384 个变量，XOR 子句总数保持不变。这里保留了原来的层末相等约束，便于单独比较消去中间量这一项改动。

```bash
# 原模式
./run.sh --mode diff --active 6 6 27 --xor-encoding binary \
  --prefix results/diff_6_6_27_binary

# 合并模式
./run.sh --mode diff --active 6 6 27 --xor-encoding inline \
  --prefix results/diff_6_6_27_inline
```

Python 接口对应 `HsiluDifferential(3, xor_encoding='inline')`。当前合并模式用于差分传播；线性掩码模式仍用 `binary`。

`--xor-max-inputs` 控制合并后每个 XOR 的输入数上限，默认 3，可设为 2…8。直接编码 k 输入 XOR 的输出等式需要 2^k 条子句：两个二输入 XOR 共 8 条，合成一个三输入 XOR 仍为 8 条；三个二输入 XOR 共 12 条，合成一个四输入 XOR 则需 16 条。因此默认不继续扩成更长的 XOR。对于论文 Hsilu 参数，默认上限 3 已经能消去上述三个中间字。

运行输出和结果 JSON 会记录 XOR 模式、删除的中间字编号及各输入长度的 XOR 数量。变量减少不代表求解一定更快，应固定同一组 S 盒子句和计数条件后比较日志。

## 线性计数扫描

`scan_linear.py` 实现参考 `main.py` 的逐组查询方式：固定第一层 i，遍历第二层 j，第三层加上界 `max_total−i−j`，逐组保存 CNF、日志和结果。例如：

```bash
../tools/sat/venv/bin/python scan_linear.py \
  --first 2 --start 3 --max-total 18 --output results/linear_scan
```

可加 `--stop-on-sat`，或 `--timeout 600` 限制每个实例。脚本不沿用原代码针对其他参数组写死的 `min2` 剪枝表。超时会单独报告，不能用超时推断无路径。

## 参考实现的兼容选项

按本次要求，默认值是 `profile='hsilu', relation='exact', topology='core'`。原参考模型可显式选：

```bash
./run.sh --profile mrm13 --relation reference --topology reference \
  --mode diff --active 6 6 27 --prefix results/old_mrm13
```

- `profile=mrm12/mrm13/mrm14/mrm15` 保存的是参考源码的参数，不应与论文 Table 9 的 Hsilu 混称。mrm14、mrm15 无对应原线性支集，只可使用 exact 关系。
- `relation=reference` 原样使用源码局部约束枚举出的支集；这些支集不等于标准 χ₅ 的精确 DDT／LAT。复现模式遇到标准表不允许的转移时，精确重量返回 null 并列出具体位置。
- `topology=reference` 复现原 `remove_outer_layer=True` 的 `L→χ→L` 三层计数语义，首层在 y0 计数。
- 为响应本次修改，新接口的默认 `remove_outer_layer=True` 只省略最后一轮，使用 `χ→L→χ→L`；需要原行为时必须显式指定 `topology='reference'`。`remove_outer_layer=False` 仍建立完整轮。
- 三个选项相互独立；精确论文模型使用默认值即可。

参考参数与支集快照在 `reference_data.json`，来源路径和 SHA-256 已记录。`extract_reference.py` 可以从本机 MRM 目录重新生成；它仅用具体整数执行提取出的局部约束方法，不依赖 STP。原代码许可证保留在 `LICENSE`。

## 验证与实验范围

```bash
../tools/sat/venv/bin/python test_hsilu_sat.py
```

测试包括精确重量所用的 DDT／LAT、独立生成器输出、手动修改 CNF 后不按预期支集拒绝求解结果、参考支集快照、论文参数的全部 320 个基向量、线性转置内积恒等式、CNF 线性电路的全部基向量，以及 Kissat 的非零 SAT、冲突 UNSAT、基数约束、重复查询和输出文件更新。

`examples/diff_nonzero.*` 和 `examples/linear_nonzero.*` 是实际运行的非零示例；它们只验证基本搜索与输出流程，不能用于判断困难计数实例的求解速度。`examples/diff_6_6_27.cnf` 只是导出的查询，没有宣称该计数实例已经求解。

论文参数提取文本和表 9 页面保存在本目录，供追溯此次实现使用的定义。

基数编码可用 `--card-encoding` 选择，默认使用 `totalizer`。六种取值为
`seqcounter`、`totalizer`、`mtotalizer`、`kmtotalizer`、`sortnetwrk`、`cardnetwrk`。
该选择统一用于每层活跃 S 盒数、每层活跃比特数，以及临时总活跃数的
等于/至少/至多约束；S 盒子句、线性层和计数语义保持原样。

```bash
# 三层统一使用 kmtotalizer。
./run.sh --mode diff --active 6 6 27 --card-encoding kmtotalizer \
  --prefix results/diff_6_6_27_kmtotalizer

# 六种编码逐个运行；每组 600 秒，输出目录自动按时间命名。
../tools/sat/venv/bin/python benchmark_cardinality.py --timeout 600

# 两层活动数 (6,5)，模型为 χ→L，记录 UNSAT/SAT/timeout。
../tools/sat/venv/bin/python benchmark_cardinality.py --active 6 5 --timeout 600

# 两层活动数 (27,2)，六个独立 Kissat 进程并行，每个求解进程仍为单线程。
../tools/sat/venv/bin/python benchmark_cardinality.py --active 27 2 --jobs 6 --timeout 600

# 可自行选择子集、重复次数和新的输出目录。
../tools/sat/venv/bin/python benchmark_cardinality.py \
  --encodings seqcounter kmtotalizer cardnetwrk --repeats 3 --timeout 600 \
  --output results/cardinality_selected_repeats
```

Python 接口对应 `HsiluDifferential(3, card_encoding='kmtotalizer')`；
`HsiluLinear` 和 `scan_linear.py --card-encoding ...` 也支持相同选择。

`benchmark_cardinality.py` 固定论文 Hsilu、差分模式、`core` 拓扑、`binary` XOR
和 Kissat 默认参数。活动数默认 `(6,6,27)`，可用 `--active` 指定，列表长度
决定轮数，例如 `--active 6 5` 使用两个计数状态 `y0,x1` 和一个 χ、一个 L。
它使用主程序当前启用的手工
S 盒子句，默认逐个启动新求解进程；`--jobs 6` 同时运行六个独立的单线程
Kissat 进程，模型分别构造，输出文件互不覆盖。不会自动选择 31/46 条版本或添加新约束。
同一轮实验的模型在进程启动时载入，并保存到 `model_snapshot/`。
`metadata.json` 记录实际 S 盒子句、模型和 Kissat 的 SHA-256、版本及实验条件。
为保护旧结果，指定的输出目录必须尚不存在。

每组保留 `.cnf`、`.log`、`.json`，SAT 时另存 `.trail.txt`。
`summary.csv` 和 `summary.json` 每组结束后更新，记录变量/子句数、建模时间、
求解进程墙钟时间、日志中的 Kissat CPU 时间、包含导出/结果解析的总时间，
以及搜出路径的精确完成重量。并行结果按完成顺序写入，`metadata.json` 记录并发数；
比较并行与串行实验时应区分墙钟时间和 CPU 时间。
超时单列为 `timeout`，不视为 UNSAT；不同编码得到的可行路径重量可以不同，
因为本实验只有活动数约束，没有最小化重量。默认每种编码运行一次，适合作为
首次对照，不能据此断言某编码在所有实例上更快。

## 双侧活动连接、变量合并与可选预处理

默认使用 `totalizer`，启用 `shared_activity=True` 与 `merge_equalities=True`；
`preprocess=False`，额外预处理需显式开启。

- 对已建模的 χ 层，同一个活动变量同时等价于输入五比特的 OR 和输出五比特的 OR。
  相比单侧连接，每盒增加 6 条子句。末端省略的 χ 仍只连接输入，不创建虚假的输出状态。
  第一轮计数位置仍标记为 y0。对置换 χ 的合法差分／线性支集，这些是冗余约束；
  若手工粘贴的其他关系允许单侧为零，新约束会排除这些点，可用 `--no-shared-activity` 关闭。
- 线性层边界使用变量别名，避免为相等的位再分配变量和添加双向蕴含子句。
  两个差分线性层共减少 640 个变量和 1280 条子句；同时启用双侧活动连接后，
  净减少 640 个变量和 512 条子句。状态读取、固定状态、比特约束仍使用相同接口。
- `--preprocess` 使用 PySAT `Processor.process(rounds=1)` 的默认技术开关。
  原始公式保存为 `.original.cnf`，交给 Kissat 的公式保存为 `.cnf`。
  SAT 后通过同一个 Processor 的 `restore()` 恢复原始变量，并检查原始 CNF 和路径约束。
  不进行手工 S 盒 CNF 与 DDT 支集的对应检查。
  `--export-only` 仍只导出原始建模 CNF；预处理和赋值恢复在实际求解时配套进行。

```bash
# 两项建模改进，预处理关闭（默认）
./run.sh --active 6 6 27 --card-encoding totalizer --prefix results/improved

# 两项改进加可选预处理
./run.sh --rounds 2 --active 6 5 --preprocess --prefix results/improved_pre

# 当前 y0 计数位置下，关闭改进的对照组
./run.sh --active 6 6 27 --no-shared-activity --no-merge-equalities \
  --prefix results/baseline_y0

# 两个实例 × 三个版本，六个独立进程并行测试
../tools/sat/venv/bin/python benchmark_improvements.py --jobs 6 --timeout 600
```

Python 构造函数具有相同参数，例如
`HsiluDifferential(3, card_encoding='totalizer', shared_activity=True, merge_equalities=True, preprocess=True)`。
每次求解重新进行预处理；临时总数约束不累积，修改后查询不复用旧的恢复栈。
结果的 `solve_stats` 记录预处理／恢复耗时及处理前后公式规模。
`solver_variables` 是 DIMACS 最大变量编号，`solver_used_variables` 是实际出现的不同变量数，
预处理后的编号可能有空洞，不能仅比较最大编号判断消去了多少变量。
基准 `solve_seconds` 是 Kissat 阶段墙钟时间，`total_seconds` 包含建模、预处理、
导出、Kissat、恢复和结果解析。`--timeout` 限制 Kissat 阶段，不包含额外预处理。

预处理接口说明：https://pysathq.github.io/docs/html/api/process.html

## 列循环移位对称

通过 `--column-symmetry` 显式启用，Python 对应 `column_symmetry=True`。
当第一层有正的活跃盒等值／下界约束时，在完整查询 CNF 末尾添加单位子句
`[a[0][0]]`，即第一轮计数状态 y0 的第 0 列必须活跃。
所有状态共同循环移位与逐列 χ、字内旋转、XOR 相容，活跃数与比特数也不变。
因此每条第一轮非零的路径都可整体移位，使某个活跃列来到第 0 列；
每个移位轨道至少保留一个代表，但一般不是唯一代表，不能据此声称 64 倍加速。

第一轮现在默认非零，因此启用对称选项时可以直接固定第 0 列活跃。
仅当显式 `--allow-zero` 且没有已知的第一层正活跃盒约束时，才改为
`a[0][j] -> a[0][0]`（j=1..63），保留全零轨道。
为防止在不具有循环对称性的查询中误用，启用时不接受固定状态、固定活跃位置
或 `additional_constraints`。直接修改底层 CNF 的调用方需自行保证循环移位不变性。
此规则不依赖重量优化；仍按既有方式在求解后计算精确重量。

```bash
# 现有两项改进加列对称约束，预处理关闭
./run.sh --rounds 2 --active 6 5 --column-symmetry --prefix results/symmetry_6_5

# 在关闭两项改进的基线上单独测试对称约束
./run.sh --rounds 2 --active 6 5 --column-symmetry \
  --no-shared-activity --no-merge-equalities --prefix results/symmetry_baseline_6_5
```

`benchmark_cardinality.py` 同样接受 `--column-symmetry`。
`metadata.json` 和结果的 `solve_stats.column_symmetry` 记录是否开启。

## 省略首轮 χ 的模糊重量搜索

独立入口 `fuzzy_weight_search.py`，不改变原 `hsilu_sat.py` 的手工 CNF 或普通计数搜索。
核心只有 `y0 → L0 → x1 → χ1 → y1 → L1 → x2`：一个 χ 层、两个线性层。

搜索估计重量为：

```
W_est = 2*A(y0) + sum_z f(HW(x1[z])) + sum_z f(HW(x2[z]))
f(0..5) = [0, 2, 3, 3, 4, 4]
```

后两轮使用 `f = 2*a + u + v`，其中 a 是已有活动标志，u 等价于五比特和 ≥2，
v 等价于五比特和 ≥4。两个阈值均双向约束，不能由求解器自由降低代价。
固定三层活动数时，`2*sum(A)` 为常数，使用所选基数编码约束剩余 u、v 的总数；
活动数未全部固定时，使用 PySAT PBEnc 的 BDD 编码处理权重为 2、1 的伪布尔上界。

此模式只按模糊重量上界找候选，并非直接优化真实重量。给定同一重量界，可枚举若干
不同的核心路径，并在 `best.json` 保存其中真实补全重量最低的候选。
每个候选均保存 CNF、Kissat 日志、JSON，以及补全后的六个状态文本图。
JSON 包含 `estimated_round_weights`、`estimated_total_weight`、`exact_round_weights`、
`completed_trail_weight`、`correction_by_round`、`completed_states`。

```bash
# 固定三层活动数，估计重量 ≤81，搜索一个候选
../tools/sat/venv/bin/python fuzzy_weight_search.py \
  --bound 81 --active 6 6 27 --column-symmetry \
  --output results/fuzzy_81_fixed

# 不固定活动数，同一估计界内最多枚举十条不同核心路径
../tools/sat/venv/bin/python fuzzy_weight_search.py \
  --bound 81 --candidates 10 --column-symmetry \
  --output results/fuzzy_81_open
```

省略 `--timeout` 时默认不限时，每次求解直到 SAT 或 UNSAT，也可以按 Ctrl+C 手动停止。
显式添加 `--timeout 600` 才将每个候选的 Kissat 求解时间限制为 600 秒；
`--candidates 10` 最多进行十次候选查询，每次分别计时。
输出目录必须尚不存在，防止覆盖旧实验。使用 `--export-only` 可仅导出估计重量查询。
默认 totalizer、双侧活动连接、变量合并、首轮非零；预处理和列对称均默认关闭，按需显式开启。

χ1 使用主程序当前手工粘贴的差分 CNF；当前为排除输入 11111 的 47 条。
省略的 χ0、χ2 用完整精确 DDT 补全：按固定 y0 为每列选择概率最大的首轮输入，
再为 x2 选择合法的末轮输出。若要允许 χ1 输入 11111，应手动切换主程序为完整 46 条。
不做手工 CNF 与完整支集的运行时对应检查；若用户的自定义 CNF 允许非法 DDT 转移，
结果会标明中间层非法点和空的真实总重量，不将该候选选作 best。

对标准 χ5，新估计仅在输入为 `00111` 的五个循环移位时低估 1；其他输入精确。
因此真实最优补全重量满足：

```
W_exact = W_est + t0 + t1 + t2
```

t0 是首轮 y0 中最优补全成本为 3 的列数；t1、t2 分别是 x1、x2 中属于
`00111, 01110, 11100, 11001, 10011` 的列数。
这只是完整 DDT 后处理，不把这些修正量加到搜索目标中。

不固定省略的端点 x0 或 y2。Python 的固定状态接口只接受 y0、x1、y1、x2。
`FuzzyHsilu(bound=81)` 可设置 `estimate_bound` 后重复查询，新的成本上界不累积到原 CNF。
`block_current_core()` 只排除已找到的具体核心；枚举后 UNSAT 指该查询与这些排除条件下
没有剩余候选，不应脱离所用 S 盒关系、活动数或枚举范围解释为全局最优性。

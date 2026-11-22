# Gaston 线性层活跃 S 盒 MILP

使用 C++17 和 Gurobi C++ API，对 Gaston 的完整 320 比特线性层建立精确 MILP。支持差分／状态传播与线性掩码传播，并分别支持固定输入活跃数量、最小化输出，以及固定输出活跃数量、最小化输入。省略两个数量参数时，要求输入非零，并最小化输入加输出活跃 S 盒列数。

固定活跃列数量很小时，可使用 [direct_enumeration](direct_enumeration/README.md) 中不依赖求解器的精确枚举程序。输出恰好两个活跃列的差分实例已完整枚举，最小输入活跃列数为 48，本机耗时为毫秒级；命令示例和证书见该子目录。

## 编译与运行

```bash
cd /home/yxb/分析/gaston_active_pattern
make
```

默认使用 `/opt/gurobi1302/linux64`。也支持 `make gaston_active_pattern`、`make SRC=gaston_active_pattern.cpp`；其他安装目录可用 `make GUROBI_HOME=/实际路径` 指定，主版本不为 13 时还需调整 `LDLIBS` 中的 Gurobi 库名。

四种用法如下，将 `3` 换成需要的数量即可：

```bash
# 差分：输入恰好 3 个活跃 S 盒，求最小输出活跃数
./gaston_active_pattern --mode difference --input-active 3

# 差分：输出恰好 3 个活跃 S 盒，求最小输入活跃数
./gaston_active_pattern --mode difference --output-active 3

# 线性掩码：输入掩码恰好 3 个活跃 S 盒，求最小输出掩码活跃数
./gaston_active_pattern --mode linear --input-active 3

# 线性掩码：输出掩码恰好 3 个活跃 S 盒，求最小输入掩码活跃数
./gaston_active_pattern --mode linear --output-active 3
```

`--input-active` 和 `--output-active` 最多选一个，数量范围为 `0..64`。显式指定时固定的是**恰好等于 K**，不是至多 K。省略 `--mode` 时默认为 `difference`。

两个数量参数都省略时，添加 `sum(input[i][j], i=0..4, j=0..63) >= 1`，目标为 `min (sum(active_input) + sum(active_output))`，即最小化输入、输出两端活跃 S 盒列数之和。默认还添加 `sum(input[i][0], i=0..4) >= 1`，要求输入第 0 列非零，以减少循环移位对称性。这里要求输入至少有一个非零比特，并没有固定输入比特重量或活跃列数。同一列位置在输入和输出都活跃时，两端各计一次；不计算两端列位置的并集。

```bash
# 输入非零，默认要求输入第 0 列活跃，求最小输入加输出活跃列数
./gaston_active_pattern

# 同样的搜索，额外要求 E=0
./gaston_active_pattern --theta-effect zero

# 只要求整个输入状态非零，关闭第 0 列代表元约束
./gaston_active_pattern --no-symmetry

# 线性掩码：输入掩码非零，求最小输入加输出掩码活跃列数
./gaston_active_pattern --mode linear
```

显式指定 `--input-active 0` 或 `--output-active 0` 时仍按原来的精确数量处理，不添加默认的输入非零约束。

完整求解日志保持开启，**不设求解时限**。按 Ctrl+C 请求终止；如果已有可行解，仍会复核并输出该解。只有 `Status: 2 (已证明最优)` 才标记为最优，否则显示当前可行解、最优下界和 MIPGap。没有可行解时，不输出解报告和解文件。一般 K 的最优性证明耗时依实例而定。

可选参数：

| 参数 | 用途 |
|---|---|
| `--prefix diff_in3` | 求解前导出 `diff_in3.lp`；有解且复核通过后导出 `diff_in3.sol` 和 `diff_in3.json` |
| `--threads 8` | 指定求解线程数；默认 0，由 Gurobi 自动选择 |
| `--theta-effect any\|zero\|nonzero` | 仅差分模式：不限制 E／要求 E=0／要求 E≠0，默认 `any` |
| `--no-symmetry` | 关闭第 0 列活跃的循环移位代表元约束 |
| `--export-only` | 仅导出 LP，不求解；需要 `--prefix` |
| `--help` | 查看帮助 |

前缀的父目录需要已存在。建议不同运行使用不同前缀；同名前缀会覆盖对应的新输出文件，没有生成新解时以前的解文件仍可能存在。普通运行不指定前缀时只输出控制台结果，不写文件。JSON 中的 64 位字以十六进制字符串保存，避免浮点解析丢失精度。

控制台和 JSON 分别报告活跃列数及比特重量，并通过 `total_active`（控制台的 `total`）报告两端活跃列数之和。JSON 使用 `objective_side` 表示目标：显式固定数量时为 `"input"` 或 `"output"`；省略两个数量参数时为 `"input+output"`，且 `objective = input_active + output_active = total_active`。此时 `fixed_side`、`fixed_count` 为 `null`，`input_bit_weight_min` 为 `1`。

## 状态与传播定义

状态由 5 个 64 位行组成。第 j 列是 `(state[0][j], ..., state[4][j])`，**该 5 比特向量非零时，计为一个活跃 S 盒**。每端最多 64 个；不是统计 320 个比特的汉明重量。

所有比特下标均模 64，字中的最低位编号为 0。`ROL(v,s)` 将位置 j 移至 `j+s mod 64`。控制台比特图从左到右显示 `j=0..63`，与通常十六进制展开的阅读方向不同。工作区其他程序若使用 plus-index 坐标，需按 `j_plus = -j mod 64` 换算。

采用原设计的 Algorithm 1、Table 6，线性层顺序为：

```text
X --ρeast--> A --θ--> B --ρwest--> Y
```

参数和字级关系：

```text
east = (0, 60, 22, 27, 4)
west = (0, 56, 31, 46, 43)
twin = (25, 32, 52, 60, 63)
r = 1, s = 18, u = 23

A[i] = ROL(X[i], east[i])
P = A[0] XOR ... XOR A[4]
Q = XOR_i ROL(A[i], twin[i])
E = P XOR ROL(P,1) XOR Q XOR ROL(Q,18)
Y[i] = ROL(A[i] XOR ROL(E,23), west[i])
```

差分模式满足 `Y = λ(X)`。线性掩码模式始终把 `input` 解释为输入掩码 α、`output` 解释为输出掩码 β，关系是 **`α = λᵀ(β)`**。固定输入掩码时求的是 `β = λ⁻ᵀ(α)`，不会把转置误作逆映射。

掩码转置的独立字级表达式为：

```text
D[i] = ROR(β[i], west[i])
P = D[0] XOR ... XOR D[4]
α[i] = ROR(D[i] XOR ROR(P,23) XOR ROR(P,24)
           XOR ROR(P,23+twin[i]) XOR ROR(P,41+twin[i]), east[i])
```

原始规范的本地副本为 [Gaston 设计论文](../gaston_ref/gaston.pdf)，公式可对照 [工作区参考实现](../MitM/src/gaston.py)。本目录的编译、运行和验证不依赖这些外部文件。

## θ-effect 核内／核外约束

MILP 入口 `gaston_active_pattern` 新增 `--theta-effect any|zero|nonzero`，用于自由选择是否要求 E 等于零。此参数目前仅支持 `--mode difference`；`--mode linear` 下显式传入任意取值（包括 `any`）都会报错。线性模式不传此参数时照常运行。

按照论文第 9 页式 (2)、第 11 页 Definition 5，θ-effect 核是 `{A : E(A)=0}`，这里的 **A 是进入 θ 的状态**，即完整线性层输入 X 经 `ρeast` 后的状态。程序中的 `effect_j` 表示上节公式中的 E 的第 j 比特，尚未做最后的 `ROL(E,23)`。

| 参数取值 | 搜索范围 | 添加的约束 |
|---|---|---|
| `any` | 不限制 E，默认行为 | 无 |
| `zero` | 核内，E=0 | `effect_j = 0`，j=0..63 |
| `nonzero` | 核外，E≠0 | `sum_j effect_j >= 1` |

`nonzero` 要求至少一个 E 比特为 1，不要求所有比特为 1。`zero` 只要求 E=0，不要求 P、Q 分别为零；后者会把搜索范围错误地缩小到 PQ-kernel。

对于 `nonzero`，另外加入输入、输出活跃列数均至少为 1 的冗余约束：E≠0 蕴含 X≠0，而 λ 可逆进一步蕴含 Y≠0。这使固定数量 K=0 的矛盾可以直接被预处理检出，不改变可行域。

例如，固定输出活跃列数为 6：

```bash
# 核内搜索：E=0
./gaston_active_pattern --mode difference --output-active 6 --theta-effect zero

# 核外搜索：E!=0
./gaston_active_pattern --mode difference --output-active 6 --theta-effect nonzero

# 不限制，等价于省略 --theta-effect
./gaston_active_pattern --mode difference --output-active 6 --theta-effect any

# 也可以固定输入数量；E 始终按 rho_east(X) 计算
./gaston_active_pattern --mode difference --input-active 6 --theta-effect zero
```

对于 `zero` 和 `nonzero`，无论固定输入还是固定输出，均使用包含 P、Q、E 的稀疏正向电路建立 `Y=λ(X)`，直接约束 E；不再为固定输出的情形另外构造逆矩阵。`any` 保留原来的模型方向选择。E 的零／非零性质在同步循环移位下不变，因此默认的第 0 列活跃约束对这些限制仍然有效。

E=0 时 θ 对 A 的作用是恒等，但完整线性层仍包括两个行移位，满足 `Y[i]=ROL(X[i],east[i]+west[i])`，并不要求 X=Y。两端活跃列数仍分别按各自的列位置计算。

求解后会从输入 X 独立重算 E，并核对选定条件。控制台显示 E 的十六进制值、比特重量和是否在核内；差分结果 JSON 增加：

```json
{
  "theta_effect_constraint": "zero",
  "theta_effect_hex": "0x0000000000000000",
  "theta_effect_weight": 0,
  "in_theta_kernel": true
}
```

以上仅展示新增字段，不代表某个具体数量问题的完整求解结果。若给定数量与核条件不相容，求解器会报告不可行；没有可行解时不会生成新的 `.sol`、`.json` 文件。

## MILP 约束与精确性

输入、输出状态每端使用 320 个二元变量，每列各有一个二元活跃变量。对任意一列的五个比特 `b[i]` 和活跃变量 `a`，加入：

```text
b[i] <= a             对 i=0..4
a <= sum_i b[i]
```

这使 `a=1` 当且仅当该列非零。不限制 E 时，差分固定输入和掩码固定输出使用 P、Q、E 或转置的共享中间量；反方向通过 GF(2) 高斯消元构造精确逆矩阵，再对其中的异或表达式进行公共子表达式提取和树分解。限制 E 时的模型方向见上节。两个数量参数都省略时，差分模式使用稀疏正向电路，线性掩码模式使用稀疏转置电路。

异或建模遵循更新后的 `milp-gurobi-cpp` 技能。**模型中的所有变量都是二元变量**；异或使用 Jeroslow 奇偶多面体不等式，不使用整数松弛变量。对于 `z = XOR(x[1],...,x[n])`，把输入和输出合并为 `v=(x[1],...,x[n],z)`，对每个奇数大小的下标集合 S 加入：

```text
sum(v[i], i in S) - sum(v[i], i not in S) <= |S| - 1
```

每个基本块完整展开奇数子集，得到该块的精确凸包：

| 输入数量 | 实现 | 辅助二元变量 | 不等式数 |
|---|---|---:|---:|
| 2 | 直接 xor2 凸包 | 0 | 4 |
| 3 | 直接 xor3 凸包 | 0 | 8 |
| 4 | xor3 + xor2 | 1 | 12 |
| 5 | 两个 xor3 | 1 | 16 |
| ≥5 | 以 xor3 为主体的平衡树；输入数为偶数时有一个末端 xor2 | 按树结构确定 | 每个 xor3 为 8 条，xor2 为 4 条 |

逆矩阵的多个输出表达式之间先执行公共子表达式提取（CSE）：统计重复出现的变量对，按出现次数选择互不重叠的变量对，用一个共享 xor2 变量替换其全部出现位置，然后重复统计。这样，多个方程会使用同一个共享变量及其凸包约束。此处的 xor2 属于技能中的 CSE 例外；提取完成后的每个长表达式仍按 xor3 树构造。原有的 P、Q、E、F、G 共享字也继续显式复用。

这些约束精确表示二元取值上的异或传播，包括异或抵消。“精确凸包”指每个基本异或块；多个共享块组成的整个模型，不据此宣称已给出全部可行状态的整体凸包。

固定输入时加 `sum(active_input)=K`，目标为 `min sum(active_output)`；固定输出时交换两端。两个数量参数都省略时，加输入 320 个比特之和至少为 1，目标为 `min (sum(active_input) + sum(active_output))`。此时还加入 `sum(active_output)>=1` 的有效冗余约束：λ 和 λᵀ 均可逆，所以输入非零蕴含输出非零。所有模型的相对和绝对 MIPGap 容许值均设置为 0。

Gaston 线性层、转置以及逆映射均与全体行同时循环移位可交换。因此，在 K>0 时，固定数量的那一端可以要求第 0 列活跃；两个数量参数都省略时，则要求输入第 0 列的五个比特之和至少为 1。任何非零可行模式都可以移位得到这样的代表元，两端活跃数量不变，因此这些约束不改变相应问题的最优值。它们减少移位对称性，但不保证每个等价类只保留一个代表元。K=0 时不加代表元约束；使用 `--no-symmetry` 可以关闭。若自行修改 LP 并额外指定绝对列位置，应关闭或移除该代表元约束。

程序用确定性样本构造初始可行解作为 MIP start；样本只帮助求解器找到解，**不会限制搜索空间或作为目标下界**。两个数量参数都省略时，分别从单活跃列输入和单活跃列输出构造候选状态对，选择两端活跃列数之和最小的可用样本；启用对称性约束时，同步移位两端以满足输入第 0 列活跃。模型本身没有固定任一端数量。设置 E 条件时会先过滤不满足条件的样本；若没有合适的样本，直接交给求解器搜索，不提交不符合约束的初始解。每次输出前，均使用独立于 MILP 约束构建过程的字级运算复核传播关系、活跃标记、固定数量或输入非零条件、目标值及 E 条件。

该结果精确描述单个线性层两端的活跃列关系。模型不包含 χ 的 DDT／LAT、多轮连接或差分概率／线性相关权重。

## 验证

```bash
make check
```

验证脚本额外需要 Python 3 和 `gurobipy`；主求解程序只需要 Gurobi C++ API。脚本从“单个输入比特流向哪些输出位置”独立构造 320×320 二元矩阵，并检查：

- 正向及转置矩阵的逆关系。
- 导出模型中全部变量均为二元变量，所有 xor2/xor3 块均包含完整且系数正确的奇数子集不等式，四输入异或采用 xor3+xor2，并检查共享子表达式确实被复用。
- 四种模式下 K=0 和 K=1 的求解结果及 JSON 证书。
- 导出的实际 MILP 对零状态、全 1 状态、边界位置、稀疏状态和稠密状态的传播与活跃标记。
- 错误目标比特与错误活跃标记的不可行性。
- 核内／核外约束的实际 LP 系数、两种固定方向中的 36 组固定状态检查，包括 P、Q 都非零但 E=0 的 row twin、PQ-kernel 中的 vortex、循环移位后的核内状态，以及单比特扰动后的核外状态。
- E 在 `rho_east(X)` 上计算，且报告的是移位 23 位之前的 E；核条件的 JSON 输出、与全零状态相容／不相容的情况，以及显式 `any` 与省略参数时导出模型的一致性。
- 线性模式使用 `--theta-effect` 时的明确报错，包括三种取值和不同参数顺序。
- 省略两个数量参数的 8 个模型：差分的三种 E 条件及线性掩码模式，分别开启／关闭对称性约束。检查输入非零约束、输入第 0 列五比特约束及目标中 128 个活跃变量的单位系数，并对 80 组固定状态检查可行性及两端活跃列数之和。
- 默认差分、关闭对称性的默认差分、核外默认差分及默认线性掩码的实际命令行求解：加载初始可行解后发送 SIGINT，检查结果报告、JSON 中的总活跃列数及最优性状态。回归测试不要求完成这些自由数量问题的全局最优性证明，也不把初始可行解当作最优结果。
- 导出模式与参数错误处理。

单活跃列在循环移位等价下只有 31 种非零取值，可以独立穷举得到精确参考最优值：

| 传播模式 | 固定条件 | 另一端最小活跃数 |
|---|---|---:|
| difference | 输入 = 1 | 18 |
| difference | 输出 = 1 | 57 |
| linear | 输入掩码 = 1 | 59 |
| linear | 输出掩码 = 1 | 17 |

这些数值分别对应四个固定数量问题，不是自由选择两端数量时的分支数。

本机已通过 `make check`，上述四组已求得最优的模型、解和 JSON 证书保存在 [examples](examples/) 目录。另已实测 SIGINT 中断后的可行解报告与证书导出。

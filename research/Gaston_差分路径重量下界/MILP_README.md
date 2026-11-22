# Gaston 三轮 lossy MILP

研究口径为标准Gaston-320的非零XOR差分路径，输入输出自由，目标为局部DDT限制重量之和的下界。χ采用《Multiple Rows Mixers and Hsilu》第8.1.1节表14的lossy模型，线性层逐位精确编码。

**最新进展（2026-10-24）**：结合已有单层传播结果与新的精确GF(2)活动集合穷尽计算，已得到三轮重量下界44；不代表最优值已知等于44。C_1=6、7、8的双核分支已由独立代数枚举排除，C_1=9解析下界为45，因此已正常停止旧MILP队列并保留日志。方法、证明依赖及复现命令见[双核分支排除44.md](双核分支排除44.md)，最终校验见[proof JSON](results/support_space_proof_20261024.json)。下文保留MILP模型和历史运行说明。

## 数学模型

记第r层χ的输入输出为a_r、b_r，实际建模拓扑为

\[
a_0\xrightarrow{\chi_{\rm lossy}}b_0
\xrightarrow{\lambda}a_1\xrightarrow{\chi_{\rm lossy}}b_1
\xrightarrow{\lambda}a_2\xrightarrow{\chi_{\rm lossy}}b_2.
\]

官方轮次最前面的λ可以消去，因为置换输入自由且λ可逆。验证程序对官方C实现的全部320个基向量核对传播，并确认λ的GF(2)秩为320。三个χ层均保留在模型中；默认模式下所有θ-effect都不受额外限制。

线性层λ=ρ_west∘θ∘ρ_east，参数为

```text
east = (0,60,22,27,4)
west = (0,56,31,46,43)
twin = (25,32,52,60,63)
r = 1, s = 18, u = 23
```

按64位左旋，A_i=ROL(b_i,east_i)，P=⊕A_i，Q=⊕ROL(A_i,twin_i)，
E=P⊕ROL(P,1)⊕Q⊕ROL(Q,18)，输出为ROL(A_i⊕ROL(E,23),west_i)。
所有状态位和异或辅助位均为二元变量。二输入异或使用4个奇偶多面体面，三输入异或使用8个面；五输入异或拆成两个三输入异或，四输入异或拆成三输入和二输入异或。旋转复用变量索引，P、Q、E由所有相关输出共享。

对每个χ列，设输入汉明重量为h、输出汉明重量为k，引入活动变量d∈{0,1}。四条支撑约束为

\[
h\ge d,\quad h\le5d,\quad k\ge d,\quad k\le5d.
\]

引入整数代价变量e∈{0,1,2}，增加两条约束

\[
h\le d+2e,\qquad e\le2d,
\]

并在目标中计入2d+e。对固定的输入输出列，最小化自动得到表14的代价：

| h | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|---:|
| 最小2d+e | 0 | 2 | 3 | 3 | 4 | 4 |

e仅用于代价编码，不是异或整数松弛变量。未达到最优时，某些可行赋值的e可能高于其最小值；结果因此同时记录求解器目标和按状态重算的标准lossy代价。

非零条件是第一层活动数至少1。默认加d_(0,0)=1，利用整体列循环平移选择一个代表；不会固定某一行或某个具体列值。所有真实路径均能经同步循环移位满足这一条件。可用--no-symmetry关闭。

默认模型没有固定活动分布、排除11111、固定核类型或设置额外目标上界。当前默认加入已完成的活跃列传播下界、核外条件下界、三轮数量组合约束，以及逐比特活动约束。核标志由实际effect精确决定，核内和核外仍全部保留。新增可选参数--kernel-both、--middle-active K、--weight-cap B用于显式分支研究；它们限制搜索范围，结论必须结合覆盖论证解释。

当前共3266个整数变量，其中3074个二元变量、192个0至2代价变量；默认11420条约束。基础版本为3072个变量、9346条约束。用--no-activity-cuts关闭传播表及组合约束，用--no-or-hull关闭逐比特活动加强；同时使用两个选项恢复基础模型。

本地结果的导入、适用范围、哈希和数量推论见[activity_profiles.json](references/activity_profiles.json)。重新导入使用python3 src/build_activity_profiles.py，然后make。新增约束及方法分析见[活跃列约束加强.md](活跃列约束加强.md)。

## 下界含义

设真正的最小路径重量为W_3，lossy模型的最优目标为L_3。完整DDT的每个合法转移都在lossy模型中，而且lossy代价不超过真实重量，因此

\[
\text{Gurobi有效目标下界}\le L_3\le W_3.
\]

求解器的incumbent只是L_3的上界。只有当它的全部χ转移通过DDT检查时，才能由其精确DDT重量给出W_3的上界。最优lossy目标也不必等于W_3。

当前采用Gurobi数值分支定界，保存模型、参数、日志及验证记录；未生成独立UNSAT证明文件。后续论文结论应区分该求解器下界、编码穷尽验证和额外的独立最优性证书。

## 编译和运行

默认Gurobi路径为/opt/gurobi1302/linux64，可用GUROBI_HOME覆盖。

```bash
make
# 以下两种单文件构建方式也支持
make SRC=gaston_3round_lossy.cpp
make gaston_3round_lossy
```

前台运行时必须使用尚不存在的输出目录：

```bash
./gaston_3round_lossy --output results/my_run --threads 16
```

程序设置Threads=16、MIPFocus=3、MIPGap=0、MIPGapAbs=0。线程参数只允许1至16。未设置TimeLimit，也未设置节点数、工作量、解数量等停止上限。MIPFocus=3优先改善下界，直到求解完成或用户主动中断。完整Gurobi控制台日志保持开启。

Ctrl+C会请求Gurobi正常中断；若已有解，程序保存当前可行解及目标、界和状态，不将中断标为最优。无可行解时不输出解值。

导出与验证：

```bash
./gaston_3round_lossy --output results/validation_new --threads 16 --export-only
python3 src/verify_lossy_model.py results/validation_new
```

已验证的源码和二进制可通过持久会话启动：

```bash
python3 src/launch_unlimited.py \
  --verification results/validation_new/verification.json \
  --threads 16
```

启动器核对验证记录中的源码和二进制哈希，建立独立进程会话，不设置外部超时；完整控制台输出保存在console.log。运行记录保存到results/latest_run.json，并保留源码快照。聊天回合结束不会主动停止该进程。

## 已完成的验证

[当前验证记录](results/activity_strengthened_validation_v2_20261023/verification.json)为PASS。验证直接读取C++导出的MPS，检查：

- 1408个异或单元的全部局部二元赋值，确认约束精确等价于异或。
- 192个χ列均具有指定的4条支撑、2条代价和10条逐比特活动约束，所有变量类型、上下界、目标系数与额外约束符合模型。
- 穷尽1024个χ差分对：lossy支撑共有962个，包含全部317个真实DDT合法对；对应最小代价全部不超过真实DDT重量。
- 每个线性层以零向量、全部320个基向量和32个确定性随机样本，共353例，对照官方C实现。
- λ满秩，以及同步循环移位与传播相容。
- 传播表、核标志、核外输出6列条件、数量选择变量及三轮组合约束与导出MPS逐项一致；独立组合数量关系复核36与44下界。

旧运行模型与当时的验证文件逐字节哈希一致。2026-10-23启动的四分支实验已于2026-10-24按用户要求停止，随后以核内专用编码及必要数量类型加强，按C_1=6、7、8顺序启动新队列。C_1=9已由解析下界45排除。详见[双核分支排除44.md](双核分支排除44.md)。

核内专用模式增加--kernel-reduced，要求同时使用--kernel-both、--middle-active 6至8、--weight-cap 43，保留活动传播约束。该模式直接通过变量索引表达旋转输出，以四变量偶校验表示完整E核，并从src/kernel_type_data.hpp加载已独立核验的必要数量类型。生成类型表使用python3 src/build_kernel_types.py；然后make。使用python3 src/kernel_campaign.py prepare验证三个模型，再使用launch启动其顺序队列。

## 输出文件

运行目录包含launch.json、console.log、source/和solve/。solve/中：

- model.lp、model.mps：完整原始模型。
- parameters.prm、settings.json：求解参数及模型口径。
- encoding.json：异或单元的变量映射，用于独立核对导出约束。
- progress.json：约每10秒记录求解时间、当前lossy目标上界、目标下界和节点数；结束后写入最终状态。
- incumbent.json：回调报告的可行候选及逐列DDT检查；结束后写入最终保留的可行解。进度中的最优目标以progress.json及控制台日志为准。
- incumbent.sol：求解结束且有解时保存完整变量赋值。
- result.json：求解结束时的状态、时间；有解时包含目标、下界及gap。

invalid_ddt_columns>0时，exact_ddt_weight记为null，不把该lossy候选表述为真实差分路径。

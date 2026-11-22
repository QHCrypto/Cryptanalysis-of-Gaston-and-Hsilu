# Gaston 核外输出 7、8、9 个活跃列：三个下界的反例证书

计算日期：2026-10-20。

**三个待验证下界全部不成立**。这里给出的结论是由显式反例严格证明的，不依赖求解器关于不可行性或最优性的判断。

| 输出恰好活跃列数 | 待验证的输入下界 | 已验证反例的输入活跃列数 | θ-effect | 找到反例的 MILP 求解时间 |
|---:|---:|---:|---|---:|
| 7 | ≥12 | 11 | `0x0000000000000801` | 93.7279 秒 |
| 8 | ≥11 | 10 | `0x0000000000000801` | 9.73786 秒 |
| 9 | ≥10 | 7 | `0x0000000000000801` | 1.86211 秒 |

这些是并行运行中对应成功实例报告的求解时间，不含建模、其他探索实例和证书复核时间，也不代表在其他机器上的保证速度。

**没有证明 11、10、7 是各自的最小输入活跃列数**。本次解决的是三个“是否至少为”的判定问题；反例足以逐项给出否定答案。

## 问题与传播定义

状态 `X=(X0,...,X4)`、`Y=(Y0,...,Y4)` 的每行是 64 位整数，最低位编号为 0。
第 j 列的 5 比特不全为 0 时，该列活跃；用 `wc` 表示活跃列数。

使用完整差分线性层 `Y = λ(X)`：

```text
east = (0,60,22,27,4)
west = (0,56,31,46,43)
twin = (25,32,52,60,63)
A_i = ROL64(X_i,east_i)
P = XOR_i A_i
Q = XOR_i ROL64(A_i,twin_i)
E = P XOR ROL64(P,1) XOR Q XOR ROL64(Q,18)
Y_i = ROL64(A_i XOR ROL64(E,23),west_i)
```

核外条件是 `E(A) != 0`，其中 `A=ρeast(X)`。本次求解的三个反例问题为：

```text
Y = λ(X), E(ρeast(X)) != 0,
wc(Y) = s, wc(X) <= 18-s,     s ∈ {7,8,9}.
```

三个实例均可满足。所有证书均有 `E=0x801`，其置位位置是 0、11；加入 θ 状态各行的是 `ROL64(E,23)=0x0000000400800000`。

## 可独立复核的反例

每组按第 0 至第 4 行给出，所有整数均为十六进制。

输出 7 列、输入 11 列：

```text
X = (0000000000800000, 0000005000000000, 0000000000001482,
     1000000000000480, 0000000040084000)
Y = (0000000400000000, 0000000001008000, 9000000000000000,
     0000000000080000, 2000000000000000)
input columns  = {1,7,10,12,14,19,23,30,36,38,60}
output columns = {15,19,24,34,60,61,63}
```

输出 8 列、输入 10 列：

```text
X = (0000000040800000, 0000000008000000, 0000000003001002,
     0000000002000000, 0000080004080000)
Y = (0000000440000000, 0000000004000000, 0000000000006000,
     0000000400010020, 0000000004002200)
input columns  = {1,12,19,23,24,25,26,27,30,43}
output columns = {5,9,13,14,16,26,30,34}
```

输出 9 列、输入 7 列：

```text
X = (0000000000800001, 0000000000000020, 0000000000000002,
     0000000000000084, 0000000040000000)
Y = (0000000400000001, 0200000004008000, 0000000000000002,
     0000000000000820, 0000000000000004)
input columns  = {0,1,2,5,7,23,30}
output columns = {0,1,2,5,11,15,26,34,57}
```

每组的 P、Q、θ 输入状态、每列 5 比特值、比特重量与文件哈希均见 `out{s}.counterexample.json`。

## 只复核结论，不重新搜索

运行以下命令只需要 Python 3 标准库，无需 Gurobi、许可证或 SAT 求解器：

```bash
cd /home/yxb/分析/gaston_active_pattern/quasi_branch_proof
python3 verify_certificates.py
```

验证脚本以精确整数运算进行以下检查：

1. 根据单个输入比特的目的位置构造完整 320×320 二进制矩阵，检查 `λ(X)=Y`。
2. 独立用字级循环移位计算同一传播及 P、Q、E，检查 `E!=0`。
3. 在 GF(2) 上求矩阵逆，确认秩为 320，并检查 `λ⁻¹(Y)=X`。
4. 重新统计每个端点的活跃列、列值与比特重量，逐项验证反例条件。
5. 检查三个完整 SAT 赋值分别满足原始精确判定 CNF 的全部 9541、9653、9753 条子句。
6. 检查 SAT 赋值中的 X、Y、E 与显式反例一致，以及证书相关文件的 SHA-256。

逆矩阵仅用于额外的交叉检查。反例有效性由直接前向传播和列数检查即可证明。
机器可读的复核结果保存为 `verification.json`。

## 搜索方法与材料

先尝试低重量 E 的受限反例搜索。`E=1+x^11` 取自此前输出 6 列核外最优见证的 θ-effect 经同步循环移位归一化后的形状。
在此 E 上，三个精确 MILP 实例分别由 4 个线程求解，并行运行。
同时试算另一种重量 2 的 E 及不固定 E 的 SAT 实例；找到三个反例后终止剩余探索。

限制 `E=0x801` 只能用来寻找反例；如果受限问题不可行，不能据此证明一般核外问题不可行。
本次找到的均是可行反例，所以这一限制不影响否定结论的严格性。

MILP 以现有程序导出的 `base.lp` 为基础，通过 C++ Gurobi API 加入输入上界及指定 E，目标设为常数 0。
异或使用原程序的二元/三元异或凸包及三元树分解，没有整数松弛变量。
求解器状态 `OPTIMAL` 在这里仅指常数目标的可行性问题求解完成，不能解释为输入列数最优。

主要文件：

| 文件 | 用途 |
|---|---|
| `results.json` | 三个问题的汇总与完整反例 |
| `out7.counterexample.json`、`out8.counterexample.json`、`out9.counterexample.json` | 显式状态证书 |
| `verify_certificates.py`、`verification.json` | 不依赖求解器的精确验证及结果 |
| `unrestricted_out{s}.cnf`、`.map.json` | 不固定 E 的精确 SAT 判定模型及变量映射 |
| `out{s}.sat.model` | 满足原始判定模型的完整 SAT 赋值 |
| `e2_gap11_out{s}.lp`、`.sol`、`.log`、`.status.json` | 成功搜索的 MILP 模型、完整解和日志 |
| `base.lp` | 未加输入上界、未固定 E 的基础 MILP |
| `check_bound.cpp`、`Makefile` | 判定问题的 C++ 驱动程序 |
| `model_source_snapshot.cpp` | 本次基础 MILP 生成器的源代码快照 |
| `make_cnf.py` | 独立 SAT 编码器 |
| `export_certificates.py` | 提取 MILP 见证并补全独立 CNF 赋值 |
| `environment.json`、`SHA256SUMS` | 软件版本、命令、文件完整性信息 |
| `exploration/` | 未完成/不再需要的探索日志及模型，不作为下界证明 |

由于三个问题都是 SAT，严格证书是具体可满足赋值。中断的 DRAT 搜索轨迹没有作为不可行性证明保留。

## 重建与重复搜索

本机使用 Gurobi 13.0.2；默认路径 `/opt/gurobi1302/linux64`，可通过 `GUROBI_HOME` 覆盖。

```bash
cd /home/yxb/分析/gaston_active_pattern/quasi_branch_proof
make
# 可选：从保存的源代码重建基础模型，不覆盖原始 base.lp。
make SRC=model_source_snapshot.cpp
./model_source_snapshot --mode difference --output-active 7 \
  --theta-effect nonzero --no-symmetry --prefix base_rebuilt --export-only

# 重跑三个成功搜索；前台日志同时保存，不设置求解时间上限。
mkdir -p reproduced
set -o pipefail
./check_bound base.lp 7 11 801 4 reproduced/out7 2>&1 | tee reproduced/out7.log &
./check_bound base.lp 8 10 801 4 reproduced/out8 2>&1 | tee reproduced/out8.log &
./check_bound base.lp 9 9 801 4 reproduced/out9 2>&1 | tee reproduced/out9.log &
wait
```

原始命令和运行环境见 `environment.json`。多线程求解的时间及首先找到的见证可能变化。

`make_cnf.py` 和 `export_certificates.py` 需要 `python-sat`；本机可使用
`/home/yxb/分析/Hsilu_trail/tools/sat/venv/bin/python`。
生成器用异或真值表子句、精确 OR 活跃标志及顺序计数器，表达 `wc(X)<=18-s`、`wc(Y)=s`。
E 不为零时具有偶数比特重量，因此加入冗余约束 `2<=wt_bits(E)<=18`。
仅在不固定 E 的模型中设置 `E[0]=1` 以去除同步循环移位对称性，没有同时增加输入或输出列 0 活跃限制。

导出 SAT 证书时固定已知的 X，求解器只补全其他变量；得到的完整赋值会再次直接代入**未增加这些固定输入单位子句的原始 CNF**进行检查。

## 对核外拟分支数的含义

按照 `/home/yxb/code/sbox/double-chi/eurocrypt_paper/main.tex` 的最小值刻画，将定义域限制为上述核外状态，且 `p=wc(X), s=wc(Y)`：

```text
output quasi = min [s + max(p,s)]
input quasi  = min [p + max(p,s)]
quasi        = max(output quasi, input quasi)
```

三组反例的输出代价都是 18。因此核外输出拟分支数至多为 18，候选值 19 已被否定。
其中 `(p,s)=(7,9)` 还给出输入拟分支数至多为 16，所以不带方向的核外拟分支数也至多为 18。
**这些都是上界，尚未证明相应数值取等号**。此前输出 1 至 6 列的穷举结果不受这些反例影响。

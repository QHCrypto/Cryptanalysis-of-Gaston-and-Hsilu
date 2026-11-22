# Gaston 小活跃列数量的精确枚举

这是独立的 C++17 枚举程序，不依赖 Gurobi。适用于固定少量输入或输出活跃 S 盒列、求另一端最小活跃列数的问题。`--input-active K`、`--output-active K` 的含义与父目录 MILP 程序一致：**恰好 K 个非零的 5 比特列**，每个活跃列允许全部 31 种非零取值。

新增 [CPU/GPU 并行程序 gaston_parallel](README_PARALLEL.md)，支持 CPU 多线程、OpenCL GPU 及两者协同，进一步按位置集合的循环移位轨道去重。已完整求出差分输出 4 列和 5 列对应的最小输入列数分别为 **34、25**。本机 GTX 1080 Ti 的 GPU 搜索耗时分别为 0.067 秒、16.606 秒。下文的 `gaston_enumerate` 仍是原串行入口。

输出 6 列且要求 E≠0 的问题另有 [严格剪枝 GPU 枚举入口](README_PARALLEL.md#核外六列的严格剪枝枚举)，使用互不重叠的零输入列对反解最后两个输出列，保留完整的最优性保证。

## 编译与运行

```bash
cd /home/yxb/分析/gaston_active_pattern/direct_enumeration
make

# 与原命令对应的差分实例：输出恰好 2 个活跃 S 盒列
./gaston_enumerate --mode difference --output-active 2 --prefix result/diff_out2

# 关闭循环移位对称性约简，枚举所有候选
./gaston_enumerate --mode difference --output-active 2 --no-symmetry --prefix result/diff_out2_full
```

其余方向和掩码模式也可直接使用：

```bash
./gaston_enumerate --mode difference --input-active 2
./gaston_enumerate --mode linear --input-active 2
./gaston_enumerate --mode linear --output-active 2
```

支持 `make SRC=gaston_enumerate.cpp`、`make gaston_enumerate`。默认编译选项包含 `-march=native`，二进制按编译机器的 CPU 优化；迁移到其他机器时重新编译即可。需要通用二进制可用：

```bash
make -B CXXFLAGS='-O3 -std=c++17 -Wall -Wextra -Wpedantic'
```

`--prefix` 会保存 JSON，并自动创建父目录。未指定前缀时只显示控制台结果。程序单线程执行；启动时显示计划枚举数量，发现更优候选时显示改进，结束时显示耗时、实际遍历数、两端活跃列数、活跃比特数、活跃列及五个 64 位字。

列值的第 i 位对应状态第 i 行；列值 `0x19` 表示行 0、3、4 的比特为 1。行字的最低位对应列 0。

## 方法与完整性

首先以 320 个单位向量构造 Gaston 的二元线性矩阵 L，再按模式选择映射：

| 模式 | 固定端 | 从固定端到另一端的映射 |
|---|---|---|
| difference | input | L |
| difference | output | L⁻¹ |
| linear | input 掩码 α | L⁻ᵀ |
| linear | output 掩码 β | Lᵀ |

逆矩阵通过 GF(2) 高斯消元得到。程序对全部 320 个基向量复核原始字级传播关系，然后预计算 `table[column][value]`：一个指定列取指定非零值时，其在另一端的完整 320 比特像。

对于两个输出活跃列 j、k，列值分别为 a、b：

```text
X = table[j][a] XOR table[k][b]
input_active = popcount(X[0] OR X[1] OR X[2] OR X[3] OR X[4])
```

内层循环直接使用五个 `uint64_t`，不为每个候选重新求逆或做稠密矩阵乘法。一般 K 列通过递归增量异或表项求值，完整枚举所有位置组合和所有非零列值。

完整搜索空间为 `C(64,K) × 31^K`。Gaston 的正向、转置和逆映射均与全体行同时循环移位可交换，且循环移位保持活跃列数。因此，K>0 时可以要求固定端的第 0 列活跃，搜索空间变为 `C(63,K-1) × 31^K`。这包含每个循环移位等价类的至少一个代表，可能保留多个代表，但不会遗漏最小值。

对 K=2：

- 完整枚举：`C(64,2) × 31² = 1,937,376` 个候选。
- 固定第 0 列活跃：`63 × 31² = 60,543` 个候选。

枚举全部完成且最终状态复核通过时，结果标记 `proven_optimal=true`。按 Ctrl+C 可提前停止并保存当前最优候选，标记 `proven_optimal=false`、返回码 130。中断后尚未遍历的候选不能据此排除。程序没有默认时间限制。

K 的参数范围为 `0..64`，但计算量按上述组合数增长，实际适合小 K。启动时的 `candidate_count_estimate` 使用浮点数表示规模；`visited` 是实际遍历数量的整数计数。

## 输出两列的结果

对 `--mode difference --output-active 2`，完整枚举得到的**最小输入活跃 S 盒列数为 48**。这是完整 320 比特线性层上的固定数量最优值。

一组最优输出为：列 0 取 `0x19`，列 25 取 `0x02`，其余列为零。因此输出有两个活跃列、四个活跃比特。对应状态为：

```text
Y = (0x0000000000000001,
     0x0000000002000000,
     0x0000000000000000,
     0x0000000000000001,
     0x0000000000000001)

X = (0x2d80051301102e64,
     0xd80051101102e652,
     0x40b994b600144c04,
     0x2285cca5b000a260,
     0x52d80051301302e6)
```

X 有 48 个活跃列、100 个活跃比特，满足 `L(X)=Y`。两种枚举方式的证书分别保存在 [result/diff_out2.json](result/diff_out2.json) 和 [result/diff_out2_full.json](result/diff_out2_full.json)。

2026-10-20 本机一次运行记录：

| 枚举方式 | 候选数 | 预计算及检查 | 枚举时间 | 合计 |
|---|---:|---:|---:|---:|
| 固定第 0 列活跃 | 60,543 | 0.000603 秒 | 0.000165 秒 | 0.000768 秒 |
| 全部位置组合 | 1,937,376 | 0.000555 秒 | 0.004864 秒 | 0.005419 秒 |

这些是程序内部计时，包含矩阵构造、求逆、基向量检查、查表预计算和枚举，不包含进程启动及最终 JSON 写盘；具体耗时随机器负载变化。

## 独立验证

```bash
make check
```

只需要 Python 标准库。验证脚本从单比特目的位置独立构造二元矩阵，使用 Python 整数消元，并独立枚举 K=0、1、2 的参考最优值。它还对四种模式分别运行无对称性约简的完整 K=2 枚举，核对最优值、候选计数、状态证书及 Ctrl+C 中断行为。

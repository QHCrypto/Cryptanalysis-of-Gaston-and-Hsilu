# Gaston 活跃列精确并行枚举

`gaston_parallel.cpp` 在原枚举方法上增加 CPU/OpenCL GPU 并行与位置集合的循环移位轨道去重。支持差分和线性掩码两种模式、固定输入或固定输出，以及 K=0..5。每个活跃列的 31 种非零取值均完整遍历。

## 编译与运行

```bash
cd /home/yxb/分析/gaston_active_pattern/direct_enumeration
make gaston_parallel

# 本机实测 GPU 单独运行最快
./gaston_parallel --mode difference --output-active 4 --backend opencl --prefix result/diff_out4_gpu
./gaston_parallel --mode difference --output-active 5 --backend opencl --prefix result/diff_out5_gpu

# CPU 多线程；--threads 默认使用可用逻辑核数
./gaston_parallel --mode difference --output-active 5 --backend cpu --threads 32 --prefix result/diff_out5_cpu

# CPU/GPU 共用动态任务队列；不同硬件上速度可能不同
./gaston_parallel --mode difference --output-active 5 --backend hybrid --threads 32 --prefix result/diff_out5_hybrid
```

本机使用 GCC 13、OpenMP 与现有 NVIDIA OpenCL 驱动，无需 CUDA Toolkit、nvcc 或 Gurobi。Makefile 链接 Linux 的 `libOpenCL.so.1`。所需三个 Khronos 头文件随代码放在 `third_party/OpenCL/`，保留 Apache-2.0 许可证及来源记录。它们来自 [KhronosGroup/OpenCL-Headers](https://github.com/KhronosGroup/OpenCL-Headers)。`make SRC=gaston_parallel.cpp` 同样可用。

`--backend cpu` 为默认后端；`--backend opencl` 选择系统枚举到的第一张可用 OpenCL GPU。K=0、1 的计算量很小，直接由 CPU 完成。`--batch-positions N` 控制 GPU 每次处理多少个位置集合，默认 64；GPU 计算分批提交。CPU/GPU 协同模式以原子队列分配互不重叠的任务。

原模式参数保持一致：

```bash
./gaston_parallel --mode difference --input-active 4 --backend opencl
./gaston_parallel --mode linear --input-active 4 --backend opencl
./gaston_parallel --mode linear --output-active 4 --backend opencl
```

`--no-symmetry` 枚举全部位置组合，适合对照验证。`--limit-positions N` 仅用于性能采样：若没有遍历全部位置集合，结果明确标记 `proven_optimal=false`，返回码 130。按 Ctrl+C 时完成正在计算的批次或位置集合后，保存已经完成部分的最优候选，同样不宣称已证明最优。程序没有默认时间限制。

## 精确性与去重方式

每个位置集合用一个 64 比特掩码表示。对于非空集合，在所有循环移位中保留掩码数值最小的一个代表；这个代表必然包含第 0 列。对代表上的各个活跃列，仍然枚举全部 `31^K` 个非零列值组合。

正向映射 L、逆映射 L⁻¹、转置 Lᵀ、逆转置 L⁻ᵀ 都与整个状态的同步循环移位可交换，且同步移位不改变任一端的活跃列数。程序逐个检查全部列值映像的这一关系。因此，任意原候选都可移位到某个保留的位置集合上，对应的非零列值组合也在搜索范围内，目标值不变。

周期性位置集合仍然完整保留所有列值赋值。这里仅对位置集合去重，没有额外对列值赋值去重，因而不会因为位置集合的旋转稳定子而漏掉候选。

位置代表数使用 Burnside 公式独立复核。对移位 s，令 `g=gcd(64,s)`，则存在 g 个长度为 `64/g` 的循环。只有 `64/g` 整除 K 时才有不变的 K 元位置集合，数量为 `C(g,K/(64/g))`。对 s=0..63 求和再除以 64，即得位置轨道数。

| K | 位置轨道数 | 每个位置的列值组合数 | 实际遍历候选数 |
|---:|---:|---:|---:|
| 2 | 32 | 961 | 30,752 |
| 3 | 651 | 29,791 | 19,393,941 |
| 4 | 9,936 | 923,521 | 9,176,104,656 |
| 5 | 119,133 | 28,629,151 | 3,410,676,646,083 |

CPU 后端将最后两个列值的像预合成 961 项缓存，再枚举其余列值。GPU 后端把每个位置集合分为 8 个工作组、每组 128 个线程，以互不重叠的步长划分前 K−2 个列值的全部赋值；每个线程完整枚举最后两个列值。每组输出最小值及对应的 31 进制赋值编码，主机归并各组结果。没有启发式剪枝或提前停止规则。

`gaston_common.hpp` 由原串行程序和并行程序共同使用，提供线性层定义、GF(2) 矩阵、字级验证及 JSON 输出。

## 2026-10-20 完整搜索结果

硬件：Intel i9-13900KF，NVIDIA GeForce GTX 1080 Ti。以下为差分模式、恰好 K 个输出活跃列、最小化输入活跃列数。

| K | 最小输入活跃列数 | GPU 搜索时间 | 预计算、GPU 初始化及搜索合计 | 状态 |
|---:|---:|---:|---:|---|
| 4 | **34** | 0.067126 秒 | 0.192454 秒 | 完整枚举 |
| 5 | **25** | 16.605703 秒 | 16.669608 秒 | 完整枚举 |

上述计时不包含进程启动、最后打印及最终 JSON/TSV 写盘；GPU 初始化时间受运行时编译缓存影响，具体耗时依机器负载而变化。

最优输出状态可用如下列赋值表示，其余列为零。列值第 i 位对应第 i 行，列编号从 0 开始。

- K=4：`{0:0x1f, 12:0x05, 33:0x0f, 36:0x12}`。对应输入有 34 个活跃列、53 个活跃比特。
- K=5：`{0:0x14, 11:0x04, 12:0x13, 23:0x0f, 41:0x13}`。对应输入有 25 个活跃列、37 个活跃比特。

完整 320 比特状态、计数和耗时见 [4 列结果](result/diff_out4_gpu.json)、[5 列结果](result/diff_out5_gpu.json)。搜索日志见 [4 列日志](result/diff_out4_gpu.log)、[5 列日志](result/diff_out5_gpu.log)。

每次指定 `--prefix` 还会保存 `PREFIX.positions.tsv`，逐项记录位置代表的掩码、该位置的最小值及对应的列值赋值编码。编码的第 d 个 31 进制数位加 1，即按列号升序排列的第 d 个活跃列的非零值。保存的是每个已完成位置集合的最小值，方便独立验证与 CPU/GPU 对照。

这两个最优值适用于本程序所定义的完整 Gaston 线性层及活跃列度量，不直接给出多轮构造的结论。

## 实测后端比较

对相同的前 512 个五列位置代表，三种后端均遍历 14,658,125,312 个候选；每个位置的最小值和对应赋值编码完全一致。

| 后端 | 搜索时间 |
|---|---:|
| OpenCL GPU | 0.088463 秒 |
| 32 线程 CPU | 0.839519 秒 |
| GPU + 32 线程 CPU | 0.114658 秒 |

因此此次完整搜索使用 GPU 后端。协同计算在本机该样本上受到 CPU 收尾等开销影响，未比 GPU 单独运行更快；其他设备可用 `--limit-positions` 重新比较。

## 验证

```bash
make check
make check-parallel
python3 verify_parallel.py --results result/diff_out4_gpu result/diff_out5_gpu
```

验证包含：

- 原串行程序的独立 Python 矩阵与小数量参考枚举全部通过。
- 并行程序在四种传播方向的 K=0..3 实例中，CPU/GPU 逐位置最小值及赋值编码完全一致；K=0..2 还与 Python 独立枚举值比较。
- 四种方向均检查关闭对称约简的 K=2 完整搜索，覆盖周期性位置集合。
- 4 列额外使用 32 线程 CPU 完整搜索，全部 9,936 个位置代表的最小值及赋值编码与 GPU 完全一致。CPU 搜索耗时 0.710534 秒。
- 5 列采样的 512 个位置代表由 CPU、GPU、协同三种方式交叉检查；正式结果的全部 119,133 个逐位置证书由独立 Python 矩阵验证传播与活跃数量。
- 对正式结果逐个检查代表确实为其轨道的最小掩码、代表互不重复、数量等于 Burnside 轨道总数，确认位置轨道覆盖完整。
- 有意截取的部分搜索及 Ctrl+C 中断均不标记为已证明最优。

逐位置状态证书本身验证该目标值可达到；每个位置不存在更小值的结论来自相应后端对全部 `31^K` 个赋值的完整枚举。

## 核外六列的严格剪枝枚举

新增 [gaston_out6_nonzero.cpp](gaston_out6_nonzero.cpp)，专门求解差分模式中输出恰好 6 个活跃列、`E(rho_east(X)) != 0` 时的最小输入活跃列数。这里的核外条件与父目录 MILP 的 `--theta-effect nonzero` 相同。该程序使用组合枚举、GF(2) 线性方程反解和 OpenCL，不依赖 Gurobi。

```bash
make gaston_out6_nonzero
./gaston_out6_nonzero --prefix result/diff_out6_nonzero

# 验证已有计算结果，并生成 1～6 列核外结果汇总
python3 verify_nonzero_results.py
```

6 列共有 1,171,552 个位置旋转轨道，每个位置代表上枚举 `31^6 = 887,503,681` 种非零列值赋值。直接遍历这些代表上的全部赋值，需要处理 1,039,756,712,482,912 个候选。新程序使用如下完备性论证减少实际计算：

1. 设要排除输入活跃列数小于 B 的状态。取 B 对互不重叠的输入列；少于 B 个活跃列最多使 B−1 对包含活跃列，因此至少有一对的两个输入列均为零。
2. 对每个六列输出位置代表，完整枚举前四列的 `31^4` 种非零赋值，保留最后两列的 10 个比特为未知量。
3. 令选定的一对输入列为零，得到关于这 10 个未知量的 10 个二元线性方程。预先选择系数矩阵满秩的列对，每对方程唯一确定最后两列的值。若其中一个输出列为零则丢弃，否则计算完整输入活跃列数并检查 E≠0。
4. 对每个前缀枚举全部 B 个选定列对。任何输入少于 B 个活跃列的解都至少会被其中一对恢复，因此不会遗漏可改善当前最优值的状态。同一个状态可能被恢复多次，不影响最小值。

程序对全部 63 种非零输出列间距都构造了 26 对互不重叠、系数矩阵满秩的输入列对，并逐一验证每个反解映射的全部 1024 种右端。10 比特线性反解表压缩为两个 32 项查表的异或，存入 GPU 共享内存。五行同步移位把最后两个输出列归一化到 `0` 和 `gap`，减少预计算表的规模。

初始阈值默认 `--cutoff 26`，即搜索输入活跃列数小于 26 的核外状态。找到更小的可行解后，后续位置自动使用更小阈值。**只有遍历全部位置轨道并且找到低于初始阈值的可行解，才报告 `proven_optimal=true`**。如果完整搜索未找到阈值以下的解，只能得出输入活跃列数至少为该阈值，程序不会把阈值当成可达到的最优值。

由于 `east[0]=west[0]=0`，差分传播恒满足 `X[0] XOR Y[0] = ROL(E,23)`，因此 GPU 用 `X[0] != Y[0]` 精确筛选 E≠0。主机对返回的每个候选还会独立重算完整传播和 E，核对活跃列数。

`PREFIX.positions.tsv` 包含所有已完成位置的掩码、本次搜索阈值、阈值以下的最小值及列值编码；没有阈值以下解时记录 `null`。`tested_candidates_with_duplicates` 是实际计算的候选次数，允许重复；`covered_output_states` 表示已完成位置上的完整赋值空间大小，其中部分赋值通过上述严格论证排除，不表示逐一计算了全部赋值。

可用以下参数做有限范围的交叉检查：

```bash
# 同一个位置：GPU 严格剪枝、CPU 严格剪枝、CPU 直接遍历 31^6 个赋值
./gaston_out6_nonzero --start-position 321320 --limit-positions 1 --prefix result/check6_gpu
./gaston_out6_nonzero --backend cpu --start-position 321320 --limit-positions 1 --prefix result/check6_cpu
./gaston_out6_nonzero --backend brute --start-position 321320 --limit-positions 1 --prefix result/check6_brute

python3 verify_nonzero_results.py --six result/check6_gpu --partial-only
```

部分搜索和 Ctrl+C 中断均不会标记全局最优。`--backend brute` 用于少量位置的对照；完整问题默认使用 GPU 剪枝后端。`--cutoff` 范围为 1..26，`--batch-positions` 默认 64。

验证脚本从独立的比特目的位置矩阵出发，检查 63 种间距下的满秩列对、位置轨道的规范性与完整覆盖、所有已记录候选的传播及核外条件，以及最终 JSON。对于 1～5 列，已有完整的无核限制枚举给出下界，而其最优状态均满足 E≠0，因此同一数值也是核外限制下的严格最优值。

### 2026-10-20 核外结果

以下均为完整 Gaston 线性层的差分传播，限定 `E(rho_east(X)) != 0`，输出恰好 K 个活跃 S 盒列。

| 输出活跃列数 K | 最小输入活跃列数 | 最优性依据 |
|---:|---:|---|
| 1 | **57** | 完整无核限制枚举的最优状态满足 E≠0 |
| 2 | **48** | 同上 |
| 3 | **39** | 同上 |
| 4 | **34** | 同上 |
| 5 | **25** | 同上 |
| 6 | **13** | 核外严格剪枝枚举覆盖全部位置轨道 |

6 列正式搜索覆盖全部 **1,171,552** 个位置轨道，实际计算 **17,040,104,227,205** 次候选（包含重复恢复的状态）。GTX 1080 Ti 的搜索耗时为 **487.308458 秒**，预计算及初始化耗时为 0.585640 秒。位置编号 0、3642、321320、500000、934719、1171551 分别通过 CPU 直接遍历全部 `31^6` 种赋值，与 CPU 剪枝及 GPU 记录交叉核对一致；还检查了阈值的严格不等号边界和 SIGINT 中断标记。

一组 6 列最优输出的非零列值为：

```text
{0:0x05, 12:0x10, 23:0x02, 30:0x01, 35:0x08, 41:0x08}
```

对应输入有 13 个活跃列，`E = 0x0000020040000000 != 0`。完整状态和结果保存在：

- [1～6 列核外结果汇总](result/nonzero_output_1_to_6.json)
- [6 列结果与状态证书](result/diff_out6_nonzero.json)
- [6 列搜索日志](result/diff_out6_nonzero.log)
- [全部位置的证明记录](result/diff_out6_nonzero.positions.tsv)
- [直接枚举交叉检查](result/nonzero6_crosschecks.json)
- [独立验证日志](result/nonzero_output_1_to_6.validation.log)

这些数值是单个完整线性层在上述核外条件与活跃列度量下的最优值。

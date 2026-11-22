# Gaston 核外普通分支数 B=16：证明与复核材料

**已完成：核外普通分支数为 16**。全部 1667 个精确实例均不可满足，其 DRAT 证明已逐个通过独立检查。汇总验证见 `proof_verification.json`，固定结果摘要见 `result.json`。

## 定义与目标

用 `wc` 表示 5×64 状态中非零的 5 比特列数。研究对象为完整差分线性层
`λ = ρwest ∘ θ ∘ ρeast`，核外普通分支数定义为

\[
B_{\neg\ker E}=
\min_{\substack{X\in\mathbb F_2^{320}\\E(\rho_{\rm east}(X))\ne0}}
\bigl(\operatorname{wc}(X)+\operatorname{wc}(\lambda(X))\bigr).
\]

本证明通过总列数为 16 的显式见证，以及对总列数至多 15 的状态的完全排除，确定其精确值。
下界证明不依赖此前的 \(13\le B\) 或各端点小重量 MILP 结果。

状态的第 j 比特对应整数的 \(2^j\)，各行下标为 0 至 4。使用参数：

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

## 上界：三个输入活跃列给出十三个输出活跃列

`upper_bound16.json` 保存下面的见证，五行依次为：

```text
X = (0000000000000000, 0000008000000100, 0000000000000000,
     0000000000000000, 0000000000000001)
Y = (0400000000002000, 1004000008000020, 0000100002000000,
     0800010000000000, 0100802000000000)
E = 0040000800000000 != 0
```

输入活跃列为 `{0,8,39}`，输出活跃列为
`{5,13,25,27,37,40,44,47,50,56,58,59,60}`。
逐比特传播与独立字级传播均可检查 \(Y=\lambda(X)\)，所以 \(B_{\neg\ker E}\le16\)。
这一见证也便于手工核对：\(X=(0,2^8+2^{39},0,0,1)\)，
故 \(P=2^{35}\)、\(Q=2^{36}\)，从而 \(E=2^{35}+2^{54}\ne0\)。

## 下界第一步：θ-effect 导出的匹配下界

假设存在核外状态满足
\(\operatorname{wc}(X)+\operatorname{wc}(Y)\le15\)。
令 \(I,O\) 分别为 X、Y 的活跃列集合，\(S=\operatorname{supp}(E)\)。

对任意 \(j\in S\) 和任意行 i，由传播定义可得

\[
 X_{i,j+23-\mathrm{east}_i}\oplus
 Y_{i,j+23+\mathrm{west}_i}=1,
\]

其中列下标均模 64。因此，这两个端点至少有一个是活跃列。

建立二分图 \(G_S\)：左侧 64 个顶点表示输入列，右侧 64 个顶点表示输出列；
每个 \(j\in S\)、每行 i 加入边

\[
\bigl(j+23-\mathrm{east}_i,\ j+23+\mathrm{west}_i\bigr).
\]

\(I\cup O\) 必须覆盖所有边。若图中存在 16 条两两没有公共端点的边，则覆盖它们至少需要 16 个活跃列，与总数至多 15 矛盾。

这里仅使用“一个顶点至多覆盖一条匹配边”的直接论证，不依赖最大匹配等于最小点覆盖的算法正确性。
最终验证器直接检查每份证书中的 16 条边确实存在且端点互异，不信任生成器给出的“最大匹配”数值。

## 下界第二步：覆盖所有 θ-effect

对 E 的 64 比特求异或：

\[
\bigoplus_jE_j
=\bigoplus_jP_j\oplus\bigoplus_jP_j
 \oplus\bigoplus_jQ_j\oplus\bigoplus_jQ_j=0.
\]

因此，非零 E 的比特重量是正偶数，只需处理重量 2、4、6 和至少 8。
同步循环移位保持传播、核外条件及两端列数，所以按 E 的循环移位轨道分类。

| E 的比特重量 | 全部轨道数 | 用 16 条匹配边直接排除 | 剩余精确判定实例 |
|---:|---:|---:|---:|
| 2 | 32 | 0 | 32 |
| 4 | 9936 | 8304 | 1632 |
| 6 | 1,171,552 | 1,171,549 | 3 |

`pruned_matchings.bin` 保存每一个被排除轨道的代表及 16 条匹配边，共 **1,179,853 份证书**。
`effects_total15.tsv` 保存剩余的 1667 个代表。
验证器同时检查：

1. 每个代表的重量正确，且为其轨道的最小整数代表。
2. 所有代表互不重复。
3. 每份匹配证书的 16 条边都存在，且左右端点分别互异。
4. 剩余代表与被排除代表的总数，等于 Burnside 公式给出的全部轨道数：

\[
N_k=\frac1{64}
\sum_{\substack{0\le r<64\\64/\gcd(64,r)\mid k}}
\binom{\gcd(64,r)}{k/(64/\gcd(64,r))}.
\]

“代表合法、互不重复、数量等于全部轨道数”保证没有漏掉任何轨道。

重量 6 的三个剩余代表为：

```text
0x0000028004040101  support={0,8,18,26,39,41}
0x0002020004040101  support={0,8,18,26,41,49}
0x0000401010040401  support={0,10,18,28,36,46}
```

对重量至少 8 的任意支持 S，先取其中任意 6 元子集 T。
如果 T 属于已排除的重量 6 轨道，其 16 条匹配边也存在于 \(G_S\)，已经得到矛盾。
否则，将 T 同步旋转到上述三个代表之一，再从 \(S\setminus T\) 取两个位置。
每个代表有 \(\binom{58}{2}=1653\) 种这样的扩展。
`extension_matchings.bin` 对全部 \(3\times1653=4959\) 种扩展都给出 16 条匹配边。
验证器逐份检查其有效性、互异性和覆盖完整性。
由于扩展图是 \(G_S\) 的子图，这就严格排除了**所有重量至少 8 的 E**。

这一步是显式组合证书，不需要依赖另一个 MILP/SAT 求解器关于高重量 E 的不可行性判断。

## 下界第三步：剩余 1667 个精确实例

对每个剩余 E，判定：

\[
Y=\lambda(X),\qquad E(\rho_{\rm east}(X))=E,\qquad
\operatorname{wc}(X)+\operatorname{wc}(Y)\le15.
\]

`total15.cnf` 是未固定 E 的基础模型。每个实例仅额外增加 64 个单位子句，指定对应的 E。
没有增加输入/输出列数下界，没有限制输入每列的非零值，也没有使用线性掩码传播。

编码中：

- 所有状态位、辅助位和活跃列标志都是布尔变量。
- 异或由禁止奇校验赋值的完整 CNF 编码表示，长异或使用三元树分解。
- 活跃列标志与五个状态位的 OR 精确等价。
- 顺序计数器表达两端活跃列总数至多 15。
- `E[0]=1` 只选择同步循环移位代表；所有保留 E 代表的第 0 位都是 1。
- `2<=wt_bits(E)<=15` 是有效冗余条件。偶数性给出非零时至少为 2；对任意一行，\(E\) 的重量至多为该行输入、输出的比特重量之和，进而至多为两端活跃列数之和。

每个实例使用 Kissat 产生二进制 DRAT 不可满足证明，随后由独立程序 `drat-trim` 检查。
只有求解器退出码为 20、检查器退出码为 0 且输出 `s VERIFIED` 的实例，才标记为完成。
每个实例均保存压缩 CNF、压缩证明、求解日志、检查日志、状态与 SHA-256。

全部 1667 个实例现已通过检查，结合前两步，排除了所有非零 E 下总活跃列数至多 15 的状态。
再结合上界见证，得到
\[
\boxed{B_{\neg\ker E}=16}.
\]

## 如何复核

进入本目录：

```bash
cd /home/yxb/分析/gaston_active_pattern/ordinary_branch_proof
```

快速复核：重新检查上界、所有匹配与扩展证书、轨道覆盖、每个精确 CNF 与基础模型的关系，以及已经检查过的 DRAT 文件的完整性与检查记录：

```bash
python3 verify_proof.py --jobs 8
```

**快速复核不会重新执行全部 DRAT 推理检查**。要从保存的原始证明重新检查所有推理步骤，运行：

```bash
python3 verify_proof.py --full --jobs 8
```

上述两条命令均不需要 Gurobi、许可证、python-sat 或重新进行 SAT 搜索。
完整复核需要本目录 `tools/drat-trim/drat-trim` 检查器和系统 gzip；检查器源代码一并保存，可用以下命令重建：

```bash
make -C tools/drat-trim
```

`SHA256SUMS` 保存交付文件的哈希。为允许再次运行验证器，它不包含可更新的 `proof_verification.json`、重放日志、文件锁和工具仓库的 Git 内部文件。
`environment.json` 保存求解器、检查器版本、源代码与二进制哈希。中断的探索记录放在 `exploration/`，不参与最终证明。

也可以分别检查：

```bash
python3 verify_matchings.py
python3 validate_encoding.py
```

`validate_encoding.py` 使用独立的比特目的位置实现，检查全部 320 个基向量、零状态和 128 个确定性随机状态，验证基础 CNF 的传播、活跃标志和覆盖约束。
基础 CNF 还由保存的生成器重新生成并进行了逐字节比较。

## 重现搜索

`make_cnf.py` 使用 python-sat 生成顺序计数器。本机 Python 环境为：
`/home/yxb/分析/Hsilu_trail/tools/sat/venv/bin/python`。

```bash
make SRC=enumerate_effects.cpp
./enumerate_effects effects_total15.tsv pruned_matchings.bin
make SRC=extend_effects.cpp
./extend_effects

/home/yxb/分析/Hsilu_trail/tools/sat/venv/bin/python \
  make_cnf.py --total 15 --prefix total15
python3 solve_effects.py --jobs 8
python3 verify_proof.py --jobs 8
```

`solve_effects.py` 支持复用已经通过检查的实例，并使用文件锁避免多个运行同时写同一个证明文件。
未通过或中断的实例不能计入下界证明。

目录中的 `ordinary_milp.cpp` 和 `effect_cover.cpp` 是探索阶段的 C++ Gurobi 程序，配有 Makefile，异或来自原程序的凸包模型；最终证明不依赖这两个程序的浮点分支定界结果。

## 结论范围

这里证明的是上述固定 Gaston 线性层、5 比特列重量、θ-effect 核外差分状态上的**普通分支数**。
它不是拟分支数，不是线性掩码分支数，也不是多轮非线性传播的活跃 S 盒下界。

# Gaston 单个核外线性层两侧 χ 的差分重量：精确最小值 32

**结论：最小值恰好为 32，34 不是有效下界**。

这里采用用户确认的边界：只考虑一个线性层和两侧的 S 盒，

\[
A\xrightarrow{\chi}X\xrightarrow{\lambda}Y\xrightarrow{\chi}B,
\qquad E(\rho_{\rm east}(X))\ne0.
\]

两个 χ 层的转移都须满足真实 DDT。目标是

\[
W=w_r(A,X)+w_r(Y,B),
\qquad
w_r(U,V)=\sum_{j=0}^{63}
-\log_2\left(\frac{\mathrm{DDT}_\chi[U_j,V_j]}{32}\right).
\]

端点 A、B 自由选择。该重量是论文的 restriction weight，即
\(-\log_2(\mathrm{EDP})\)。它不等于四个状态的活跃列数之和，也不声称
两轮真实路径概率一定等于各轮概率的乘积。

## 下界证明

记 \(c(U)\) 为状态 U 的非零 5 比特列数。

1. χ 是 5 比特置换。因此有效 χ 差分转移的输入列非零，当且仅当输出列非零。
2. 对任意非零输入差分 \(a\)，χ 的 DDT 满足
   \(\mathrm{DDT}_\chi[a,b]\le8\)。故每个活跃列的差分重量至少为
   \(-\log_2(8/32)=2\)。
3. 已完成的核外普通分支数证明给出
   \[
   c(X)+c(\lambda(X))\ge16
   \quad\text{whenever }E(\rho_{\rm east}(X))\ne0.
   \]

所以，对任意满足条件的路径，

\[
\begin{aligned}
W
&=w_r(A,X)+w_r(Y,B)\\
&\ge2c(X)+2c(Y)\\
&\ge2\cdot16=\boxed{32}.
\end{aligned}
\]

第一层可按输出 X 的活跃列计数，是因为 χ 置换性保证
\(c(A)=c(X)\)。这一点不要求第一层的所有输出差分具有相同的最小逆向重量。

所依赖的分支数定理、组合匹配证书和 1667 个已独立检查的 DRAT 证明在
[ordinary_branch_proof](../ordinary_branch_proof/README.md)。
本证明没有再引入新的 MILP/SAT 最优性假设。

## 达到 32 的显式路径

行编号从 0 到 4，列 j 为 64 位字的第 j 个最低有效位。取下表中的 X、Y：

| 行 | X | Y |
|---|---|---|
| 0 | `0x0000000000000000` | `0x0400000000002000` |
| 1 | `0x0000008000000100` | `0x1004000008000020` |
| 2 | `0x0000000000000000` | `0x0000100002000000` |
| 3 | `0x0000000000000000` | `0x0800010000000000` |
| 4 | `0x0000000000000001` | `0x0100802000000000` |

逐位及按 64 位字的两种实现均验证了

\[
Y=\lambda(X),\qquad
E(\rho_{\rm east}(X))=\mathtt{0x0040000800000000}\ne0.
\]

X 的活跃列为 \(\{0,8,39\}\)，共 3 列；
Y 的活跃列为
\(\{5,13,25,27,37,40,44,47,50,56,58,59,60\}\)，共 13 列。
两端的每个活跃列都恰好只有一个非零比特。

令 \(e_k\) 为 5 比特单位向量。对所有 k，
\[
\mathrm{DDT}_\chi[e_k,e_k]=8.
\]
因此可以取 **A=X、B=Y**，得到有效差分路径
\[
X\xrightarrow[\;6\;]{\chi}
X\xrightarrow{\lambda}
Y\xrightarrow[\;26\;]{\chi}Y.
\]
这里表示 χ 的输入差分与输出差分相同，不表示具体状态被 χ 固定。

总重量为
\[
W=3\cdot2+13\cdot2=6+26=\boxed{32}.
\]
下界与构造相等，故精确最小值为 32。其期望差分概率为 \(2^{-32}\)。

单位差分的 DDT 条目还可直接手算。采用
\[
\chi_i(z)=z_i+(1+z_{i+1})z_{i+2},
\]
下标模 5。输入差分为 \(e_k\) 时，导数的三个可能非零分量分别为
\[
D_{e_k}\chi_k(z)=1,\quad
D_{e_k}\chi_{k-1}(z)=z_{k+1},\quad
D_{e_k}\chi_{k-2}(z)=1+z_{k-1}.
\]
输出差分等于 \(e_k\) 恰当且仅当
\(z_{k+1}=0,\ z_{k-1}=1\)。其余 3 比特自由，正好有 8 个基点。

## 复核

仅需 Python 3 标准库，无需 Gurobi、PySAT 或新的求解：

```bash
cd /home/yxb/分析/gaston_active_pattern/two_round_outside_proof
python3 verify.py
```

脚本会完整生成 χ 的 DDT（32 个差分 × 32 个基点），检查：

- χ 的置换性、最大非平凡 DDT 值 8、零/非零活动保持；
- 所有 5 种单位差分的恒等差分转移权重为 2；
- 给定 X、Y 的线性传播、核外条件、活跃列数；
- 两侧所有 128 个局部差分转移及其重量，得出总重量 32；
- 已证明分支数定理的结果摘要 SHA-256。

生成的 [result.json](result.json) 还逐列列出每个活跃 χ 转移的
8 个有效基点；[chi_ddt.json](chi_ddt.json) 保存完整 DDT。

如需同时检查旧分支数证明的全部组合证书、CNF、证明文件哈希和已保存的
DRAT 检查结果：

```bash
python3 verify.py --check-branch --jobs 8
```

如需重新运行独立 DRAT 检查器，检查全部 1667 个已保存证明：

```bash
python3 verify.py --full-branch --jobs 8
```

默认命令只引用已认证的分支数定理并验证摘要哈希；它不会重新检查
全部旧证明。`--check-branch` 也不会重放 DRAT；
`--full-branch` 才会重放，且无需重新进行 SAT 搜索。

为完整对接官方轮函数 \(\chi\circ\iota\circ\lambda\)，脚本额外计算
\(Z=\lambda^{-1}(X)\)，得到
\[
Z\xrightarrow{\lambda}X\xrightarrow{\chi}X
\xrightarrow{\lambda}Y\xrightarrow{\chi}Y.
\]
这个额外的前导线性层不参与用户要求的核心边界；它的 θ-effect 也非零，
为 `0x743b9ca3a0eb61a6`。因此这份见证亦可扩展为官方轮边界的两轮路径。
轮常数在差分中抵消，不改变上述局部 DDT 重量。

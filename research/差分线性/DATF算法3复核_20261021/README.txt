原论文 DATF Algorithm 3 复算，2026-10-21

先看 datf_algorithm3_report.html；精确数值、状态与全部尝试见 results.json。
所有输出数值为相关度 rho；论文 Bias 输出 epsilon=rho/2。
算法 3 消去全部 beta，不消去 alpha 定义；剩余变量按均匀独立求和。
完成仅表示这个模型下求和完成，不自动证明真实置换相关度。
resource_limit 的 correlation/bias 必须为 null，绝不能替换为 0。

复现（先进入本目录）：
python3 -m unittest -v test_algorithm3.py
python3 build_manifest.py
python3 worker.py --manifest manifest.json --index 1 --output rerun/G3a.json --seconds 180 --max-terms 5000000 --max-products 1000000000 --max-nodes 2000000
python3 worker.py --manifest manifest.json --index 13 --output rerun/K06_prefix3_for_R5.json

worker.py --input-shannon 是同一均匀形式变量和的精确求值顺序调整，非增强 DATF 模型。
令 F=ExpandBeta(e)，有 E[(-1)^F]=2^-k sum_a E[(-1)^ExpandBeta(e|X=a)]。
Alpha 仍为自由均匀形式变量，未按真实置换条件分布重新建模。

统计 product_pairs 是累计申请展开的单项式乘积对数，最后一次超预算请求尚未执行。
每个进程地址空间上限 4 GiB。运行时间、项数、节点数上限均随尝试归档。
旧报告、HNW26、采样均为历史数据。本次只重算指定保留端点，不重做全部搜索。

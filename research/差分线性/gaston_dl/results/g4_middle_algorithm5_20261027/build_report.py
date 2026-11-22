"""Render the middle-only Algorithm 3/5 comparison from archived exact values."""
import html
import json
from fractions import Fraction
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote, unquote

H=Path(__file__).resolve().parent
ROOT=H.parents[2]
data=json.loads((H/'results.json').read_text())
validation=json.loads((H/'validation.json').read_text())
assert validation['status']=='PASS'
def esc(x):return html.escape(str(x))
def number(x):
    q=Fraction(int(x['numerator']),int(x['denominator']))
    if not q:return '0'
    sign='−' if q<0 else '+';q=abs(q)
    if q==1:return sign+'1'
    return sign+str(q.numerator)+'/'+str(q.denominator)
def link(p,label):return f'<a href="{quote(str(p),safe="/")}">{esc(label)}</a>'
def coords(ws):return '{'+', '.join(f'({i},{j})' for i,w in enumerate(ws) for j in range(64) if int(w,16)>>j&1)+'}'
tr=[];details=[]
for r in data['rows']:
    code=r['code'];a=r['algorithm3'];b=r['algorithm5'];s=r['sampling']
    assert a['correlation']==b['correlation']
    tr.append(f'''<tr><td><a href="#{code}">{code}</a></td><td>{'+'.join(map(str,r['split']))}</td>
      <td>{number(r['p'])}</td><td>{number(r['q_squared'])}</td>
      <td><b>{number(a['correlation'])}</b><small>{a['elapsed_seconds']*1000:.3f} ms</small></td>
      <td><b>{number(b['correlation'])}</b><small>{b['elapsed_seconds']*1000:.3f} ms</small></td>
      <td>{number(b['splice_single_trail'])}<small>2<sup>{b['splice_single_trail']['log2_abs']:.0f}</sup></small></td>
      <td>{s['rho']:+.12f}<small>2<sup>26</sup> 对；已有采样</small></td></tr>''')
    stats=a['statistics'];bs=r['fixed_linear_trail_all_first_DDT_branches']['branches']
    branchtext='；'.join(f"χ 输出列差分 {x['first_output_pattern']}：c<sub>m</sub>={number(x['middle_crosscheck'])}" for x in bs)
    details.append(f'''<section id="{code}"><h2>{code}：原始中间边界</h2>
      <p>整体输入差分：<code>{esc(coords(r['candidate']['input_difference']))}</code><br>
      整体输出掩码：<code>{esc(coords(r['candidate']['output_mask']))}</code></p>
      <div class="endpoints"><div><b>中间段输入差分</b><pre>{chr(10).join(r['middle_job']['delta_core'])}</pre></div>
      <div><b>中间段输出掩码</b><pre>{chr(10).join(r['middle_job']['lambda_core'])}</pre></div></div>
      <p>中间段含 {r['middle_job']['rounds']} 个 χ 层，轮常数起始索引为 {r['middle_job']['start_round_index']}。
      输入差分 {r['middle_bit_weights']['input']} 比特，输出掩码 {r['middle_bit_weights']['output']} 比特。
      完整置换仍省略首个线性层；差分前段后的线性层已计入中间输入差分，线性尾段前的线性层由转置映射计入中间输出掩码。</p>
      <p>算法 3 完全展开 β 后：{stats['expanded_terms']} 个单项式、{stats['expanded_variables']} 个变量、次数 {stats['expanded_degree']}。
      建立形式状态并核验 16 组实际输入共 {r['construction_and_16_dag_checks_seconds']*1000:.3f} ms；表中求值时间另计，不包括进程启动。</p>
      <p>固定同一条线性尾段，检查首轮全部 4 个 DDT 分支（各 p=1/4）：{branchtext}。
      分支求和后的拼接值为 {number(r['fixed_linear_trail_all_first_DDT_branches']['splice_sum'])}。
      此项使用独立局部/HNW 计算器核对，不包括其他线性路径，也不是整体相关度的精确求和。</p>
      <p>已有完整四轮采样的 95% 区间：[{s['ci95'][0]:+.12f}, {s['ci95'][1]:+.12f}]。
      历史完整四轮算法 3 状态：<b>资源上限、未完成</b>；原因为 <code>{esc(r['previous_full_algorithm3']['reason'])}</code>。</p>
      <p>{' · '.join([link(H/'raw'/f'{code}_algorithm3.json','算法 3 原始结果'),link(H/'raw'/f'{code}_algorithm5.json','算法 5 原始结果与逐步记录'),link(H/'raw'/f'{code}_dag.json','完整形式状态 DAG')])}</p></section>''')

page='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>G4a / G4b 中间连接段：DATF 算法 3 与算法 5</title><style>
:root{color-scheme:light;--ink:#193043;--muted:#596e7b;--line:#d9e3e9;--blue:#075b84}
*{box-sizing:border-box}body{margin:0;background:#f1f5f7;color:var(--ink);font:16px/1.75 system-ui,"Noto Sans CJK SC",sans-serif}
main{max-width:1320px;margin:36px auto;padding:30px 38px;background:white;border:1px solid var(--line);border-radius:12px}
h1{font-size:28px;line-height:1.4;margin:0 0 10px}h2{font-size:21px;margin-top:30px}p{margin:14px 0}a{color:var(--blue);text-underline-offset:3px}
.meta,small{color:var(--muted)}small{display:block;font-size:12px;line-height:1.7}.callout{padding:15px 20px;background:#eaf5ef;border-left:4px solid #28815d}
.note{padding:13px 18px;background:#fff6df;border-left:4px solid #b98821}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px;min-width:950px}
th,td{border:1px solid var(--line);padding:12px 10px;text-align:left;vertical-align:top}th{background:#edf3f7}td b{font-size:18px}td:first-child{font-weight:700}
code,pre{font-family:ui-monospace,monospace;font-size:13px;overflow-wrap:anywhere}pre{background:#f5f7f9;padding:14px;margin:8px 0;line-height:1.65;overflow:auto}
.endpoints{display:grid;grid-template-columns:1fr 1fr;gap:24px}section{border-top:1px solid var(--line);margin-top:28px}.foot{font-size:13px;color:var(--muted)}
@media(max-width:720px){main{margin:0;padding:20px;border-radius:0}h1{font-size:24px}.endpoints{grid-template-columns:1fr;gap:8px}}
</style></head><body><main><h1>G4a / G4b 中间连接段复算</h1>
<p class="meta">2026-10-27 · 原始 DATF · JOC 算法 3 与折中算法 5 · 相关度 ρ = 2ε</p>
<div class="callout"><b>两个中间段都能计算。</b>算法 3 与算法 5 一致：G4a 为 +1/4，G4b 为 +1。
此前“未完成”指完整四轮端点，不是这里的两轮／一轮连接段。本次没有重算完整四轮，也没有新增采样。</div>
<h2>计算与拼接结果</h2><p>表中 DATF 两列仅表示中间段 c<sub>m</sub>。拼接列按原来保留的差分路径和线性路径计算 p·c<sub>m</sub>·q²。</p>
<div class="scroll"><table><thead><tr><th>路径</th><th>差分+中间+线性</th><th>p</th><th>q²</th><th>DATF(JOC算法3)<small>中间段</small></th>
<th>DATF(JOC算法5)<small>中间段</small></th><th>单路径拼接估计</th><th>完整四轮实验相关度</th></tr></thead><tbody>'''+''.join(tr)+'''</tbody></table></div>
<p class="note"><b>计算完成不等于证明整体拼接准确。</b>G4a 的拼接值 1/256 位于已有采样区间内；G4b 的单路径拼接值 1/1024 明显低于采样值（约 1/249）。
将算法 3 换为算法 5 没有改变中间值，也没有消除 G4b 的差异。其他线性路径贡献、拼接独立性和条件分布的影响仍需另行研究。</p>
<h2>为什么这些中间段很快</h2><p>G4a 虽然中间输入差分有 21 个活跃比特，但输出掩码仅 2 比特；算法 3 完全展开 β 后只有 21 个单项式、9 个变量，而且是二次式。
G4b 的输入差分／输出掩码分别有 21／40 个活跃比特，但单轮 χ 的目标差分奇偶多项式化简后恒为 0，因而局部相关度严格为 +1。
端点的活跃比特数不能单独决定求值难度。</p>
<h2>算法 5 的实现口径</h2><p>依论文第 19 页算法 5，第 2 行初始化输入分布；第 6 行按拓扑顺序用 BiasP 为每个 α、β 定义计算边缘分布；
第 12 行用这些分布计算非孤立部分；第 13 行代回孤立变量。BiasP 使用精确有理数、变量分离与加权 Shannon 分解，均匀二次式使用等价的二次型求值。</p>
<p>没有递归调用整个算法 5 来替代第 6 行，也没有强行把所有 α 设为均匀。状态构建沿用作者参考代码的单变量保留与相同／互补 β 复用规则，
转换边界沿用项目已有的 χ 前边界。与算法 3 相比，算法 5 额外假设过渡变量的独立性；所得数值是模型估计。</p>
<p>本例 G4a 新引入的 320 个 α 和 73 个 β 的边缘分布都算得均匀，最终目标没有孤立变量，故算法 5 的非均匀分布处理没有带来数值改变。
G4b 仅一个 χ 层，不需要引入过渡变量。算法 5 不是每个短例都更快：它会预计算边缘分布，而算法 3 在 G4a 中只需展开两个相关 β。</p>
'''+''.join(details)+'''<section><h2>复现与核验</h2>
<p>已通过 140 个小规模非均匀 BiasP 穷举对照、逐层边缘分布与非均匀 α 专项检查、32 组实际置换 DAG 检查、64 组算法 3 形式展开检查，
以及差分 DDT、线性 LAT、线性层边界的核验。G4b 用精确局部 DLCT 交叉核对，G4a 用独立的两轮 HNW 零值掩码模型计算器核对。</p>
<p>在本目录运行 <code>python3 run.py</code>，再运行 <code>python3 build_report.py</code>。每项求值上限为 120 秒、100 万项、50 万节点；资源上限始终记为未完成。</p>
<p>'''+ ' · '.join([link(H/'results.json','机器可读结果'),link(H/'validation.json','核验记录'),link(H/'manifest.json','中间段任务'),link(H/'algorithm5.py','算法 5 实现'),link(H/'run.py','复算脚本'),link(H/'source_manifest.json','来源与校验和'),link(ROOT/'路径分段汇总_20261027/gaston_hsilu_segmented_results.html','整体与分段汇总'),link(ROOT/'文献/Differential-Linear Cryptanalysis from an Algebraic Perspective.pdf','DATF 原论文')])+'''</p>
<p class="foot">所有端点与轮常数窗口保持原始匹配记录；完整四轮算法 3 的未完成状态保持原样。中间段的均匀输入估计不包含差分路径条件下的真实状态分布。</p>
</section></main></body></html>'''
(H/'g4_middle_algorithm5_report.html').write_text(page)
class Check(HTMLParser):
    def __init__(self):super().__init__();self.links=[];self.ids=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag=='a':self.links.append(d.get('href',''))
        if 'id' in d:self.ids.append(d['id'])
c=Check();c.feed(page)
for href in c.links:
    if href.startswith('/'):assert Path(unquote(href)).exists(),href
    if href.startswith('#'):assert href[1:] in c.ids
assert 'DATF(JOC算法3)' in page and 'DATF(JOC算法5)' in page
(H/'report_validation.json').write_text(json.dumps(dict(status='PASS',local_links_checked=len(c.links),report_rows=2),indent=2)+'\n')
print(H/'g4_middle_algorithm5_report.html')

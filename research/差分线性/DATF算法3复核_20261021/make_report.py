import hashlib
import html
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

OUT=Path(__file__).resolve().parent
manifest=json.loads((OUT/'manifest.json').read_text())
catalog=json.loads((OUT.parent/'研究汇总_20261019/results_catalog.json').read_text())
attempts={}
for folder in ['representatives','representatives_extended','shannon','batch','extensions_extended']:
    for p in sorted((OUT/folder).glob('*.json')):
        if '.anf.' in p.name:continue
        d=json.loads(p.read_text())
        assert d['status'] in ('completed','resource_limit'),(str(p),d['status'])
        if d['status']=='resource_limit':assert d['correlation'] is None and d['bias'] is None
        d['result_file']=str(p.relative_to(OUT))
        attempts.setdefault(d['id'],[]).append(d)

rows=[]
for j in manifest['jobs']:
    all_attempts=attempts[j['id']]
    complete=[a for a in all_attempts if a['status']=='completed']
    if complete:
        values={Fraction(int(a['correlation']['numerator']),int(a['correlation']['denominator'])) for a in complete}
        assert len(values)==1,(j['id'],values)
        selected=complete[-1]
    else:
        # Prefer the enlarged literal Algorithm-3 attempt for structural stats.
        selected=max(all_attempts,key=lambda a:(a.get('evaluation_order','full beta expansion then Bias').startswith('full'),a['limits']['max_terms'],a['elapsed_seconds']))
    rows.append({**j,'status':selected['status'],'algorithm3':selected['correlation'],
                 'bias':selected['bias'],'selected_attempt':selected['result_file'],
                 'statistics':selected.get('statistics',{}),'reason':selected.get('reason'),
                 'elapsed_seconds':selected['elapsed_seconds'],
                 'attempts':[{k:v for k,v in a.items() if k not in ['sampling','previous_hnw','previous_datf']} for a in all_attempts]})

counts=Counter(x['status'] for x in rows)
summary={'date':'2026-10-21','method':'Original paper Algorithm 3: all beta eliminated, uniform independent remaining variables',
         'previous_method':'custom isolated-variable estimator with recursively computed nonuniform beta marginals; NOT original Algorithm 3; not a literal implementation of original Algorithm 5',
         'correlation_equals_twice_bias':True,'first_L_omitted':True,
         'Hsilu_mask_boundary':'core output, final L omitted as in the original archived endpoints',
         'scope':manifest['scope'],'excluded':manifest['excludes'],'counts':dict(counts),'rows':rows,
         'algorithm3_is_a_model':'Alpha/input variables at different cuts are treated jointly uniform independent; no claim that the actual permutation correlation is exact unless a separate balance/structural proof exists',
         'source_sha256':manifest['source_sha256'],
         'implementation_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.py')}}
method_files=[OUT.parent/'文献/Differential-Linear Cryptanalysis from an Algebraic Perspective.pdf',
              OUT.parent/'文献/DATF/4rascon.sage',
              OUT.parent/'gaston_dl/tools/datf.py',OUT.parent/'gaston_dl/tools/datf_k06.py',
              OUT.parent/'gaston_dl/tools/datf_quadratic.py',OUT.parent/'hsilu_dl/tools/hsilu.py',
              OUT.parent/'gaston_dl/tools/linear_seed_scan.py',OUT.parent/'gaston_dl/tools/fixed_tail_search.py']
summary['method_source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in method_files}
(OUT/'results.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')

def esc(x):return html.escape(str(x),quote=True)
def val(x):
    if x is None:return '未算出'
    if isinstance(x,str):return esc(x)
    n,d=int(x['numerator']),int(x['denominator'])
    if not n:return '0'
    if d==1:return f'{n:+d}'
    if abs(n)==1 and d&(d-1)==0:return ('−' if n<0 else '+')+f'2<sup>−{d.bit_length()-1}</sup>'
    if len(str(n))+len(str(d))>70:
        log=x.get('log2_abs',math.log2(abs(n))-math.log2(d))
        return ('−' if n<0 else '+')+f'2<sup>{log:.6f}</sup>（近似记法）'
    return f'{n:+d}/{d}'
def coords(words):
    bits=[f'({i},{j})' for i,w in enumerate(words) for j in range(64) if int(w,16)>>j&1]
    if len(bits)>18:return f'{len(bits)} 比特；十六进制：<code>'+esc(' '.join(words))+'</code>'
    return '{'+', '.join(bits)+'}'
def table(head,body):
    return '<div class="tablewrap"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(c)+'</td>' for c in row)+'</tr>' for row in body)+'</tbody></table></div>'
def status(x):return '<span class="ok">已完成</span>' if x['status']=='completed' else '<span class="limit">资源限制，未得数值</span>'
def sample(x):
    s=x.get('sampling')
    if not s:return '无对应采样'
    text=f'{s["rho"]:+.12f}（2<sup>{int(math.log2(s["samples"]))}</sup> 对）'
    if x['id']=='H5':
        ref=next(c for c in catalog['endpoints'] if c['code']=='H5')['sampling_reference']
        text+='<br>独立参考实现：'+f'{ref["rho"]:+.12f}；同为 2<sup>28</sup> 对。尚未独立确认。'
    return text

mainrows=rows[:13]
changed=[r for r in rows if r['status']=='completed' and isinstance(r.get('previous_datf'),dict)
         and Fraction(int(r['previous_datf']['numerator']),int(r['previous_datf']['denominator'])) != Fraction(int(r['algorithm3']['numerator']),int(r['algorithm3']['denominator']))]
parts=[]
parts.append(f'<div class="notice"><b>结论：</b>旧实现不是原论文算法 3。本次对 {len(rows)} 个不同端点／常数窗口进行算法 3 复算，{counts["completed"]} 项完成，{counts["resource_limit"]} 项在本次资源预算内未算出。未算出不等于 0，也不是已证明无法计算。</div>')
parts.append('<h2 id="method">1　算法编号与实现纠正</h2>')
parts.append('<p>依据本地论文第 17–19 页及作者参考代码 <code>4rascon.sage</code> 核对。旧报告“普通 DATF”主要指项目中的 <code>OrdinaryDATF.estimate</code>：先分离孤立变量，只代入这些变量，并把其余部分按独立边缘分布求和；β 边缘值又递归调用该估计器。它接近算法 5 的加权孤立变量思路，但边缘更新并非原文算法 5 第 6 行规定的逐层 BiasP，不能直接标为算法 5。C++ 两轮中间段实现也采用孤立变量因子法；单轮精确 DLCT 和结构恒等式属于另行可严格求解的特例。</p>')
parts.append(table(['方法','代入规则','分布／近似位置'],[
 ['原论文算法 3','所有 β 代入到完全消失；不展开 α 定义','最后一次 Bias；将剩余原始输入及跨轮 α 视为联合均匀独立'],
 ['原论文算法 4','逐步代入孤立变量','分步因子相乘，各次 Bias 视变量均匀独立'],
 ['原论文算法 5','逐步代入孤立变量；前向 BiasP 更新边缘','保留非均匀边缘，但仍有独立性及跨因子近似'],
 ['此前项目实现','加权孤立变量法；递归计算 β 边缘','非算法 3，也非算法 5 原伪代码的严格复现'],
 ['本次实现','完整消去 β，精确求形式变量均匀和','算法 3 的模型；求和用 Separate、Shannon、二次消元和精确仿射消元']]))
parts.append('<p>本报告统一记录相关度 <b>ρ = E[(−1)<sup>e</sup>] = 2ε</b>；论文伪代码返回偏差 ε。本次没有引入增强 DATF 的转移规则或近似。G4a 还尝试了一个与算法 3 完全等价的求值顺序：对已存在的独立形式输入变量先做精确 Shannon 分解，再在每个分支消去全部 β，最后按 2<sup>−k</sup> 加权。该尝试也未完成，保留独立日志。</p>')
parts.append('<p>算法 3 保留更多 β 依赖，但不能无条件称为实际置换的“精确相关度”。原文定理 10 使用独立性假设；定理 11 另要求 X→α* 的联合映射平衡。形式变量空间里的精确求和不自动证明这一条件。G2／H3 的目标多项式恒为 0，则可沿用已有结构结论。</p>')
parts.append('<h2 id="main">2　原汇总报告 13 个代表端点</h2>')
parts.append('<p>轮数、输入差分、输出掩码、常数窗口均沿用原归档。仍省略第一个线性层；Hsilu 使用原报告的核心末端掩码。采样值仅从旧数据引入，本次没有重新采样，也没有重跑 HNW26。</p>')
parts.append(table(['编号','轮数','旧快速 DATF','算法 3 相关度','算法 3 状态','已有采样'],[[r['id'],r['rounds'],val(r['previous_datf']),val(r['algorithm3']),status(r),sample(r)] for r in mainrows]))
parts.append('<p>G3a 的新值 −3/16 = −0.1875 与已有采样 −0.187284469604 相容；旧值 −5/32 = −0.15625 不相容。另对其展开后的 20 个形式变量独立遍历全部 2<sup>20</sup> 赋值，得到 N₀=425,984、N₁=622,592，直接确认算法 3 模型值为 −3/16。这个枚举是对形式变量模型的复核。</p>')
parts.append('<h2 id="changed">3　已完成且改变旧数值的项目</h2>')
parts.append(table(['端点／前缀','轮数','旧快速 DATF','算法 3','耗时'],[[esc(r['id']),r['rounds'],val(r['previous_datf']),val(r['algorithm3']),f'{r["elapsed_seconds"]:.3f} s'] for r in changed]))
parts.append('<p>K06 两个三轮前缀均恢复 −2<sup>−23</sup>；这与早期受限贡献的数值相同，但本次计算的对象是同一端点的完整算法 3 模型多项式。对应五、六轮完整端点仍未算出。两个向前扩展的两轮前缀得到 +2<sup>−16</sup>，与先前独立精确前缀计算一致；其五、六轮端点也未算出。不能再用旧快速模型的 0 作为这些前缀真实相关度为 0 的证据。</p>')
parts.append('<h2 id="cost">4　计算规模与未完成的原因</h2>')
parts.append('<p>批次初筛每项上限为 30 秒、100 万单个多项式项、10 万 Bias 节点、4 GiB 地址空间；13 个代表端点另加大到 180 秒、500 万项、200 万节点。K06 完整端点和 G4a 延长另有扩大预算的尝试。上限按每次调用独立记录；并行执行时耗时是实测墙钟时间，不能简单相加当作串行 CPU 时间。乘法预算按待展开的单项式对数预检，预算失败的那次乘法未必实际执行。</p>')
cost=[]
for r in mainrows[3:]:
    s=r['statistics'];cost.append([r['id'],s.get('expanded_terms','展开未完成'),s.get('expanded_variables','—'),s.get('max_component_variables','—'),s.get('expanded_degree','—'),s.get('bias_nodes','尚未进入'),f'{r["elapsed_seconds"]:.3f} s',esc(r.get('reason') or '完成')])
parts.append(table(['端点','完全展开项数','剩余变量','最大连通分量变量','次数','Bias 节点','代表尝试耗时','状态原因'],cost))
parts.append('<p>关键量是 β 消去后的单项式规模和变量依赖图，而非端点各自是否超过 3 个活跃比特。单比特输入、两比特输出的四轮端点也可能展开到数百万项；K06 六比特输入的三轮前缀反而可在本模型内迅速完成。资源超限只说明当前实现与预算未得到数值，不构成计算复杂度下界。</p>')
parts.append('<h2 id="all">5　全部复算记录</h2>')
parts.append(f'<p>{len(rows)} 个去重端点来自 {sum(len(r["aliases"]) for r in rows)} 条归档引用，覆盖原表的 13 个代表、Hsilu 四五轮全部 32 个保留候选、三个三轮结构掩码、26 个 HNW26 对照端点、176 个 L2 池保留端点、K06 及重点前后扩展。相同端点若常数窗口不同，分别计算。没有把历史 24,515 次筛选记录或整个搜索空间重新跑一遍。</p>')
parts.append('<div class="controls"><input id="search" placeholder="检索编号、算法、别名" aria-label="检索记录"><select id="statusfilter" aria-label="状态筛选"><option value="">所有状态</option><option value="completed">已完成</option><option value="resource_limit">未算出</option></select><span id="visible"></span></div>')
parts.append('<div class="tablewrap"><table id="allrows"><thead><tr><th>端点 / 轮数</th><th>Δ 与 Λ（核心边界）</th><th>旧快速 DATF</th><th>算法 3</th><th>记录</th></tr></thead><tbody>')
for r in rows:
    aliases=' / '.join(r['aliases'])
    links=' · '.join(f'<a href="{esc(a["result_file"])}">尝试 {i+1}</a>' for i,a in enumerate(r['attempts']))
    parts.append(f'<tr data-status="{r["status"]}" data-search="{esc(aliases+" "+r["algorithm"])}"><td><b>{esc(r["id"])}</b><br>{r["algorithm"]} / {r["rounds"]} 轮 / start={r["start_round_index"]}<details><summary>别名</summary>{esc(aliases)}</details></td><td>Δ={coords(r["delta_core"])}<br>Λ={coords(r["lambda_core"])}</td><td>{val(r["previous_datf"])}</td><td>{val(r["algorithm3"])}<br>{status(r)}</td><td>{links}</td></tr>')
parts.append('</tbody></table></div>')
parts.append('<h2 id="validation">6　验证与复现</h2>')
parts.append('<p>这是算法 3 数学流程的 Python 移植，不是直接执行作者针对 Ascon／Serpent 的 Sage 算例。本机没有 Sage；作者算例的状态初始化也不同。Gaston／Hsilu 的轮函数和归档端点由本地已验证构造提供，新的求值器不调用旧的 estimate 或 marginal。已完成数值通过下述独立小例枚举和代数恒等检查复核。</p>')
parts.append('<p>已用 350 个随机小布尔多项式对照完整枚举检验 Bias 后端；40 组含共享 β 和自由 α 的随机 DAG 对照完整枚举，核对完整 β 代入与等价 Shannon 顺序；另外验证了资源限制不会返回 0。每个端点的代数 DAG 对照 4 组完整整数置换状态，完成的直接展开结果再对照 16 组形式变量赋值。G3a 还有独立的 2<sup>20</sup> 真值表复核。</p>')
parts.append('<p>脚本、端点、原始输入文件 SHA-256、每次尝试的限制／统计／状态和完成项的展开多项式均保存在本目录。<a href="results.json">汇总 JSON</a> · <a href="manifest.json">端点与来源</a> · <a href="algorithm3.py">算法 3 实现</a> · <a href="test_algorithm3.py">测试</a> · <a href="G3a_independent_exhaustive_check.json">G3a 独立枚举</a> · <a href="README.txt">复现说明</a></p>')
parts.append('<p>文献：<a href="../文献/Differential-Linear%20Cryptanalysis%20from%20an%20Algebraic%20Perspective.pdf#page=17">原论文第 17–19 页</a>；<a href="../文献/DATF/4rascon.sage">作者 4rascon.sage</a>。旧实现：<a href="../gaston_dl/tools/datf.py">datf.py</a>。<a href="../研究汇总_20261019/gaston_hsilu_results.html">原研究汇总</a>保留旧模型结果供对照，相关算法身份以本次核对为准。</p>')

css='''*{box-sizing:border-box}body{margin:0;background:#f3f5f8;color:#1e293b;font:16px/1.8 "Noto Sans CJK SC",sans-serif}main{max-width:1280px;margin:24px auto;background:#fff;padding:36px 44px 64px}h1{font-size:30px;line-height:1.4;margin:10px 0}h2{font-size:22px;margin-top:36px;scroll-margin-top:16px}p{margin:12px 0}.meta{color:#64748b}.notice{background:#eef4fc;border-left:4px solid #2563a5;padding:16px 20px;margin:22px 0}.tablewrap{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.65;margin:14px 0 22px}th,td{border:1px solid #d5dde7;padding:10px 12px;text-align:left;vertical-align:top}th{background:#e8eef7}tbody tr:nth-child(even){background:#f8fafc}td{overflow-wrap:anywhere}#allrows td:first-child{min-width:180px}#allrows td:nth-child(2){min-width:280px;max-width:430px}.ok{color:#156746}.limit{color:#8a4d08}a{color:#175f9a}code{font-size:13px;overflow-wrap:anywhere}nav{display:flex;gap:8px 20px;flex-wrap:wrap;margin:18px 0}.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:15px 0}input,select{padding:8px 12px;font:inherit;border:1px solid #b7c6d7;border-radius:4px}details{margin-top:8px}summary{cursor:pointer}sup{line-height:0}@media(max-width:760px){main{margin:0;padding:20px}h1{font-size:24px}th,td{padding:7px}}@media print{body{background:#fff}main{padding:0;margin:0}nav,.controls,#allrows{display:none}h2{break-after:avoid}tr{break-inside:avoid}thead{display:table-header-group}}'''
js='''const s=document.querySelector('#search'),f=document.querySelector('#statusfilter'),rows=[...document.querySelectorAll('#allrows tbody tr')];function filter(){let n=0;const q=s.value.toLowerCase();for(const r of rows){const yes=(!f.value||r.dataset.status===f.value)&&r.dataset.search.toLowerCase().includes(q);r.hidden=!yes;if(yes)n++;}document.querySelector('#visible').textContent=n+' / '+rows.length+' 项';}s.addEventListener('input',filter);f.addEventListener('change',filter);filter();'''
out='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DATF 算法 3 复算与旧实现审计</title><style>'+css+'</style></head><body><main><p class="meta">Gaston · Hsilu / 2026-10-21 / 首个线性层仍省略</p><h1>DATF 算法 3 复算与旧实现审计</h1><nav><a href="#method">算法身份</a><a href="#main">代表端点</a><a href="#changed">改变的结果</a><a href="#cost">计算规模</a><a href="#all">全部记录</a><a href="#validation">验证与复现</a></nav>'+''.join(parts)+'</main><script>'+js+'</script></body></html>'
(OUT/'datf_algorithm3_report.html').write_text(out)
(OUT/'README.txt').write_text('''原论文 DATF Algorithm 3 复算，2026-10-21

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
''')
print(json.dumps({'unique_endpoints':len(rows),'counts':dict(counts),'changed_completed':[r['id'] for r in changed],'html':str(OUT/'datf_algorithm3_report.html')},ensure_ascii=False))

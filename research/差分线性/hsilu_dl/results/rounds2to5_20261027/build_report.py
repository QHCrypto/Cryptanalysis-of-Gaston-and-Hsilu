"""Build the Hsilu 2–5 round audit and candidate report from archived evidence."""
import csv,json,math,html,hashlib,os
from pathlib import Path
from fractions import Fraction
from urllib.parse import quote
from search_new import H,ROOT,BASE,rec,save
from hsilu import inverse_function,transpose,words,encode

def read(p):return json.loads(p.read_text())
def val(d):return Fraction(int(d['numerator']),int(d['denominator']))
def decode(a):return tuple(int(w,16) for w in a)
def sample(d):
 d=dict(d)
 for k in ['samples','n0','n1']:d[k]=int(d[k])
 d['rho']=(d['n0']-d['n1'])/d['samples'];d['log2_abs']=math.log2(abs(d['rho'])) if d['rho'] else None
 return d
manifest=read(H/'manifest.json');hnw=read(H/'hnw_results.json');research=read(ROOT/'hsilu_dl/results/research_summary.json')
historical=read(ROOT/'DATF算法3复核_20261021/results.json')['rows']
segmented={r['code']:r for r in read(ROOT/'路径分段汇总_20261027/results.json')['rows'] if r['algorithm']=='Hsilu'}
initial={r['id']:sample(r) for r in csv.DictReader((BASE/'initial_samples.csv').open())}
discovery={r['id']:sample(r) for r in csv.DictReader((BASE/'pool_samples.csv').open())}
confirmed={r['original_id']:r['sampling'] for r in segmented.values() if r.get('original_id') and r.get('experimental_status')=='clear_signal'}
invT=inverse_function(transpose)
rows=[]
def full_record(R,d,m):
 for r in historical:
  if r['algorithm']=='Hsilu' and r['rounds']==R and r['delta_core']==d and r['lambda_core']==m:
   return dict(status=r['status'],reason=r.get('reason'),elapsed_seconds=r['elapsed_seconds'],source=str(ROOT/'DATF算法3复核_20261021'/r['selected_attempt']),statistics=r.get('statistics'))
 return None

def add(code,R,d,m,split,jobid,q2,sampling=None,original=None,note='',whole=None):
 raw=read(H/'raw'/f'{jobid}.json');assert raw['status']=='completed'
 c=val(raw['correlation']);q2=Fraction(q2);score=c*q2
 native=words(invT(decode(m)));assert transpose(decode(native))==decode(m)
 pair=hnw.get(jobid)
 row=dict(code=code,rounds=R,input_difference=d,core_output_mask=m,full_output_mask=native,split=split,
  evaluated_part='whole' if split==[0,R,0] else 'prefix/middle',datf_algorithm3=raw['correlation'],datf_seconds=raw['elapsed_seconds'],
  datf_source=str(H/'raw'/f'{jobid}.json'),hnw26={v:r['correlation'] for v,r in pair.items()} if pair else None,
  p=rec(Fraction(1)),q_squared=rec(q2),reported_correlation=rec(score),sampling=sampling,original_id=original,note=note,
  whole_attempt=whole or (full_record(R,d,m) if R>=4 else dict(status='completed',source=str(H/'raw'/f'{jobid}.json'))))
 rows.append(row);return row
for j in manifest['jobs']:
 if j['scope']=='whole':
  id=j['id'];R=j['rounds'];code={'R2_p1_c1_m1':'H2','R3_p8_c29_m1':'H3'}.get(id,id)
  add(code,R,j['delta_core'],j['lambda_core'],[0,R,0],id,1,initial[id],id,'结构性恒等于 +1' if id in ['R2_p1_c1_m1','R3_p8_c29_m1'] else '整体模型值；保留历史采样对照')
labels={'R4_pool_000026':'H4a','R4_pool_000027':'H4b','R4_pool_000000':'H4c','R4_pool_000001':'H4d'}
for c in research['four_candidates']:
 id=c['id'];row=add(labels.get(id,id),4,c['input_difference']['words'],c['core_output_mask']['words'],[0,2,2],manifest['candidate_to_job'][id],Fraction(1,2**(2*c['linear_weight'])),confirmed.get(id,discovery[id]),id,
 '独立 2^28 对确认' if id in confirmed else '2^26 对发现批次；未列入四个独立确认代表')
 row['historical_exact_prefix2']=c['exact_prefix_correlation'];assert val(row['datf_algorithm3'])==Fraction(c['exact_prefix_correlation'])
old5_primary=read(ROOT/'hsilu_dl/results/H5_sampling_2p28_20261021/results.json')['sampling']
old5_reference=read(ROOT/'hsilu_dl/results/H5_sampling_2p28_20261021/reference_results.json')['sampling']
for c in research['five_candidates']:
 id=c['id'];row=add('H5-old' if id=='R5_pool_000000' else id,5,c['input_difference']['words'],c['core_output_mask']['words'],[0,3,2],manifest['candidate_to_job'][id],Fraction(1,2**(2*c['linear_weight'])),old5_primary if id=='R5_pool_000000' else None,id,'旧候选：原算法 3 的三轮前缀为 0；不保留旧自定义 DATF 非零评分')
 if id=='R5_pool_000000':row['independent_reference_sampling']=old5_reference
optimized=read(H/'search_023.json')['top'][:4]
for i,c in enumerate(optimized,1):
 id=c['id'];whole=read(H/'raw'/f'{id}_whole5.json')
 row=add('H5-new'+str(i),5,c['input_difference'],c['output_mask'],[0,2,3],id+'_prefix2',val(c['q2']),None,id,'新分段候选；未实验确认，非五轮整体精确相关度',dict(status=whole['status'],reason=whole.get('reason'),elapsed_seconds=whole['elapsed_seconds'],statistics=whole['statistics'],source=str(H/'raw'/f'{id}_whole5.json')))
 assert val(row['datf_algorithm3'])==val(c['exact_prefix2'])==Fraction(1,2)
 row['trail']=c
searches={k:read(H/f'search_{k}.json') for k in ['221','122','113','023']}
validation=read(H/'validation.json');validation['native_mask_conversion']='passed for all '+str(len(rows))+' report rows';save(H/'validation.json',validation)
sources=[ROOT/'hsilu_dl/tools/hsilu.py',ROOT/'hsilu_dl/src/main.cpp',ROOT/'hsilu_dl/src/hnw_endpoint.cpp',ROOT/'hsilu_dl/results/research_summary.json',BASE/'linear_pool.json',ROOT/'DATF算法3复核_20261021/algorithm3.py',ROOT/'DATF算法3复核_20261021/worker.py',ROOT/'DATF算法3复核_20261021/results.json',ROOT/'路径分段汇总_20261027/results.json',ROOT/'hsilu_dl/results/H5_sampling_2p28_20261021/results.json',ROOT/'hsilu_dl/results/H5_sampling_2p28_20261021/reference_results.json']+list(H.glob('*.py'))
provenance=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]
doc=dict(date='2026-10-27',algorithm='Hsilu',quantity='rho = E[(-1)^(Lambda dot (F(X) xor F(X xor Delta)))]; epsilon = rho/2',
 policy='Use completed full Algorithm 3 where feasible; otherwise evaluate middle/prefix with original Algorithm 3 and splice a single linear trail.',
 boundary='Core: first chi input to last chi output; every internal L retained. Native output mask = L^{-T}(core mask). First RC absorbed into uniform input translation; last R round constants.',
 rows=rows,searches={k:{a:b for a,b in d.items() if a not in ['top','perhead']} for k,d in searches.items()},validation=validation,
 limitations=['Algorithm 3 expands all beta definitions, retaining the original DATF independent-uniform alpha model.',
 'A completed model zero is not a proof of zero for the actual full permutation.',
 'Resource-limit states have no numeric value and are never treated as zero.',
 'All five-round nonzero values here are single-tail spliced estimates; omitted trails and boundary dependencies may reinforce or cancel them.',
 'Two-round prefix exact checks validate the short prefix only, not the five-round endpoint.',
 'Search spaces are bounded families, not an exhaustive search of all Hsilu five-round endpoints.'],source_sha256=provenance)
save(H/'results.json',doc)

E=html.escape
def link(path,label):return '<a href="'+quote(os.path.relpath(path,H),safe='/._-')+'">'+E(label)+'</a>'
def bits(ws):return '{'+', '.join('('+str(i)+','+str(j)+')' for i,w in enumerate(decode(ws)) for j in range(64) if w>>j&1)+'}'
def qr(d):
 q=val(d)
 if not q:return '0'
 if q.denominator==1:return ('+' if q>0 else '−')+str(abs(q.numerator))
 sign='+' if q>0 else '−'
 if abs(q.numerator)==1 and q.denominator&(q.denominator-1)==0:
  return sign+'2<sup>−'+str(q.denominator.bit_length()-1)+'</sup>'
 return sign+str(abs(q.numerator))+'/'+str(q.denominator)
def sampling(r):
 s=r.get('sampling')
 if not s:return '<span class="muted">未验证</span>'
 log=s.get('log2_abs');rho=s['rho'];t=f'{rho:+.10g}<br><small>log₂|ρ̂| = {log:.6f}；2<sup>{int(math.log2(s["samples"]))}</sup> 对</small>'
 if s.get('ci95'):t+='<br><small>95% CI ['+', '.join(f'{x:.9g}' for x in s['ci95'])+']</small>'
 return t

def table(rs,showq=True):
 hd=['标识 / 轮数','输入差分 Δ','核心输出掩码 Λ','分段 (d,m,ℓ)','DATF(JOC算法3)','HNW26（同一段）']
 if showq:hd+=['线性因子 q²','拼接估计 ρ̃']
 hd+=['已有采样','判定与整体可计算性']
 body=[]
 for r in rs:
  h=r['hnw26'];hv=qr(h['full']) if h else '未计算'
  if h and val(h['full'])!=val(h['zero']):hv+=' / zero: '+qr(h['zero'])
  cells=[E(r['code'])+'<br><small>'+str(r['rounds'])+' 轮</small>',E(bits(r['input_difference'])),E(bits(r['core_output_mask'])),str(tuple(r['split'])),qr(r['datf_algorithm3'])+'<br><small>'+('整体' if r['evaluated_part']=='whole' else str(r['split'][1])+' 轮前缀')+' · '+link(r['datf_source'],'原始记录')+'</small>',hv]
  if showq:cells+=[qr(r['q_squared']),qr(r['reported_correlation'])]
  w=r['whole_attempt'];ws='整体算法 3 完成' if w and w['status']=='completed' else ('整体算法 3：资源限制' if w else '整体状态未归档')
  if w:ws=link(w['source'],ws)
  cells +=[sampling(r),E(r['note'])+'<br><small>'+ws+'</small>']
  body.append('<tr>'+''.join('<td>'+c+'</td>' for c in cells)+'</tr>')
 return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+x+'</th>' for x in hd)+'</tr></thead><tbody>'+''.join(body)+'</tbody></table></div>'
short=[r for r in rows if r['rounds']<=3];four=[next(r for r in rows if r['code']==c) for c in ['H4a','H4b','H4c','H4d']];old=next(r for r in rows if r['code']=='H5-old');fresh=[r for r in rows if r['code'].startswith('H5-new')]
parts=[]
parts.append('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Hsilu 2–5 轮差分线性结果 · 算法 3 复核与新五轮候选</title><style>
:root{--ink:#14273d;--muted:#53677b;--line:#d8e2ec;--accent:#075b78;--bg:#f4f7fb}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.65 system-ui,-apple-system,"Noto Sans CJK SC","Microsoft YaHei",sans-serif}main{max-width:1560px;margin:auto;padding:34px 32px 65px}header{border-bottom:3px solid var(--accent);padding:10px 0 24px}h1{font-size:32px;line-height:1.3;margin:10px 0}h2{font-size:24px;margin:35px 0 15px}h3{font-size:19px}a{color:#006c91;text-underline-offset:3px}p{max-width:1180px}nav{display:flex;flex-wrap:wrap;gap:18px;margin-top:20px}.eyebrow{letter-spacing:.12em;font-size:13px;color:var(--accent)}.muted,small{color:var(--muted)}small{font-size:12px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:24px 0}.card{background:white;border:1px solid var(--line);border-radius:10px;padding:18px}.card strong{display:block;font-size:26px;color:var(--accent)}.callout{border-left:4px solid #158197;background:#e8f4f6;padding:16px 20px;margin:20px 0}.warning{background:#fff5e5;border-color:#b87b20}.scroll{overflow:auto;background:white;border:1px solid var(--line);border-radius:8px;margin:14px 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:12px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{background:#eaf0f6;white-space:nowrap;color:#173c56}tbody tr:nth-child(even){background:#f8fbfd}td:first-child{font-weight:650;white-space:nowrap}code,pre{font-family:ui-monospace,"DejaVu Sans Mono",monospace}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#eef3f7;padding:16px;border-radius:8px;font-size:12px}details{background:white;border:1px solid var(--line);padding:16px;border-radius:8px;margin:16px 0}summary{cursor:pointer;font-weight:650}.formula{font:20px/1.7 Georgia,"Noto Serif CJK SC",serif;background:white;border:1px solid var(--line);padding:18px;text-align:center}.badge{font-size:12px;border:1px solid #6e939e;border-radius:12px;padding:2px 10px}.trail{display:flex;align-items:center;gap:10px;flex-wrap:wrap;padding:20px;background:white;border:1px solid var(--line);border-radius:8px}.block{padding:14px 22px;background:#e4f1f4;border-radius:8px;text-align:center}.tail{background:#edf0fb}.foot{margin-top:35px;border-top:1px solid var(--line);padding-top:20px;font-size:13px}@media(max-width:800px){main{padding:18px 14px}h1{font-size:26px}.cards{grid-template-columns:1fr 1fr}th,td{padding:9px 8px}.formula{font-size:17px}}@media print{body{background:white}main{max-width:none;padding:0}.scroll{overflow:visible}table{font-size:9px}th,td{padding:5px}details{break-inside:avoid}nav{display:none}}
</style><main><header><div class="eyebrow">HSILU · RESEARCH RECORD · 2026-10-27</div><h1>2–5 轮差分线性结果<br>DATF 算法 3 复核与新五轮分段候选</h1><p>整体算法 3 能完成时采用整体结果；达到资源上限时采用连接段的算法 3 值与单条线性尾段拼接。采样数据保留原批次与原样本量。旧自定义 DATF 非零评分不作为 JOC 算法 3 的结果。</p><nav><a href="#definition">边界与计算口径</a><a href="#short">2–3 轮</a><a href="#four">4 轮</a><a href="#five">5 轮</a><a href="#search">搜索范围</a><a href="#audit">复核与数据</a></nav></header>
<div class="cards"><div class="card">2 轮 / 整体<strong>ρ = +1</strong><small>算法 3 + 结构证明 + 采样</small></div><div class="card">3 轮 / 整体<strong>ρ = +1</strong><small>另有 ±3/4 的整体结果</small></div><div class="card">4 轮 / 已确认<strong>log₂|ρ̂| ≈ −9.014</strong><small>H4a · 独立 2²⁸ 对采样</small></div><div class="card">新 5 轮 / 分段<strong>ρ̃ = +2<sup>−63</sup></strong><small>4 个端点 · 尚无实验确认</small></div></div>
<section id="definition"><h2>1. 边界、分段与数值含义</h2><p>轮数按 χ 层数计算。沿用已有 Hsilu 实验：Δ 位于首个 χ 输入，Λ 位于最后一个 χ 输出，内部所有 L 层均保留。原生轮顺序为 RC → χ → L；首个 RC 对均匀输入可用平移吸收，末尾 L 用外部掩码 Λ<sub>native</sub> = L<sup>−T</sup>Λ 处理。因此不是另行删除一个内部线性层。每个 R 轮结果使用 12 轮常数的最后 R 个；四轮从索引 8 开始，五轮从索引 7 开始。</p><div class="formula">ρ = E[ (−1)<sup>Λ·(F(X) ⊕ F(X ⊕ Δ))</sup> ]， &nbsp; ε = ρ/2<br>分段 (d,m,ℓ)： &nbsp; ρ̃ = p · c<sub>m</sub> · q²， &nbsp; R = d + m + ℓ</div><p>本报告主要采用 (0,2,2)、(0,3,2)、(0,2,3)，即 p = 1，把前缀整体交给算法 3，再接一条线性尾段。表头 <b>DATF(JOC算法3)</b> 显示所标明整体或连接段的有符号相关度；“拼接估计”才是整个 R 轮的模型评分。q 为一条线性路径的相关度，表中的 q² 已经平方。</p><div class="callout warning"><b>证据层级：</b>算法 3 完整消去 β，但仍有原 DATF 模型对剩余 α 的独立均匀假设；完成计算不自动意味着实际置换的严格精确值。短前缀另做精确枚举，可确认该前缀本身。分段还引入边界独立性与单条线性路径主导假设。资源限制表示“未求出数值”；模型结果 0 也不能直接升级为真实整体相关度为 0。</div></section>''')
parts.append('<section id="short"><h2>2. 两轮与三轮：使用整体算法 3</h2><p>本次重新计算下表 8 个端点，整体算法 3 均完成。H2 与 H3 的 +1 还有既有的逐点零导数结构证明；其余采样结果用于对照模型值。H3 的输出还可选 {(1,29)} 或 {(0,29),(1,29)}，已有结构证明同样给出 +1，详见原研究记录。</p>'+table(short,False)+'</section>')
parts.append('<section id="four"><h2>3. 四轮：整体受限，改用两轮前缀＋两轮线性尾段</h2><p>四个独立采样确认的代表 H4a–H4d，历史整体算法 3 在 180 秒预算内未完成。本次重算它们的两轮前缀，算法 3、HNW full/zero 与已有精确前缀结果一致。下表中线性尾段均满足 |q| = 2<sup>−5</sup>，故 q² = 2<sup>−10</sup>。</p>'+table(four)+'<div class="callout warning"><b>拼接值与采样值有可见差异。</b>H4a/H4b 的分段估计是 +2<sup>−11</sup>，采样约为 +2<sup>−9.01</sup> 与 +2<sup>−9.03</sup>，约为单条路径估计的 4 倍；H4c/H4d 的 +2<sup>−10</sup> 估计与采样相近。这说明不能因为中间段算得准确，就把拼接值当作整体精确值。其他路径贡献与段间依赖的各自作用尚未分离。</div><details><summary>展开全部 16 个四轮保留候选</summary>'+table([r for r in rows if r['rounds']==4])+'</details></section>')
parts.append('<section id="five"><h2>4. 五轮：旧 H5 的零结果与四个新端点</h2><h3>旧 H5：三轮前缀重新计算仍为 0</h3><p>不只旧 H5：旧池保留的 16 个五轮代表，对应 8 个不同的“输入差分＋三轮前缀输出掩码”；本次原始算法 3 全部完成并得到 0，HNW full/zero 也均为 0。因此旧 (0,3,2) 分段不再保留旧自定义 DATF 的非零评分。这里判定的是该模型与该分段，真实五轮相关度仍没有被严格判零。</p>'+table([old]))
parts.append('<p>旧 H5 两个独立的 2<sup>28</sup> 对批次分别为 <b>'+f'{old5_primary["rho"]:+.12f}'+'</b> 与 <b>'+f'{old5_reference["rho"]:+.12f}'+'</b>；后者的 95% 区间为 ['+', '.join(f'{x:+.12f}' for x in old5_reference['ci95'])+']，包含 0，且符号相反。故继续标为“未独立确认”，不以第一批单独作为有效信号。</p>')
parts.append('''<h3>新候选：两轮前缀＋三轮线性尾段</h3><div class="trail"><div class="block">首两轮 χ<br>算法 3 / 精确枚举：c<sub>m</sub> = +1/2</div><span>→</span><div class="block tail">向前补入的一轮线性 χ<br>|q₁| = 2<sup>−26</sup></div><span>→</span><div class="block tail">原两轮线性尾段 H0064<br>|q₂| = 2<sup>−3</sup> · 2<sup>−2</sup></div></div><div class="formula">(d,m,ℓ) = (0,2,3)， &nbsp; |q| = 2<sup>−31</sup>， &nbsp; q² = 2<sup>−62</sup><br>ρ̃ = (+1/2) × 2<sup>−62</sup> = +2<sup>−63</sup></div><p>这四个输入差分都与旧 H5 不同，输出掩码均为 {(0,10),(2,23)}。两轮前缀并非只取一个 DDT 分支：它对全部首轮取值贡献求和；算法 3、HNW full/zero 与独立精确两轮枚举均得到 +1/2。表中每个完整五轮端点都已尝试整体算法 3，但触及项数或单项式乘积数上限，故按要求采用分段。</p>'''+table(fresh))
parts.append('<div class="callout warning"><b>这次得到的是新的非零分段候选，尚不是已验证的五轮区分器。</b>单条线性尾段之外仍可能有增强或抵消。若将 |ρ̃| = 2<sup>−63</sup> 暂作真实相关度，常用独立样本尺度约为 |ρ̃|<sup>−2</sup> = 2<sup>126</sup> 对（忽略置信度和成功率常数），本机无法直接采样确认；这也不构成完整模式的复杂度结论。</div>')
for r in fresh:
 c=r['trail'];p=c['prepended_chi'];w=r['whole_attempt']
 detail={'input_difference':r['input_difference'],'prefix_output_mask':c['middle_mask'],'native_output_mask':r['full_output_mask'],'prepended_linear_chi':p,'old_first_linear_chi':c['head']['first_chi'],'old_last_linear_chi':c['last_chi'],'same_sign_first_DDT_output_patterns':c['matched_branches'],'full_Algorithm3_status':{k:v for k,v in w.items() if k!='statistics'}}
 parts.append('<details><summary>'+E(r['code'])+'：连接掩码、三层线性路径与整体限制</summary><p>掩码的十六进制数组按行 0–4 排列，最低有效位为列 0。末层原生掩码已校验 LᵀΛ<sub>native</sub> = Λ。</p><pre>'+E(json.dumps(detail,ensure_ascii=False,indent=2))+'</pre></details>')
parts.append('<details><summary>展开全部 16 个旧五轮候选的算法 3 前缀复核</summary>'+table([r for r in rows if r['rounds']==5 and not r['code'].startswith('H5-new')])+'</details></section>')
parts.append('''<section id="search"><h2>5. 这次如何搜索其他五轮候选</h2><p>采用已有线性池并显式保留连接约束，没有使用 MILP。首先尝试两种现有池匹配；未得到非零值后，将原两轮线性尾段向前增加一个 χ。对 χ 的每个输出掩码，其非零 LAT 输入掩码集合可表示为仿射子空间，而且非零 LAT 绝对值相同。程序逐一校验此性质，将中间 χ 的非零导数相关条件转化为 GF(2) 线性方程组，直接求可连接掩码。</p><p>第一阶段得到 1+1+3 的 ±2<sup>−64</sup> 候选；随后将前两轮合并并做精确求和，再求“多个首轮 DDT 分支同号”的联立方程，得到四个 +2<sup>−63</sup> 候选。该优化只改变连接掩码与输入差分，三轮线性尾段的平方权重保持 62。</p><div class="scroll"><table><thead><tr><th>搜索家族</th><th>具体范围</th><th>本次结果</th><th>证据范围</th></tr></thead><tbody>''')
searchrows=[('2+2+1','10 个已认证两轮差分前缀 × 64 旋转 × 20 个单比特输出一轮线性路径',f'{searches["221"]["tested"]:,} 次；非零 0','两轮 HNW zero 模型筛选；不是整体无解证明'),('1+2+2','64 列 × 31 个首轮输出差分模式 × 384 个已有两轮线性头',f'{searches["122"]["tested"]:,} 次；非零 0','同上；首轮输出先分解枚举'),('1+1+3','384 个线性头逐个向前延长；每头枚举 64×31 个差分输出，求 LAT 支持内所有变量的可行性',f'{searches["113"]["tested"]:,} 个可行性问题；{searches["113"]["nonzero"]:,} 个可行匹配；所取自由变量解产生 {searches["113"]["distinct_endpoints"]:,} 个不同外部端点','每个可行系统只取一个自由变量解；并未枚举所有线性尾段'),('0+2+3 复核','保留前 128 个不同端点，把首两轮全部贡献求和','127 个所选前缀非零，1 个为 0；算法 3 与精确枚举逐一一致','仍仅检查所选线性尾段'),('0+2+3 改进','固定最小平方权重 62 的 H0064 延长族；320 个单比特输入，枚举四个 DDT 输出的非空子集并约束同号',f'{searches["023"]["tests"]:,} 个方程组；{searches["023"]["found"]} 个输入取得可行解；四个最好为 +2⁻⁶³','每个可行系统取一个自由变量解；不是全部五轮端点的最优性证明')]
for x in searchrows:parts.append('<tr>'+''.join('<td>'+E(v)+'</td>' for v in x)+'</tr>')
parts.append('''</tbody></table></div><p>循环旋转只用于相对位置搜索；涉及多轮实际计算时均使用对应轮常数和固定坐标重新求值。本报告不把固定轮常数下任意端点的 64 个旋转自动认作等相关度。</p><p>后续更有价值的改进是扩大三轮线性尾段池，降低当前平方权重 62，并对短前缀累加贡献继续复核。当前完整扫描仅覆盖旧的三类稀疏线性桥及其 384 个头；增加活跃列数可能得到不同的三轮代价，尚未搜索。</p></section>''')
parts.append('<section id="audit"><h2>6. 复核记录与可重现材料</h2><ul><li>2–3 轮整体算法 3：8 个端点；四轮前缀：8 个独立任务覆盖 16 个候选；旧五轮前缀：8 个独立任务覆盖 16 个候选。</li><li>新候选两轮前缀：128 次算法 3 与独立精确枚举对照；改进后四个最好前缀再次核对。每个算法 3 任务验证 4 个实际状态的 DAG，并校验 16 个形式变量赋值的展开一致性。</li><li>新保留线性路径：173 条逐层核对非零 LAT、平方因子和转置连接关系。</li><li>快速两轮筛选器：120 个独立 DDT/DLCT 枚举对照、20 个仿射导数映射检查、2 个原生 HNW 控制案例；连接前缀与短轮整体另做 28×2 个 HNW full/zero 任务。</li><li>新五轮整体尝试：7 个端点均记录资源限制。默认上限为 1,000,000 个多项式项、100,000,000 个单项式乘积、200,000 个 Bias 节点、30 秒和 4 GiB 地址空间；本次新五轮主要先触及项数或乘积数限制。</li><li>没有为极弱新五轮相关度运行会被噪声淹没的采样；没有改写旧采样数据。</li></ul>')
parts.append('<p>'+ ' · '.join(link(H/name,label) for name,label in [('results.json','完整结果 JSON'),('validation.json','复核清单'),('middle2_validation.json','快速筛选器交叉验证'),('search_113.json','仿射匹配结果'),('search_023.json','同号贡献优化'),('manifest.json','既有候选任务'),('optimized_manifest.json','改进五轮任务'),('build_report.py','报告生成脚本')])+'</p>')
parts.append('<details><summary>重现入口与来源 SHA-256</summary><p>在本目录运行 python3 evaluate.py、python3 search_new.py、python3 search_113.py、python3 evaluate_new.py、python3 improve_prefix.py。改进候选算法 3 可按 optimized_manifest.json 的 index 0–7 调用原 worker.py；最后运行 finish_checks.py 与 build_report.py。既有历史完整四轮资源限制由来源文件复用，本次未重复消耗相同预算。</p><pre>'+E(json.dumps(provenance,ensure_ascii=False,indent=2))+'</pre></details></section><footer class="foot">生成于 2026-10-27 · 所有数值均可追溯到本地 JSON、原始脚本与历史采样记录。相关度、拼接估计、实验观察和严格结构结论分别标注。</footer></main></html>')
(H/'hsilu_rounds2to5_report.html').write_text(''.join(parts))
print('report rows',len(rows),'new best',[r['reported_correlation'] for r in fresh]);print(H/'hsilu_rounds2to5_report.html')

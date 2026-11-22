"""Build the two representative-result tables without using legacy DATF values."""
import hashlib
import html
import json
import math
from fractions import Fraction
from pathlib import Path
from urllib.parse import quote

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
AUDIT=ROOT/'DATF算法3复核_20261021'
def read(p): return json.loads(p.read_text())
def dump(p,x): p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def frac(x): return Fraction(int(x['numerator']),int(x['denominator']))
def rec(x):
    return dict(numerator=str(x.numerator),denominator=str(x.denominator),
        log2_abs=math.log2(abs(x.numerator))-math.log2(x.denominator) if x else None)
def esc(x): return html.escape(str(x))
def url(p): return quote(str(p),safe='/')
def link(p,text): return f'<a href="{url(p)}">{esc(text)}</a>'
def number(x,model_zero=False):
    q=frac(x) if isinstance(x,dict) else Fraction(x)
    if not q:return '0'+('<small>模型值</small>' if model_zero else '')
    s='−' if q<0 else '+';a=abs(q)
    if a==1:return s+'1'
    if a.numerator==1 and a.denominator&(a.denominator-1)==0:
        return f'{s}2<sup>−{a.denominator.bit_length()-1}</sup>'
    return f'{s}{a.numerator}/{a.denominator}'
def bits(ws):
    return '{'+', '.join(f'({i},{j})' for i,w in enumerate(ws) for j in range(64) if int(w,16)>>j&1)+'}'
def pow2(n):
    assert n>0 and n&(n-1)==0
    return f'2<sup>{n.bit_length()-1}</sup>'

catalog=read(ROOT/'研究汇总_20261019/results_catalog.json')
source_endpoints={e['code']:e for e in catalog['endpoints']}
audits={r['id']:r for r in read(AUDIT/'results.json')['rows']}
cuts=read(OUT/'cuts.json')
jobs={j['id']:j for j in read(OUT/'middle_manifest.json')['jobs']}
validation=read(OUT/'validation.json')
assert validation['status']=='PASS'
order=['G2','G3a','G3b','G4a','G4b','G4c','G4d','H3','H4a','H4b','H4c','H4d','H5']
rows=[]
discovery={
 'G2':'单列输入差分、单比特输出掩码扫描。',
 'G3a':'单列输入差分、单比特输出掩码扫描。',
 'G3b':'差分池与线性池经一轮中间段匹配，原搜索分段 1+1+1。',
 'G4a':'差分池与线性池经两轮中间段匹配，原搜索分段 1+2+1。',
 'G4b':'使用新增 L2_2to2_pool 与差分池匹配，原搜索分段 1+1+2。',
 'G4c':'固定论文中的 D3_chi2_to_3 两轮线性尾段，枚举单列输入差分；保留两个同号连接分支。',
 'G4d':'差分池与固定 D3_chi2_to_3 线性尾段匹配，保留一个差分分支。',
 'H3':'结构分析与扫描得到，另有相关度 1 的结构证据。',
 'H4a':'固定两轮线性尾段，搜索两轮前缀，再用独立采样复核。',
 'H4b':'固定两轮线性尾段，搜索两轮前缀，再用独立采样复核。',
 'H4c':'固定两轮线性尾段，搜索两轮前缀，再用独立采样复核。',
 'H4d':'固定两轮线性尾段，搜索两轮前缀，再用独立采样复核。',
 'H5':'固定两轮线性尾段，搜索三轮前缀；旧搜索使用的估计器不是 JOC 算法 3。'
}
for code in order:
    e=source_endpoints[code]
    r={k:e[k] for k in ['code','algorithm','rounds','delta_core','lambda_core']}
    r.update(start_round_index=12-e['rounds'],discovery=discovery[code],
        sampling=e.get('sampling'),sampling_reference=e.get('sampling_reference'),
        lambda_full=e.get('lambda_full'),original_id=e.get('original_id'))
    if code not in cuts:
        a=audits[code]
        assert a['algorithm']==e['algorithm'] and a['rounds']==e['rounds']
        assert a['start_round_index']==r['start_round_index']
        assert a['delta_core']==e['delta_core'] and a['lambda_core']==e['lambda_core']
        r.update(evaluation_scope='whole',evaluation_split=[0,e['rounds'],0],
            datf_algorithm3=dict(status=a['status'],correlation=a['algorithm3'],
                source=str(AUDIT/a['selected_attempt']),reason=a.get('reason'),
                elapsed_seconds=a['elapsed_seconds'],statistics=a['statistics'],
                attempts=[dict(file=str(AUDIT/t['result_file']),status=t['status'],
                    reason=t.get('reason'),elapsed_seconds=t['elapsed_seconds']) for t in a['attempts']]),
            hnw26=dict(status='completed',correlation=e['hnw'],
                variant=e.get('hnw_variant','full / zero according to archived evaluation'),
                source=str(ROOT/'研究汇总_20261019/results_catalog.json')),
            experimental_status='clear_signal' if code not in ['G2','H3'] else 'unit_correlation_observed')
    else:
        cut=cuts[code];branches=[];datf_total=Fraction(0);hnw_total=Fraction(0)
        for b in cut['branches']:
            j=jobs[b['job_id']];dp=OUT/'raw'/f"{j['id']}_DATF.json";d=read(dp)
            hp=OUT/'raw'/f"{j['id']}_HNW{'' if j['rounds']==1 else '_full'}.json";h=read(hp)
            assert d['status']==h['status']=='completed'
            probability=Fraction(b['probability'])
            branches.append(dict(b,middle_job=j,probability=rec(probability),
                datf_algorithm3=dict(status='completed',correlation=d['correlation'],source=str(dp),
                    elapsed_seconds=d['elapsed_seconds'],statistics=d['statistics']),
                hnw26=dict(status='completed',correlation=h['correlation'],source=str(hp),
                    elapsed_seconds=h['elapsed_seconds'])))
            datf_total+=probability*frac(d['correlation']);hnw_total+=probability*frac(h['correlation'])
        q2=Fraction(1,2**cut['linear_squared_weight'])
        r.update(evaluation_scope='middle_only',evaluation_split=cut['split'],
            linear_squared_correlation=rec(q2),tail_id=cut['tail_id'],
            branches=branches,aggregation=cut['aggregation'],cut_source=str(ROOT/cut['source']),
            datf_spliced_correlation=rec(datf_total*q2),hnw_spliced_correlation=rec(hnw_total*q2),
            experimental_status='not_independently_confirmed' if code=='H5' else 'not_measured')
        assert sum(cut['split'])==e['rounds']
        for b in branches:
            assert b['middle_job']['start_round_index']==r['start_round_index']+cut['split'][0]
        assert datf_total==hnw_total
    rows.append(r)

assert len(rows)==13 and sum(r['evaluation_scope']=='whole' for r in rows)==10
assert sum(r.get('datf_algorithm3',{}).get('status')=='resource_limit' for r in rows)==6
assert frac(next(r for r in rows if r['code']=='G3a')['datf_algorithm3']['correlation'])==Fraction(-3,16)
assert frac(next(r for r in rows if r['code']=='G4c')['datf_spliced_correlation'])==Fraction(1,2**39)
assert frac(next(r for r in rows if r['code']=='G4d')['datf_spliced_correlation'])==Fraction(-1,2**40)
assert frac(next(r for r in rows if r['code']=='H5')['datf_spliced_correlation'])==0
for r in rows:
    for sample in [r.get('sampling'),r.get('sampling_reference')]:
        if sample:
            assert sample['n0']+sample['n1']==sample['samples']
            assert sample['rho']==(sample['n0']-sample['n1'])/sample['samples']

data=dict(title='Gaston 与 Hsilu：整体与分段评估结果',date='2026-10-27',
    datf_column_name='DATF(JOC算法3)',quantity='signed correlation rho = 2 epsilon',
    policy='experimentally resolved representatives: whole; unconfirmed representatives: middle only',
    new_full_endpoint_evaluations=0,new_middle_algorithm3_jobs=4,
    boundaries={'Gaston':'first L omitted; first chi input to last chi output',
                'Hsilu':'first chi input to last chi output; last L absorbed into the external output mask'},
    model_limits=['Algorithm 3 expands every beta; remaining formal alpha variables are assumed independent uniform.',
                  'Spliced p*c_m*q^2 is a model estimate; it is not an exact endpoint correlation or a proven bound.',
                  'A resource limit has no numerical correlation; a model zero is not proof of actual zero.'],rows=rows)
dump(OUT/'results.json',data)

def sample_html(r):
    s=r.get('sampling')
    if not s:return '<span class="muted">未采样验证</span>'
    def v(x):return f'{x["rho"]:+.12f}' if abs(x['rho'])!=1 else f'{x["rho"]:+.0f}'
    text=f'<b>{v(s)}</b><small>{pow2(s["samples"])} 对</small>'
    if r.get('sampling_reference'):
        t=r['sampling_reference']
        text=f'核心：{v(s)}<br>参考：{v(t)}<small>各 {pow2(s["samples"])} 对；尚未独立确认</small>'
    return text

def table(algorithm):
    tr=[]
    for r in rows:
        if r['algorithm']!=algorithm:continue
        code=r['code'];sp='+'.join(map(str,r['evaluation_split']))
        ends=f'<span class="endpoint">Δ={esc(bits(r["delta_core"]))}</span><span class="endpoint">Λ={esc(bits(r["lambda_core"]))}</span>'
        if r['evaluation_scope']=='whole':
            scope=f'<span class="pill whole">整体</span><small>{sp}</small>'
            pq='<span class="muted">不使用拼接公式</span>'
            d=r['datf_algorithm3']
            dc=number(d['correlation'],True) if d['status']=='completed' else '<span class="incomplete">资源上限<br>未完成</span>'
            dc=link(d['source'],'原始记录')+'<br>'+dc
            hc=number(r['hnw26']['correlation'],True)
            score='<span class="muted">—</span>'
        else:
            scope=f'<span class="pill middle">仅中间段</span><small>{sp}</small>'
            p='两支，各 p=2<sup>−2</sup>' if code=='G4c' else ('p=2<sup>−2</sup>' if code=='G4d' else 'p=1（无差分外段）')
            pq=p+f'<small>q²={number(r["linear_squared_correlation"]).lstrip("+")}</small>'
            def midcell(field):
                values=[b[field]['correlation'] for b in r['branches']]
                assert all(frac(x)==frac(values[0]) for x in values)
                s='c<sub>m</sub>='+number(values[0],True)
                if len(values)>1:s+='<small>两支分别均为此值</small>'
                else:s+='<small>中间段相关度</small>'
                return s
            dc=midcell('datf_algorithm3');hc=midcell('hnw26')
            score=number(r['datf_spliced_correlation'],True)+'<small>两方法拼接值一致</small>'
        tr.append(f'<tr data-code="{code}"><td><a href="#{code}"><b>{code}</b></a><small>{r["rounds"]} 轮</small></td><td>{ends}</td><td>{scope}</td><td>{pq}</td><td>{dc}</td><td>{hc}</td><td>{score}</td><td>{sample_html(r)}</td></tr>')
    return '<div class="tablewrap"><table class="main"><thead><tr><th>编号 / 轮数</th><th>整体端点（core 坐标）</th><th>评估区间<br>(r<sub>d</sub>,r<sub>m</sub>,r<sub>l</sub>)</th><th>外段因子</th><th>DATF(JOC算法3)</th><th>HNW26</th><th>拼接整体估计 ρ̂</th><th>已有实验相关度</th></tr></thead><tbody>'+''.join(tr)+'</tbody></table></div>'

def details(r):
    code=r['code'];parts=[f'<details id="{code}"><summary>{code} · {r["algorithm"]} {r["rounds"]} 轮：端点、区间和证据</summary>',
        f'<p>{esc(r["discovery"])}</p>',f'<p>整体起始轮索引：{r["start_round_index"]}（从 0 开始）；每个十六进制字按行 0–4 排列，位号从最低位 0 起算。</p>',
        '<div class="wordgrid"><b>整体输入差分 Δ</b><pre>'+ '\n'.join(r['delta_core'])+'</pre><b>整体输出掩码 Λ</b><pre>'+ '\n'.join(r['lambda_core'])+'</pre></div>']
    if r.get('lambda_full'):
        parts+=['<p>Hsilu 含最后线性层的等价输出掩码：</p><pre>'+'\n'.join(r['lambda_full'])+'</pre>']
    if r['evaluation_scope']=='whole':
        d=r['datf_algorithm3']
        parts+=['<p>本次复用已有整条路径评估记录：'+link(d['source'],'JOC 算法 3 原始记录')+'；'+link(ROOT/'研究汇总_20261019/results_catalog.json','HNW26 与采样归档')+'。</p>']
        if d['status']!='completed':
            parts+=['<p class="note">算法 3 未产生数值，保留“未完成”。选定归档尝试的终止原因：<code>'+esc(d['reason'])+'</code>。已有尝试记录如下。</p><ul>']
            for a in d['attempts']:
                parts+=['<li>'+link(a['file'],Path(a['file']).parent.name+'/'+Path(a['file']).name)+f'：{a["elapsed_seconds"]:.3f} 秒，'+esc(a['reason'])+'</li>']
            parts+=['</ul>']
        else:
            parts+=[f'<p>算法 3：{number(d["correlation"])}；归档耗时 {d["elapsed_seconds"]:.3f} 秒。</p>']
    else:
        parts+=['<p>本次仅评估以下中间段。尾段：<code>'+esc(r['tail_id'])+'</code>；'+link(r['cut_source'],'原搜索和切分记录')+'。</p>']
        for i,b in enumerate(r['branches']):
            j=b['middle_job'];d=b['datf_algorithm3'];h=b['hnw26']
            parts+=[f'<h4>中间段 {i+1}：{j["rounds"]} 轮 χ，起始索引 {j["start_round_index"]}</h4>',
                f'<p>p={number(b["probability"]).lstrip("+")}；DATF(JOC算法3) c<sub>m</sub>={number(d["correlation"])}（{d["elapsed_seconds"]:.4f} 秒）；HNW26 c<sub>m</sub>={number(h["correlation"])}（{h["elapsed_seconds"]:.4f} 秒）。</p>',
                '<div class="wordgrid"><b>中间输入差分</b><pre>'+'\n'.join(j['delta_core'])+'</pre><b>中间输出掩码</b><pre>'+'\n'.join(j['lambda_core'])+'</pre></div>',
                '<p>'+link(d['source'],'DATF 原始结果')+' · '+link(h['source'],'HNW26 原始结果')+'</p>']
        if code in ['G4c','G4d']:
            a=validation['first_round_branch_audits'][code]
            parts+=['<p>全部第一轮 DDT 分支核对（仅计算各分支的一轮中间段）：</p><ul>']
            for b in a['branches']:
                parts+=[f'<li>首轮输出列模式 {b["first_output_pattern"]}：p={number(b["probability"]).lstrip("+")}，c<sub>m</sub>={number(b["middle_correlation"])}。</li>']
            parts+=['</ul><p>Σ p<sub>i</sub>c<sub>m,i</sub>='+number(a['signed_sum'])+'；再乘同一线性尾段 q²=2<sup>−38</sup>。这些因子相乘是分段模型估计。</p>']
        else:
            parts+=['<p class="note">H5 沿用原搜索的 0+3+2：没有独立归档的差分外段，前三轮整体作为中间段。此次三轮中间段完整展开了 72 个 β 定义，得到 1,228 项、143 个形式变量、四次多项式；Bias 消元得到模型值 0。输入变量先做精确 Shannon 分解的另一求值顺序也得到 0。HNW26 full 与 zero 均覆盖全部 8 个首轮 DDT 分支，结果也为 0。</p>',
                '<p>'+link(OUT/'raw/H5_prefix3_DATF_shannon.json','DATF 第二求值顺序复核')+' · '+link(OUT/'raw/H5_prefix3_HNW_zero.json','HNW26 zero 复核')+'</p>',
                '<p>旧前缀估计 −2<sup>−48.959588</sup>、旧拼接值 −2<sup>−58.959588</sup>来自此前的自定义估计器，不能作为 DATF(JOC算法3) 的结果。新分段值 0 只针对当前切分与固定尾段，不能据此证明真实五轮端点相关度为 0。</p>']
    for label,s in [('核心采样',r.get('sampling')),('参考实现采样',r.get('sampling_reference'))]:
        if not s:continue
        ci=s.get('ci95')
        cis=f'[{ci[0]:+.12f}, {ci[1]:+.12f}]' if ci else '未单独归档（H3 另有结构证据）'
        parts+=[f'<p>{label}：N={pow2(s["samples"])}，n₀={s["n0"]:,}，n₁={s["n1"]:,}，ρ̂={s["rho"]:+.12f}；归档 95% 区间 {cis}；seed={esc(s["seed"])}。</p>']
    parts+=['</details>']
    return ''.join(parts)

css='''
:root{--ink:#18334a;--muted:#607386;--line:#d9e4ec;--blue:#175a8c;--bg:#f3f7fa}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.65 system-ui,"Noto Sans CJK SC","Microsoft YaHei",sans-serif}main{max-width:1560px;margin:0 auto;padding:36px 28px 70px}header{padding:28px 32px;background:#12354e;color:white;border-radius:16px}header p{max-width:1000px;color:#deebf3}h1{font-size:29px;line-height:1.35;margin:8px 0 18px}h2{font-size:24px;margin:34px 0 13px}h3{font-size:18px}h4{margin-bottom:8px}a{color:var(--blue);text-decoration:none}a:hover{text-decoration:underline}.tag{font-size:13px;letter-spacing:1px;color:#aad1ea}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:15px;margin:20px 0}.card{background:white;border:1px solid var(--line);padding:20px;border-radius:12px}.card strong{font-size:22px;display:block}.note{padding:15px 18px;border-left:4px solid #d39939;background:#fff8e9;border-radius:4px}.formula{font-size:20px;text-align:center;margin:20px 0;padding:18px;background:white;border:1px solid var(--line);border-radius:10px}.tablewrap{overflow-x:auto;background:white;border:1px solid var(--line);border-radius:12px}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}.main{min-width:1240px}th{background:#e7eff5;color:#163b57;text-align:left;position:sticky;top:0;z-index:1;font-size:13px}th,td{padding:15px 12px;border-bottom:1px solid var(--line);vertical-align:top}td{font-size:14px}tr:last-child td{border-bottom:0}.main th:nth-child(2){width:22%}.main th:nth-child(8){width:19%}.main td:nth-child(5),.main td:nth-child(6),.main td:nth-child(7){font-size:16px}small{display:block;font-size:12px;color:var(--muted);line-height:1.5;margin-top:5px}.endpoint{display:block;font:12.5px/1.7 ui-monospace,monospace;overflow-wrap:anywhere;margin-bottom:5px}.pill{display:inline-block;font-size:12px;padding:2px 9px;border-radius:20px;white-space:nowrap}.whole{background:#dfefe7;color:#245b40}.middle{background:#e2ecfa;color:#284e8b}.incomplete{color:#925f18;font-size:14px}.muted{color:var(--muted);font-size:12px}details{background:white;border:1px solid var(--line);border-radius:10px;padding:0 20px;margin:12px 0;scroll-margin-top:16px}summary{font-weight:650;padding:18px 0;cursor:pointer}details p{max-width:1100px}.wordgrid{display:grid;grid-template-columns:170px 1fr;gap:10px;align-items:start}pre{background:#f0f5f8;border-radius:6px;padding:12px;margin:0;overflow:auto;font:13px/1.6 ui-monospace,monospace}code{font-family:ui-monospace,monospace;overflow-wrap:anywhere}.links{line-height:2.2}.foot{color:var(--muted);font-size:13px;margin-top:28px}@media(max-width:700px){main{padding:16px 10px}header{padding:22px 18px}.cards{grid-template-columns:1fr}h1{font-size:25px}.wordgrid{grid-template-columns:1fr}}@media print{body{background:white}main{padding:0}header{color:#12354e;background:#edf3f7}.cards{display:none}.tablewrap{overflow:visible}.main{min-width:0}th,td{font-size:10px;padding:7px}.main td:nth-child(n){font-size:10px}.endpoint{font-size:9px}h2{break-after:avoid}tr{break-inside:avoid}details{break-inside:avoid}}
'''
page='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston 与 Hsilu：整体与分段评估结果</title><style>'''+css+'''</style></head><body><main>
<header><div class="tag">研究结果 · 2026-10-27 · 13 条代表端点</div><h1>Gaston 与 Hsilu<br>差分线性结果：整体评估与分段评估</h1><p>已得到可靠采样信号的代表路径保留整条路径评估；尚未得到可靠实验确认的代表路径，仅评估中间段，并用外段因子计算拼接评分。DATF 一律采用 JOC 算法 3 的口径。</p></header>
<div class="cards"><div class="card"><strong>10 条整体评估</strong>Gaston 2–4 轮、Hsilu 3–4 轮；复用已完成的计算尝试与采样。</div><div class="card"><strong>3 条分段评估</strong>G4c、G4d、H5；本次补算 4 个中间段任务，不重算这些端点的整条置换。</div><div class="card"><strong>6 项算法 3 未完成</strong>保留资源上限状态，既不填 0，也不沿用旧估计器的数值。</div></div>
<p>本表数值是带符号相关度 <b>ρ=2Pr[Λ·(F(X)⊕F(X⊕Δ))=0]−1</b>，不是偏差 ε；两者关系为 ρ=2ε。整体行中 DATF 与 HNW26 列是整体 ρ；“仅中间段”行中这两列是 c<sub>m</sub>，右侧另列拼接整体估计。</p>
<div class="formula">单分支：ρ̂ ≈ p · c<sub>m</sub> · q²　　同一尾段多分支：ρ̂ ≈ (Σ<sub>i</sub> p<sub>i</sub>c<sub>m,i</sub>) · q²</div>
<p>轮数按 χ 层数计，R=r<sub>d</sub>+r<sub>m</sub>+r<sub>l</sub>。整体评估标为 (0,R,0)，这表示本次评估范围，不改变原来的搜索分段。p 为外部差分路径概率，q 为所保留线性尾路径的相关度；q² 使用尾路径的 LAT 乘积平方。</p>
<p class="note"><b>证据口径：</b>拼接值是模型估计，不是完整端点的精确相关度或上下界；HNW26 的 0 只表示模型给出 0。若实验已观察到非零信号，则保留这种模型与实验的差异。JOC 算法 3 完整代回 β，仍对剩余形式 α 变量采用联合独立均匀假设。</p>
<h2>Gaston</h2>'''+table('Gaston')+'''
<p><b>G4c：</b>两个第一轮差分事件各有 p=1/4，中间段分别为 +1；Σp<sub>i</sub>c<sub>m,i</sub>=1/2，乘 q²=2<sup>−38</sup> 得 +2<sup>−39</sup>。其余两条首轮 DDT 分支的中间值为 0。这里不能把 p=1/2 表述为一条差分路径的概率。</p>
<p><b>G4d：</b>选定首轮差分事件 p=1/4，中间段 −1，q²=2<sup>−38</sup>，得到 −2<sup>−40</sup>；其余三条首轮 DDT 分支的中间值为 0。</p>
<h2>Hsilu</h2>'''+table('Hsilu')+'''
<p><b>H5：</b>原搜索是三轮前缀接固定两轮线性尾段，因此本次沿用 <b>0+3+2</b>：p=1，q²=2<sup>−10</sup>。仅三轮前缀的 JOC 算法 3、HNW26 full 与 zero 都为 0，当前切分下的拼接值也为 0。两批 2<sup>28</sup> 对采样符号不同，参考实现那一批的置信区间包含 0，仍标为“尚未独立确认”。</p>
<h2>边界与计算记录</h2><p>所有端点沿用历史结果的 core 坐标，范围从第一个 χ 的输入到最后一个 χ 的输出。Gaston 省略首个线性层；Hsilu 的原生轮顺序不同，其最后线性层由输出掩码的逆转置变换吸收。未更改已有输入差分、输出掩码或轮常数窗口。Hsilu 原始 worker 的通用字段 first_L_omitted 沿用了 Gaston 的命名，本报告以明确的 core 边界及具体置换实现为准。</p>
<p>一轮中间段采用 HNW26 双线性 LAT 算子的闭式形式（精确局部 DLCT 乘积），避免对稠密输入差分做无必要的首层 DDT 笛卡尔积枚举；该公式已对所有 32×32 个局部输入差分/输出掩码核对。H5 三轮中间段使用现有 HNW26 可执行程序，分别运行 full 和 zero，完整枚举 8 个首层 DDT 分支。所有新中间段都完成了 JOC 算法 3 的全部 β 展开；H5 还以精确输入 Shannon 分解的另一求值顺序交叉复核。</p>
'''+''.join(details(r) for r in rows)+'''
<h2>复现与来源</h2><div class="links">'''+ ' · '.join([
 link(OUT/'results.json','整理后的机器可读结果'),link(OUT/'middle_manifest.json','中间段端点与轮常数窗口'),
 link(OUT/'validation.json','分支与尾段核对记录'),link(OUT/'prepare.py','准备脚本'),
 link(OUT/'evaluate_middle.py','中间段评估脚本'),link(OUT/'build_report.py','报告生成脚本'),
 link(ROOT/'研究汇总_20261019/gaston_hsilu_results.html','历史汇总'),
 link(AUDIT/'datf_algorithm3_report.html','历史 JOC 算法 3 复核')])+'''</div>
<p>复现顺序：在本报告目录运行 <code>python3 prepare.py</code>，再运行 <code>python3 evaluate_middle.py</code>，最后运行 <code>python3 build_report.py</code>。加 <code>--reuse</code> 可复用本次已归档的中间段原始结果。历史整体评估与采样文件保持原样。</p>
<p class="foot">本次没有新增采样，也没有重算未确认候选的完整置换。计算耗时为各程序记录的计算区间，不包括进程启动和建表的全部时间。完整来源与校验和见 source_manifest.json。</p>
</main></body></html>'''
(OUT/'gaston_hsilu_segmented_results.html').write_text(page)

sources={ROOT/'研究汇总_20261019/results_catalog.json',AUDIT/'results.json',
    ROOT/'gaston_dl/results/fixed_paper_tail/search.json',ROOT/'gaston_dl/results/datf_pools/paper_full24/search_results.json',
    ROOT/'hsilu_dl/results/research_summary.json',AUDIT/'algorithm3.py',AUDIT/'worker.py',
    ROOT/'gaston_dl/tools/datf_k06.py',ROOT/'hsilu_dl/tools/hsilu.py',ROOT/'hsilu_dl/src/hnw_endpoint.cpp',
    ROOT/'hsilu_dl/build/hnw_endpoint',ROOT/'gaston_dl/tools/fixed_tail_search.py',
    OUT/'prepare.py',OUT/'evaluate_middle.py',OUT/'build_report.py',OUT/'cuts.json',OUT/'middle_manifest.json'}
sources.update(OUT.glob('raw/*.json'))
for r in rows:
    if r['evaluation_scope']=='whole':sources.add(Path(r['datf_algorithm3']['source']))
dump(OUT/'source_manifest.json',dict(files=[dict(path=str(p),bytes=p.stat().st_size,
    sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(sources)]))

# Validate links, table row identities, column labels, and resource-limit semantics.
from html.parser import HTMLParser
class Check(HTMLParser):
    def __init__(self):super().__init__();self.codes=[];self.links=[];self.tables=0;self.ids=[]
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if tag=='tr' and 'data-code' in a:self.codes.append(a['data-code'])
        if tag=='a':self.links.append(a.get('href',''))
        if 'id' in a:self.ids.append(a['id'])
        if tag=='table':self.tables+=1
checker=Check();checker.feed(page)
assert checker.codes==order and checker.tables==2
assert len(checker.ids)==len(set(checker.ids))
from urllib.parse import unquote
for href in checker.links:
    if href.startswith('/'):assert Path(unquote(href)).exists(),href
    if href.startswith('#'):assert href[1:] in checker.ids,href
assert page.count('<th>DATF(JOC算法3)</th>')==2
dump(OUT/'report_validation.json',dict(status='PASS',main_tables=2,representatives=13,
    whole_evaluation=10,middle_only_evaluation=3,whole_algorithm3_completed=4,
    whole_algorithm3_resource_limited=6,new_middle_algorithm3_completed=4,
    local_links_checked=len(checker.links),legacy_datf_values_used_as_algorithm3=False))
print('PASS: two tables, 13 representatives; 4 completed new middle jobs; links and arithmetic checked.')

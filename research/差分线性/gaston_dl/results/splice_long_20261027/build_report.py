import json,sys,math,gzip,html,hashlib
from pathlib import Path
from fractions import Fraction
from urllib.parse import quote,unquote
from html.parser import HTMLParser
H=Path(__file__).resolve().parent;G=H.parents[1]
sys.path.insert(0,str(H));from prepare import read,save,state,ws,rec,fraction
sys.path.insert(0,str(G/'tools'))
from fixed_tail_search import S,column,linear,rotate

def ddt(a,b):return math.prod(Fraction(sum(S[x]^S[x^column(a,j)]==column(b,j) for x in range(32)),32) for j in range(64))
def esc(s):return html.escape(str(s))
def link(p,s):return f'<a href="{quote(str(p),safe="/")}">{esc(s)}</a>'
def fmt(x):
    q=fraction(x) if isinstance(x,dict) else Fraction(x)
    if not q:return '0'
    sign='−' if q<0 else '+';q=abs(q)
    if q==1:return sign+'1'
    if q.numerator==1 and q.denominator&(q.denominator-1)==0:return f'{sign}2<sup>−{q.denominator.bit_length()-1}</sup>'
    return sign+str(q)
def bits(a):return '{'+', '.join(f'({i},{j})' for i,w in enumerate(a) for j in range(64) if int(w,16)>>j&1)+'}'

def main():
    data=read(H/'candidates.json');rows=data['rows'];jobs={j['id']:j for j in read(H/'manifest.json')['jobs']}
    meta=read(G/'results/datf_pools/pools/pools.json')
    wanted={r['extra']['original_candidate']['differential_id'] for r in rows if r['family']=='pool_single'}
    pool={}
    for key in ['differential_one','differential_pq_two']:
        with gzip.open(G/'results/datf_pools/pools'/meta['pools'][key]['file'],'rt') as f:
            for line in f:
                x=json.loads(line)
                if x['id'] in wanted:pool[x['id']]=x
    assert set(pool)==wanted
    branch_checks=0
    for r in rows:
        # Recompute every archived outer DDT product and all difference joins.
        if r['family']=='pool_single':
            c=r['extra']['original_candidate'];d=pool[c['differential_id']];rot=c['differential_rotation']
            def ds(k):return rotate(state(d[k]),rot)
            assert ds('input')==state(r['input_difference'])
            assert ds('output')==state(c['differential_output'])
            if d['rounds']==1:p=ddt(ds('input'),ds('output'))
            else:
                assert linear(ds('first_output'))==ds('second_input')
                p=ddt(ds('input'),ds('first_output'))*ddt(ds('second_input'),ds('output'))
            assert linear(ds('output'))==state(jobs[r['branches'][0]['job_id']]['delta_core'])
            assert p==fraction(r['branches'][0]['probability']);branch_checks+=1
        elif r['family']=='single_column_cluster':
            for b,s in zip(r['branches'],r['extra']['candidate']['prefix']['branches']):
                dout=state(s['first_chi_output_difference'])
                assert ddt(state(r['input_difference']),dout)==fraction(b['probability'])
                assert linear(dout)==state(jobs[b['job_id']]['delta_core']);branch_checks+=1
        elif r['family']=='pq_cluster':
            e=r['extra'];delta=state(r['input_difference']);first=state(e['first_output']);second=state(e['second_input'])
            assert linear(first)==second
            for b,c in zip(r['branches'],e['connections']):
                assert ddt(delta,first)*ddt(second,state(c['b2']))==fraction(b['probability'])
                assert linear(state(c['b2']))==state(jobs[b['job_id']]['delta_core']);branch_checks+=1
        else:
            c=r['extra']['candidate'];rd=r['split'][0];delta=state(r['input_difference'])
            if rd==0:p=Fraction(1);last=delta
            elif rd==1:p=ddt(delta,state(c['first_output']));last=linear(state(c['first_output']))
            else:
                assert linear(state(c['first_output']))==state(c['second_input'])
                p=ddt(delta,state(c['first_output']))*ddt(state(c['second_input']),state(c['second_output']));last=linear(state(c['second_output']))
            assert p==fraction(r['branches'][0]['probability'])
            assert last==state(jobs[r['branches'][0]['job_id']]['delta_core']);branch_checks+=1
        for method,suffix in [('datf_algorithm3','datf'),('hnw26','hnw')]:
            total=Fraction(0)
            for b in r['branches']:
                p=H/'raw'/f"{b['job_id']}_{suffix}.json";raw=read(p)
                assert raw['status']=='completed'
                b[method]=dict(correlation=raw['correlation'],elapsed_seconds=raw['elapsed_seconds'],source=str(p))
                if suffix=='hnw':
                    assert any(v is True for k,v in raw['certificate'].items() if k.startswith('full_equals_zero_certified'))
                total+=fraction(b['probability'])*fraction(raw['correlation'])
            r[method+'_score']=rec(total*fraction(r['q_squared']))
        assert r['datf_algorithm3_score']==r['hnw26_score']
        r['heuristic_sample_exponent']=-2*r['datf_algorithm3_score']['log2_abs'] if fraction(r['datf_algorithm3_score']) else None
    best=[max((r for r in rows if r['rounds']==R),key=lambda r:abs(fraction(r['datf_algorithm3_score']))) for R in range(5,10)]
    assert [r['datf_algorithm3_score']['log2_abs'] for r in best]==[-59,-83,-116,-140,-164]
    # Records with no archived long-round match: absence is scoped to the old search.
    uncovered=[]
    for path,key in [('datf_pools/linear_two_match/search_results.json','rows'),('datf_pools/search_final/results.json','summaries')]:
        for row in read(G/'results'/path)[key]:
            if row.get('best') is not None:continue
            if row.get('tail') not in ['L1_pool','L2_2to2_pool']:continue
            uncovered.append(dict(source=str(G/'results'/path),family=row['family'],middle=row['middle_rounds'],tail=row['tail'],
                archived_exhaustive=row.get('exhaustive',False),status='no nonzero candidate archived under the previous search model and scope'))
    result=dict(date='2026-10-27',first_L_omitted=True,scope=data['scope'],
        datf_column_name='DATF(JOC算法3)',middle_only=True,new_full_endpoint_jobs=0,
        linear_inventory=data['linear_inventory'],short_linear_pools={k:v for k,v in meta['pools'].items() if k.startswith('linear')},
        no_archived_long_match=uncovered,rows=rows,best_per_round=[r['id'] for r in best],
        all_middle_methods_agree=True,all_selected_two_chi_full_zero_equivalences_certified=True,
        heuristic_samples_note='N ~ |rho|^-2 is only a model-based scale, not a proven attack complexity; R9 exceeds unique input-pair capacity.')
    save(H/'results.json',result)
    save(H/'validation.json',dict(status='PASS',candidates=len(rows),middle_jobs=len(jobs),outer_DDT_branches_checked=branch_checks,
        tail_LAT_and_transpose_checks=len(data['linear_inventory']),completed_DATF_middle_jobs=len(jobs),completed_HNW_middle_jobs=len(jobs),
        max_rounds=9,full_endpoint_jobs=0,best_log2_abs=[r['datf_algorithm3_score']['log2_abs'] for r in best]))

    def tr(r,compact=False):
        n=len(r['branches']);p=r['branches'][0]['probability'];cms=[b['datf_algorithm3']['correlation'] for b in r['branches']]
        assert all(c==cms[0] for c in cms)
        prob=fmt(p).lstrip('+')+(f'<small>{n} 个同端点分支，各自概率</small>' if n>1 else '')
        cm=fmt(cms[0])+(f'<small>{n} 支均如此</small>' if n>1 else '')
        return f'<tr><td>{r["rounds"]}</td><td>{"+".join(map(str,r["split"]))}</td><td><code>{r["tail_id"]}</code></td><td>{prob}</td><td>{cm}</td><td>{cm}</td><td>{fmt(r["q_squared"]).lstrip("+")}</td><td><strong>{fmt(r["datf_algorithm3_score"])}</strong></td><td>'+('2<sup>'+str(round(r['heuristic_sample_exponent']))+'</sup>')+'</td></tr>'
    header='<tr><th>总轮数</th><th>差分+中间+线性</th><th>线性尾段</th><th>p</th><th>DATF(JOC算法3)<small>中间相关度 cₘ</small></th><th>HNW26<small>中间相关度 cₘ</small></th><th>q²</th><th>拼接估计 ρ̂</th><th>形式样本尺度 |ρ̂|⁻²</th></tr>'
    main_table='<div class="scroll"><table><thead>'+header+'</thead><tbody>'+''.join(tr(r) for r in best)+'</tbody></table></div>'
    inv='<div class="scroll"><table><thead><tr><th>尾段</th><th>线性轮数</th><th>q²</th></tr></thead><tbody>'+''.join(f'<tr><td>{t["id"]}</td><td>{t["rounds"]}</td><td>{fmt(t["q_squared"]).lstrip("+")}</td></tr>' for t in data['linear_inventory'])+'</tbody></table></div>'
    det=[]
    for r in best:
        det+=[f'<details id="R{r["rounds"]}"><summary>{r["rounds"]} 轮代表：{esc(r["id"])}</summary>',
            f'<p>整体从索引 {r["start_round_index"]} 开始，省略第一个 L；本次只算中间 {r["split"][1]} 个 χ。所保留的线性路径固定为 {esc(r["tail_id"])}。</p>',
            '<p>整体输入差分 Δ：'+esc(bits(r['input_difference']))+'</p><p>整体输出掩码 Λ：'+esc(bits(r['output_mask']))+'</p>',
            '<pre>Δ = '+'; '.join(r['input_difference'])+'\nΛ = '+'; '.join(r['output_mask'])+'</pre>']
        for i,b in enumerate(r['branches']):
            j=jobs[b['job_id']]
            det+=[f'<p>分支 {i+1}：中间起始索引 {j["start_round_index"]}；p={fmt(b["probability"])}；cₘ={fmt(b["datf_algorithm3"]["correlation"])}。</p>',
                '<pre>中间 Δ = '+'; '.join(j['delta_core'])+'\n中间 Λ = '+'; '.join(j['lambda_core'])+'</pre>',
                '<p>'+link(Path(b['datf_algorithm3']['source']),'DATF(JOC算法3) 原始记录')+' · '+link(Path(b['hnw26']['source']),'HNW26 原始记录')+'</p>']
        det+=['<p>'+link(G/'results'/r['source'],'历史候选来源')+'</p></details>']
    controls=[r for r in rows if r['id'] in ['cluster1_D3_chi2_to_4','cluster2_D3_chi2_to_3','D2_M1_D3_chi2_to_4','D2_M2_D3_chi3_to_4','cluster1_D2_repeat_4','cluster2_D2_repeat_6']]
    text='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston 5–9 轮分段拼接复核</title><style>
body{margin:0;background:#f4f7fa;color:#173348;font:15px/1.7 system-ui,"Noto Sans CJK SC",sans-serif}main{max-width:1440px;margin:auto;padding:30px 24px 60px}header{background:#153f5b;color:white;border-radius:14px;padding:26px 30px}h1{margin:5px 0 15px;font-size:29px}h2{margin-top:30px;font-size:23px}a{color:#176c9e}table{width:100%;border-collapse:collapse;background:white;font-variant-numeric:tabular-nums}td,th{padding:12px 10px;border-bottom:1px solid #dce5ed;vertical-align:top;text-align:left}th{background:#e7eff5;font-size:13px}small{display:block;color:#60778b;font-size:12px}code{font-size:12px;overflow-wrap:anywhere}.scroll{overflow:auto;border:1px solid #d6e2eb;border-radius:10px}.note{padding:15px 20px;border-left:4px solid #b77c22;background:#fff7e7}.formula{background:white;padding:18px;text-align:center;border-radius:8px;font-size:20px}details{margin:12px 0;background:white;border:1px solid #d6e2eb;border-radius:9px;padding:0 18px}summary{padding:17px 0;cursor:pointer;font-weight:650}pre{white-space:pre-wrap;word-break:break-all;background:#eff4f8;padding:12px;border-radius:8px;font:13px/1.8 monospace}.muted{color:#5f7688}.links{line-height:2.1}@media(max-width:700px){main{padding:16px 10px}header{padding:20px 18px}table{min-width:980px}h1{font-size:25px}}
</style></head><body><main><header><div>2026-10-27 · 已有路径池复核</div><h1>Gaston：使用分段拼接能延长到多少轮？</h1><p>现有路径池已有 5–9 轮非零拼接候选。本次复核 49 条代表候选、76 个不同的一轮或两轮中间段。DATF(JOC算法3) 与 HNW26 的中间段结果全部一致；没有计算这些长轮端点的完整 DATF/HNW26，也没有新增实验采样。</p></header>
<p>范围为历史各搜索分组的最优归档候选，以及已归档的同端点分支簇。这是对已有候选的重新核对，不是对整个搜索空间的重新穷举，也不宣称下表为全局最优。所有状态均为 5×64=320 比特，省略第一个线性层，末轮停在 χ 输出，坐标按 (行,最低位起算的比特位置)。</p>
<div class="formula">ρ̂ ≈ p · cₘ · q²　或　ρ̂ ≈ (Σᵢ pᵢcₘ,ᵢ) · q²</div>
<h2>本次复核候选中，每个轮数的较大拼接评分</h2>'''+main_table+'''
<p>5 轮使用 1+2+2：2<sup>−4</sup>×2<sup>−13</sup>×2<sup>−42</sup>=2<sup>−59</sup>。6 轮使用 K06 的两个同端点分支，各自 p=2<sup>−24</sup>、cₘ=−1，故总分为 −2<sup>−23</sup>×2<sup>−60</sup>=−2<sup>−83</sup>；单个 K06 分支则为 −2<sup>−84</sup>。这些和式只包含归档保留的分支，不是完整相关度的全部贡献之和。</p>
<p>7–9 轮使用同一类 1+2+n 连接：p=2<sup>−2</sup>，cₘ=2<sup>−18</sup>，D2_repeat_n 的 q²=2<sup>−24n</sup>。因此 ρ̂=2<sup>−(20+24n)</sup>，n=4,5,6 分别给出 7、8、9 轮。中间段在各自轮常数窗口下都已重新计算。</p>
<p class="note"><b>九轮的含义需要限定。</b>2<sup>−164</sup> 是非零分段模型评分，不能直接称为有效九轮区分器。常用形式样本尺度为 2<sup>328</sup>；对于固定非零输入差分，320 比特状态只有 2<sup>319</sup> 个不同的无序输入对，该尺度已超出可用的不重复输入对数量。还需考虑理想置换的背景相关度波动，不能把这个公式当作已证明的区分复杂度。5–8 轮的表中尺度也只是拼接模型的推算，尚无实验确认。</p>
<h2>与原有一轮中间段组合对照</h2><div class="scroll"><table><thead>'''+header+'''</thead><tbody>'''+''.join(tr(r) for r in controls)+'''</tbody></table></div>
<p>只看 1+1+n，会遗漏两轮中间段提供的较好评分；而给差分前缀多加一轮也未必更好，要同时比较 p、cₘ 和 q²。当前 6 轮受益于 6→6 的 K06 前缀，7–9 轮则由一轮差分配两轮中间段的组合取得更大的保留评分。</p>
<h2>已有线性池及局限</h2><p>单轮线性池有 12,730 个规范代表；L2_2to2_pool 有 16,000 个规范代表，其中最小线性路径权重为 4（q²=2<sup>−8</sup>）。论文导入尾段共有下列 11 个，最长 6 轮。</p>'''+inv+'''
<p>L2_2to2_pool 的已归档结果中，D1+M2+L2 的五轮组合以及 D2+M2+L2 的六轮组合在当时的有限批次中没有非零匹配；D2+M1+L2 的五轮组合在当时指定差分池的完整枚举中没有匹配。L1 的 D2+M2+L1 五轮组合也未归档非零候选。这里没有把旧有限搜索的“未发现”升级成全局不存在，也没有用这些零记录替代 JOC 算法 3 的新搜索。</p>
<p>如果继续改进，应先扩大 L2_2to2_pool 配两轮中间段的候选覆盖，争取低代价五轮连接；同时为更长轮数建立比 D2_repeat_n 更低权重的 3–6 轮线性尾段。重复 D2 尾段每增加一轮，q² 就再损失 24 个二进制指数，因此单纯加长该尾段很快削弱评分。这里提出的是下一步方向，本次未运行新的全池搜索。</p>
<h2>中间段核对方法</h2><p>DATF(JOC算法3) 使用此前的原始算法 3 实现，全部代回 β，再对剩余独立均匀形式变量求期望。76 个任务均完成。HNW26 使用零值掩码模型，保留首层 DDT 完整求和；一轮中间段为精确局部 DLCT 乘积。两轮中间段利用 χ 导数的仿射性，将首层 DDT 的求和转化为线性方程约束下的二次符号求和，不增加新的截断或独立性假设。</p>
<p>该加速求和与现有原生 HNW26 程序做了 8 组交叉核对，包含 0、1/2、1/64 等结果。对本次全部两轮中间段，还逐项验证了第一个 L 后所有非零值掩码的局部分量均被消去，所以这些端点的 HNW26 full 与 zero 在该局部模型中相同。两种模型相同不等于真实整体相关度已得到验证。</p>
<h2>代表路径的完整端点与中间边界</h2>'''+''.join(det)+'''
<h2>原始数据与复现</h2><div class="links">'''+ ' · '.join(link(H/f,s) for f,s in [('results.json','全部 49 条候选及评估'),('manifest.json','76 个中间段任务'),('validation.json','外段与拼接核对'),('hnw_validation.json','HNW26 原生交叉验证'),('prepare.py','候选提取'),('run_datf.py','算法 3 批处理'),('hnw_middle.py','HNW26 局部求值'),('build_report.py','报告生成')])+'''</div><p>复现顺序：prepare.py → run_datf.py → hnw_middle.py → build_report.py。结果保留单条路径与分支簇的区别，任何拼接评分都不作为完整端点的精确值或上下界。</p></main></body></html>'''
    page=H/'gaston_long_splice_report.html';page.write_text(text)
    class Check(HTMLParser):
        def __init__(self):super().__init__();self.links=[]
        def handle_starttag(self,tag,attrs):
            if tag=='a':self.links.append(dict(attrs).get('href',''))
    c=Check();c.feed(text)
    for href in c.links:
        if href.startswith('/'):assert Path(unquote(href)).exists(),href
    sources=[G/'results/datf_pools/pools/pools.json',G/'results/fixed_paper_tail/search.json',G/'results/pq_kernel_six_all/search.json',
        G/'results/datf_pools/paper_full24/search_results.json',G/'results/fixed_tail_extension/results.json',
        G/'results/datf_pools/linear_two_match/search_results.json',G/'results/datf_pools/search_final/results.json',
        G.parent/'DATF算法3复核_20261021/algorithm3.py',G.parent/'DATF算法3复核_20261021/worker.py',G/'src/hnw_endpoint.cpp']
    sources+=list(H.glob('*.py'))+list((H/'raw').glob('*_datf.json'))+list((H/'raw').glob('*_hnw.json'))
    save(H/'source_manifest.json',dict(files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]))
    print('PASS',len(rows),'candidates',len(jobs),'middle jobs,',branch_checks,'outer DDT checks; best scores:')
    for r in best:print(r['rounds'],r['id'],r['datf_algorithm3_score'])
if __name__=='__main__':main()

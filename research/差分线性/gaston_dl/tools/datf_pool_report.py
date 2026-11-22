#!/usr/bin/env python3
"""Generate a Chinese report of explicit trail-pool matching and validation."""
import argparse
import csv
import hashlib
import html
import json
from fractions import Fraction
from pathlib import Path


def esc(x):
    return html.escape(str(x))


def number(x):
    if x is None:
        return '未计算'
    n, d = int(x['numerator']), int(x['denominator'])
    if not n:
        return '0'
    sign = '−' if n < 0 else '+'
    w = d.bit_length()-1
    if abs(n) == 1:
        return sign+'1' if not w else f'{sign}2<sup>−{w}</sup>'
    return f'{sign}{abs(n)}/2<sup>{w}</sup>'


def bits(words):
    return ', '.join(f'({i},{j})' for i, w in enumerate(words) for j in range(64) if int(w,16)>>j&1)


def state_block(label, words):
    return '<div><b>'+esc(label)+'</b><pre>'+'\n'.join(f'x[{i}] = 0x{w}' for i,w in enumerate(words))+'</pre></div>'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    args=p.parse_args()
    root=args.results
    result=json.loads((root/'results.json').read_text())
    validation=json.loads((root/'validation.json').read_text())
    endpoint_file=root/'endpoints/results.json'
    endpoints=json.loads(endpoint_file.read_text())['results'] if endpoint_file.exists() else []
    endpoint_lookup={r['id']:r for r in endpoints}
    experiments=[]
    for filename in ('experiments.csv','confirmation.csv'):
        if (root/filename).exists():
            for row in csv.DictReader((root/filename).open()):
                experiments.append({**row,'file':filename})
    rows=[];details=[];flat=[]
    for i,s in enumerate(result['summaries']):
        b=s['best'];middle=s['middle_rounds'];family=s['family'];tail=s['tail']
        tag=f'{family}_M{middle}_{tail}'
        e=endpoint_lookup.get(tag)
        scope='穷尽' if s['exhaustive'] else '最佳评分已证' if s['optimal_score_certified'] else '批次搜索'
        wD=b['r_d'] if b else int(family[-1])
        rl=(b['r_l'] if b else 1 if tail=='L1_pool' else 2 if tail=='L2_2to2_pool' else next(l['rounds'] for l in result['pool_metadata']['paper_tails'] if l['id']==tail))
        total=wD+middle+rl
        split=f'{wD} + {middle} + {rl}'
        values=([f"{int(-b['differential_probability']['log2_abs'])}",
                 number(b['middle_DATF_correlation']),f"{int(-b['linear_squared_correlation']['log2_abs'])}",number(b['score'])] if b else ['—','未找到','—','—'])
        rows.append(f'<tr data-family="{family}" data-middle="{middle}" data-nonzero="{int(bool(b))}"><td>{total}<small>{split}</small></td><td>{esc(tail)}</td>'+''.join('<td>'+v+'</td>' for v in values)+f'<td>{s["nonzero_evaluated_pairs"]:,}<small>/ {s["evaluated_pairs"]:,}</small></td><td><span class="badge {"bounded" if scope=="批次搜索" else "complete"}">{scope}</span></td><td>{number(e["correlation"]) if e else "未复评"}</td><td>'+ (f'<a href="#candidate-{i}">边界</a>' if b else '—')+'</td></tr>')
        flat.append({'family':family,'r_m':middle,'tail':tail,'R':total,'split':split,
                     'evaluated_pairs':s['evaluated_pairs'],'nonzero_evaluated_pairs':s['nonzero_evaluated_pairs'],
                     'candidate_pair_scope':s['candidate_pair_scope'],'status':s['scope_status'],
                     'best_score_numerator':b['score']['numerator'] if b else '',
                     'best_score_denominator':b['score']['denominator'] if b else '',
                     'best_score_weight':b['score_weight'] if b else '',
                     'differential_id':b['differential_id'] if b else '',
                     'linear_id':b['linear_id'] if b else '',
                     'differential_rotation':b['differential_rotation'] if b else '',
                     'full_DATF_numerator':e['correlation']['numerator'] if e and e['correlation'] else '',
                     'full_DATF_denominator':e['correlation']['denominator'] if e and e['correlation'] else ''})
        if b:
            details.append(f'<details id="candidate-{i}"><summary>{split} · {esc(tail)} · S = {number(b["score"])}</summary><p>差分 {esc(b["differential_id"])}，共同左移 {b["differential_rotation"]} 位；线性 {esc(b["linear_id"])}。完整区间采用末 {total} 轮，零起始轮号 {b["start_round_index"]}。此处展示排名第一条；其余已保存候选见 results.json。</p><div class="states">'+state_block('完整输入差分 Δ',b['input_difference'])+state_block('完整输出掩码 Λ',b['output_mask'])+state_block('中间段输入差分 L(δD,out)',b['middle_input_difference'])+state_block('中间段输出掩码 Lᵀ(αL,in)',b['middle_output_mask'])+'</div><p>输入活跃位置：'+esc(bits(b['input_difference']))+'<br>输出活跃位置：'+esc(bits(b['output_mask']))+'</p></details>')
    with (root/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    erows=[]
    for r in endpoints:
        b=r['matching_winner']
        erows.append(f'<tr><td>{b["r_d"]}+{b["r_m"]}+{b["r_l"]}</td><td>{esc(b["linear_id"] if b["linear_id"].startswith("L1") else r["id"].split("_",2)[2])}</td><td>{number(b["score"])}</td><td>{number(r["correlation"])}</td><td>{esc(r["status"])}</td></tr>')
    srows=[]
    for r in experiments:
        srows.append(f'<tr><td>{esc(r["id"])}</td><td>{int(r["samples"]):,}</td><td>{float(r["rho"]):.9f}</td><td>[{float(r["ci95_low"]):.9f}, {float(r["ci95_high"]):.9f}]</td><td>{"是" if r["model_inside_ci95"]=="1" else "否"}</td></tr>')
    main_candidate=next(s['best'] for s in result['summaries'] if s['family']=='D1' and s['middle_rounds']==2 and s['tail']=='L1_pool')
    pool_rows=[]
    for key,label in [('differential_one','单轮差分'),('differential_pq_two','两轮 PQ 六比特差分'),('linear_one','单轮线性'),('linear_two_2to2','两轮 2→2 线性')]:
        if key not in result['pool_metadata']['pools']:
            continue
        a=result['pool_metadata']['pools'][key]
        pool_rows.append(f'<tr><td>{label}</td><td>{a["representatives"]:,}</td><td>{a["expanded_records"]:,}</td></tr>')
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston · DATF 路径池匹配</title>
<style>
:root{color-scheme:light;--ink:#172e39;--muted:#5c6e76;--line:#d8e2e4;--accent:#0f786e;--paper:#fff;--bg:#f3f6f5}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.7 system-ui,"Noto Sans CJK SC",sans-serif}main{max-width:1400px;margin:auto;padding:42px 34px 80px}header{border-bottom:2px solid var(--accent);padding-bottom:26px}h1{font-size:36px;line-height:1.3;margin:8px 0 16px}h2{font-size:24px;margin:0 0 14px}h3{font-size:18px}p{margin:10px 0 16px}a{color:var(--accent)}.eyebrow{font-size:13px;font-weight:700;letter-spacing:.12em;color:var(--accent)}.subtitle{max-width:1080px;color:var(--muted)}section{margin-top:28px;background:var(--paper);border:1px solid var(--line);border-radius:12px;padding:26px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-top:24px}.card{border:1px solid var(--line);border-radius:10px;padding:18px;background:white}.card strong{display:block;font-size:27px;line-height:1.4}.card span{color:var(--muted);font-size:14px}.callout{background:#e8f4f0;border-left:4px solid var(--accent);padding:18px 22px;border-radius:6px;margin:18px 0}.caution{background:#fff5e5;border-color:#b78029}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:12px 10px;border-bottom:1px solid var(--line);vertical-align:top}th{background:#eaf0ef;font-weight:650;white-space:nowrap}small{display:block;color:var(--muted);font-size:12px}sup{font-size:.75em}.scroll{overflow:auto}tbody tr:hover{background:#f5faf8}.badge{display:inline-block;white-space:nowrap;padding:2px 8px;border-radius:5px;font-size:12px}.complete{background:#e1f2eb;color:#12654d}.bounded{background:#fff0d5;color:#8a5716}.filters{display:flex;gap:16px;align-items:center;flex-wrap:wrap;padding:8px 0 18px}select{padding:7px;border:1px solid var(--line);border-radius:5px;font:inherit}button{font:inherit}label{white-space:nowrap}.formula{font:19px/1.7 ui-monospace,monospace;padding:15px;background:#edf3f2;border-radius:6px}pre{font-size:12px;line-height:1.7;overflow:auto;background:#f2f5f5;padding:14px;border-radius:6px}code{font-size:.9em}details{border:1px solid var(--line);border-radius:8px;margin:12px 0;padding:12px 16px}summary{font-weight:600;cursor:pointer}.states{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}li{margin:8px 0}.muted{color:var(--muted)}.links{display:flex;gap:22px;flex-wrap:wrap}.status{font-size:13px;color:var(--muted)}@media(max-width:800px){main{padding:24px 16px}.cards{grid-template-columns:repeat(2,1fr)}.states{grid-template-columns:1fr}section{padding:18px}h1{font-size:28px}}@media print{.filters{display:none}main{max-width:none;padding:0}section{break-inside:avoid}details{break-inside:avoid}.scroll{overflow:visible}}
</style><main><header><div class="eyebrow">GASTON / ORDINARY DATF / 2026-10-16</div><h1>差分池 × 线性池：一轮与两轮中间段匹配</h1><p class="subtitle">完整保存两端路径池，按普通 DATF 计算中间相关度。只省略整个区间的首个线性层；保留其后的线性层和参考轮常数。以下严格区分单条拼接评分、完整端点模型估计和实际置换采样。</p></header>
<div class="cards"><div class="card"><span>单轮端点比特上限</span><strong>2</strong><span>差分输入 / 线性输出</span></div><div class="card"><span>差分路径代表</span><strong>@D_COUNT@</strong><span>@D1_COUNT@ 单轮 + @D2_COUNT@ 两轮</span></div><div class="card"><span>线性路径代表</span><strong>@L_COUNT@</strong><span>@L1_COUNT@ 单轮 + 11 论文尾段</span></div><div class="card"><span>两轮中间段搜索规模</span><strong>@BATCH@ / 类</strong><span>所选差分代表覆盖全部移位</span></div></div>
<section><h2>值得继续研究的四轮候选</h2><div class="callout"><strong>1 轮差分 + 2 轮中间段 + 1 轮线性</strong><p>输入差分活跃位置：@FEATURE_INPUT@；输出掩码活跃位置：@FEATURE_OUTPUT@。坐标按 (行,位) 表示，均从 0 开始。</p><div class="formula">p<sub>D</sub> = @FEATURE_P@，ρ̂<sub>M</sub> = @FEATURE_M@，c<sub>L</sub><sup>2</sup> = @FEATURE_L@ ⇒ S = @FEATURE_S@</div><p>@FEATURE_ESTIMATE@</p></div><p>@FEATURE_CONCLUSION@</p>'''
    page+='<div class="states">'+state_block('输入差分 Δ',main_candidate['input_difference'])+state_block('输出掩码 Λ',main_candidate['output_mask'])+'</div></section>'
    page+='''<section><h2>评分定义与轮边界</h2><div class="formula">S = p<sub>D</sub> · ρ̂<sub>M</sub> · c<sub>L</sub><sup>2</sup>，R = r<sub>D</sub> + r<sub>M</sub> + r<sub>L</sub></div><p>p<sub>D</sub> 是一条差分路径的概率，c<sub>L</sub> 是一条线性路径的相关度乘积；ρ̂<sub>M</sub> 是把中间输入视为均匀变量后，普通 DATF 对中间段的估计。每一条评分只对应一条差分路径和同一线性路径对，未对其他差分路径、线性路径或交叉项求和。</p><p>中间段入口差分为 δ<sub>M</sub> = L(δ<sub>D,out</sub>)，出口掩码为 γ<sub>M</sub> = Lᵀ(α<sub>L,in</sub>)。这两个边界线性层都已计入。轮数按 χ 层计算，整体采用参考实现末 R 轮，零起始轮号为 12−R。</p><div class="callout caution">S 是筛选候选的单条拼接评分，不是完整端点相关度，也不是其严格上下界。中间输入的均匀性、线性尾段的处理以及普通 DATF 内部独立性都涉及模型假设。相关度 ρ 与偏差 ε 的关系为 ε = ρ/2。</div></section>'''
    page+='<section><h2>48 组匹配结果</h2><p>w<sub>D</sub> = −log₂p<sub>D</sub>，w<sub>L²</sub> = −log₂c<sub>L</sub>²。非零数只统计实际已计算的路径对。“未找到”在批次搜索中不表示整个池没有匹配。完整 DATF 列只复评本行排名第一的端点。</p><div class="filters"><label>差分段 <select id="family"><option value="">全部</option><option value="D1">1 轮</option><option value="D2">2 轮 PQ 6→6</option></select></label><label>中间段 <select id="middle"><option value="">全部</option><option value="1">1 轮</option><option value="2">2 轮</option></select></label><label><input type="checkbox" id="nonzero"> 只看非零匹配</label><span class="status" id="count"></span></div><div class="scroll"><table id="matches"><thead><tr><th>R / 分段</th><th>线性池 / 尾段</th><th>wD</th><th>ρ̂M</th><th>wL²</th><th>最佳 S</th><th>非零 / 已计算</th><th>覆盖状态</th><th>完整 DATF</th><th>详情</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table></div><p class="muted">固定尾段的单轮中间段已穷尽；单轮线性池的单轮中间段已逐层排除全部更低评分权重，故最佳评分已证，但没有枚举所有同分者。两轮中间段均为当前批次的结果。</p></section>'
    page+='<section><h2>完整端点复评与采样</h2><p>另选八组匹配的第一名，从完整均匀输入传播到最终掩码，使用原有普通 DATF 实现重新计算；不固定差分路径事件，也不再次乘上线性尾段因子。每项都检查了 32 组符号 DAG / 直接置换一致性，并独立审计模型因子的算术。</p><div class="scroll"><table><thead><tr><th>分段</th><th>线性路径</th><th>拼接 S</th><th>完整 DATF</th><th>计算状态</th></tr></thead><tbody>'+''.join(erows)+'</tbody></table></div><h3>完整置换采样</h3><div class="scroll"><table><thead><tr><th>候选</th><th>样本对数</th><th>经验相关度</th><th>95% Wilson 区间</th><th>完整 DATF 在区间内</th></tr></thead><tbody>'+''.join(srows)+'</tbody></table></div><p class="muted">轮常数使用参考实现，首个 L 省略。前三组各使用 2²⁴ 对输入，最后一组用独立种子作 2²⁶ 对确认。计数式 SplitMix64 提供可复现的伪随机仿真；区间逐项给出，未作多重比较校正。区间包含 0 表示当前样本未分辨出非零信号，不证明真实值等于 0。</p></section>'
    page+='<section><h2>路径池、覆盖范围与循环移位</h2><div class="scroll"><table><thead><tr><th>池</th><th>循环移位代表数</th><th>全部移位展开数</th></tr></thead><tbody>'+''.join(pool_rows)+'</tbody></table></div><ul><li>上限 4 对应 434,223,120 个非零端点，尚未展开 DDT/LAT 分支。因此按授权降为上限 2：51,360 个非零端点。单轮池包含上限内全部端点以及全部非零局部转移。</li><li>两轮差分池包含全部 11 个成对 P=Q=0 的六比特类别；每类保留唯一最大概率的首轮前驱，枚举全部 4⁶ 个第二轮输出。该明确范围内每条两轮差分路径概率为 2⁻²⁴。没有扩展为任意 θ-effect kernel 状态。</li><li>线性论文池包含 D2_repeat_2 至 D2_repeat_6，以及 D3_chi2_to_3、D3_chi2_to_4、D3_chi3_to_4、D3_chi1_to_2、D3_chi1_to_3、D3_chi1_to_4。</li><li>两轮中间段先选每类差分池最多 @BATCH@ 个代表，全部循环移位：D1 为 @M2_D1@ 条，D2 为 @M2_D2@ 条。D2 按 11 个 PQ 类轮流取样；优先取旧单轮连接及中间输入较稀疏的候选。配单轮线性池时，本批取权重和入口活跃数排序的前 @L1_BATCH@ 个代表，每类预算 @L1_BUDGET@ 对。</li><li>只展开差分侧的相对循环移位即可；线性侧保持规范代表。存在周期 32 的单轮路径时只展开其不同移位，避免重复。这里说的是五行共同循环移位，不是每行独立移位。</li><li>一、两轮中间段的本实现中，轮常数最多影响记录的符号相位，因此普通 DATF 的单条评分绝对值在共同移位下不变。此约简不证明固定轮常数下完整置换相关度不变，也不适用于未加证明的带符号路径和。</li></ul><p>完整笛卡尔积每种中间长度超过 470 亿对。当前把单轮中间的较低代价搜索做全，并把两轮中间的覆盖边界明确保存；没有把未计算候选记作模型零。</p></section>'
    page+='<section><h2>验证与文件</h2><p>新增验证全部通过：240 个小规模二次多项式的完全求和，全部 1,024 个单 χ 差分/掩码组合，64 个完整 Lχ 符号状态；25,440 条单轮池转移；'+str(validation['compiled_vs_original_python_evaluations'])+' 项 C++ / 原始 Python 普通 DATF 对照，其中 '+str(validation['nonzero_oracle_cases'])+' 项非零；1,536 项共同移位与轮常数检查。全部 11 个论文尾段的旧 PQ 单轮中间连接计数也已复现。</p><p>加速只替换普通 DATF 内部 Bias 的求和后端：一、两轮情形化为均匀二次多项式，用精确消元计算。同样使用原孤立变量规则与独立性假设，没有加入增强 DATF 的划分、回溯或转换规则。</p><div class="links"><a href="results.json">全部匹配数据</a><a href="summary.csv">汇总 CSV</a><a href="validation.json">验证记录</a><a href="endpoints/results.json">完整端点复评</a><a href="experiments.csv">第一批采样</a><a href="confirmation.csv">独立确认采样</a><a href="../pools/pools.json">路径池清单</a><a href="reproducibility.json">SHA-256 清单</a></div></section><section><h2>已找到候选的具体边界</h2><p>每组展示一条最佳代表。点击展开；全部已保存的同组候选与分段权重在 results.json。</p>'+''.join(details)+'</section><footer class="muted"><p>研究对象：Gaston 减轮核心 · 只省略第一个 L · 原始普通 DATF · 结果范围见各表覆盖状态</p></footer></main>'
    page+='''<script>
function filterRows(){let family=document.querySelector('#family').value,middle=document.querySelector('#middle').value,nz=document.querySelector('#nonzero').checked,n=0;document.querySelectorAll('#matches tbody tr').forEach(r=>{let show=(!family||r.dataset.family===family)&&(!middle||r.dataset.middle===middle)&&(!nz||r.dataset.nonzero==='1');r.hidden=!show;if(show)n++});document.querySelector('#count').textContent=`显示 ${n} / 48 组`}
document.querySelectorAll('.filters select,.filters input').forEach(x=>x.addEventListener('change',filterRows));filterRows();
document.querySelectorAll('a[href^="#candidate-"]').forEach(a=>a.addEventListener('click',()=>{document.querySelector(a.getAttribute('href')).open=true}));
</script></html>'''
    featured_full=endpoint_lookup.get('D1_M2_L1_pool')
    estimates='完整四轮的普通 DATF 输出为 '+(number(featured_full['correlation']) if featured_full else '未计算')+'。'
    conclusion='完整端点的模型估计和采样结果需要单独验证；不能从减轮实验直接推得完整 Gaston 的结论。'
    confirmation=next((r for r in experiments if r['id']=='D1_M2_L1_pool' and r['file']=='confirmation.csv'),None)
    if confirmation:
        lo,hi=float(confirmation['ci95_low']),float(confirmation['ci95_high'])
        estimates+=f"独立确认采样使用 {int(confirmation['samples']):,} 对输入，得到 ρ̂MC = {float(confirmation['rho']):.9f}，95% 区间 [{lo:.9f}, {hi:.9f}]。"
        if lo>0 or hi<0:
            conclusion='此采样区间排除 0，为该四轮端点提供非零相关度的实验支持。'
        if confirmation['model_inside_ci95']=='0':
            conclusion+='完整 DATF 估计落在该区间之外，显示出本端点上的定量估计偏差。'
        conclusion+='采样不能确定相关度的精确值，也不能从减轮实验直接推得完整 Gaston 的结论。'
    ps=result['pool_metadata']['pools']
    replacements={'FEATURE_INPUT':bits(main_candidate['input_difference']),
                  'FEATURE_OUTPUT':bits(main_candidate['output_mask']),
                  'FEATURE_P':number(main_candidate['differential_probability']),
                  'FEATURE_M':number(main_candidate['middle_DATF_correlation']),
                  'FEATURE_L':number(main_candidate['linear_squared_correlation']),
                  'FEATURE_S':number(main_candidate['score']),
                  'FEATURE_ESTIMATE':estimates,'FEATURE_CONCLUSION':conclusion,
                  'D_COUNT':f"{ps['differential_one']['representatives']+ps['differential_pq_two']['representatives']:,}",
                  'D1_COUNT':f"{ps['differential_one']['representatives']:,}",
                  'D2_COUNT':f"{ps['differential_pq_two']['representatives']:,}",
                  'L_COUNT':f"{ps['linear_one']['representatives']+ps.get('linear_two_2to2',{}).get('representatives',0)+len(result['pool_metadata']['paper_tails']):,}",
                  'L1_COUNT':f"{ps['linear_one']['representatives']:,}",
                  'BATCH':str(result['middle2_batch']['representatives_per_differential_family']),
                  'L1_BUDGET':f"{result['linear1_pair_budget_per_family_and_middle']:,}",
                  'L1_BATCH':str(result['linear1_middle2_representatives'])}
    for family in ('D1','D2'):
        replacements['M2_'+family]=f"{next(s['evaluated_pairs'] for s in result['summaries'] if s['family']==family and s['middle_rounds']==2 and s['tail']=='D2_repeat_2'):,}"
    for key,value in replacements.items():page=page.replace('@'+key+'@',value)
    page=page.replace('48 组', f'{len(rows)} 组')
    if 'linear_two_2to2' in ps:
        count=ps['linear_two_2to2']['representatives']
        page=page.replace('单轮 + 11 论文尾段',f'单轮 + {count:,} 条两轮 2→2 + 11 论文尾段')
        page=page.replace('</ul><p>完整笛卡尔积',
            '<li>两轮 2→2 线性池包含 10 个内部桥接类别及全部 χ 补全；新增 16,000 个代表。'
            '单轮中间段按权重逐层搜索，两轮中间段使用独立预算；覆盖以各行状态为准。</li></ul><p>完整笛卡尔积')
    (root/'datf_pool_matches.html').write_text(page)
    project=Path(__file__).resolve().parents[1]
    source_paths=[project/name for name in ('tools/datf_pool_build.py','tools/datf_pool_search.py','tools/datf_pool_endpoint_check.py','tools/datf_pool_report.py','include/datf_middle.hpp','src/datf_pool_match.cpp','src/datf_middle_cli.cpp','src/datf_pool_verify.cpp','tests/test_datf_middle.cpp','tests/test_datf_pools.py','tools/datf.py','tools/datf_k06.py','tools/datf_fixed_tail.py','include/gaston.hpp','include/experiment.hpp','Makefile','README.md')]
    paths=source_paths+[f for f in root.rglob('*') if f.is_file() and f.name!='reproducibility.json']
    paths += [f for f in (root.parent/'pools').iterdir() if f.is_file()]
    manifest={'scope':'source, persisted pools, matching results, endpoint estimates and sampling files',
              'files':[{'path':str(f.resolve()),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in paths]}
    (root/'reproducibility.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(root/'datf_pool_matches.html')


if __name__=='__main__':main()

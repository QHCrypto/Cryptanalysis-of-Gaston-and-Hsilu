#!/usr/bin/env python3
"""Merge the pool report and completed paper-tail report without rerunning searches."""
import argparse
import csv
import hashlib
import html
import json
import os
from collections import Counter
from pathlib import Path

from datf_pool_report import bits, number, state_block

ROOT = Path(__file__).resolve().parents[1]
LABELS = {
    'found_nonzero': '找到非零',
    'all_endpoints_model_zero': '全部候选端点为模型零',
    'no_nonzero_matches': '无非零中间匹配',
    'first_nonzero_only': '首名非零；仅复评首名',
    'first_zero_only': '首名为零；其余未复评',
    'batch_no_match': '批次内未找到匹配',
    'not_evaluated': '尚未复评',
}
CSS = '''
:root{color-scheme:light;--ink:#203438;--muted:#647877;--line:#d7e2de;--accent:#12675d;--bg:#f2f5f1}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:25px}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.75 system-ui,"Noto Sans CJK SC",sans-serif}main{max-width:1540px;margin:auto;padding:42px 32px 70px}header{padding:8px 0 24px;border-bottom:2px solid var(--accent)}h1{font-size:34px;line-height:1.45;margin:8px 0 14px}h2{font-size:23px;line-height:1.5;margin:0 0 12px}h3{font-size:18px;margin:10px 0}p{margin:10px 0 16px}.eyebrow{font-size:12px;letter-spacing:.12em;font-weight:700;color:var(--accent)}.subtitle,.muted,small{color:var(--muted)}.subtitle{max-width:1120px}.nav{display:flex;flex-wrap:wrap;gap:12px 24px;margin:18px 0 0}a{color:var(--accent);text-underline-offset:3px}section{background:white;border:1px solid var(--line);border-radius:12px;padding:26px;margin-top:24px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-top:24px}.card{background:white;border:1px solid var(--line);border-radius:10px;padding:18px}.card strong{display:block;font-size:29px;line-height:1.5}.card span{color:var(--muted);font-size:13px}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}.feature{padding:20px;background:#f0f7f2;border-radius:8px;border:1px solid #d7e7da}.feature.paper{background:#f2f5fa;border-color:#dce3ed}.tag,.badge{display:inline-block;font-size:12px;border-radius:4px;padding:2px 7px;white-space:nowrap}.tag{background:#dbeadf;color:#345c42}.badge.good{background:#e3f1e8;color:#216044}.badge.pending{background:#fff0d6;color:#805c21}.badge.neutral{background:#edf0f2;color:#566873}.formula{padding:12px 15px;background:#edf2f0;border-radius:6px;font-size:17px;overflow:auto}.callout{border-left:4px solid #bb9450;padding:14px 18px;background:#fff7e8;margin:18px 0;border-radius:4px}.filters{display:flex;flex-wrap:wrap;gap:12px 18px;align-items:center;margin:18px 0}.filters label{font-size:13px;white-space:nowrap}select,input[type=search],button{font:inherit;padding:6px 9px;border:1px solid #c8d5cf;border-radius:5px;background:white;color:inherit}input[type=search]{width:215px}button{cursor:pointer}button:hover{background:#edf4ef}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:10px 9px;border-bottom:1px solid var(--line);vertical-align:top;text-align:left}th{background:#edf3ee;white-space:nowrap;font-weight:650}tbody tr:hover{background:#f6faf7}small{display:block;font-size:11px;line-height:1.6}.scroll{overflow:auto;border:1px solid var(--line);border-radius:6px}#matches{min-width:1400px}.nowrap{white-space:nowrap}.numeric{font-variant-numeric:tabular-nums}.tail{font-family:ui-monospace,monospace;font-size:12px}.empty{text-align:center;padding:20px;color:var(--muted)}details{border:1px solid var(--line);padding:12px 16px;margin:12px 0;border-radius:7px}summary{cursor:pointer;font-weight:650}details:target{border-color:var(--accent)}.states{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}pre{font:12px/1.7 ui-monospace,monospace;overflow:auto;background:#f2f5f3;padding:12px;border-radius:5px}.links{display:flex;gap:10px 22px;flex-wrap:wrap}li{margin:7px 0}.scope{font-size:12px}footer{font-size:13px;padding:24px 2px}.legend{font-size:13px}.report-count{font-size:13px;color:var(--muted)}[hidden]{display:none!important}@media(max-width:850px){main{padding:22px 14px}.cards{grid-template-columns:repeat(2,1fr)}.grid,.states{grid-template-columns:1fr}section{padding:18px}h1{font-size:27px}.filters{align-items:flex-start}}@media print{.filters,.nav,button{display:none}main{max-width:none;padding:0}body{background:white}section{border-radius:0}.scroll{overflow:visible}#matches{min-width:0;font-size:9px}th,td{padding:4px}.feature{break-inside:avoid}}
'''
SCRIPT = '''
const selects=[...document.querySelectorAll('.filters select')];
const search=document.querySelector('#query');
function filterRows(){
  const filters=Object.fromEntries(selects.map(s=>[s.id,s.value]));
  const query=search.value.trim().toLowerCase();let count=0;
  document.querySelectorAll('#matches tbody tr').forEach(row=>{
    const show=Object.entries(filters).every(([key,value])=>!value||row.dataset[key]===value)
      &&(!query||row.textContent.toLowerCase().includes(query));
    row.hidden=!show;if(show)count++;
  });
  document.querySelector('#shown').textContent=`显示 ${count} / 48 行`;
  document.querySelector('#empty').hidden=count!==0;
}
selects.forEach(s=>s.addEventListener('change',filterRows));search.addEventListener('input',filterRows);
document.querySelector('#reset').addEventListener('click',()=>{selects.forEach(s=>s.value='');search.value='';filterRows()});
function openHash(){const target=document.getElementById(decodeURIComponent(location.hash.slice(1)));if(target&&target.tagName==='DETAILS')target.open=true}
document.querySelectorAll('a[href^="#detail-"]').forEach(a=>a.addEventListener('click',()=>{const target=document.getElementById(a.hash.slice(1));if(target)target.open=true}));
window.addEventListener('hashchange',openHash);openHash();filterRows();
'''


def esc(value):
    return html.escape(str(value), quote=True)


def load(path):
    return json.loads(path.read_text())


def nonzero(value):
    return value is not None and int(value['numerator']) != 0


def normalize(old, review, endpoints):
    current = {r['id']: r for r in review['rows']}
    rows = []
    for s in old['summaries']:
        rid = f'{s["family"]}_M{s["middle_rounds"]}_{s["tail"]}'
        if s['tail'] != 'L1_pool':
            r = current[rid]
            first = r['first_ranked_endpoint_result']
            winner = r['first_nonzero']
            row = dict(id=rid, family=r['family'], middle=r['r_m'], tail=r['tail'],
                       r_d=r['r_d'], r_m=r['r_m'], r_l=r['r_l'], R=r['R'], kind='paper',
                       source='paper_full24', coverage='exhaustive', coverage_label='已穷尽',
                       evaluated_pairs=r['evaluated_pairs'], nonzero_matches=r['nonzero_matches'],
                       best=r['best'], selected=winner['candidate'] if winner else None,
                       first_correlation=first['correlation'] if first else None,
                       selected_correlation=r.get('full_DATF') if winner else None,
                       selected_bias=r.get('bias') if winner else None,
                       selected_rank=winner['first_match_rank'] if winner else None,
                       tested=r['evaluated_unique_endpoints'], unique_endpoints=r['unique_endpoints'],
                       status=r['status'], matches_file=r['matches_file'],
                       ranked_endpoints_file=r['ranked_endpoints_file'])
        else:
            b = s['best']; e = endpoints.get(rid)
            rho = e['correlation'] if e else None
            status = ('first_nonzero_only' if nonzero(rho) else 'first_zero_only') if e else ('batch_no_match' if not b else 'not_evaluated')
            coverage = 'certified' if s['optimal_score_certified'] else 'bounded'
            row = dict(id=rid, family=s['family'], middle=s['middle_rounds'], tail=s['tail'],
                       r_d=int(s['family'][-1]), r_m=s['middle_rounds'], r_l=1,
                       R=int(s['family'][-1])+s['middle_rounds']+1, kind='linear1', source='search_final',
                       coverage=coverage, coverage_label='最佳拼接评分已证' if coverage=='certified' else '批次搜索',
                       evaluated_pairs=s['evaluated_pairs'], nonzero_matches=s['nonzero_evaluated_pairs'],
                       candidate_pair_scope=s['candidate_pair_scope'], best=b,
                       selected=b if nonzero(rho) else None, first_correlation=rho,
                       selected_correlation=rho if nonzero(rho) else None,
                       selected_bias=e['bias'] if nonzero(rho) else None,
                       selected_rank=1 if nonzero(rho) else None, tested=int(e is not None),
                       unique_endpoints=None, status=status)
        row['split'] = f'{row["r_d"]}+{row["r_m"]}+{row["r_l"]}'
        row['nonzero'] = 'yes' if nonzero(row['selected_correlation']) else 'no'
        rows.append(row)
    assert len(rows) == len({r['id'] for r in rows}) == 48
    assert sum(r['kind']=='paper' for r in rows) == 44
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pool-results', type=Path, default=ROOT/'results/datf_pools/search_final')
    parser.add_argument('--paper-results', type=Path, default=ROOT/'results/datf_pools/paper_full24')
    parser.add_argument('--output', type=Path, default=ROOT/'results/datf_pools/merged')
    args = parser.parse_args()
    oldroot, paperroot, out = args.pool_results.resolve(), args.paper_results.resolve(), args.output.resolve()
    assert out not in (oldroot, paperroot), 'Write the merged report into a separate directory.'
    out.mkdir(parents=True, exist_ok=True)
    inputs = [oldroot/'results.json', oldroot/'endpoints/results.json', oldroot/'validation.json',
              oldroot/'experiments.csv', oldroot/'confirmation.csv', oldroot/'datf_pool_matches.html',
              paperroot/'search_results.json', paperroot/'reassessment/results.json',
              paperroot/'validation.json', paperroot/'paper_full_datf.html',
              paperroot/'backend_validation.json', paperroot/'original_R9_controls.json',
              paperroot/'supplemental_factor_audits.json']
    old, search, review = load(inputs[0]), load(paperroot/'search_results.json'), load(paperroot/'reassessment/results.json')
    oldcheck, newcheck = load(oldroot/'validation.json'), load(paperroot/'validation.json')
    endpoint_records = load(oldroot/'endpoints/results.json')['results']
    endpoints = {r['id']:r for r in endpoint_records}
    assert old['first_L_omitted'] and search['first_L_omitted'] and not old['enhanced_DATF']
    assert old['pool_metadata']==search['pool_metadata']
    assert oldcheck['status']==newcheck['status']=='PASS' and review['complete']
    rows = normalize(old, review, endpoints)
    lookup = {r['id']:r for r in rows}
    experiments = []
    for name in ('experiments.csv', 'confirmation.csv'):
        with (oldroot/name).open() as f:
            experiments.extend({**r, 'source_file':name} for r in csv.DictReader(f))
    for e in experiments:
        b = lookup[e['id']]['best']
        assert e['input_words'].split(';')==b['input_difference']
        assert e['mask_words'].split(';')==b['output_mask']
    def href(path): return esc(os.path.relpath(path, out))
    def link(path, text): return f'<a href="{href(path)}">{esc(text)}</a>'
    def detail_link(rid, text): return f'<a href="#detail-{esc(rid)}">{esc(text)}</a>'
    sources = [{'path':str(f), 'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in inputs]
    merged = {'title':'Gaston 差分线性路径池综合报告', 'date':'2026-10-16',
              'merge_rule':'44 paper-tail rows replaced by completed paper_full24 results; 4 L1 rows and historical samples retained from search_final',
              'first_L_omitted':True, 'enhanced_DATF':False, 'rows':rows,
              'experiments':experiments, 'sources':sources}
    (out/'merged_results.json').write_text(json.dumps(merged, ensure_ascii=False, indent=2)+'\n')
    table, details = [], []
    for r in rows:
        b = r['best']; selected = r['selected']; shown = selected or b
        coverage_style = 'pending' if r['coverage']=='bounded' else 'good'
        result_style = 'good' if r['nonzero']=='yes' else 'pending' if r['kind']=='linear1' else 'neutral'
        tested = (f'{r["tested"]:,} / {r["unique_endpoints"]:,}' if r['kind']=='paper'
                  else '1；仅首名' if r['tested'] else '0；未复评')
        full = number(r['first_correlation']) if b else '—'
        selected_value = number(r['selected_correlation']) if selected else '—'
        if selected: selected_value += f'<small>匹配排名 {r["selected_rank"]}</small>'
        attrs = ' '.join(f'data-{key}="{esc(r[key])}"' for key in ('family','middle','kind','coverage','nonzero'))
        cells = [f'<b>{r["R"]} 轮</b><small>{r["split"]}</small>',
                 f'<span class="tail">{esc(r["tail"])}</span><small>{"论文尾段 · 全量更新" if r["kind"]=="paper" else "单轮线性池 · 保留原结果"}</small>',
                 f'<span class="badge {coverage_style}">{r["coverage_label"]}</span>',
                 f'{r["nonzero_matches"]:,}<small>/ {r["evaluated_pairs"]:,} 已计算</small>',
                 number(b['score']) if b else '—', full, selected_value, tested,
                 f'<span class="badge {result_style}">{LABELS[r["status"]]}</span>',
                 detail_link(r['id'],'展开边界') if b else '—']
        table.append(f'<tr {attrs}>'+''.join('<td>'+c+'</td>' for c in cells)+'</tr>')
        if shown:
            if r['kind']=='paper':
                note = ('本行按 |S| 排序复评，展示首个完整 DATF 非零端点。' if selected else
                        '本行所有非零匹配对应的不同端点均已复评为模型零；下面保留拼接首名的边界。')
                refs = link(paperroot/r['matches_file'],'全部非零匹配')+' · '+link(paperroot/r['ranked_endpoints_file'],'按排名去重的端点')
            else:
                note = '此处展示原报告的拼接首名，只复评了该端点；没有继续检查本行其余候选。'
                refs = link(oldroot/'results.json','原匹配候选')+' · '+link(oldroot/'endpoints/results.json','原完整端点复评')
            formula = ('p<sub>D</sub> = '+number(shown['differential_probability'])+'，ρ̂<sub>M</sub> = '+number(shown['middle_DATF_correlation'])+
                       '，c<sub>L</sub>² = '+number(shown['linear_squared_correlation'])+'，S = '+number(shown['score']))
            correlation = r['selected_correlation'] or r['first_correlation']
            detail = f'<details id="detail-{r["id"]}"><summary>{r["split"]} · {esc(r["tail"])} · {LABELS[r["status"]]}</summary><p>{note}</p>'
            detail += f'<p>差分路径 <code>{esc(shown["differential_id"])}</code>；线性路径 <code>{esc(shown["linear_id"])}</code>；差分共同左移 {shown["differential_rotation"]} 位；完整区间零起始轮号 {shown["start_round_index"]}。</p><div class="formula">{formula}</div>'
            detail += '<p>本端点完整普通 DATF 相关度：'+number(correlation)+'。'+('偏差 ε̂ = '+number(r['selected_bias'])+'。' if selected else '')+'</p><div class="states">'
            for label,key in [('完整输入差分 Δ','input_difference'),('完整输出掩码 Λ','output_mask'),('中间段输入差分 L(δD,out)','middle_input_difference'),('中间段输出掩码 Lᵀ(αL,in)','middle_output_mask')]:
                detail += state_block(label,shown[key])
            detail += '</div><p class="scope">输入活跃位置：'+esc(bits(shown['input_difference']))+'<br>输出活跃位置：'+esc(bits(shown['output_mask']))+'</p><p>'+refs+'</p></details>'
            details.append(detail)
    feature = lookup['D1_M2_L1_pool']; paper_feature = lookup['D1_M1_D3_chi2_to_3']
    confirmation = next(e for e in experiments if e['source_file']=='confirmation.csv')
    feature_blocks = f'''<div class="grid"><article class="feature"><span class="tag">原报告 · 已有采样</span><h3>1 + 2 + 1：单轮线性池四轮候选</h3>
<p>拼接评分 S = {number(feature['best']['score'])}；完整普通 DATF ρ̂ = {number(feature['first_correlation'])}。</p>
<div class="formula">确认采样 ρ̂<sub>MC</sub> = {float(confirmation['rho']):.9f}</div><p>{int(confirmation['samples']):,} 对输入；95% 区间 [{float(confirmation['ci95_low']):.9f}, {float(confirmation['ci95_high']):.9f}]。区间排除 0，为该端点的非零相关度提供实验支持，但完整 DATF 估计落在区间外。</p>
<p class="scope">该行仍属批次搜索，不能声称已找到池内最优完整端点。</p>{detail_link(feature['id'],'查看完整边界')}</article>
<article class="feature paper"><span class="tag">全量更新 · 模型结果</span><h3>1 + 1 + 2：固定论文尾段四轮候选</h3><p>尾段 <code>D3_chi2_to_3</code>；首名拼接评分 S = {number(paper_feature['best']['score'])}。</p>
<div class="formula">完整普通 DATF ρ̂ = {number(paper_feature['selected_correlation'])}</div><p>偏差 ε̂ = {number(paper_feature['selected_bias'])}；|ρ̂| ≈ 2<sup>−40.830</sup>。在 44 行固定论文尾段组合中，只有这一行按既定停止规则找到了非零。</p><p class="scope">该行第 1 名即非零后停止，未进一步优化完整 DATF 相关度；尚无本端点的采样结果。</p>{detail_link(paper_feature['id'],'查看完整边界')}</article></div>'''
    sampling_rows = []
    for e in experiments:
        sampling_rows.append('<tr>'+''.join('<td>'+c+'</td>' for c in [detail_link(e['id'],e['id']),
            '独立确认' if e['source_file']=='confirmation.csv' else '第一批', f'{int(e["samples"]):,}',
            number(lookup[e['id']]['first_correlation']), f'{float(e["rho"]):.9f}',
            f'[{float(e["ci95_low"]):.9f}, {float(e["ci95_high"]):.9f}]',
            '是' if e['model_inside_ci95']=='1' else '否',esc(e['seed'])])+'</tr>')
    pool_rows = []
    for key,label in [('differential_one','D1：单轮差分'),('differential_pq_two','D2：两轮 PQ 六比特差分'),('linear_one','L1：单轮线性')]:
        p=old['pool_metadata']['pools'][key]
        pool_rows.append(f'<tr><td>{label}</td><td>{p["representatives"]:,}</td><td>{p["expanded_records"]:,}</td></tr>')
    tail_rows = ''.join(f'<tr><td class="tail">{t["id"]}</td><td>{t["rounds"]}</td><td>{t["weight_squared"]}</td></tr>' for t in old['pool_metadata']['paper_tails'])
    history_rows = ''.join('<tr>'+''.join('<td>'+c+'</td>' for c in [esc(e['id']),number(e['matching_winner']['score']),number(e['correlation']),e['factor_audit']['status'],
                            '仍用于单轮线性池首名列' if e['id'].endswith('L1_pool') else '该类结果已更新至全量复评'])+'</tr>' for e in endpoint_records)
    counts = Counter(r['status'] for r in rows if r['kind']=='paper')
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="合并 Gaston 路径池匹配、固定论文尾段全量 DATF 复评与已有采样的研究报告。"><title>Gaston · 差分线性路径池综合报告</title><style>{CSS}</style></head><body><main>
<header><div class="eyebrow">GASTON / ORDINARY DATF / 2026-10-16</div><h1>差分线性路径池综合报告</h1><p class="subtitle">两类差分池 × 单轮线性池及 11 个论文尾段。中间段分别为 1 轮、2 轮；统一呈现候选匹配、完整端点 DATF 复评和已有置换采样。轮数按 χ 层计，仍只省略整个区间的第一个线性层。</p><nav class="nav" aria-label="报告导航"><a href="#highlights">四轮候选</a><a href="#all-results">48 行总表</a><a href="#sampling">采样结果</a><a href="#scope">池与覆盖范围</a><a href="#method">方法与核验</a><a href="#boundaries">具体边界</a><a href="#sources">来源文件</a></nav></header>
<div class="cards"><div class="card"><span>合并后组合</span><strong>48 行</strong><span>44 行论文尾段 + 4 行单轮线性池</span></div><div class="card"><span>论文尾段覆盖</span><strong>44 / 44</strong><span>匹配穷尽；按规则完成逐行复评</span></div><div class="card"><span>论文尾段完整端点复评</span><strong>{review['completed_endpoint_jobs']:,}</strong><span>不同端点；重复路径共用结果</span></div><div class="card"><span>有非零完整 DATF 代表的组合</span><strong>{sum(r['nonzero']=='yes' for r in rows)} 行</strong><span>单轮线性池 2 行 + 论文尾段 1 行</span></div></div>
<section id="highlights"><h2>两项四轮候选：采样与模型估计</h2>{feature_blocks}<p class="legend">三轮的 1+1+1 候选也保留在总表中：完整 DATF 为 −2<sup>−4</sup>，已有采样与该估计相容。上面的两项四轮结果分别来自不同线性池，不能用“论文尾段只有一行非零”概括整个合并报告。</p></section>
<section id="all-results"><h2>48 行统一结果</h2><p>论文尾段的 44 行采用全量更新结果，单轮线性池的 4 行采用原报告结果。首名 S 与首名完整 DATF 对应同一端点；“获选完整 DATF”展示按本行规则确认的非零候选。</p>
<div class="filters"><label>差分段 <select id="family"><option value="">全部</option><option value="D1">D1 · 1 轮</option><option value="D2">D2 · 2 轮</option></select></label><label>中间段 <select id="middle"><option value="">全部</option><option value="1">1 轮</option><option value="2">2 轮</option></select></label><label>线性池 <select id="kind"><option value="">全部</option><option value="linear1">单轮线性池</option><option value="paper">论文尾段</option></select></label><label>匹配覆盖 <select id="coverage"><option value="">全部</option><option value="exhaustive">已穷尽</option><option value="certified">最佳拼接评分已证</option><option value="bounded">批次搜索</option></select></label><label>完整 DATF <select id="nonzero"><option value="">全部</option><option value="yes">有非零代表</option><option value="no">其余结果</option></select></label><label>搜索 <input id="query" type="search" placeholder="尾段名称、状态或分段"></label><button id="reset" type="button">重置</button><span id="shown" class="report-count" aria-live="polite"></span></div>
<div class="scroll"><table id="matches"><thead><tr><th>R / 分段</th><th>线性池 / 尾段</th><th>匹配覆盖</th><th>非零 / 已计算路径对</th><th>首名 S</th><th>首名完整 DATF</th><th>获选完整 DATF</th><th>已复评 / 不同端点</th><th>复评结论</th><th>边界</th></tr></thead><tbody>{''.join(table)}</tbody></table><div id="empty" class="empty" hidden>当前筛选条件下没有结果。</div></div>
<p class="legend">“最佳拼接评分已证”只针对指定池内的单条 |S| 目标，不表示全部同分路径都已枚举，更不表示完整 DATF 最优。“首名为零”只评价这一端点；“全部候选端点为模型零”才表示本行全部非零匹配对应的不同端点均已检查。“批次内未找到”不表示整个池没有匹配。</p></section>
<section id="sampling"><h2>保留原报告的完整置换采样</h2><p>以下 4 条实验记录原样合并，没有新增采样。它们均对应单轮线性池的候选，与论文尾段的模型复评分开解释。</p><div class="scroll"><table id="sampling-table"><thead><tr><th>候选</th><th>批次</th><th>样本对数</th><th>完整 DATF</th><th>经验相关度</th><th>95% Wilson 区间</th><th>DATF 在区间内</th><th>随机种子</th></tr></thead><tbody>{''.join(sampling_rows)}</tbody></table></div><p class="legend">首批各 2²⁴ 对，独立确认使用 2²⁶ 对。均采用参考轮常数，首个 L 省略；计数式 SplitMix64 提供可复现的伪随机仿真。区间逐项给出，未作多重比较校正；包含 0 不证明真实相关度为零。采样结论仅适用于相应减轮端点。</p></section>
<section id="scope"><h2>路径池、移位与覆盖范围</h2><div class="grid"><div><h3>差分池与单轮线性池</h3><div class="scroll"><table><thead><tr><th>池</th><th>循环移位代表</th><th>展开全部移位</th></tr></thead><tbody>{''.join(pool_rows)}</tbody></table></div></div><div><h3>固定的 11 个论文尾段</h3><div class="scroll"><table><thead><tr><th>尾段</th><th>线性轮数</th><th>wL² = −log₂ cL²</th></tr></thead><tbody>{tail_rows}</tbody></table></div></div></div>
<ul><li>单轮池的外端重量上限为 2：差分限制首个 χ 的输入差分，线性限制末个 χ 的输出掩码。上限 4 对应 434,223,120 个非零端点，因此先采用上限 2 的 51,360 个非零端点，并展开非零局部转移。</li><li>两轮差分池限定于既定 11 类成对 P=Q=0 的六比特模式：每类保留唯一最大概率首轮前驱，枚举全部 4⁶ 个第二轮输出；每条路径概率 2⁻²⁴。这里的“完整池”指这一明确范围。</li><li>固定论文尾段的两种中间长度均已穷尽：{search['total_logical_pairs']:,} 个逻辑路径对，{search['total_nonzero_matches']:,} 条非零匹配全部保存，没有只留前 32 条的截断。</li><li>论文尾段复评结果为 {counts['found_nonzero']} 行找到非零、{counts['all_endpoints_model_zero']} 行全部候选端点为模型零、{counts['no_nonzero_matches']} 行无非零中间匹配。匹配共导出 {search['unique_full_endpoint_jobs']:,} 个不同完整端点；按首个非零停止规则实际计算 {review['completed_endpoint_jobs']:,} 个，其余 3,158 个属于已找到非零的同一行。</li><li>单轮线性池：中间 1 轮的两行已有最佳 |S| 证明；中间 2 轮的两行仍为批次搜索，差分侧每类取 {old['linear1_middle2_representatives']} 个代表展开移位，每类预算 {old['linear1_pair_budget_per_family_and_middle']:,} 对。原报告只复评已找到候选的首名，没有进行遇零继续的全候选复评。</li><li>只展开差分侧共同循环移位，线性侧保留规范代表。周期 32 的代表只展开 32 个不同移位，其余为 64 个。一、两轮中间段的单条 |S| 在该模型下具有共同移位不变性；这不证明固定轮常数下完整相关度或带符号路径和也具有同样对称性。</li></ul></section>
<section id="method"><h2>方法、轮边界与核验</h2><div class="formula">S = p<sub>D</sub> · ρ̂<sub>M</sub> · c<sub>L</sub>²；R = r<sub>D</sub> + r<sub>M</sub> + r<sub>L</sub>；ε = ρ / 2</div><p>中间入口差分 δ<sub>M</sub> = L(δ<sub>D,out</sub>)，出口掩码 γ<sub>M</sub> = Lᵀ(α<sub>L,in</sub>)。整体采用参考实现末 R 轮，零起始轮号 12−R，只省略整个区间第一个 L；后续线性层和参考轮常数保留。</p><p>完整 DATF 从均匀的 320 比特输入直接传播到最终输出掩码，不限定差分路径事件，也不再次乘尾段相关度平方。论文尾段按 |S| 递减顺序复评，首名为零则继续；跨路径的相同完整端点只计算一次。确认首个非零前，所有更早的不同端点均必须已经算为零。</p><div class="callout">单条拼接评分 S 不是完整端点相关度，也不是其严格上下界。普通 DATF 的零和非零都是包含独立性假设的模型估计；本报告保留已有采样用于对照。首个非零停止不等于完整 DATF 全局最优，减轮结果也不能直接表述为完整构造的结论。</div>
<p>论文尾段全量匹配使用 24 个独立 C++ 进程，耗时 {search['search_seconds']:.1f} 秒；完整端点复评使用 24 个 Python 进程槽位，耗时 {review['elapsed_seconds_this_run']:.1f} 秒。最终核验检查全部 {newcheck['all_nonzero_matches_checked']:,} 条非零匹配及 {newcheck['full_endpoint_artifacts_checked']:,} 个端点，所有因子审计与逐行排名前缀均通过。</p><p>后端仍采用普通 DATF 的转换规则、孤立变量循环与独立性假设，只改用精确二次求和。核对包括 200 个小规模完全求和、200 个原 Shannon Bias、28 个已有端点及 2 个九轮原后端对照；另有 3 个触及初始因子审计枚举上限的端点，已按不相交变量分量精确求和补核并用原后端复算通过。原池验证还覆盖 25,440 条单轮转移、122 项 C++ / 原 Python 对照及 1,536 项共同移位与轮常数检查。</p>
<details><summary>保留原报告的 8 项端点复评记录</summary><p>这是历史检查记录，不作为额外组合重复计数；论文尾段以总表中的最新覆盖状态为准。</p><div class="scroll"><table id="historical-table"><thead><tr><th>原记录</th><th>原拼接 S</th><th>完整 DATF</th><th>因子审计</th><th>合并后的用途</th></tr></thead><tbody>{history_rows}</tbody></table></div></details></section>
<section id="boundaries"><h2>候选的完整端点与分段边界</h2><p>点击总表中的“展开边界”可直达对应记录。每项保留具体差分、输出掩码、轮号、旋转量及分段概率；没有匹配的行不生成虚构边界。</p>{''.join(details)}</section>
<section id="sources"><h2>报告来源与原始记录</h2><p>合并规则：以“差分池 + 中间轮数 + 线性尾段”为键，用最新的 44 行论文尾段结果替换旧版本同类行，保留 4 行单轮线性池结果及已有采样。本次仅整理报告，没有重新搜索或采样。</p><h3>综合结果</h3><div class="links">{link(out/'merged_results.json','48 行合并数据与来源 SHA-256')}{link(out/'merge_validation.json','报告合并核对记录')}{link(Path(__file__),'报告生成脚本')}</div><h3>原路径池报告</h3><div class="links">{link(oldroot/'datf_pool_matches.html','原报告 HTML')}{link(oldroot/'results.json','原匹配数据')}{link(oldroot/'endpoints/results.json','原端点复评')}{link(oldroot/'experiments.csv','首批采样')}{link(oldroot/'confirmation.csv','独立确认采样')}{link(oldroot/'validation.json','原池验证')}{link(oldroot.parent/'pools/pools.json','路径池清单')}{link(oldroot/'reproducibility.json','原来源清单')}</div><h3>论文尾段全量报告</h3><div class="links">{link(paperroot/'paper_full_datf.html','全量报告 HTML')}{link(paperroot/'search_results.json','全量搜索数据')}{link(paperroot/'reassessment/results.json','完整 DATF 复评')}{link(paperroot/'validation.json','最终完成核验')}{link(paperroot/'original_R9_controls.json','九轮原后端对照')}{link(paperroot/'supplemental_factor_audits.json','3 项补核记录')}{link(paperroot/'reproducibility.json','全量 SHA-256 清单')}</div></section>
<footer>Gaston 减轮核心 · 普通 DATF · 省略首个 L · 合并日期：2026-10-16 · 主报告的文字、样式和筛选脚本均内嵌，可离线打开；原始数据链接相对于当前结果目录。</footer></main><script>{SCRIPT}</script></body></html>'''
    (out/'datf_merged.html').write_text(page)
    print(out/'datf_merged.html')


if __name__ == '__main__':
    main()

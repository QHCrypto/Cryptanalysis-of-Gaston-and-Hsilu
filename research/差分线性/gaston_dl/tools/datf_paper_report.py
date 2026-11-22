#!/usr/bin/env python3
"""Render the exhaustive paper-tail search and ranked full-DATF reassessment."""
import argparse
import csv
import html
import json
from pathlib import Path


def esc(x):return html.escape(str(x))


def value(r):
    if not r:return '未计算'
    n,d=int(r['numerator']),int(r['denominator'])
    if not n:return '0'
    sign='−' if n<0 else '+';w=d.bit_length()-1
    if abs(n)==1:return sign+('1' if w==0 else f'2<sup>−{w}</sup>')
    return sign+str(abs(n))+f'/2<sup>{w}</sup>'


def state(label,words):
    return '<div><b>'+esc(label)+'</b><pre>'+'\n'.join(f'x[{i}] = 0x{w}' for i,w in enumerate(words))+'</pre></div>'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True)
    args=p.parse_args();root=args.results
    search=json.loads((root/'search_results.json').read_text())
    reassessment=json.loads((root/'reassessment/results.json').read_text())
    statuses={'found_nonzero':'找到非零','all_endpoints_model_zero':'全部端点为模型零',
              'no_nonzero_matches':'无非零中间匹配','running':'复评中'}
    counts={s:sum(r['status']==s for r in reassessment['rows']) for s in statuses}
    rows=[];details=[];flat=[]
    for i,r in enumerate(reassessment['rows']):
        first=r['first_ranked_endpoint_result'];winner=r['first_nonzero'];best=r['best']
        full_first=first['correlation'] if first else None
        selected=winner['candidate'] if winner else None
        summary=statuses[r['status']];split=f'{r["r_d"]}+{r["r_m"]}+{r["r_l"]}'
        rows.append(f'<tr data-family="{r["family"]}" data-middle="{r["r_m"]}" data-status="{r["status"]}"><td>{r["R"]}<small>{split}</small></td><td>{esc(r["tail"])}</td><td>{r["evaluated_pairs"]:,}</td><td>{r["nonzero_matches"]:,}</td><td>{r["evaluated_unique_endpoints"]:,}<small>/ {r["unique_endpoints"]:,}</small></td><td>{value(best["score"]) if best else "—"}</td><td>{value(full_first) if best else "—"}</td><td>{winner["first_match_rank"] if winner else "—"}</td><td>{value(r.get("full_DATF")) if winner else "—"}</td><td><span class="badge {r["status"]}">{summary}</span></td><td>'+(f'<a href="#row-{i}">边界</a>' if best else '—')+'</td></tr>')
        flat.append({'row':r['id'],'R':r['R'],'split':split,'tail':r['tail'],'evaluated_pairs':r['evaluated_pairs'],
                     'nonzero_matches':r['nonzero_matches'],'unique_endpoints':r['unique_endpoints'],
                     'evaluated_unique_endpoints':r['evaluated_unique_endpoints'],'status':r['status'],
                     'first_splice_weight':best['score_weight'] if best else '',
                     'first_full_DATF_numerator':full_first['numerator'] if full_first else '',
                     'first_full_DATF_denominator':full_first['denominator'] if full_first else '',
                     'selected_match_rank':winner['first_match_rank'] if winner else '',
                     'selected_splice_weight':selected['score_weight'] if selected else '',
                     'selected_full_DATF_numerator':r['full_DATF']['numerator'] if winner else '',
                     'selected_full_DATF_denominator':r['full_DATF']['denominator'] if winner else ''})
        if best:
            b=selected or best
            details.append(f'<details id="row-{i}"><summary>{split} · {esc(r["tail"])} · {summary}</summary><p>'+('展示按评分顺序复评后首个非零端点。' if selected else '展示原拼接首名；本行尚无复评非零端点。')+f'差分路径 {esc(b["differential_id"])}，旋转 {b["differential_rotation"]} 位；完整区间零起始轮号 {b["start_round_index"]}。</p><div class="states">'+state('完整输入差分 Δ',b['input_difference'])+state('完整输出掩码 Λ',b['output_mask'])+state('中间段输入差分',b['middle_input_difference'])+state('中间段输出掩码',b['middle_output_mask'])+'</div><p>拼接评分 S = '+value(b['score'])+'。'+('完整普通 DATF ρ̂ = '+value(r['full_DATF'])+'，ε̂ = '+value(r['bias'])+'。' if selected else '')+'</p><p><a href="'+esc(r['matches_file'])+'">全部非零匹配</a> · <a href="'+esc(r['ranked_endpoints_file'])+'">按排名去重的完整端点</a></p></details>')
    with (root/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    complete=reassessment['complete'] and all(r['status']!='running' for r in reassessment['rows'])
    audit_path=root/'validation.json'
    audited=audit_path.exists() and json.loads(audit_path.read_text()).get('status')=='PASS'
    banner=('全量搜索、逐行复评与完成核验均已通过' if audited else '全量搜索与逐行复评已完成；等待最终核验') if complete else '全量搜索已完成；完整端点复评仍在运行'
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston · 论文尾段全量匹配与 DATF 复评</title><style>
*{box-sizing:border-box}body{margin:0;background:#f2f6f5;color:#1b303b;font:16px/1.7 system-ui,"Noto Sans CJK SC",sans-serif}main{max-width:1520px;margin:auto;padding:36px 28px 72px}header{border-bottom:2px solid #187568;padding-bottom:20px}h1{font-size:32px;line-height:1.4;margin:8px 0}h2{font-size:23px;margin:0 0 12px}.kicker{color:#187568;font-size:13px;font-weight:700;letter-spacing:.1em}section{background:white;border:1px solid #d4e1de;border-radius:10px;padding:24px;margin:24px 0}.banner{padding:14px 20px;background:#e0f0e9;border-radius:6px;font-weight:650;margin:18px 0}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:24px 0}.card{background:white;border:1px solid #d4e1de;padding:16px;border-radius:8px}.card strong{display:block;font-size:28px}.card span,small,.muted{color:#5f7278;font-size:13px}small{display:block}a{color:#187568}table{width:100%;border-collapse:collapse;font-size:13px}td,th{padding:10px 9px;text-align:left;vertical-align:top;border-bottom:1px solid #dce5e3}th{background:#eaf1ef;white-space:nowrap}.scroll{overflow:auto}.badge{display:inline-block;white-space:nowrap;border-radius:4px;padding:2px 6px;font-size:12px;background:#edf0f0}.found_nonzero{background:#def1e5;color:#145c3d}.running{background:#fff0cf;color:#875e11}.all_endpoints_model_zero{color:#665967}.filters{display:flex;gap:18px;flex-wrap:wrap;align-items:center;margin:18px 0}select{font:inherit;padding:5px;border:1px solid #ccd9d6;border-radius:4px}details{border:1px solid #d4e1de;padding:12px 16px;margin:12px 0;border-radius:6px}summary{font-weight:650;cursor:pointer}.states{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}pre{overflow:auto;background:#f1f5f4;padding:12px;border-radius:5px;font-size:12px;line-height:1.6}.callout{border-left:4px solid #c49a46;background:#fff6e6;padding:14px 20px}li{margin:8px 0}@media(max-width:800px){main{padding:20px 12px}.cards{grid-template-columns:repeat(2,1fr)}.states{grid-template-columns:1fr}section{padding:16px}h1{font-size:25px}}
</style><main><header><div class="kicker">GASTON · ORDINARY DATF · 24 PROCESSES</div><h1>两类完整差分池 × 11 个论文尾段</h1><p>中间段分别取 1 轮、2 轮；按单条拼接评分排序，再复评完整端点，直到每行找到非零或耗尽全部候选。仅省略整个区间的第一个线性层，保留后续线性层及参考轮常数。</p></header>'''
    page+=f'<div class="banner">{banner}</div><div class="cards"><div class="card"><span>全部路径对</span><strong>{search["total_logical_pairs"]:,}</strong><span>44 行均已穷尽</span></div><div class="card"><span>保留的非零匹配</span><strong>{search["total_nonzero_matches"]:,}</strong><span>没有前 32 条截断</span></div><div class="card"><span>已复评的不同完整端点</span><strong>{reassessment["completed_endpoint_jobs"]:,}</strong><span>重复端点跨行复用</span></div><div class="card"><span>已找到完整 DATF 非零的行</span><strong>{counts["found_nonzero"]} / 44</strong><span>模型零 {counts["all_endpoints_model_zero"]} 行，未连接 {counts["no_nonzero_matches"]} 行</span></div></div>'
    page+='<section><h2>如何阅读结果</h2><p>评分 S = pD · ρ̂M · cL²，其中 ρ̂M 是普通 DATF 对中间段的估计。完整端点复评直接从均匀的 320 比特输入计算到最终掩码，不限定差分路径，也不再次乘尾段相关度平方。ρ̂ 是相关度，偏差 ε̂ = ρ̂/2。</p><ul><li>每行按 |S| 递减排列；同分按原差分池顺序及旋转量排序。先评第一名，为零就继续。表中“成功路径排名”是原始非零匹配列表中的排名。</li><li>同一轮区间、输入差分和输出掩码只复评一次；对应的多条路径共享结果。并行任务可以提前计算后续候选，但只有全部更早端点已经证实为模型零，才接受该行首个非零结果。</li><li>“全部端点为模型零”表示本行全部非零匹配对应的不同端点都已复评完；“无非零中间匹配”表示中间段匹配阶段已穷尽且没有连接，两者不能混同。</li><li>两轮 PQ 差分池仍限定为此前 11 类成对 P=Q=0 的六比特模式，保留唯一最大概率首轮前驱和全部 4⁶ 个第二轮输出。单轮差分池输入重量 ≤2。没有加入单轮线性池。</li></ul><div class="callout">完整 DATF 的非零和零都属于普通方法的模型估计。本次没有对真实置换采样，不能把模型零表述为实际相关度为零。各行在首个非零处停止，因此本表也不声称找到了完整 DATF 相关度的全局最优值。</div></section>'
    page+='<section><h2>全量结果与逐行复评</h2><div class="filters"><label>差分段 <select id="family"><option value="">全部</option><option>D1</option><option>D2</option></select></label><label>中间段 <select id="middle"><option value="">全部</option><option>1</option><option>2</option></select></label><label>状态 <select id="status"><option value="">全部</option>'+''.join(f'<option value="{k}">{v}</option>' for k,v in statuses.items())+'</select></label><span id="shown" class="muted"></span></div><div class="scroll"><table id="results"><thead><tr><th>R / 分段</th><th>论文尾段</th><th>已匹配路径对</th><th>非零匹配</th><th>已复评 / 不同端点</th><th>首名 S</th><th>首名完整 DATF</th><th>成功路径排名</th><th>获选完整 DATF</th><th>状态</th><th>详情</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table></div></section>'
    page+=f'<section><h2>运行与验证口径</h2><p>匹配使用 24 个独立 C++ 进程，完整运行 {search["search_seconds"]:.1f} 秒；复评使用 24 个 Python 进程槽位，本次累计运行 {reassessment["elapsed_seconds_this_run"]:.1f} 秒。复评算术后端将均匀二次 Bias 用精确消元求和，非均匀情形仍作精确条件加权；普通 DATF 的变量转换、孤立变量循环和独立性假设均保持不变。</p><p>分片覆盖、44 行计数、全体非零匹配去重、旧 22 行单轮中间结果均已核对。精确求和后端通过 200 个小规模完全求和、200 个原始 Shannon Bias 对照，以及 28 个已有普通 DATF 端点对照。每个复评端点还保存符号 DAG、逐步因子和独立因子审计。</p><p><a href="search_results.json">全量搜索汇总</a> · <a href="summary.csv">表格 CSV</a> · <a href="reassessment/results.json">完整复评汇总</a> · <a href="configuration.json">24 分片配置</a> · <a href="reassessment/configuration.json">复评配置</a></p></section>'
    controls=root/'original_R9_controls.json'
    if controls.exists():
        c=json.loads(controls.read_text())
        page+='<section><h2>长轮数交叉核对</h2><p>另外使用原始、未加速的普通 DATF 对两类差分的九轮端点各复算一例，两项结果均与本次后端一致。<a href="original_R9_controls.json">查看原始九轮对照记录</a>。</p></section>' if len(c)==2 and all(r['status']=='PASS' for r in c) else ''
    supplement=root/'supplemental_factor_audits.json'
    if supplement.exists():
        count=len(json.loads(supplement.read_text()))
        page+=f'<section><h2>附加因子审计</h2><p>{count} 个端点的初始审计触及有偏变量枚举上限，已按互不相交的变量分量分别精确求和完成审计，并用原始 Shannon 后端复算完整端点，结果全部一致。原受限记录仍保留。<a href="supplemental_factor_audits.json">查看补核记录</a>。</p></section>'
    if complete and audited:
        page+='<section><h2>最终完成核验</h2><p><a href="validation.json">逐行排名前缀与全部已计算端点的核验记录</a> · <a href="reproducibility.json">源码和结果的 SHA-256 清单</a></p></section>'
    page+='<section><h2>具体边界与原始候选</h2>'+''.join(details)+'</section><p class="muted">状态为“复评中”的行尚未完成，不能据此给出不存在非零端点的结论。</p></main><script>function update(){let f=document.querySelector("#family").value,m=document.querySelector("#middle").value,s=document.querySelector("#status").value,n=0;document.querySelectorAll("#results tbody tr").forEach(r=>{r.hidden=!!((f&&r.dataset.family!==f)||(m&&r.dataset.middle!==m)||(s&&r.dataset.status!==s));if(!r.hidden)n++});document.querySelector("#shown").textContent=`显示 ${n} / 44 行`}document.querySelectorAll("select").forEach(s=>s.addEventListener("change",update));update();document.querySelectorAll("a[href^=\"#row-\"]").forEach(a=>a.addEventListener("click",()=>document.querySelector(a.getAttribute("href")).open=true));</script></html>'
    (root/'paper_full_datf.html').write_text(page)
    print(root/'paper_full_datf.html')


if __name__=='__main__':main()

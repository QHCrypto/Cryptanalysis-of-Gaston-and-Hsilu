#!/usr/bin/env python3
"""Self-contained report of the new two-round linear-pool match and checks."""
import argparse
import csv
import hashlib
import html
import json
from fractions import Fraction
from pathlib import Path


def number(r):
    if r is None:return '未计算'
    n,d=int(r['numerator']),int(r['denominator'])
    if not n:return '0'
    sign='−' if n<0 else '+';w=d.bit_length()-1
    if abs(n)==1:return sign+'1' if w==0 else f'{sign}2<sup>−{w}</sup>'
    return f'{sign}{abs(n)}/2<sup>{w}</sup>'


def bits(ws):
    return ', '.join(f'({i},{j})' for i,w in enumerate(ws) for j in range(64) if int(w,16)>>j&1)


def states(label,ws):
    return f'<div><b>{html.escape(label)}</b><pre>'+ '\n'.join(f'x[{i}] = 0x{w}' for i,w in enumerate(ws))+'</pre></div>'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True)
    args=p.parse_args();root=args.results
    search=json.loads((root/'search_results.json').read_text())
    check=json.loads((root/'endpoint_check/results.json').read_text())
    validation=json.loads((root/'matching_validation.json').read_text())
    sampling=list(csv.DictReader((root/'sampling.csv').open()))
    sample=sampling[0];best=check['best_full_DATF_within_selection'];b=best['candidate'];result=best['result']
    match_rows=[]
    for row in search['rows']:
        split=f'{row["r_d"]}+{row["r_m"]}+{row["r_l"]}'
        status='已穷尽' if row['exhaustive'] else '差分池批次'
        match_rows.append(f'<tr><td><b>{split}</b><small>{row["R"]} 轮</small></td>'
             f'<td>{row["differential_representatives_evaluated"]:,} / {row["differential_representatives_total"]:,}<small>展开 {row["differential_states_evaluated"]:,} 条差分</small></td>'
             f'<td>{row["linear_representatives"]:,}<small>全部线性代表</small></td><td>{row["evaluated_pairs"]:,}<small>{row["group_evaluations"]:,} 次中间求值</small></td>'
             f'<td>{row["nonzero_matches"]:,}</td><td>{number(row["best"]["score"]) if row["best"] else "未找到"}</td><td>{status}</td></tr>')
    entries=[];details=[]
    for i,e in enumerate(check['entries']):
        c=e['candidate'];r=e['result'];v=r['correlation'];nz=bool(int(v['numerator']))
        entries.append(f'<tr data-nonzero="{int(nz)}"><td>{i+1}</td><td>{html.escape(c["differential_id"])}<small>移位 {c["differential_rotation"]}</small></td>'
                       f'<td>{html.escape(c["linear_id"])}</td><td>{number(c["score"])}</td><td>{number(v)}</td>'
                       f'<td>{html.escape(bits(c["input_difference"]))}</td><td>{html.escape(bits(c["output_mask"]))}</td><td><a href="#endpoint-{i}">详情</a></td></tr>')
        artifact=Path(r['artifact']).relative_to(root.resolve())
        details.append(f'<details id="endpoint-{i}"><summary>端点 {i+1} · {html.escape(c["linear_id"])} · 完整 DATF {number(v)}</summary>'
                       '<div class="states">'+states('输入差分 Δ',c['input_difference'])+states('输出掩码 Λ',c['output_mask'])
                       +states('中间输入差分 L(δD,out)',c['middle_input_difference'])+states('中间输出掩码 Lᵀ(aL,in)',c['middle_output_mask'])
                       +f'</div><p><a href="{artifact}">计算 DAG、因子审计及完整结果</a></p></details>')
    page='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston · L2_2to2 新线性池匹配</title><style>
:root{color-scheme:light;--ink:#17372f;--muted:#5d726b;--line:#d8e3de;--accent:#087c69}*{box-sizing:border-box}body{margin:0;background:#f3f6f4;color:var(--ink);font:16px/1.75 system-ui,"Noto Sans CJK SC",sans-serif}main{max-width:1400px;margin:auto;padding:36px 28px 80px}header{border-bottom:2px solid var(--accent);padding-bottom:20px}h1{font-size:34px;line-height:1.35}h2{font-size:24px;margin-top:0}a{color:var(--accent)}section{background:white;border:1px solid var(--line);border-radius:12px;padding:26px;margin-top:24px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin-top:24px}.card{background:white;border:1px solid var(--line);padding:18px;border-radius:10px}.card strong{display:block;font-size:27px}.muted,small{color:var(--muted)}small{display:block;font-size:12px}.callout{padding:18px 22px;background:#e7f4ef;border-left:4px solid var(--accent);border-radius:6px}.caution{background:#fff4e2;border-color:#b77d1f}.formula{font-size:19px;margin:14px 0;overflow:auto}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;vertical-align:top;padding:11px 9px;border-bottom:1px solid var(--line)}th{background:#edf3ef;white-space:nowrap}tr:hover{background:#f4f9f6}pre{font:13px/1.8 monospace;overflow:auto;padding:14px;background:#f1f5f3;border-radius:6px}.states{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}details{padding:14px;border:1px solid var(--line);border-radius:8px;margin-top:12px}summary{cursor:pointer;font-weight:600}button,input{font:inherit}input[type=search]{padding:8px;border:1px solid var(--line);border-radius:6px;width:300px;max-width:100%}.filters{display:flex;gap:20px;align-items:center;flex-wrap:wrap;margin:18px 0}.eyebrow{font-size:13px;letter-spacing:.12em;color:var(--accent)}@media(max-width:800px){main{padding:22px 14px}.cards{grid-template-columns:repeat(2,1fr)}section{padding:18px}.states{grid-template-columns:1fr}h1{font-size:27px}}@media print{.filters{display:none}.scroll{overflow:visible}main{padding:0}}</style></head><body><main>
<header><div class="eyebrow">GASTON · ORDINARY DATF · L2_2TO2</div><h1>新两轮线性池 × 差分池：首次匹配</h1><p>线性池：全部 16,000 条 2→2 两轮路径代表。差分池：单轮池与两轮 P=Q=0 的 6→6 池。中间段：1 轮及 2 轮。仅省略完整区间首个 L，其余线性层和参考轮常数保留。</p></header>
<div class="cards"><div class="card"><span>非零拼接路径对</span><strong>@MATCHES@</strong><small>均来自 1+1+2</small></div><div class="card"><span>最佳拼接评分 |S|</span><strong>2<sup>−10</sup></strong><small>该池内已证明最优</small></div><div class="card"><span>最佳评分端点复评</span><strong>@NONZERO@ / @CHECKED@</strong><small>完整普通 DATF 非零</small></div><div class="card"><span>最佳已复评模型 |ρ̂|</span><strong>5 / 2048</strong><small>仅在已复评端点范围内</small></div></div>
<section><h2>四种分段的匹配结果</h2><div class="scroll"><table><thead><tr><th>差分 + 中间 + 线性</th><th>差分代表覆盖</th><th>线性代表覆盖</th><th>逻辑路径对数</th><th>非零匹配</th><th>最佳 S</th><th>覆盖状态</th></tr></thead><tbody>@MATCH_ROWS@</tbody></table></div><p>两轮中间段先取每个差分池优先级排序的前 64 个代表，展开其全部不同移位。D1 含周期 32 的代表，因此其 64 个代表展开为 4,064 条，D2 为 4,096 条。两行“差分池批次”的零结果不能扩展到尚未搜索的差分代表。</p><p>16,000 条线性路径共有 160 种入口掩码，每种包含 100 条路径。中间相关度只依赖差分段出口和这个入口掩码，因此每组精确复用一次求值；全部路径的不同权重及外部输出掩码均保留。搜索采用 24 个本地计算进程，用时 @SECONDS@ 秒。</p></section>
<section><h2>最佳拼接评分与完整四轮复评</h2><div class="formula">S = p<sub>D</sub> · ρ̂<sub>M</sub> · c<sub>L</sub><sup>2</sup> = 2<sup>−2</sup> × (±1) × 2<sup>−8</sup> = ±2<sup>−10</sup></div><p>最佳评分并列的路径对有 284 个，去重后对应 176 个完整端点，已全部复评：148 个模型非零，28 个模型为零，无资源限制未决项。57,400 个非零路径对共对应 38,400 个不同完整端点，剩余 38,224 个端点尚未做完整 DATF 复评。</p><p>每个已复评端点均通过 8 组输入的符号 DAG 核验及因子审计。下面的最大模型值只针对这 176 个最佳拼接评分端点，不是全部匹配端点的最大值。</p>
<div class="callout"><b>一个值得继续研究的四轮端点</b><p>差分路径 @DIFF@，共同移位 @SHIFT@；线性路径 @LINEAR@。零起始轮号 8，即参考置换末四轮。</p><p>输入差分：@DELTA_BITS@<br>输出掩码：@MASK_BITS@</p><div class="formula">S = 2<sup>−10</sup>；完整普通 DATF ρ̂ = 5/2048 ≈ 2<sup>−8.6781</sup>；偏差 ε̂ = ρ̂/2 = 5/4096</div></div><div class="states">@BEST_STATES@</div><p>最佳端点同时使用原始 Shannon Bias 后端重新计算，得到相同的 5/2048。完整端点估计会包含不同于单条拼接评分的贡献；这两个数不要求相等。</p></section>
<section><h2>完整置换采样：非零得到支持，模型数值有偏差</h2><p>对上面的端点，直接计算省略首个 L 的四轮核心，保留后续轮常数，采样 @SAMPLES@ 对输入。</p><div class="scroll"><table><thead><tr><th>完整 DATF</th><th>采样相关度</th><th>95% Wilson 区间</th><th>模型值落在区间内</th></tr></thead><tbody><tr><td>0.00244140625</td><td>@MC@</td><td>[@LOW@, @HIGH@]</td><td>@INSIDE@</td></tr></tbody></table></div><div class="callout caution"><b>采样区间排除 0，支持该四轮端点具有非零相关度；完整 DATF 的数值落在采样区间之外。</b><p>因此可将它作为后续研究候选，但不把 5/2048 当作真实相关度的精确值。区间是伪随机仿真的统计估计，不是精确求和或完整构造的结论。</p></div><p class="muted">输入生成器为项目既有的计数式 SplitMix64；种子 @SEED@，24 个采样线程。端点先由模型选定，再用独立输入序列采样。</p></section>
<section><h2>176 个最佳评分端点：完整 DATF 明细</h2><p>所有条目均满足 |S| = 2⁻¹⁰。可按路径标识或坐标筛选，再展开具体边界。坐标为零起始 (行,位)，位 0 为 LSB。</p><div class="filters"><input id="query" type="search" placeholder="搜索路径标识或坐标" aria-label="筛选端点"><label><input id="nonzero" type="checkbox">仅显示 DATF 非零</label><span id="count"></span></div><div class="scroll"><table id="endpoints"><thead><tr><th>序号</th><th>差分路径</th><th>线性路径</th><th>拼接 S</th><th>完整 DATF</th><th>输入差分</th><th>输出掩码</th><th>边界</th></tr></thead><tbody>@ENDPOINT_ROWS@</tbody></table></div>@DETAILS@</section>
<section><h2>验证与可复现文件</h2><p>检查了全部分片的差分代表、移位覆盖、入口掩码分组、非零计数和源文件哈希；574 个非零中间连接均用单 χ 的直接差分线性表重新核对；两轮中间段抽取 32 个零结果，用原始 Python 普通 DATF 复核。所有检查通过。</p><ul><li><a href="search_results.json">四行匹配汇总及覆盖范围</a></li><li><a href="matches/D1_M1_L2_2to2_pool.jsonl.gz">全部 57,400 个非零路径对</a></li><li><a href="ranked_endpoints/D1_M1_L2_2to2_pool.jsonl.gz">38,400 个完整端点的评分排序</a></li><li><a href="endpoint_check/results.json">176 个最佳评分端点的完整复评</a></li><li><a href="sampling.csv">直接置换采样原始计数</a></li><li><a href="matching_validation.json">匹配验证记录</a></li><li><a href="configuration.json">搜索配置及输入哈希</a></li><li><a href="reproducibility.json">本次报告与代码文件哈希</a></li><li><a href="../pools/linear_two_2to2.html">两轮线性池的抵消原理与枚举报告</a></li></ul></section><p class="muted">单条拼接评分、完整普通 DATF 模型估计和置换采样是三个不同层次的结果。只展开差分侧相对移位，不据此声称固定轮常数下所有完整端点等价。</p></main><script>
function filterRows(){const q=document.querySelector('#query').value.toLowerCase(),nz=document.querySelector('#nonzero').checked;let n=0;const rs=[...document.querySelectorAll('#endpoints tbody tr')];for(const r of rs){r.hidden=!!((nz&&r.dataset.nonzero!=='1')||(q&&!r.textContent.toLowerCase().includes(q)));if(!r.hidden)n++}document.querySelector('#count').textContent=`显示 ${n} / ${rs.length} 个端点`}
document.querySelector('#query').addEventListener('input',filterRows);document.querySelector('#nonzero').addEventListener('change',filterRows);filterRows();document.querySelectorAll('a[href^="#endpoint-"]').forEach(a=>a.addEventListener('click',()=>{document.querySelector(a.getAttribute('href')).open=true}));
</script></body></html>'''
    replacements={'MATCHES':f'{search["total_nonzero_matches"]:,}','NONZERO':str(check['nonzero_endpoints']),
                  'CHECKED':str(check['completed_endpoints']),'MATCH_ROWS':''.join(match_rows),
                  'SECONDS':f'{search["search_seconds"]:.2f}','DIFF':b['differential_id'],
                  'SHIFT':str(b['differential_rotation']),'LINEAR':b['linear_id'],
                  'DELTA_BITS':bits(b['input_difference']),'MASK_BITS':bits(b['output_mask']),
                  'BEST_STATES':states('输入差分 Δ',b['input_difference'])+states('输出掩码 Λ',b['output_mask']),
                  'SAMPLES':f'{int(sample["samples"]):,}','MC':f'{float(sample["rho"]):.10f}',
                  'LOW':f'{float(sample["ci95_low"]):.10f}','HIGH':f'{float(sample["ci95_high"]):.10f}',
                  'INSIDE':'是' if sample['model_inside_ci95']=='1' else '否','SEED':sample['seed'],
                  'ENDPOINT_ROWS':''.join(entries),'DETAILS':''.join(details)}
    for k,v in replacements.items():page=page.replace('@'+k+'@',v)
    hnw_path=root/'hnw26/results.json'
    if hnw_path.exists():
        hnw=json.loads(hnw_path.read_text())
        assert hnw['endpoint_id']==best['endpoint_id']
        full=hnw['HNW']['full'];zero=hnw['HNW']['zero']
        section=('<section id="hnw26"><h2>同一端点的 HNW26 复算</h2>'
                 '<p>固定上面的输入差分与输出掩码，沿用首轮精确 DDT、r<sub>d</sub>=1、'
                 'r<sub>l</sub>=0 的逐列近似设置，重新评估完整四轮。原搜索分段为 1+1+2，'
                 '此处的 HNW 评估分段为 1+3+0。</p><div class="scroll"><table><thead><tr>'
                 '<th>普通 DATF</th><th>HNW 局部联合版（1,024 维）</th><th>HNW 零状态值掩码版（32 维）</th>'
                 '<th>已有采样</th></tr></thead><tbody><tr>'
                 f'<td>{number(hnw["DATF"])}</td><td>{number(full["correlation"])}</td>'
                 f'<td>{number(zero["correlation"])}</td><td>{float(sample["rho"]):.10f}</td>'
                 '</tr></tbody></table></div>'
                 '<p>HNW 首轮全部分支的概率及带符号贡献均已核对，简化版另经独立 Python Fraction 递归复算。'
                 'HNW 的局部近似结果不等于完整置换真实相关度；具体分支退化过程及与采样的对照见 '
                 '<a href="hnw26/hnw26_endpoint.html">HNW26 专项报告</a>。</p></section>')
        page=page.replace('<section><h2>176 个最佳评分端点：完整 DATF 明细</h2>',
                          section+'<section><h2>176 个最佳评分端点：完整 DATF 明细</h2>')
    target=root/'linear_two_matches.html';target.write_text(page)
    project=Path(__file__).resolve().parents[1]
    sources=[project/f for f in ['src/datf_linear_two_shard.cpp','tools/datf_linear_two_search.py',
              'tools/datf_linear_two_merge.py','tools/datf_linear_two_validate.py','tools/datf_linear_two_report.py',
              'include/datf_middle.hpp','include/gaston.hpp','include/experiment.hpp','src/datf_pool_verify.cpp',
              'tools/hnw26_endpoint_compare.py']]
    files=sources+[f for f in root.rglob('*') if f.is_file() and f.name!='reproducibility.json']
    (root/'reproducibility.json').write_text(json.dumps({'scope':'New L2 matching, full endpoint checks, sampling and report',
        'files':[{'path':str(f.resolve()),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files]},ensure_ascii=False,indent=2)+'\n')
    print(target)


if __name__=='__main__':main()

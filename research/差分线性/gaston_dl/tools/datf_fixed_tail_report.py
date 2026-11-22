#!/usr/bin/env python3
"""Chinese comparison report for ordinary DATF on saved 1+1+n endpoints."""
import argparse
import csv
import hashlib
import html
import json
import os
import platform
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def h(value):
    return html.escape(str(value))


def number(record):
    if record is None:
        return '未计算完成'
    n,d = int(record['numerator']),int(record['denominator'])
    if not n:
        return '0'
    if d==1:
        return str(n)
    if d & (d-1)==0:
        k=d.bit_length()-1
        return (f'2^-{k}' if n==1 else f'{n}/2^{k}')
    return str(Fraction(n,d))


def table(headers, rows):
    def row(cells, tag):
        return '<tr>'+''.join(f'<{tag}>{h(c)}</{tag}>' for c in cells)+'</tr>'
    return '<div class="scroll"><table><thead>'+row(headers,'th')+'</thead><tbody>'+''.join(row(r,'td') for r in rows)+'</tbody></table></div>'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    args=parser.parse_args()
    folder=args.results.resolve()
    data=json.loads((folder/'results.json').read_text())
    validation=json.loads((folder/'validation.json').read_text())
    assert validation['status']=='PASS'
    jobs=data['results']
    full=[r for r in jobs if r['kind']=='full_endpoint']
    prefixes={(r['tail_id'],r['candidate_rank']):r for r in jobs if r['kind']=='prefix_control'}
    rows=[]; best=[]; summaries=[]; details=[]; traces=[]
    for r in full:
        p=prefixes[r['tail_id'],r['candidate_rank']]
        logarithm=r['correlation']['log2_abs'] if r['correlation'] else None
        row=[r['rounds'],f'1 + 1 + {r["r_l"]}',r['tail_id'],r['candidate_rank'],
             number(p['exact_prefix_correlation']),number(p['correlation']),
             number(r['previous_splice_score']),number(r['correlation']),
             '—' if logarithm is None else f'{logarithm:.6f}']
        rows.append(row)
        if r['was_best_by_round']:
            best.append(row)
        traces.append([r['tail_id'],r['rounds'],
                       ' × '.join(number({'numerator':t['factor_numerator'],'denominator':t['factor_denominator']}) for t in r['trace']),
                       r['engine_statistics'].get('bias_nodes',0),r['status'],r.get('factor_audit',{}).get('status','未复核')])
        summaries.append({'R':r['rounds'],'rd':1,'rm':1,'rl':r['r_l'],'tail_id':r['tail_id'],
                          'candidate_rank':r['candidate_rank'],'was_best_by_round':r['was_best_by_round'],
                          'prefix_exact_rho':number(p['exact_prefix_correlation']),
                          'prefix_DATF_rho':number(p['correlation']),
                          'prior_splice_score':number(r['previous_splice_score']),
                          'full_DATF_rho':number(r['correlation']),'full_DATF_log2_abs':logarithm,
                          'full_DATF_bias':number(r['bias']),'true_full_correlation':'unknown',
                          'start_round_index':r['start_round_index'],
                          'remaining_constants':';'.join(r['remaining_constants']),
                          'input_words':';'.join(r['input_difference']['words']),
                          'join_words':';'.join(r['join_mask']['words']),
                          'output_words':';'.join(r['output_mask']['words'])})
        details.append(f'<details><summary>{h(r["tail_id"])}，候选 {r["candidate_rank"]}：端点和常数</summary>'
                       +table(['边界','五行 64 位字，bit 0 为最低位'],[
                           ['第一个 χ 输入差分 Δ','\n'.join(r['input_difference']['words'])],
                           ['两轮前缀输出掩码 γ','\n'.join(r['join_mask']['words'])],
                           ['最终输出掩码 λ','\n'.join(r['output_mask']['words'])]])
                       +f'<p>起始轮号 s={r["start_round_index"]}（从 0 起），首个 χ 之后的常数：{h(", ".join(r["remaining_constants"]))}。</p>'
                       +f'<p><a href="{h(p["id"])}.json.gz">两轮前缀完整代数 DAG</a> · '
                       +f'<a href="{h(r["id"])}.json.gz">完整端点代数 DAG</a>（gzip JSON）</p></details>')
    with (folder/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summaries[0]))
        writer.writeheader();writer.writerows(summaries)
    headers=['总轮数','分段','尾段','候选','前缀精确 ρ','前缀 DATF ρ̂','原拼接评分','整体 DATF ρ̂','log₂|ρ̂|']
    nonzero=[r for r in full if r['correlation'] and int(r['correlation']['numerator'])]
    zeros=[r for r in full if r['correlation'] and not int(r['correlation']['numerator'])]
    page=f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston：1+1+n 候选的普通 DATF 相关度</title>
<style>
body{{margin:0;background:#f3f6f8;color:#243b49;font:16px/1.8 system-ui,"Noto Sans CJK SC",sans-serif}}
main{{max-width:1450px;margin:28px auto;padding:36px;background:white;border-radius:12px}}
h1{{font-size:29px;line-height:1.5}}h2{{font-size:22px;margin-top:30px}}p{{max-width:1150px}}
.muted{{color:#617783}}.note{{background:#eaf5f3;border-left:4px solid #268071;padding:14px 20px}}
.scope{{background:#fff6e8;border-left:4px solid #ad803d;padding:14px 20px}}
.formula{{background:#eff3f7;padding:16px;font-family:ui-monospace,monospace;overflow:auto}}
.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;margin:18px 0;font-size:14px;font-variant-numeric:tabular-nums}}
th,td{{padding:10px 12px;border-bottom:1px solid #dce5ea;white-space:pre;text-align:left;vertical-align:middle}}
th{{background:#edf3f6;border-top:2px solid #395364}}tbody tr:last-child td{{border-bottom:2px solid #395364}}
a{{color:#086d84}}summary{{cursor:pointer;font-weight:650}}details{{margin:20px 0}}code{{background:#edf2f5;padding:2px 5px}}
@media(max-width:700px){{main{{margin:0;padding:18px}}h1{{font-size:24px}}}}
</style></head><body><main>
<h1>Gaston：1 + 1 + n 候选的普通 DATF 评估</h1>
<p class="muted">差分段 1 轮 + 中间连接 1 轮 + 固定线性尾段 n 轮 · n=2,…,6 · 仅省略首个 L</p>
<p class="note">本次评估 {len(full)} 个完整端点和 {len(prefixes)} 个两轮前缀对照。
完整端点中 {len(nonzero)} 项普通 DATF 估计非零，{len(zeros)} 项为模型零。
原四轮最佳端点得到 <strong>ρ̂ = 7/2⁴³ ≈ 2<sup>−40.192645</sup></strong>，偏差 ε̂=ρ̂/2=7/2⁴⁴。
所有前缀对照均与已有精确值一致；完整端点的实际相关度尚未验证。</p>
<h2>此前汇总表的五个最佳端点</h2>
{table(headers,best)}
<p>“此前最佳”按原先的两轮前缀精确相关度乘固定线性路径相关度平方排序。
本次没有重新搜索端点，也没有声称上表是 DATF 在全部候选中的最优结果。</p>
<h2>全部已评估尾段</h2>
{table(headers,rows)}
<p>从 <a href="../fixed_paper_tail/table13_fixed_tail.html">原固定论文尾段结果</a> 中选取每个非空尾段的前 {data['top_per_tail']} 个候选。
未评估的空候选尾段：{h(', '.join(data['skipped_empty_tails']))}。这些尾段此前没有保存可拼接端点，不能记作本次 DATF 零结果。
原结果保存的其他排名候选也不在本次默认评估范围内。</p>
<h2>四轮非零结果的含义</h2>
<p>四轮 D3_chi2_to_3 的输入差分仍为第 4 行第 1 位，即 Δ=(0,0,0,0,0x2)。
其两轮前缀相关度精确为 1/2，固定线性尾段平方权重为 38，旧拼接评分因此为 2⁻³⁹。
本次直接对整个四轮输入差分与最终输出掩码构造普通 DATF，得到独立的整体模型估计。</p>
<div class="formula">ρ̂₄ = 2⁻²⁰ × (7/2²²) × 2⁻¹ = 7/2⁴³<br>ρ̂₄ / 2⁻³⁹ = 7/16 = 0.4375</div>
<p>上式的三个因子对应原 Algorithm 1 的多项式剥离步骤，不分别等同于三个轮段的权重。
DATF 估计绝对值为原评分的 43.75%，负对数权重增加约 1.192645 比特。
这个比较不证明 DATF 更接近真实值，也不把旧评分当作上下界。</p>
<h2>计算方法与验证</h2>
<p>复用 <a href="../datf_k06/datf_k06.html">K06 评估</a> 中同一套普通 DATF 引擎，
采用原始孤立变量 Algorithm 1，未启用增强 DATF。输入是全部均匀的 320 比特状态，未施加差分路径事件条件。
完整端点直接评估至最后一个 χ 的输出，不把前缀 DATF 值再乘线性权重。</p>
<p>R 轮边界为 χ → (L,iota,χ) → … → (L,iota,χ)，含 R 个 χ、R−1 个 L。
首个常数通过均匀输入平移吸收，后续常数沿用原候选：s=12−R。
两轮前缀对照使用所属完整候选的常数，不独立改用末两轮常数。</p>
{table(['尾段','R','估计器因子','Bias 节点数','计算状态','独立因子审计'],traces)}
<p>计算采用精确有理数。每项上限 {data['resource_limit']['max_bias_nodes_per_job']:,} 个 Bias 节点，超限返回未知，不返回零。
独立复核使用二次消元，并在存在非均匀边际时显式加权，或构造保持全部非均匀坐标不变的平衡平移。
该平移只核对独立模型因子为零，不参与 DATF 估计，也不意味着真实输入空间存在相同平移。</p>
<p>验证通过：重新精确计算 {validation['exact_prefix_controls_recomputed']} 个两轮前缀；
构造时每项检查 32 组完整状态及差分，共 {sum(r['dag_verified_samples'] for r in jobs)} 组；
从压缩文件恢复 DAG 后再检查 {validation['serialized_DAG_parity_checks']} 组最终奇偶；
完成 {validation['factor_audits_completed']} 项独立因子审计。
另外通过小规模完整真值表核对平衡证书，并验证资源不足状态不会被写成零。</p>
<p class="scope"><strong>结论范围</strong>：两轮前缀的 1/2 或 1 有既有精确计算作为依据。
完整四轮的 7/2⁴³ 是普通 DATF 的非零模型估计；其他端点的 0 是普通 DATF 的模型零。
它们都不是完整置换的实测值或精确值，不能据此宣布五至八轮实际相关度为零。</p>
<h2>端点与复现数据</h2>
{''.join(details)}
<p><a href="results.json">结果、常数、全部估计轨迹与审计</a> · <a href="summary.csv">汇总 CSV</a> ·
<a href="validation.json">验证记录</a> · <a href="reproducibility.json">SHA-256 清单</a> ·
<a href="../../README.md">项目说明与复现命令</a></p>
</main></body></html>'''
    (folder/'datf_one_one_tail.html').write_text(page)
    paths=[*folder.glob('*.json'),*folder.glob('*.json.gz'),folder/'summary.csv',folder/'datf_one_one_tail.html',
           ROOT/'README.md',ROOT/'Makefile',Path(data['source']['path']),
           *[ROOT/'tools'/n for n in ('datf.py','datf_k06.py','datf_audit.py','datf_fixed_tail.py',
                                      'datf_fixed_tail_report.py','fixed_tail_search.py','linear_seed_scan.py')],
           ROOT/'tests/test_datf_fixed_tail.py']
    manifest={'python':platform.python_version(),'method':data['method'],
              'sha256':{os.path.relpath(p,folder):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in paths if p.name!='reproducibility.json'}}
    (folder/'reproducibility.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(folder/'datf_one_one_tail.html')


if __name__=='__main__':
    main()

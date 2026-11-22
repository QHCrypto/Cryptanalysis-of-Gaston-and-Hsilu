#!/usr/bin/env python3
"""Generate the ordinary DATF K06 report, with explicit model/true distinction."""
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
        return '未完成'
    n, d = int(record['numerator']), int(record['denominator'])
    if not n:
        return '0'
    if abs(n) & (abs(n) - 1) == 0 and d & (d - 1) == 0:
        exponent = abs(n).bit_length() - d.bit_length()
        return ('−' if n < 0 else '') + f'2^{exponent}'
    return str(Fraction(n, d))


def table(headers, rows):
    def row(cells, tag):
        return '<tr>' + ''.join(f'<{tag}>{h(c)}</{tag}>' for c in cells) + '</tr>'
    return ('<div class="scroll"><table><thead>' + row(headers, 'th')
            + '</thead><tbody>' + ''.join(row(cells, 'td') for cells in rows)
            + '</tbody></table></div>')


def job_name(job):
    if job['rounds'] == 3:
        return f'完整三轮前缀（接 {job["parent_total_rounds"]} 轮候选）'
    return f'完整 {job["rounds"]} 轮端点'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    folder = args.results.resolve()
    data = json.loads((folder / 'results.json').read_text())
    audits = json.loads((folder / 'zero_audit.json').read_text())
    validation = json.loads((folder / 'validation.json').read_text())
    assert validation['status'] == 'PASS'
    jobs = data['results']
    rows, traces, endpoint_sections, summary = [], [], [], []
    for job in jobs:
        factors = [number({'numerator': t['factor_numerator'], 'denominator': t['factor_denominator']})
                   for t in job['trace']]
        rows.append([job_name(job), job['start_round_index'],
                     ', '.join(job['round_constants_after_first_chi']),
                     number(job['previous_restricted_score']), number(job['correlation']),
                     number(job['bias']), '尚未确定'])
        traces.append([job_name(job), ' × '.join(factors), job['engine_statistics']['bias_nodes'],
                       f'{job["elapsed_seconds"]:.3f} s', job['status']])
        endpoint_sections.append('<details><summary>' + h(job_name(job)) + '：精确端点与中间规模</summary>'
            + table(['边界', '五个 64 位行字（x[0] 至 x[4]，bit 0 为最低位）'],
                    [['第一个 χ 的输入差分 Δ', '\n'.join(job['input_difference']['words'])],
                     ['最后一个 χ 的输出掩码 λ', '\n'.join(job['output_mask']['words'])]])
            + table(['已计算 χ 层数', 'DAG 变量', 'α 变量', 'β 变量', '差分总项数', '每比特最多项数'],
                    [[s[k] for k in ('chi_layers', 'dag_variables', 'alpha_variables', 'beta_variables',
                                     'difference_terms', 'max_bit_difference_terms')] for s in job['stages']])
            + f'<p><a href="{h(job["id"])}.json">下载本项结果、目标多项式和全部变量定义</a></p></details>')
        summary.append({'id': job['id'], 'rounds': job['rounds'], 'start_round_index': job['start_round_index'],
                        'omitted_first_L': True, 'input_words': ';'.join(job['input_difference']['words']),
                        'output_mask_words': ';'.join(job['output_mask']['words']),
                        'prior_restricted_contribution': number(job['previous_restricted_score']),
                        'ordinary_DATF_rho': number(job['correlation']),
                        'ordinary_DATF_epsilon': number(job['bias']),
                        'true_correlation_status': 'unknown', 'status': job['status']})
    with (folder / 'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    audit_rows = []
    for audit in audits:
        for factor in audit['factor_audits']:
            branches = factor['conditioned_branches']
            audit_rows.append([audit['id'], factor['step'], factor['variables'], factor['terms'],
                               len(factor['biased_variable_marginals']), number(factor['independent_quadratic_correlation']),
                               '；'.join('h={' + ','.join(map(str, b['balance_certificate']['shift_variable_ids'])) + '}'
                                         if b['balance_certificate'] else '非零因子'
                                         for b in branches)])
    result_table = table(['评估对象', '起始轮号 s（从 0 起）', '首个 χ 之后的常数',
                          '原限定贡献', '普通 DATF：ρ̂', '普通 DATF：ε̂', '真实相关度'], rows)
    trace_table = table(['对象', '原算法依次计算的因子', 'Bias 递归节点', '计算耗时', '执行状态'], traces)
    audit_table = table(['对象', '步骤', '变量数', '项数', '非均匀边际数', '独立二次消元结果',
                         '零因子平衡平移证书（变量编号）'], audit_rows)
    page = f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Gaston K06：普通 DATF 三轮前缀与五、六轮评估</title>
<style>
body{{margin:0;background:#f3f6f8;color:#243b49;font:16px/1.8 system-ui,"Noto Sans CJK SC",sans-serif}}
main{{max-width:1400px;margin:28px auto;padding:36px;background:#fff;border-radius:12px}}
h1{{font-size:29px;line-height:1.5}}h2{{font-size:22px;margin-top:32px}}p{{max-width:1150px}}
.muted{{color:#617783}}.note{{background:#edf5f8;border-left:4px solid #327b98;padding:14px 20px}}
.scope{{background:#fff6e8;border-left:4px solid #ad803d;padding:14px 20px}}
.formula{{padding:16px;background:#eff3f7;font-family:ui-monospace,monospace;overflow:auto}}
.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;margin:18px 0;font-size:14px;font-variant-numeric:tabular-nums}}
th,td{{padding:10px 12px;border-bottom:1px solid #dce5ea;text-align:left;vertical-align:middle;white-space:pre}}
th{{background:#edf3f6;border-top:2px solid #395364}}tbody tr:last-child td{{border-bottom:2px solid #395364}}
tr:hover td{{background:#f6f9fb}}a{{color:#086d84}}code{{background:#edf2f5;padding:2px 5px}}
summary{{cursor:pointer;font-weight:650}}details{{margin:20px 0}}li{{margin:8px 0}}pre{{overflow:auto}}
@media(max-width:700px){{main{{margin:0;padding:18px}}h1{{font-size:24px}}}}
</style></head><body><main>
<h1>Gaston K06：普通 DATF 评估</h1>
<p class="muted">K06 / 22 · 完整三轮前缀 + 完整五轮、六轮端点 · 仅省略第一个线性层</p>
<p class="note"><strong>四项普通 DATF 计算均完成，模型相关度均为 0</strong>。
计算使用精确有理数，未触及每项 {data['resource_limit']['max_bias_nodes_per_job']:,} 个 Bias 节点的资源上限；
独立二次多项式消元得到同样结果。它们不是采样结果，也不是实际置换相关度为零的证明。</p>
<h2>结果表</h2>
{result_table}
<p>ρ = E<sub>X</sub>[ (−1)<sup>λ·(F(X)⊕F(X⊕Δ))</sup> ]，ε = Pr[λ·(F(X)⊕F(X⊕Δ))=0]−1/2 = ρ/2。
帽号表示普通 DATF 模型估计。原限定贡献列沿用此前相关度口径，不能作为真实相关度或其上下界。</p>
<p>三轮前缀分两项评估，是因为五轮与六轮候选的起始轮号不同。两项使用相同 Δ 和 γ，
但前缀常数分别为 (0x78, 0x69) 和 (0x87, 0x78)。五、六轮整体计算直接采用对应最终输出掩码，
没有再将三轮 DATF 结果乘以固定线性路径的平方。</p>
<h2>轮边界与输入分布</h2>
<div class="formula">X₁ = χ(X)<br>Xᵣ₊₁ = χ(L(Xᵣ) ⊕ RC[s+r])，r = 1,…,R−1<br>F(X) = Xᴿ，X 均匀遍历全部 320 比特</div>
<p>输入位于第一个 χ 之前，输出位于最后一个 χ 之后。R 轮含 R 个 χ、R−1 个 L。
首个常数通过均匀输入平移吸收；后续常数全部保留。三轮前缀没有限制第一轮输出 B1 或第二轮输出 B2，
也没有只选择此前匹配的两条分支。</p>
<h2>普通 DATF 的实现范围</h2>
<p>程序同时保存基准状态 a 和差分 d，利用布尔环中的 χ 关系：</p>
<div class="formula">a′ᵢ = aᵢ ⊕ aᵢ₊₂ ⊕ aᵢ₊₁aᵢ₊₂<br>
d′ᵢ = dᵢ ⊕ dᵢ₊₂ ⊕ aᵢ₊₁dᵢ₊₂ ⊕ dᵢ₊₁aᵢ₊₂ ⊕ dᵢ₊₁dᵢ₊₂</div>
<p>每次在后续 χ 之前完成 L 和常数加法，然后引入 α、β 过渡变量。一个表达式含两个以上变量时进行替换；
零变量、单变量表达式保留，相同或互补的 β 表达式复用变量。保存的全部定义构成精确代数 DAG。</p>
<p>估计采用原始 Algorithm 1：把当前多项式拆成非孤立部分 e* 和孤立变量之和 eₗ，
计算 e* 在独立边际分布下的相关度并乘入累计结果；同时回代 eₗ 中所有过渡变量，再重复。
输入变量和 α 的边际取均匀，β 的边际用同一原算法递归估计。
这里会假设变量及先后剥离的因子独立。没有采用增强 DATF 的转换阈值、最新变量回溯、输入空间划分或增强 Bias 预处理。</p>
<p>转换规则参考 <a href="https://github.com/QHCrypto/DATF">原作者公开实现</a>，固定提交
<code>e5759a5fbdad660ead7d0058914c5a46b79e0f0f</code>；估计循环核对
<a href="https://eprint.iacr.org/2025/1632">2025/1632</a> 第 2.4 节重述的原始 Algorithm 1。
这是面向 Gaston 的独立实现，不是对作者 Ascon 专用脚本末尾手工 β 展开过程的逐行移植。
参考文件及 MIT 许可证保存在 <a href="../../vendor/datf_reference/README.md">vendor/datf_reference</a>。</p>
<h2>零因子的来源及独立复核</h2>
{trace_table}
<p>上表的因子按原算法剥离次序列出，不分别对应差分段、中间段和线性段。
三轮前缀先得到 2⁻²²，回代孤立变量后得到含 131 个变量、187 项的二次多项式，
其独立均匀模型下相关度严格为 0。独立消元还找到了平衡平移 h：f(z⊕h)=f(z)⊕1。</p>
<p>该零因子包含 13 个第一层 β 和 11 个原始输入坐标。它们共同作为 24 个独立坐标参与模型计算，
但其实际仿射秩仅为 12，因此有 12 条线性关系被独立性假设忽略。这说明模型确实丢失了依赖，
不能由此认定实际相关度也是零；也尚未证明这些关系就是偏差估计失效的唯一原因。</p>
<p>复核程序对每个二次因子采用独立的消元恒等式
E<sub>u,v</sub>(−1)<sup>uv+uA+vB+C</sup> = ½(−1)<sup>C+AB</sup>。
五轮零因子内有一个非均匀边际变量，先分别固定其值，再按精确边际概率加权；两支均获得平衡证书。
该复核只核对估计器内部计算，不改变 DATF 或恢复被忽略的相关性。</p>
<details><summary>展开每个因子的独立审计与平衡证书</summary>{audit_table}</details>
<p class="scope"><strong>本次结论：普通 DATF 在当前转换边界与这些端点上未给出非零估计</strong>。
先前两条负路径贡献之和 −2⁻²³，五、六轮限定拼接贡献 −2⁻⁶¹、−2⁻⁸³，仍是此前限定范围的数值。
它们与这次 0 不构成数值矛盾。当前结果既不能证明其他贡献恰好抵消，也不能证明原限定贡献接近真实相关度。
真实三轮、五轮、六轮相关度均仍待确定。</p>
<h2>验证与可复现文件</h2>
<ul>
<li>全部 32×32=1024 个单 χ 差分/掩码组合，与直接遍历 32 个输入的真值表一致。</li>
<li>全部 320 个 L 的基向量，符号传播与整数实现一致。</li>
<li>40 组加权 Bias、80 组均匀二次消元、40 组加权二次消元，均与完整枚举一致。</li>
<li>文献比较例复现：原 Algorithm 1 估计 ρ̂=0，而完全展开的精确值为 ρ=1/4。
另一个既有两轮正对照得到 ρ̂=1/2。</li>
<li>四项各 64 组完整状态与差分检查，共 256 组，精确 DAG 与直接置换一致；
四项全部因子经过独立代数审计。</li>
</ul>
{''.join(endpoint_sections)}
<p><a href="results.json">全部结果与估计过程</a> · <a href="summary.csv">汇总 CSV</a> ·
<a href="zero_audit.json">独立代数审计</a> · <a href="validation.json">测试记录</a> ·
<a href="reproducibility.json">来源与 SHA-256</a> · <a href="../../README.md">项目说明和复现命令</a> ·
<a href="../pq_kernel_six_all/table13_pq_kernel.html">此前限定路径结果</a></p>
</main></body></html>'''
    (folder / 'datf_k06.html').write_text(page, encoding='utf-8')
    files = [*sorted(folder.glob('*.json')), folder / 'summary.csv', folder / 'datf_k06.html',
             ROOT / 'README.md', ROOT / 'Makefile', ROOT / 'include/gaston.hpp',
             *[ROOT / 'tools' / name for name in ('datf.py', 'datf_k06.py', 'datf_audit.py',
                                                  'datf_report.py', 'linear_seed_scan.py', 'fixed_tail_search.py')],
             ROOT / 'tests/test_datf.py', Path(data['source']['path']),
             *sorted((ROOT / 'vendor/datf_reference').glob('*'))]
    manifest = {'python': platform.python_version(), 'platform': platform.platform(),
                'method': data['method'], 'true_correlation': 'unknown',
                'sources': {'original_DATF': 'https://github.com/QHCrypto/DATF',
                            'original_Algorithm_1_recap': 'https://eprint.iacr.org/2025/1632#page=10'},
                'sha256': {os.path.relpath(p, folder): hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in files if p.is_file() and p.name != 'reproducibility.json'}}
    (folder / 'reproducibility.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(folder / 'datf_k06.html')


if __name__ == '__main__':
    main()

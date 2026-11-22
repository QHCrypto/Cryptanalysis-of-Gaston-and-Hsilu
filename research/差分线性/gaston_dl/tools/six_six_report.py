#!/usr/bin/env python3
"""Report the restricted 2+1+rl experiment without claiming full correlation."""
import argparse
import csv
import html
import json
from pathlib import Path


def h(x):
    return html.escape(str(x))


def row(cells):
    return '<tr>' + ''.join(f'<td>{h(c)}</td>' for c in cells) + '</tr>'


def table(headers, rows):
    return '<div class="scroll"><table><thead>' + row(headers) + '</thead><tbody>' + ''.join(rows) + '</tbody></table></div>'


def signed_power(value):
    if value is None:
        return '未找到非零和'
    if not value['numerator']:
        return '0（限定分支抵消）'
    return ('−' if value['numerator'] < 0 else '+') + f'2^{int(value["log2_abs"])}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    root = args.results.resolve()
    data = json.loads((root / 'search.json').read_text())
    with (root / 'conditional_experiments.csv').open() as f:
        experiments = list(csv.DictReader(f))
    main_rows, summary = [], []
    for r in data['reports']:
        if r['best_restricted_splice_score']:
            cells = [r['R'], '2 + 1 + ' + str(r['r_l']), '6 → 6', r['id'], r['best_rotation'],
                     '2^-12 × 2^-12', '−1', r['tail_w_squared'], signed_power(r['best_restricted_splice_score'])]
            main_rows.append(row(cells))
            summary.append({'R': r['R'], 'rd': 2, 'rm': 1, 'rl': r['r_l'], 'tail': r['id'], 'rotation': r['best_rotation'],
                            'first_chi_active': 6, 'second_chi_active': 6, 'differential_weight': 24,
                            'middle_correlation': -1, 'tail_squared_weight': r['tail_w_squared'],
                            'restricted_score_sign': -1, 'log2_abs_restricted_score': r['best_restricted_splice_score']['log2_abs']})
    with (root / 'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)
    selected = table(['R', 'rd + rm + rl', '差分活跃数', '固定尾段', '移位', '精确差分路径概率', '中间 g', 'wL²', '限定拼接贡献'], main_rows)
    scans = table(['尾段', 'R', 'wL²', '可连接路径数', '非零分支和的移位', '最佳限定贡献'],
                  [row([r['id'], r['R'], r['tail_w_squared'], r['feasible_trails'],
                        ', '.join(map(str, r['noncancelling_rotations'])) or '—', signed_power(r['best_restricted_splice_score'])]) for r in data['reports']])
    chosen = next(r for r in data['reports'] if r['id'] == 'D3_chi2_to_3')
    group = data['groups'][chosen['connection_group']]
    best = group['rotations'][24]
    c = best['connections'][0]
    stages = [('第一轮 χ 输入 Δ₀', best['d0']), ('第一轮 χ 输出 B₁', best['b1']),
              ('第二轮 χ 输入 A₂ = L(B₁)', best['a2']), ('第二轮 χ 输出 B₂', c['b2']),
              ('中间 χ 输入差分 A₃ = L(B₂)', c['a3']), ('中间 χ 输出掩码 γ', group['gamma']),
              ('五轮最终输出掩码 Λ', chosen['output_mask'])]
    state_rows = table(['边界', '活跃列数', 'x[0] … x[4]'],
                       [row([name, len(s['active_columns']), '\n'.join(s['words'])]) for name, s in stages])
    empirical = []
    for e in experiments:
        assert e['bad_first'] == e['bad_middle'] == '0'
        rid, shift = e['id'].rsplit('_rot', 1)
        r = next(r for r in data['reports'] if r['id'] == rid)
        gr = data['groups'][r['connection_group']]
        s = gr['rotations'][int(shift)]
        for field, expected in [('input_words', s['d0']), ('first_output_words', s['b1']),
                                ('second_output_words', s['connections'][0]['b2']), ('join_words', gr['gamma'])]:
            assert [int(w, 16) for w in e[field].split(';')] == [int(w, 16) for w in expected['words']]
        empirical.append(row([e['id'], f'{int(e["conditional_samples"]):,}', e['second_event_hits'],
                              f'{float(e["p_second_estimate"]):.9f}',
                              f'[{float(e["ci95_low"]):.9f}, {float(e["ci95_high"]):.9f}]',
                              f'{e["positive_hits"]} / {e["negative_hits"]}']))
    measurements = table(['候选', '首轮条件样本数', '第二轮命中数', '第二轮条件概率实测', '逐项 95% 区间', '中间符号 + / −'], empirical)
    cancellations = table(['入口', '共同移位', '可连接路径数', '正 / 负路径数', '限定前缀贡献之和'],
                          [row([key, s['rotation'], len(s['connections']),
                                f'{sum(c["middle_correlation"] == 1 for c in s["connections"])} / {sum(c["middle_correlation"] == -1 for c in s["connections"])}',
                                signed_power(s['restricted_signed_prefix_sum'])])
                           for key, gr in data['groups'].items() for s in gr['rotations'] if s['connections']])
    content = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston：6→6 两轮差分与固定线性尾段的拼接</title>
<style>body{{margin:0;background:#f4f6f8;color:#203440;font:16px/1.75 system-ui,"Noto Sans CJK SC",sans-serif}}main{{max-width:1300px;background:white;margin:30px auto;padding:36px;border-radius:12px}}h1{{font-size:29px;line-height:1.4}}h2{{font-size:22px;margin-top:30px}}.muted{{color:#566b78}}.note{{padding:14px 20px;background:#edf7f7;border-left:4px solid #247f83}}.limit{{padding:14px 20px;background:#fff3df;border-left:4px solid #b47b20}}.flow{{text-align:center;padding:20px;background:#f1f5f7;font-size:18px}}.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums;margin:18px 0}}td{{border-bottom:1px solid #dce4e7;text-align:center;vertical-align:middle;padding:11px 12px;white-space:pre}}thead{{background:#edf3f5;font-weight:650;border-top:2px solid #284957}}tbody tr:last-child{{border-bottom:2px solid #284957}}a{{color:#096d86}}details{{margin-top:24px}}summary{{cursor:pointer;font-weight:650}}code{{background:#eef3f5;padding:2px 5px}}@media(max-width:700px){{main{{margin:0;padding:20px}}h1{{font-size:24px}}}}</style>
<main><h1>Gaston：6→6 两轮差分接上固定线性尾段</h1>
<p class="muted">2 轮差分 + 1 轮中间连接 + 原有线性尾段 · 首个 L 仍然省略</p>
<p class="note">已找到 5 轮与 6 轮可拼接候选。两轮差分各 6 个活跃 χ，差分路径概率为 2⁻²⁴，中间相关度为 −1。
原先的两轮、三轮线性尾段保持不变。下面的数值是限定路径贡献，不是完整置换相关度。</p>
<div class="flow">χ₁：6 活跃，p₁=2⁻¹² → L → χ₂：6 活跃，p₂=2⁻¹² → L → χ₃：g=−1 → 固定线性尾段</div>
{selected}
<p>五轮评分：2⁻¹² × 2⁻¹² × (−1) × 2⁻³⁸ = −2⁻⁶²。六轮评分：2⁻¹² × 2⁻¹² × (−1) × 2⁻⁶⁰ = −2⁻⁸⁴。
原先 1+1+2 的总轮数为 4；增加一轮差分后，总轮数相应变为 5。</p>
<h2>搜索范围</h2>
<p>使用已有 <code>diff_in6.json</code> 的 L 输入/输出模式，验证两侧都为 6 个活跃列。
穷举其全部 64 个共同移位。首轮逐列最大 DDT 概率的前驱在本模式下唯一；
第二轮每个活跃 χ 有 4 种可达输出差分，穷举 4⁶=4096 个组合，再连接到已有的全部 11 个尾段候选。
四个不同入口共检查 1,048,576 个组合，MITM 匹配与独立直接枚举全部一致。</p>
<p>这没有穷尽所有 6→6 传播模式。对没有命中的尾段，只能表述为本次模式及移位范围内没有找到连接。</p>
{scans}
<h2>为什么两轮路径概率可以精确写成 2⁻²⁴</h2>
<p>第一轮 6 个局部 DDT 事件各有概率 1/4，因此 p₁=2⁻¹²。第二轮差分导数产生 12 个独立线性条件。
将这些条件通过 L 拉回后，它们在首轮非活跃列对应的变量上仍有秩 12。
这些列在给定首轮 DDT 事件后仍独立均匀，因此第二轮所选差分的条件概率精确为 2⁻¹²，
从而所选两轮差分路径概率精确为 2⁻²⁴。常数只改变仿射偏移，不改变秩。</p>
<p>对于本 χ，非平凡 DDT 概率最大为 1/4，所以在“每轮 6 个活跃 χ”的单路径 DDT 乘积模型中，
权重 24 已达到逐 S 盒代价给出的下界。这不等于所有多轮真实概率或整体区分器的最优性。</p>
<h2>符号抵消检查</h2>
<p>首选入口为 <code>join_2</code>（D.3 第 2 个 χ 的入口，供 2–3 与 2–4 尾段共用）。
共同移位 7 位时有 4 条连接、两正两负，等权贡献抵消；移位 24 位时只有 1 条连接，贡献为负且未被同组其他可连接分支抵消。
<code>join_3</code> 为 D.3 第 3 个 χ 的入口，其 32 条连接按各输入差分分组后均抵消。</p>
{cancellations}
<p class="limit"><strong>结论范围：</strong>上述带符号求和只包含固定首轮输出差分事件与中间 DLCT 为 ±1 的分支。
它没有计算完整三轮前缀相关度；DLCT 为 0 的分支在条件分布下不能未经核验就全部丢弃。
完整 R 轮还包含其他差分事件、其他线性路径与交叉项，所以 −2⁻⁶²、−2⁻⁸⁴ 是限定拼接贡献，不是实测相关度或上下界。</p>
<h2>条件采样核对</h2>
<p>先在首轮指定 DDT 事件内均匀取样，再观察第二轮差分是否命中。
每组使用 2²² 个条件样本，理论命中数为 1024，第二轮理论条件概率为 1/4096≈0.000244141。
所有命中样本再经过第三轮检查，奇偶符号均与 g=−1 一致。
这是差分路径与中间连接验证，不是完整五轮、六轮相关度测量。</p>
{measurements}
<details><summary>查看移位 24 的完整路径与五轮端点</summary>
<p>状态按 x[0]…x[4] 列出，bit 0 为 LSB。前两个 χ 输入活跃列均为 6；中间 χ 输入差分有 51 个活跃列，连接掩码有 11 个活跃列。
中间差分较稠密仍能满足零权重连接，说明这里约束的是掩码下的奇偶关系。</p>{state_rows}</details>
<h2>数据与复现</h2><p><a href="search.json">完整搜索、条件秩与逐轮路径</a> · <a href="summary.csv">汇总 CSV</a> ·
<a href="conditional_experiments.csv">条件采样原始计数</a> · <a href="../fixed_paper_tail/table13_fixed_tail.html">原先一轮差分的结果</a></p>
</main></html>'''
    path = root / 'table13_six_six.html'
    path.write_text(content, encoding='utf-8')
    print(path)


if __name__ == '__main__':
    main()

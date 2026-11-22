#!/usr/bin/env python3
"""Build a self-contained Chinese result page for all six-bit paired PQ states."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

from pq_kernel_six_search import decode, encode, rotate


def h(value):
    return html.escape(str(value))


def power(value):
    if value is None:
        return '未找到非零和'
    if not value['numerator']:
        return '0'
    return ('−' if value['numerator'] < 0 else '+') + f'2^{value["log2_abs"]:g}'


def table(headers, rows):
    def row(cells, kind):
        return '<tr>' + ''.join(f'<{kind}>{h(c)}</{kind}>' for c in cells) + '</tr>'
    return ('<div class="scroll"><table><thead>' + row(headers, 'th') + '</thead><tbody>'
            + ''.join(row(r, 'td') for r in rows) + '</tbody></table></div>')


def csv_write(path, records):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def pairs_string(pairs):
    return '; '.join(f'{p["column"]}:({p["rows"][0]},{p["rows"][1]})' for p in pairs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    root = args.results.resolve()
    data = json.loads((root / 'search.json').read_text())
    pattern_data = json.loads((root / 'patterns.json').read_text())
    patterns = pattern_data['orbits']
    selected = [r for r in data['reports'] if r['id'] in data['best_by_round']]
    selected.sort(key=lambda r: r['R'])
    rows, summary = [], []
    for r in selected:
        b = r['best']
        old = '−2^-62 → −2^-61，绝对值 ×2' if r['R'] == 5 else ('−2^-84 → −2^-83，绝对值 ×2' if r['R'] == 6 else '原先未找到；新增连接')
        rows.append([r['R'], f'2 + 1 + {r["r_l"]}', r['id'], f'{b["orbit"]} / {b["rotation"]}',
                     f'{b["positive"]} / {b["negative"]}', power(b['restricted_signed_prefix_sum']),
                     r['tail_w_squared'], power(b['restricted_splice_score']), old])
        summary.append({'R': r['R'], 'rd': 2, 'rm': 1, 'rl': r['r_l'], 'tail': r['id'],
                        'orbit': b['orbit'], 'rotation': b['rotation'], 'positive': b['positive'], 'negative': b['negative'],
                        'p_single_trail': '2^-24', 'log2_abs_restricted_prefix_sum': b['restricted_signed_prefix_sum']['log2_abs'],
                        'tail_squared_weight': r['tail_w_squared'],
                        'score_sign': 1 if b['restricted_splice_score']['numerator'] > 0 else -1,
                        'log2_abs_restricted_score': b['restricted_splice_score']['log2_abs'],
                        'input_words': ';'.join(b['d0']['words']), 'output_mask_words': ';'.join(r['output_mask']['words'])})
    csv_write(root / 'summary.csv', summary)
    selected_table = table(['总轮数', 'rd + rm + rl', '固定线性尾段', '类别 / 移位', '正 / 负路径',
                            '限定前缀贡献', 'wL²', '限定拼接贡献', '相对旧模式'], rows)
    orbit_rows, orbit_csv, rotated_csv = [], [], []
    for p in patterns:
        counts = []
        for gid, gr in data['groups'].items():
            scans = [s for s in gr['scans'] if s['orbit'] == p['id']]
            nonzero = [s for s in scans if s['positive'] != s['negative']]
            if nonzero:
                best = max(nonzero, key=lambda s: abs(s['positive'] - s['negative']))
                counts.append(f'{gid}: {power(best["restricted_signed_prefix_sum"])} @ {best["rotation"]}')
        kind = 'vortex ' + str(p['vortex_rows']) if p['vortex_rows'] else '额外的成对 PQ 模式'
        orbit_rows.append([p['id'] + ('（旧模式）' if p['is_previous_orbit'] else ''), kind,
                           pairs_string(p['P_pairs']), pairs_string(p['Q_pairs']), '\n'.join(counts) or '无非零限定和'])
        orbit_csv.append({'id': p['id'], 'kind': kind, 'previous': p['is_previous_orbit'],
                          'P_pairs': pairs_string(p['P_pairs']), 'Q_pairs': pairs_string(p['Q_pairs']),
                          'theta_input_words': ';'.join(p['theta_input']['words']),
                          'first_output_words': ';'.join(p['first_output']['words']),
                          'second_input_words': ';'.join(p['second_input']['words'])})
        for shift in range(64):
            rotated_csv.append({'orbit': p['id'], 'rotation': shift,
                                'theta_input_words': ';'.join(encode(rotate(decode(p['theta_input']), shift))['words']),
                                'first_output_words': ';'.join(encode(rotate(decode(p['first_output']), shift))['words']),
                                'second_input_words': ';'.join(encode(rotate(decode(p['second_input']), shift))['words'])})
    csv_write(root / 'orbits.csv', orbit_csv)
    csv_write(root / 'all_704_states.csv', rotated_csv)
    orbit_table = table(['类别', '结构 / 行号', 'P 列：行对', 'Q 列：行对', '各入口最佳限定前缀贡献 @ 移位'], orbit_rows)
    connection_rows, connection_csv = [], []
    for gid, gr in data['groups'].items():
        for s in gr['scans']:
            if not s['connections']:
                continue
            connection_rows.append([gid, s['orbit'], s['rotation'], s['positive'], s['negative'],
                                    power(s['restricted_signed_prefix_sum'])])
            for i, c in enumerate(s['connections']):
                connection_csv.append({'group': gid, 'orbit': s['orbit'], 'rotation': s['rotation'], 'index': i,
                                       'sign': c['middle_correlation'], 'p_single_trail': '2^-24',
                                       'b2_words': ';'.join(c['b2']['words']), 'a3_words': ';'.join(c['a3']['words'])})
    csv_write(root / 'connections.csv', connection_csv)
    connections = table(['线性入口', '类别', '移位', '正路径', '负路径', '同一输入的限定贡献和'], connection_rows)
    tails = table(['尾段', '总轮数', '入口', 'wL²', '路径数（跨输入）', '非零和的输入数', '最佳类别 / 移位', '最佳限定贡献'],
                  [[r['id'], r['R'], r['connection_group'], r['tail_w_squared'], r['feasible_trails'],
                    r['noncancelling_input_count'], f'{r["best"]["orbit"]} / {r["best"]["rotation"]}' if r['best'] else '—',
                    power(r['best']['restricted_splice_score']) if r['best'] else '未找到非零和'] for r in data['reports']])
    endpoint_sections = []
    for r in selected:
        b, gr = r['best'], data['groups'][r['connection_group']]
        scan = next(s for s in gr['scans'] if (s['orbit'], s['rotation']) == (b['orbit'], b['rotation']))
        stages = [('χ1 输入 Δ0', b['d0']), ('χ1 输出 B1', b['b1']), ('χ2 输入 A2=L(B1)', b['a2'])]
        for i, c in enumerate(scan['connections']):
            stages.extend([(f'χ2 输出 B2[{i}], g={c["middle_correlation"]}', c['b2']),
                           (f'χ3 输入 A3[{i}]=L(B2[{i}])', c['a3'])])
        stages.extend([('中间 χ 输出掩码 γ', gr['gamma']), ('最终输出掩码 Λ', r['output_mask'])])
        endpoint_sections.append(f'<details><summary>{r["R"]} 轮：{h(r["id"])}，{h(b["orbit"])} / {b["rotation"]}，全部边界</summary>'
                                 + table(['边界', '活跃列数', '五行十六进制（x[0] 至 x[4]）'],
                                         [[name, len(s['active_columns']), '\n'.join(s['words'])] for name, s in stages]) + '</details>')
    with (root / 'conditional_experiments.csv').open() as f:
        experiments = list(csv.DictReader(f))
    job_map = json.loads((root / 'conditional_job_map.json').read_text())
    empirical_rows = []
    for e in experiments:
        assert e['bad_first'] == e['bad_middle'] == '0'
        job = next(j for j in job_map if j['id'] == e['id'])
        report = next(r for r in data['reports'] if r['id'] == job['report'])
        gr = data['groups'][report['connection_group']]
        scan = next(s for s in gr['scans'] if (s['orbit'], s['rotation']) == (job['orbit'], job['rotation']))
        c = scan['connections'][job['connection_index']]
        for key, state in [('input_words', report['best']['d0']), ('first_output_words', report['best']['b1']),
                           ('second_output_words', c['b2']), ('join_words', gr['gamma'])]:
            assert tuple(int(w, 16) for w in e[key].split(';')) == decode(state)
        empirical_rows.append([e['id'], f'{int(e["conditional_samples"]):,}', e['second_event_hits'],
                               f'{float(e["p_second_estimate"]):.9f}',
                               f'[{float(e["ci95_low"]):.9f}, {float(e["ci95_high"]):.9f}]',
                               f'{e["positive_hits"]} / {e["negative_hits"]}'])
    empirical = table(['采样候选', '首轮条件样本数', '第二轮命中', '条件概率实测', '逐项 95% 区间', '中间 + / −'], empirical_rows)
    content = f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston：全部 11 类六比特 P=Q=0 模式及差分线性拼接</title>
<style>
body{{margin:0;background:#f3f6f8;color:#243b49;font:16px/1.8 system-ui,"Noto Sans CJK SC",sans-serif}}
main{{max-width:1400px;margin:28px auto;background:white;padding:36px;border-radius:12px}}
h1{{font-size:29px;line-height:1.45}}h2{{font-size:22px;margin-top:32px}}p{{max-width:1100px}}
.note{{background:#eaf5f3;border-left:4px solid #268071;padding:14px 20px}}
.scope{{background:#fff5e4;border-left:4px solid #b58437;padding:14px 20px}}
.flow{{background:#f0f4f8;padding:16px;text-align:center}}.scroll{{overflow-x:auto}}
table{{width:100%;border-collapse:collapse;margin:18px 0;font-size:14px;font-variant-numeric:tabular-nums}}
th,td{{padding:10px 12px;border-bottom:1px solid #dce5ea;white-space:pre;text-align:center;vertical-align:middle}}
th{{background:#edf3f6;border-top:2px solid #395364}}tbody tr:last-child td{{border-bottom:2px solid #395364}}
tr:hover td{{background:#f6f9fb}}a{{color:#086d84}}code{{background:#edf2f5;padding:2px 5px}}
details{{margin:18px 0}}summary{{cursor:pointer;font-weight:650}}.muted{{color:#617783}}
@media(max-width:700px){{main{{margin:0;padding:18px}}h1{{font-size:24px}}}}
</style></head><body><main>
<h1>Gaston：全部六比特成对 P=Q=0 模式</h1>
<p class="muted">差分段 2 轮 + 中间连接 1 轮 + 固定论文线性尾段 · 仅省略整个区间首个 L</p>
<p class="note"><strong>完整枚举得到 11 类、704 个不同移位状态。</strong>
旧模式是 K09 的一个代表，共 64 个移位；本次补齐另外 10 类、640 个状态。
五轮与六轮的最佳限定拼接贡献分别改进为 −2⁻⁶¹、−2⁻⁸³；另找到七至九轮的连接。</p>
<div class="flow">χ₁：6 活跃，p₁=2⁻¹² → L → χ₂：6 活跃，p₂=2⁻¹² → L → χ₃：g=±1 → 固定线性尾段</div>
<h2>Table 13 形式的结果表</h2>
{selected_table}
<p>每条两轮差分路径的精确概率均为 2⁻²⁴。K06 / 22 有两条负贡献，和为 −2⁻²³；
K03 / 52 有四条正贡献，和为 +2⁻²²。表中的改善来自同一输入下的带符号叠加。</p>
<p class="scope"><strong>数值含义：</strong>只聚合固定首轮输出 B1 事件中，中间规范化 DLCT 为 ±1 的分支，
随后乘以一条固定线性路径的相关度平方。其他 B1 分支、条件分布下 DLCT 为 0 的分支，以及其他线性路径和交叉项尚未聚合。
这些数值是限定贡献，不是完整三轮前缀或完整 R 轮相关度，也不是它们的上下界。</p>
<h2>枚举范围与完整性</h2>
<p>状态 A 位于 ρeast 之后、θ 之前。A 恰有 6 个非零比特，分布于三列、每列两个，故 P=0。
逐行左移 t=(25,32,52,60,63) 后同样是三列、每列两个，故 Q=0，继而 E=0、θ(A)=A。
这里专门研究用户指定的这一子集，不把它等同于定义 5 的全部 E=0 情形。</p>
<p>把一个比特 (i,p) 看成 P 列 p 与 Q 列 p+tᵢ 之间的一条边，得到 3+3 个顶点、每个顶点度为 2 的简单二分图。
tᵢ 互异，因此没有重边；六条边只能组成一个六边形。固定一个 P 列为 0，枚举前五条边的行号、由闭合条件确定第六条边，
即可穷尽所有模式。独立的“三个 P 列对，其 Q 签名异或为零”枚举得到相同结果：66 个有向六边形表示、33 个 P=0 列锚定表示，去重后均为 11 类。</p>
<p>只按五行共同循环移位去重，不交换行号。规范代表是五个无符号 64 位行字按字典序最小的移位。
每类恰有 64 个不同状态。转换回完整 L 的边界：B1[i]=ROR(A[i],east[i])，A2[i]=ROL(A[i],west[i])。
全部 11 类都验证了 L(B1)=A2，且 B1、A2 各有六个活跃列，所以对应前两轮各六个活跃 S 盒。</p>
<h2>11 类的具体 P、Q 列对</h2>
<p>行号、列号均从 0 开始，最低位为 bit 0；“28:(1,3)”表示第 28 列的第 1、3 行各有一个比特。
其中 10 类由三行 vortex 公式产生：Aᵢ=Xᵗʲ+Xᵗᵏ，Aⱼ=Xᵗᵏ+Xᵗⁱ，Aₖ=Xᵗⁱ+Xᵗʲ；K05 是额外的一类，行比特数为 2、1、3。</p>
{orbit_table}
<p>K05 的规范代表为 A1=X⁰+X²⁸，A2=X³⁶，A3=X⁰+X²⁸+X³⁶，A0=A4=0。
它满足同样的成对 P=Q=0 条件，但本次固定入口池中未找到连接。
旧文件 <code>diff_in6.json</code> 对应 K09 / 22；旧文件再移位 24 的候选即本报告 K09 / 46。</p>
<h2>固定线性尾段与全部匹配结果</h2>
<p>11 个既有尾段只有四个不同入口 γ：join_0 为 D.2；join_1 为 D.3 从第 1 个 χ 起；
join_2 为 D.3 从第 2 个 χ 起；join_3 为 D.3 从第 3 个 χ 起。实际使用 γ=Lᵀα，因此没有额外省略线性尾段前的 L。</p>
<p>每个状态的第二轮输出共有 4⁶=4096 种。全部范围等价于检查 11×64×4×4096=11,534,336 个组合，
通过 MITM 精确匹配完成；四个入口分别有 4、0、8、113 条连接。不同输入间的路径计数不能直接相加为相关度。</p>
{tails}
<details><summary>按同一输入分组的全部非空连接，包括符号抵消</summary>{connections}</details>
<h2>概率与独立验证</h2>
<p>首轮每个单比特输出的最大概率前驱唯一，且 Δ0=B1，故 p₁=(1/4)⁶=2⁻¹²。
第二轮导数有 12 个独立条件，拉回后在不受首轮事件限制的输出列上仍为秩 12。
这在 11 个代表上逐一核对，并由共同移位保持，因此任意可行第二轮输出的条件概率精确为 2⁻¹²，路径概率为 2⁻²⁴。
其他首轮前驱具有同一活跃列支撑，条件秩证明仍适用，而 p₁ 不会更大，因此不能提高固定 B1 下的这项限定贡献。</p>
<p>搜索对所有连接逐列穷举 32 个输入，确认中间导数恒为 ±1；核对了每条固定尾段的 LAT、转置关系与轮常数。
回归程序另用正向 L 核对 11 类的条件秩，以直接 4⁶ 枚举核对全部非空输入组及每类移位 0，并恢复旧模式所有移位的结果。</p>
<p>下列四组各取 2²² 个首轮条件样本。第二轮理论命中数为 1024，理论条件概率为 1/4096。
所有命中样本的中间符号都符合精确结果；该实验只验证路径和连接，未测量完整减轮相关度。</p>
{empirical}
<h2>候选完整边界</h2>{''.join(endpoint_sections)}
<h2>数据与复现</h2>
<p><a href="patterns.json">11 类代表及完备性信息</a> · <a href="all_704_states.csv">全部 704 个状态</a> ·
<a href="orbits.csv">11 类简表</a> · <a href="summary.csv">结果表</a> · <a href="connections.csv">全部 125 条连接</a> ·
<a href="search.json">逐尾段完整结果</a> · <a href="conditional_experiments.csv">条件采样原始计数</a> ·
<a href="reproducibility.json">来源与文件校验和</a></p>
<p>复现命令与定义见项目 README。搜索仅依赖 Python 标准库，采样复用现有 C++17 程序。</p>
</main></body></html>'''
    (root / 'table13_pq_kernel.html').write_text(content, encoding='utf-8')
    project = Path(__file__).resolve().parents[1]
    sources = [project / 'tools/pq_kernel_six_search.py', Path(__file__), project / 'tests/test_pq_kernel_six.py',
               project / 'tools/six_six_prefix_search.py', project / 'src/conditional_trail_verify.cpp',
               project.parents[1] / 'gaston_ref/gaston.pdf', project.parents[1] / 'gaston_ref/gaston_slides.pdf']
    files = sorted(p for p in root.iterdir() if p.is_file() and p.name != 'reproducibility.json')
    manifest = {'sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files + sources},
                'sample_count_per_job': 4194304, 'seed': 2026101541, 'threads': 4,
                'scope': data['scope'], 'limits': data['limits']}
    (root / 'reproducibility.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(root / 'table13_pq_kernel.html')


if __name__ == '__main__':
    main()

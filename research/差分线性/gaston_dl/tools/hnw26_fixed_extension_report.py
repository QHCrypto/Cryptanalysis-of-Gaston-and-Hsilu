#!/usr/bin/env python3
"""Audit cross-mode equality and render the five/six-round HNW26 report."""
import argparse
import gzip
import hashlib
import html
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from hnw26_fixed_extension import ROOT, frac, save


def num(value):
    q = frac(value) if isinstance(value, dict) else Fraction(value)
    if not q:
        return '0'
    if abs(q.numerator) == 1 and q.denominator & (q.denominator - 1) == 0:
        exponent = q.denominator.bit_length() - 1
        if not exponent:
            return '+1' if q > 0 else '−1'
        return ('+' if q > 0 else '−') + f'2<sup>−{exponent}</sup>'
    return html.escape(str(q))


def state(title, value):
    words = value['words']
    return '<div><h3>' + html.escape(title) + '</h3><pre>' + '\n'.join(
        f'x{i}: {int(w, 16):016x}' for i, w in enumerate(words)) + '</pre></div>'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=ROOT / 'results/fixed_tail_extension/hnw26')
    args = p.parse_args()
    root = args.source.resolve()
    result = json.loads((root / 'results.json').read_text())
    config = json.loads((root / 'configuration.json').read_text())
    validation = json.loads((root / 'validation.json').read_text())
    shortcut = json.loads((root / 'shortcut_validation.json').read_text())
    assert result['status'] == 'completed' and validation['status'] == shortcut['status'] == 'PASS'
    for source in config['sources']:
        assert hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() == source['sha256']
    branch_checks = 0
    for row in result['rows']:
        full, zero = row['HNW']['full'], row['HNW']['zero']
        for a, b in zip(full['shards'], zero['shards']):
            assert (a['branch_begin'], a['branch_end']) == (b['branch_begin'], b['branch_end'])
            with gzip.open(a['artifact'], 'rt') as f:
                aa = json.load(f)['branches']
            with gzip.open(b['artifact'], 'rt') as f:
                bb = json.load(f)['branches']
            assert aa == bb, 'Full and zero models differ branch by branch'
            branch_checks += len(aa)
        assert frac(full['prefix_correlation']) == frac(row['exact_prefix']['correlation'])
        assert full['nonzero_prefix_branch_correlations'] == 1
        assert full['nonzero_branch_correlations'] == 0
        assert sum(full['absorbed_histogram'].values()) == 65536
    audit = {'status': 'PASS', 'full_vs_zero_branch_equality_checks': branch_checks,
             'all_source_hashes_unchanged': True, 'prefix_matches_exact_for_both_intervals': True,
             'shortcut_regression': shortcut}
    save(root / 'cross_validation.json', audit)
    table_rows = []
    endpoint_sections = []
    raw_links = []
    c = config['rows'][0]['candidate']
    for row, selected in zip(result['rows'], config['rows']):
        assert row['row_id'] == selected['id']
        full, zero = row['HNW']['full'], row['HNW']['zero']
        weight = selected['candidate']['single_characteristic_weight']
        tail_weight = selected['tail_squared_weight']
        splice = frac(full['prefix_correlation']) * Fraction(1, 2**tail_weight)
        assert splice == Fraction(1, 2**weight)
        table_rows.append(f'<tr><td><strong>{row["rounds"]} 轮</strong></td><td>1 + 1 + {selected["r_l"]}</td>'
                          f'<td>{html.escape(row["tail_id"])}</td><td>{num(splice)}</td>'
                          f'<td>{num(row["DATF"])}</td><td>{num(full["correlation"])}</td>'
                          f'<td>{num(zero["correlation"])}</td></tr>')
        br = full['nonzero_prefix_branches'][0]
        assert frac(br['prefix_conditional_model_correlation']) == 1
        assert frac(br['probability']) == Fraction(1, 65536)
        endpoint_sections.append(f'<section id="r{row["rounds"]}"><h2>{row["rounds"]} 轮端点 · {html.escape(row["tail_id"])}</h2>'
            f'<p>原路径分段：1 轮差分 + 1 轮连接 + {selected["r_l"]} 轮固定线性尾段。'
            f'完整 HNW 计算分段：r<sub>d</sub> = 1、r<sub>m</sub> = {row["rounds"] - 1}、r<sub>l</sub> = 0。'
            f'参考轮区间为 {row["start_round_index"] + 1}–12（此处轮号从 1 开始）。</p>'
            f'<div class="states">{state("输入差分 Δ", selected["candidate"]["input_difference"])}'
            f'{state("最终输出掩码 Λ", selected["output_mask"])}</div>'
            f'<p>线性尾段单路径 c<sub>L</sub>² = 2<sup>−{tail_weight}</sup>。'
            f'若用两轮前缀 {num(full["prefix_correlation"])} 与该固定单条线性路径继续拼接，得到 {num(splice)}；'
            '这个值仍依赖尾段拼接近似，未成为完整置换相关度的证明。</p>'
            f'<p>首轮编号为 {br["branch_index"]} 的分支是前缀读出唯一非零项：概率 {num(br["probability"])}, '
            f'条件模型相关度 {num(br["prefix_conditional_model_correlation"])}。'
            f'该分支在继续向后传播时，到第 {br["absorbed_by_round"]} 轮也完全退化为模型零。</p>'
            f'<p class="muted">端点 ID：<code>{row["endpoint_id"]}</code></p></section>')
        for mode in ('full', 'zero'):
            for shard in row['HNW'][mode]['shards']:
                rel = Path(shard['artifact']).relative_to(root).as_posix()
                raw_links.append(f'<li><a href="{rel}">R{row["rounds"]} {mode} · 分支 '
                                 f'{shard["branch_begin"]}–{shard["branch_end"] - 1}</a></li>')
    hist = result['rows'][0]['HNW']['full']['absorbed_histogram']
    assert all(r['HNW'][m]['absorbed_histogram'] == hist for r in result['rows'] for m in ('full', 'zero'))
    h3, h4 = hist.get('3', 0), hist.get('4', 0)
    p3 = h3 / 65536 * 100
    prefix_rows = ''.join(f'<tr><td>{r["rounds"]} 轮候选的前两轮</td><td>0</td>'
                         f'<td>{num(r["HNW"]["full"]["prefix_correlation"])}</td>'
                         f'<td>{num(r["HNW"]["zero"]["prefix_correlation"])}</td>'
                         f'<td>{num(r["exact_prefix"]["correlation"])}</td></tr>' for r in result['rows'])
    page = f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston 五、六轮固定尾段 · HNW26 相关度评估</title><style>
:root{{--ink:#16324b;--muted:#566b7c;--line:#d9e3e9;--blue:#17668b;--paper:#fff;--bg:#f2f5f7}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.85 system-ui,-apple-system,sans-serif}}
main{{max-width:1220px;margin:auto;padding:40px 28px 80px}}h1{{font-size:34px;line-height:1.35;letter-spacing:-.03em;margin:12px 0 18px}}h2{{font-size:23px;line-height:1.5;margin:0 0 16px}}h3{{font-size:16px}}
section{{background:var(--paper);border:1px solid var(--line);border-radius:12px;padding:26px;margin-top:24px}}
a{{color:var(--blue)}}.kicker{{font-size:13px;letter-spacing:.12em;color:var(--muted)}}.lead{{max-width:960px;font-size:18px}}
.cards{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px;margin-top:26px}}.card{{background:white;border:1px solid var(--line);border-top:4px solid var(--blue);border-radius:8px;padding:20px}}.card span{{display:block;color:var(--muted);font-size:14px}}.card strong{{display:block;font-size:36px;line-height:1.6}}.card small{{font-size:14px}}
.callout{{background:#fff5e2;border-left:4px solid #b77a20;padding:15px 20px;margin-top:20px}}.note{{background:#eaf3f8;border-left:4px solid var(--blue);padding:15px 20px}}
.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:13px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}th{{background:#edf3f6;white-space:nowrap}}sup{{line-height:0}}.states{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}}pre{{background:#f0f4f7;padding:16px;font:14px/1.9 ui-monospace,monospace;overflow:auto;border-radius:6px}}code{{overflow-wrap:anywhere;font-size:13px}}.muted{{color:var(--muted);font-size:14px}}.formula{{overflow:auto;background:#f3f6f8;padding:16px;text-align:center;font-family:serif;font-size:19px}}.bar{{height:30px;display:flex;border-radius:6px;overflow:hidden;background:#b17924}}.bar .r3{{width:{p3:.8f}%;background:#377d9d}}.legend{{display:flex;gap:22px;flex-wrap:wrap}}.legend i{{display:inline-block;width:12px;height:12px;margin-right:8px;background:#377d9d}}.legend .r4{{background:#b17924}}details{{margin-top:18px}}summary{{cursor:pointer;color:var(--blue)}}.files{{columns:3;font-size:13px;padding-left:20px}}.files li{{break-inside:avoid;margin:5px 0}}nav{{display:flex;gap:18px;flex-wrap:wrap;font-size:14px}}
@media(max-width:760px){{main{{padding:24px 14px 60px}}h1{{font-size:28px}}section{{padding:18px}}.cards,.states{{grid-template-columns:1fr}}.card strong{{font-size:30px}}.files{{columns:1}}}}
</style></head><body><main>
<div class="kicker">GASTON · FIXED LINEAR TAIL · HNW26</div>
<h1>五轮、六轮候选的 HNW26 相关度评估</h1>
<p class="lead">对同一八比特差分前缀接 D3_chi1_to_3 与 D3_chi1_to_4 的两个端点，完整 HNW26 估计均为 0。两轮前缀的 HNW26 结果为 +2<sup>−16</sup>，与此前的严格计算一致。</p>
<nav><a href="#comparison">完整结果</a><a href="#prefix">前缀对照</a><a href="#zero">零值来源</a><a href="#method">方法与范围</a><a href="#files">复核文件</a><a href="../fixed_tail_extension.html">原搜索报告</a></nav>
<div class="cards"><div class="card"><span>完整五轮 · HNW26</span><strong>0</strong><small>联合版与零状态值掩码版一致</small></div><div class="card"><span>完整六轮 · HNW26</span><strong>0</strong><small>每种设置覆盖 65,536 个分支</small></div><div class="card"><span>完整两轮前缀</span><strong>+2<sup>−16</sup></strong><small>HNW26 与严格前缀计算一致</small></div></div>
<div class="callout"><strong>结论范围：</strong>五、六轮的 0 是 HNW26 局部列近似内的零；不能据此证明真实完整相关度等于 0。固定尾段的 +2<sup>−76</sup>、+2<sup>−98</sup> 也仍然是拼接估计。当前严格确认的是两轮前缀的 +2<sup>−16</sup>。</div>
<section id="comparison"><h2>完整端点与固定尾段拼接值</h2><div class="scroll"><table><thead><tr><th>总轮数</th><th>原路径分段</th><th>固定尾段</th><th>固定尾段拼接值<br>依赖尾段近似</th><th>普通 DATF<br>完整端点</th><th>HNW26 联合版<br>完整端点</th><th>HNW26 零值版<br>完整端点</th></tr></thead><tbody>{''.join(table_rows)}</tbody></table></div>
<p>完整 HNW 计算固定输入差分 Δ 与最终输出掩码 Λ，精确枚举首个 χ 的所有 DDT 分支，然后对剩余全部轮次做局部相关向量传播。原来的线性尾段不再作为一个外乘系数：此处没有再次乘入 p<sub>D</sub> 或 c<sub>L</sub>²。</p><p>固定尾段拼接值则使用两轮前缀相关度 × 单条线性尾段相关度平方。两种评估保留的信息与采用的近似不同，数值不同并不矛盾。</p></section>
<section id="prefix"><h2>前两轮：HNW26 恢复了普通 DATF 丢失的非零值</h2><div class="scroll"><table><thead><tr><th>区间</th><th>普通 DATF</th><th>HNW26 联合版</th><th>HNW26 零值版</th><th>严格前缀结果</th></tr></thead><tbody>{prefix_rows}</tbody></table></div>
<p>八个活跃输入列各有 4 个等概率首轮 DDT 输出，因此共有 4<sup>8</sup> = 65,536 个分支，每个分支概率 2<sup>−16</sup>。对于连接处掩码 γ，仅首轮输出差分 b = Δ 的分支给出条件相关度 +1，其余分支在 HNW 的前缀读出下均为 0，最终得到 +2<sup>−16</sup>。</p>
<div class="states">{state('公共输入差分 Δ', c['input_difference'])}{state('两轮前缀输出掩码 γ', c['join_mask'])}</div>
<p>坐标约定：5 行 × 64 位，行号与位号从 0 开始，最低有效位为位 0。γ 是固定线性尾段首个输入掩码经 L<sup>T</sup> 回传后的掩码。前缀评估终止于第二个 χ 的输出，其后的 L 已计入 γ。</p>
<p>相关度 ρ = Pr[输出差分奇偶为 0] − Pr[输出差分奇偶为 1]；偏差 ε = ρ/2。因此前缀的偏差为 +2<sup>−17</sup>。完整五、六轮的两种 HNW 模型偏差均为 0。</p></section>
<section id="zero"><h2>完整模型为什么变成 0</h2><p>所有四次计算（两个端点 × 两个版本）得到相同的局部相关向量退化统计：</p><div class="bar" role="img" aria-label="65378 个分支到第三轮完全退化，158 个分支到第四轮完全退化"><div class="r3"></div></div><div class="legend"><p><i></i>第 3 轮：{h3:,} 个分支（{p3:.4f}%）</p><p><i class="r4"></i>第 4 轮：{h4:,} 个分支（{100-p3:.4f}%）</p></div>
<p>“完全退化”表示模型保存的每一列局部相关向量，只剩零掩码坐标为 1，所有非零掩码坐标都严格为 0；后续模型更新无法从这个状态恢复非零值。因此，到五轮或六轮的任何非零输出掩码读出均为模型零。</p>
<p>此次最终求和中没有非零条件贡献：65,536 个条件模型相关度全部为 0。这不是正负非零分支相消，也不是浮点精度、下溢或人为小值阈值造成的零。</p>
<div class="note">HNW26 在每列内部保留信息，并在跨列传播及多列读出时使用乘积近似。即使每一列的边缘分布都变为均匀，真实状态仍可能存在跨列联合关系。本结果定位了模型信息消失的位置，不能将“局部量为零”直接升级为“真实完整端点相关度为零”。</div></section>
{''.join(endpoint_sections)}
<section id="method"><h2>计算设置与数学范围</h2>
<p>保持此前的轮约定：仅省略完整区间开头的线性层 L；R 轮包含 R 个 χ 与 R−1 个 L。首轮常数可吸收到均匀首轮输入中，后续使用参考实现的轮常数。五轮从参考第 8 轮开始，六轮从第 7 轮开始，均到第 12 轮结束。</p>
<p>联合版保存每列 1,024 个联合状态／差分 Fourier 坐标；零状态值掩码版只保存 32 个差分坐标。联合版保留局部条件状态值、差分及轮常数相位，但没有保存全部跨列联合分布。</p>
<div class="formula">ρ̂<sub>HNW</sub>(Δ, Λ) = Σ<sub>b</sub> Pr[χ(X) ⊕ χ(X ⊕ Δ) = b] · ∏<sub>j</sub> q<sub>R,j</sub><sup>(b)</sup>(0, Λ<sub>j</sub>)</div>
<p>分支概率由 DDT 精确取得，局部 χ 更新采用精确 Fourier 变换或 LAT 平方更新；误差来源是跨列乘积近似，不是有理数运算。计算使用 GMP 任意精度整数与 2 的幂分母，没有浮点判零。</p>
<p>为加快运行，对全部首轮分支证明了首个 L 输出的 1,984 个非零局部状态值掩码均被严格零因子消去。该证书允许缓存第二轮局部初始化，不改变 HNW 模型。另用原更新过程对 40 个组合做了逐分支对照，覆盖两种模式、稠密前缀掩码、单比特输出及 10 个分支位置。</p>
<p>后续研究应把这两个候选保留为“前缀严格非零、完整五／六轮尚未确认”。若继续区分真实零与近似丢失信息，需要针对第 3–4 轮保留相关列组，或对固定尾段的贡献及可能抵消项做更细的带符号求和；单纯把本次 0 当作排除证据并不充分。</p></section>
<section id="files"><h2>完成情况与复核文件</h2><p>已完成 256 个批次、262,144 次分支评估。每个端点与每个模式的分支概率和均严格为 1；独立重建了每个分支的 DDT 概率与枚举序号，检查了加权和、偏差换算与完整覆盖。联合版和零值版的 131,072 对分支记录完全一致。计算及逐批验证用时约 {result['elapsed_seconds']:.1f} 秒。</p>
<p><a href="results.json">完整汇总 JSON</a> · <a href="configuration.json">端点、设置与源文件哈希</a> · <a href="validation.json">全分支核验</a> · <a href="cross_validation.json">两模式交叉核验</a> · <a href="shortcut_validation.json">优化前后回归检查</a> · <a href="backend_validation.json">算术后端自测</a> · <a href="../exact_prefix_audit.json">此前严格前缀计算</a></p>
<details><summary>查看 256 个原始分支文件（JSON gzip）</summary><ul class="files">{''.join(raw_links)}</ul></details>
<p class="muted">完整端点的 HNW26 结果均为近似模型结论。报告中的“两轮严格值”引用已有精确计算，本次以 HNW26 独立对照；未进行完整五／六轮穷举或新的置换采样。</p></section>
</main></body></html>'''
    target = root / 'hnw26_fixed_extension.html'
    target.write_text(page)
    save(root / 'report_provenance.json', {'script': str(Path(__file__).resolve()),
         'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'report_sha256': hashlib.sha256(target.read_bytes()).hexdigest(), 'status': 'completed'})
    print(target)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Render the fixed-paper-tail search, keeping prefix validation separate."""
import argparse
import csv
import html
import json
from pathlib import Path


def h(x):
    return html.escape(str(x))


def tr(values):
    return "<tr>" + "".join(f"<td>{h(x)}</td>" for x in values) + "</tr>"


def table(headers, rows):
    return '<div class="scroll"><table><thead>' + tr(headers) + '</thead><tbody>' + ''.join(rows) + '</tbody></table></div>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    root = args.results.resolve()
    data = json.loads((root / "search.json").read_text())
    patterns = json.loads((root / "pattern_prefix_search.json").read_text())
    with (root / "prefix_experiments.csv").open() as f:
        experiments = list(csv.DictReader(f))
    by_id = {r["id"]: r for r in data["reports"]}
    # Refuse a report that joins unrelated empirical endpoints to a candidate.
    for e in experiments:
        r = by_id[e["id"]]
        assert e["rounds"] == "2" and int(e["start"]) == r["start_round_index"]
        assert [int(x, 16) for x in e["mask_words"].split(";")] == [int(x, 16) for x in r["join_mask"]["words"]]
        if r["top_candidates"]:
            c = r["top_candidates"][0]
            assert [int(x, 16) for x in e["input_words"].split(";")] == [int(x, 16) for x in c["input_difference"]["words"]]
            assert int(e["exact_numerator"]) == c["prefix"]["numerator"]
    selected = table(["总轮数 R", "差分 rd", "中间 rm", "线性 rl", "固定尾段", "wL²", "精确前缀 ρG", "log₂|拼接评分|", "完整 R 轮实测"],
                     [tr([r["R"], 1, 1, r["r_l"], r["tail_id"], r["tail_w_squared"],
                          f'{r["prefix_numerator"]}/{r["prefix_denominator"]}', int(r["log2_abs_score"]), "未测量"])
                      for r in data["best_by_R"]])
    all_rows = []
    for r in data["reports"]:
        score = int(r["top_candidates"][0]["log2_abs_score"]) if r["top_candidates"] else "无非零前缀"
        all_rows.append(tr([r["id"], r["R"], r["tail_w_squared"], len(r["join_mask"]["active_columns"]),
                            r["nonzero_inputs"], score]))
    scans = table(["固定尾段", "R", "wL²", "连接掩码活跃列", "非零输入 / 1984", "最佳 log₂|评分|"], all_rows)
    measured = table(["尾段代表", "前缀轮数", "精确 ρG", "实测 ρG", "逐项 95% Wilson 区间", "样本对数"],
                     [tr([e["id"], 2, f'{int(e["exact_numerator"])/int(e["exact_denominator"]):.6f}',
                          f'{float(e["rho"]):.9f}', f'[{float(e["ci95_low"]):.9f}, {float(e["ci95_high"]):.9f}]',
                          f'{int(e["samples"]):,}']) for e in experiments])
    endpoints = []
    for r in data["best_by_R"]:
        body = [f'<h3>R = {r["R"]} · {h(r["tail_id"])}</h3>']
        for name, key in (("输入差分 Δ", "input_words"), ("前缀输出 / 尾段输入掩码 γ", "join_words"), ("最终输出掩码 Λ", "output_words")):
            body.append(f'<p>{name}</p><pre>{h(chr(10).join(r[key].split(";")))}</pre>')
        endpoints.append(''.join(body))
    links = ' · '.join(f'<a href="{name}">{label}</a>' for name, label in [
        ("summary.csv", "汇总 CSV"), ("search.json", "所有候选与逐轮掩码"),
        ("imported_trails.json", "论文导入与核验"), ("prefix_experiments.csv", "前缀采样计数"),
        ("pattern_prefix_search.json", "六列 / 七列试验")])
    content = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston：固定论文线性尾段，搜索差分前缀</title>
<style>
body{{margin:0;background:#f4f6f8;color:#192e3b;font:16px/1.75 system-ui,"Noto Sans CJK SC",sans-serif}}
main{{max-width:1250px;margin:32px auto;padding:36px;background:white;border-radius:12px}}
h1{{font-size:29px;line-height:1.4;margin:0 0 10px}}h2{{font-size:22px;margin-top:34px}}h3{{font-size:18px}}
.sub{{color:#536674}}.note{{padding:14px 20px;background:#edf6f8;border-left:4px solid #237d8a}}
.limits{{padding:14px 20px;background:#fff3df;border-left:4px solid #b07a20}}
.flow{{font:19px/1.8 ui-monospace,monospace;text-align:center;padding:18px;background:#f6f8fa}}
.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;font-size:14px;font-variant-numeric:tabular-nums;margin:18px 0}}
td{{text-align:center;padding:12px 9px;border-bottom:1px solid #dce3e8;white-space:nowrap}}thead{{font-weight:650;background:#f2f6f8;border-top:2px solid #254555}}tbody tr:last-child{{border-bottom:2px solid #254555}}
code,pre{{background:#f2f5f7;border-radius:6px}}code{{padding:2px 5px}}pre{{font-size:14px;padding:14px;overflow-x:auto}}
details{{margin-top:25px}}summary{{cursor:pointer;font-weight:650}}a{{color:#086f89}}footer{{color:#536674;font-size:13px;margin-top:30px}}
@media(max-width:700px){{main{{margin:0;padding:20px}}h1{{font-size:24px}}}}
</style><main>
<h1>Gaston：固定论文线性尾段，搜索前面的差分</h1>
<p class="sub">ePrint 2024/1867 · 附录 D.2 / D.3 · 首个线性层省略 · 本轮实际搜索结果</p>
<p class="note">已逐层核对论文路径，再固定尾段搜索前缀。下表保留当前候选池中每个总轮数的最高评分。
前缀相关度通过整数枚举精确计算；完整置换相关度尚未验证。</p>
<div class="flow">[χ]<sub>差分 1 轮</sub> → [L → ι → χ]<sub>中间 1 轮</sub> → [(L → ι → χ)<sup>rl</sup>]<sub>固定线性尾段</sub></div>
<p>R = rd + rm + rl；只省略整个区间的第一个 L。尾段最前面的 L 保留，连接掩码为 γ = Lᵀα，
α 是论文所列首个 χ 的输入掩码。后续常数按参考实现的末 R 轮保留。</p>
{selected}
<p class="limits"><strong>评分含义：</strong>S = ρG(Δ,γ) × c<sub>trail</sub>²。
ρG 是两轮前缀的精确差分线性相关度；c<sub>trail</sub> 是固定尾段各 χ 的 LAT 相关度乘积。
该值保留同一尾路径配对的一项贡献，不是完整 R 轮的实测相关度，也不是其上下界。
其他尾路径、线性包及交叉项可能改变或抵消它。</p>
<h2>为什么 4 轮评分是 2⁻³⁹</h2>
<p>取附录 D.3 的第 2–3 个 χ，线性平方相关度权重为 18+20=38。
输入差分只置第 4 行第 1 位（行、位均从 0 编号，位 0 是 LSB）。首个 χ 的两个输出差分分支各有概率 1/4，
分别为同列模式 0x10 与 0x1c；它们在连接处的规范化 DLCT 均为 +1，其余分支贡献为零。
因此 ρG=1/4+1/4=1/2，S=2⁻¹×2⁻³⁸=2⁻³⁹。仅取其中一条差分分支时评分为 2⁻⁴⁰。</p>
<h2>所有固定尾段的比较</h2>
<p><code>D2_repeat_l</code> 是附录 D.2 的迭代路径延长至 l 轮；<code>D3_chia_to_b</code>
表示附录 D.3 的第 a–b 个 χ 连同中间线性层。每个尾段搜索全部 64×31=1984 个单活跃列输入差分。</p>
{scans}
<p>附录 D.3 完整 4 轮路径虽然 wL²=82，但连接掩码有 61 个活跃列。对这个固定掩码，
单列范围 1984 个候选、双列范围 1,937,376 个候选的两轮前缀相关度都精确为零。
双列结果由逐输出差分分支的条件掩码检查得出，与常数取值无关。
这只限制当前前缀与连接掩码，不能推出完整 6 轮置换的相关度为零。</p>
<p>另外，把已有六列、七列 L 传播模式作为首个 χ 输出差分种子，分别选择一个局部最大 DDT 概率的输入前驱，
并尝试全部 64 个共同移位，共 {len(patterns['candidates'])} 个输入候选；精确前缀仍全部为零。
此试验未枚举全部多列差分，也未枚举所有并列最优前驱。</p>
<h2>前缀采样核对</h2>
<p>这部分只测量两轮前缀 G。4 组代表各使用 2²² 对均匀状态，保留对应常数。
全部精确值落在各自的 95% Wilson 区间内；区间未作多重比较校正。</p>
{measured}
<p>没有用小样本声称验证 2⁻³⁹ 等完整评分。若完整相关度确实为 2⁻³⁹，直接采样的数量级也需约 2⁷⁸ 对；
当前 2²² 对实验的用途是检验前缀计算与位序。</p>
<details><summary>查看最佳候选的完整 320 位端点</summary><p>各状态依次列出 x[0]…x[4] 的 64 位十六进制字。</p>{''.join(endpoints)}</details>
<h2>下一步</h2>
<p>优先扩充论文路径的截取尾段与其他入口，按“精确前缀代价 + 尾段平方相关度权重”比较；
围绕 4、5 轮候选保留多个相关尾路径，评估线性包和交叉项。
完整 wL²=82 路径需要更丰富的多列差分或不同的中间连接，再与 wL²=96 的可拼接路径比较。</p>
<h2>复现数据</h2><p>{links}</p>
<footer>源论文：Symmetric Twin Column Parity Mixers and their Applications，附录 D.2（印刷页 30）、D.3（印刷页 31）。
位序、逐 χ 掩码、权重、轮常数、DDT 分支、采样种子和原始计数均保存在数据文件中。
本表独立于原先 r_l=0 的 Table 13 基线。</footer>
</main></html>'''
    output = root / "table13_fixed_tail.html"
    output.write_text(content, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()

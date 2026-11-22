#!/usr/bin/env python3
"""Create a self-contained HTML table and summary CSV from saved raw results."""
import argparse
import csv
import html
import json
from pathlib import Path


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def endpoint(col, pattern):
    col, pattern = int(col), int(pattern)
    if pattern.bit_count() == 1:
        return f"({pattern.bit_length() - 1}, {col})"
    bits = ", ".join(str(i) for i in range(5) if pattern & (1 << i))
    return f"第 {col} 列，行 {{{bits}}}（0x{pattern:02x}）"


def number(value, digits=6):
    return f"{float(value):.{digits}f}" if value else "—"


def power(n):
    n = int(n)
    if n and n & (n - 1) == 0:
        return f"2^{n.bit_length() - 1}"
    return str(n)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", type=Path, required=True)
    parser.add_argument("--verification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    meta = json.loads((args.search / "search.json").read_text())
    verification = json.loads((args.verification / "verification.json").read_text())
    candidates = rows(args.search / "top.csv")
    experiments = rows(args.verification / "experiments.csv")
    best = {int(row["rounds"]): row for row in candidates if row["rank"] == "1"}
    observed = {}
    for row in experiments:
        if row["rank"] == "1" and row["rotation"] == "0":
            observed[int(row["rounds"])] = row
    # Join by full endpoints, not only round number. Reject unrelated measurements.
    for r, exp in observed.items():
        if r not in best:
            raise ValueError(f"experiment for R={r} has no selected model candidate")
        for field in ("input_col", "input_pattern", "output_col", "output_pattern"):
            if exp[field] != best[r][field]:
                raise ValueError(f"experiment endpoint mismatch at R={r}: {field}")
        if exp["model_rho"] != best[r]["model_rho"]:
            raise ValueError(f"experiment model score mismatch at R={r}")

    summary = []
    for r in meta["rounds"]:
        c, exp = best.get(r), observed.get(r)
        item = {"rounds": r, "status": "no_nonzero_model_candidate"}
        if c:
            item.update({field: c[field] for field in (
                "input_col", "input_pattern", "output_col", "output_pattern", "input_words", "output_words", "model_rho", "model_log2_abs")})
            item["input_endpoint"] = endpoint(c["input_col"], c["input_pattern"])
            item["output_endpoint"] = endpoint(c["output_col"], c["output_pattern"])
            item["status"] = "model_candidate"
        if exp:
            item.update({field: exp[field] for field in (
                "experimental_rho", "experimental_log2_abs", "samples", "n0", "n1", "ci95_low", "ci95_high", "significant", "start_round", "constants", "seed")})
            item["status"] = "experimentally_detected" if exp["significant"] == "1" else "not_significant_at_this_sample_size"
        summary.append(item)

    args.output.mkdir(parents=True, exist_ok=True)
    fields = ["rounds", "input_endpoint", "output_endpoint", "model_log2_abs", "experimental_log2_abs", "model_rho", "experimental_rho", "samples", "ci95_low", "ci95_high", "status", "input_col", "input_pattern", "output_col", "output_pattern", "input_words", "output_words", "n0", "n1", "significant", "start_round", "constants", "seed"]
    with (args.output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)

    main_rows = []
    for row in summary:
        if row["status"] == "no_nonzero_model_candidate":
            cells = [row["rounds"], "—", "—", "模型未检出", "—", "—"]
        else:
            exp_log = number(row.get("experimental_log2_abs", "")) if row.get("significant") == "1" else "未验证或未显著检出"
            cells = [row["rounds"], row["input_endpoint"], row["output_endpoint"], number(row["model_log2_abs"]), exp_log,
                     power(row["samples"]) if "samples" in row else "—"]
        main_rows.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in cells) + "</tr>")

    experiment_rows = []
    for row in experiments:
        cells = [row["rounds"], row["rank"], row["rotation"], endpoint(row["input_col"], row["input_pattern"]),
                 endpoint(row["output_col"], row["output_pattern"]), number(row["model_rho"]), number(row["experimental_rho"]),
                 f'[{number(row["ci95_low"])}, {number(row["ci95_high"])}]', power(row["samples"])]
        experiment_rows.append("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in cells) + "</tr>")
    scope = "输入与输出均为单比特" if meta["scope"] == "bit" else "输入与输出各有一个活跃 χ 实例（每端均允许 31 种非零列模式）"
    constant_note = ("保留后续参考轮常数" if verification["constants"] == "reference" else "后续轮常数全部置零")
    start_note = ("对 R 轮采用参考实现的末 R 轮，起始轮号为 12−R（从 0 编号）" if verification["start_round"] == "last"
                  else f"起始轮号固定为 {verification['start_round']}（从 0 编号）")
    underflow_note = "<p class='warning'>本次搜索检测到浮点下溢；不能将数值零视为模型的严格零值。</p>" if meta["underflow_detected"] else ""
    source_links = [
        ("模型配置", args.search / "search.json"), ("完整端点评分", args.search / "scores.csv"),
        ("候选列表", args.search / "top.csv"), ("采样原始计数", args.verification / "experiments.csv"),
        ("采样配置", args.verification / "verification.json"), ("汇总 CSV", args.output / "summary.csv")]
    links = " · ".join(f'<a href="{html.escape(str(p.resolve()), quote=True)}">{title}</a>' for title, p in source_links)
    content = f"""<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston 差分线性区分器：Table 13 型结果</title>
<style>
:root{{color-scheme:light}} body{{margin:0;background:#f4f6f8;color:#182631;font:16px/1.75 system-ui,"Noto Sans CJK SC",sans-serif}}
main{{max-width:1180px;margin:36px auto;padding:40px;background:white;border-radius:14px;box-shadow:0 4px 24px #14273810}}
h1{{font-size:30px;margin:0 0 8px}} h2{{font-size:21px;margin-top:32px}} .muted{{color:#526575}} .scope{{border-left:4px solid #277f91;padding:8px 18px;background:#eff8fa}}
.scroll{{overflow-x:auto}} table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:15px;margin:18px 0}}
th,td{{padding:12px 14px;text-align:center;border-bottom:1px solid #dce4e8;white-space:nowrap}} th{{border-top:2px solid #234250;background:#f7fafb}} tbody tr:last-child td{{border-bottom:2px solid #234250}}
code{{font-size:.94em;background:#eef3f5;padding:2px 5px;border-radius:4px}} a{{color:#116d82}} .warning{{padding:12px 16px;background:#fff4dd;border-left:4px solid #b87920}} .formula{{text-align:center;font-size:19px;margin:24px 0}}
details{{margin-top:28px}} summary{{cursor:pointer;font-weight:600}} footer{{font-size:13px;color:#526575;margin-top:36px}} @media(max-width:700px){{main{{margin:0;padding:22px;border-radius:0}}h1{{font-size:24px}}}}
</style><main>
<h1>Gaston 差分线性区分器</h1>
<p class="muted">Table 13 型结果 · 首个线性层省略 · r<sub>d</sub>=1，r<sub>l</sub>=0</p>
<p class="scope">搜索范围：{scope}。每个轮数覆盖 {meta['endpoint_pairs_per_round']:,} 对端点；
{'利用模型中的共同循环移位等价性' if meta['rotation_reduction'] else '直接枚举所有列位置'}，实际计算 {meta['evaluated_representatives_per_round']:,} 个代表。
最优性仅针对本次端点范围和近似模型。</p>
<p class="formula">χ → L → χ → ··· → L → χ &nbsp;（R 个 χ，R−1 个 L）</p>
<p>输入从第一个 χ 的输入边界开始，输出位于最后一个 χ 之后。上图省略加常数符号；实验{constant_note}。
{start_note}。第一个 χ 之前的常数通过均匀输入变量的平移吸收。</p>
<div class="scroll"><table><thead><tr><th>轮数 R</th><th>输入差分 Δ</th><th>输出掩码 Λ</th><th>逐轮近似<br>log₂|ρ|</th><th>实验测量<br>log₂|ρ|</th><th>样本对数 N</th></tr></thead>
<tbody>{''.join(main_rows)}</tbody></table></div>
<p>单比特端点 <code>(i,j)</code> 表示第 i 个 64 位字的第 j 位；两者均从 0 编号，第 0 位是最低位。
多比特列模式另行列出其行集合。实验列对应每轮模型排名第一的同一端点，不根据验证样本重新选择候选。</p>
{underflow_note}
<h2>模型含义与适用范围</h2>
<p>实现采用 HNW 的零状态值掩码优化：第一轮按 5 比特 χ 的 DDT 精确加权，后续用 LAT 相关度的平方更新 32 维局部向量，
通过 L 的转置传播掩码，并用各列相关度的乘积近似跨列依赖。完整固定常数置换的相关度并未被精确枚举。</p>
<p>两轮的单活跃端点结果还经过了线性满秩检查和独立局部精确枚举。三轮的模型值由采样检验。
对于本配置的四轮及后续轮，局部乘积模型丢失了继续传播所需的相关信息：测试确认第四轮前所有相关的非平凡线性乘积均含精确零因子，
并非浮点下溢。“模型未检出”不表示真实置换不存在非零差分线性相关度，也不是完整构造的界限结论。</p>
<p>每个实验使用 {power(verification['samples_per_candidate'])} 个样本对；置信区间是逐项 95% Wilson 区间（映射到有符号相关度），不提供多重比较的同时覆盖保证。
共同循环移位只用于近似模型的枚举约简；固定轮常数下的真实相关度没有被假定具有同样的对称性。</p>
<details><summary>查看全部 {len(experiments)} 组采样核对</summary><div class="scroll"><table><thead><tr>
<th>R</th><th>模型排名</th><th>移位</th><th>Δ</th><th>Λ</th><th>模型 ρ</th><th>实验 ρ</th><th>95% 区间</th><th>N</th>
</tr></thead><tbody>{''.join(experiment_rows)}</tbody></table></div></details>
<h2>复现数据</h2><p>{links}</p>
<p>方法来源：<a href="https://eprint.iacr.org/2026/358">HNW，Round-Based Approximation of (Higher-Order) Differential-Linear Correlation</a>；
对照了<a href="https://github.com/hukaisdu/Round-based-approximation/tree/783e8edc560b848d0f2ba63f11838aecd23da40f">作者公开实现</a>的零状态值掩码优化流程。</p>
<footer>模型：{html.escape(meta['model'])} · 浮点：{meta['numeric_type']} / {meta['mantissa_bits']} 位尾数 · 数据由程序输出，表格由原始 CSV 自动生成。</footer>
</main></html>"""
    (args.output / "table13.html").write_text(content, encoding="utf-8")
    print(args.output.resolve() / "table13.html")


if __name__ == "__main__":
    main()

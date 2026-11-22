#!/usr/bin/env python3
"""Validate saved evidence and generate the Chinese Hsilu DL research report."""
import csv
import hashlib
import html
import json
import math
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path
from hsilu import transpose, words
from pool_search import completions, decode, local_score

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'results'
BASE = RESULTS/'baseline'
esc = html.escape


def read(path):
    return json.loads(path.read_text())


def stats(path):
    result = {}
    for source in csv.DictReader(path.open()):
        row = {k: (v if k == 'id' else float(v) if k == 'rho' else int(v))
               for k, v in source.items()}
        n, n0 = row['samples'], row['n0']
        assert n0 + row['n1'] == n
        rho = (n0-row['n1'])/n
        assert rho == row['rho']
        z, p = 1.959963984540054, n0/n
        den = 1+z*z/n
        center = (p+z*z/(2*n))/den
        radius = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
        row.update(log2_abs=math.log2(abs(rho)) if rho else None,
                   ci95=[2*(center-radius)-1, 2*(center+radius)-1])
        assert row['id'] not in result
        result[row['id']] = row
    m = len(result)
    for row in result.values():
        row['joint_hoeffding_threshold'] = math.sqrt(2*math.log(2*m/.05)/row['samples'])
        row['joint_significant'] = abs(row['rho']) > row['joint_hoeffding_threshold']
    return result


def validate_pool():
    pool = read(BASE/'linear_pool.json')
    total = 0
    for bridge in pool['bridges']:
        b1, a2 = decode(bridge['first_output']), decode(bridge['second_input'])
        assert transpose(a2) == b1
        assert {decode(x['output']) for x in bridge['last_chi_choices']} == set(completions(a2, False))
        for last in bridge['last_chi_choices']:
            assert last == local_score(a2, decode(last['output']))
        heads = [h for h in pool['heads'] if h['bridge_id'] == bridge['id']]
        assert {decode(h['first_chi']['input']) for h in heads} == set(completions(b1, True))
        for head in heads:
            a1 = decode(head['first_chi']['input'])
            assert decode(head['join_mask']) == transpose(a1)
            assert head['first_chi'] == local_score(a1, b1)
        total += len(heads)*len(bridge['last_chi_choices'])
    assert len(pool['heads']) == 384 and total == 38400
    heads = {h['id']: h for h in pool['heads']}
    unique, nonzero = set(), 0
    for row in csv.DictReader((BASE/'prefix_matches.csv').open()):
        nonzero += 1
        assert int(row['prefix_denominator']) == 32
        assert 0 < abs(int(row['prefix_numerator'])) <= 32
        head = heads[row['head_id']]
        for last in pool['bridges'][head['bridge_id']]['last_chi_choices']:
            unique.add((int(row['input_col']), int(row['input_pattern']), tuple(last['output']['words'])))
    saved = read(BASE/'pool_candidates.json')
    assert nonzero == saved['prefix_nonzero_matches'] == 11805
    assert len(unique) == saved['unique_endpoints'] == 102700
    five = [json.loads(s) for s in (BASE/'five_prefixes.jsonl').read_text().splitlines()]
    assert {(x['input_col'], x['input_pattern']) for x in five} == {(c, p) for c in range(64) for p in range(1, 32)}
    assert sum(x['nonzero'] for x in five) == 6005
    assert sum(x['resource_limited'] for x in five) == 0
    retained = set()
    for x in five:
        assert len(x['top_heads']) <= 4
        for chosen in x['top_heads']:
            head = heads[chosen['head_id']]
            tails = pool['bridges'][head['bridge_id']]['last_chi_choices']
            best = min(t['weight'] for t in tails)
            for t in tails:
                if t['weight'] == best:
                    retained.add((x['input_col'], x['input_pattern'], tuple(t['output']['words'])))
    assert len(retained) == 1000
    return {'heads':384, 'two_chi_trails':total, 'four_prefix_evaluations':384*1984,
            'four_nonzero_prefixes':nonzero, 'four_unique_endpoints':len(unique),
            'five_prefix_evaluations':len(five)*384, 'five_nonzero_prefix_DATF':6005,
            'five_resource_limited':0, 'five_retained_unique_endpoints':len(retained)}


def validate_jobs(path, selected, full=False):
    jobs = {}
    for line in path.read_text().splitlines():
        ident, r, d, m = line.split()
        row = selected[ident]
        assert int(r) == row['rounds']
        assert d.split(';') == row['input_difference']['words']
        assert m.split(';') == row['full_output_mask' if full else 'core_output_mask']['words']
        jobs[ident] = row
    return set(jobs)


def frac(record):
    return Fraction(int(record['numerator']), int(record['denominator']))


def bitlabel(record):
    return '{' + ', '.join(f'({i},{j})' for i, j in record['bits']) + '}'


def main():
    audit = read(RESULTS/'structure_audit.json')
    assert audit['status'] == 'PASS'
    coverage = validate_pool()
    four = {r['id']:r for r in read(BASE/'pool_selection.json')}
    five = {r['id']:r for r in read(BASE/'five_selection.json')}
    datf4 = {r['id']:r for r in read(BASE/'pool_datf.json')['results']}
    datf5 = {r['id']:r for r in read(BASE/'five_datf.json')['results']}
    initial = stats(BASE/'initial_samples.csv')
    discovery, confirm = stats(BASE/'pool_samples.csv'), stats(BASE/'confirmation.csv')
    samples5, reference = stats(BASE/'five_samples.csv'), stats(BASE/'reference_samples.csv')
    assert set(four) == set(datf4) == set(discovery) == validate_jobs(BASE/'pool_jobs.txt', four)
    assert set(five) == set(datf5) == set(samples5) == validate_jobs(BASE/'five_jobs.txt', five)
    assert set(confirm) == validate_jobs(BASE/'confirmation_jobs.txt', four)
    assert set(reference) == validate_jobs(BASE/'reference_jobs.txt', four, True)
    assert len(four) == len(five) == 16 and len(confirm) == 4 and len(reference) == 2
    assert {v['seed'] for v in discovery.values()}.isdisjoint({v['seed'] for v in confirm.values()})
    hnws = read(RESULTS/'hnw/summary.json')['results']
    hnw = {}
    for row in hnws:
        raw = read(ROOT/row['path'])
        assert raw['correlation'] == row['correlation']
        hnw[row['id'], row['mode']] = str(frac(row['correlation']))
    for family, estimates in [(four,datf4),(five,datf5)]:
        for ident, row in family.items():
            assert transpose(decode(row['full_output_mask'])) == decode(row['core_output_mask'])
            estimate = estimates[ident]
            assert estimate['status'] == 'completed'
            assert estimate['input_difference'] == row['input_difference']
            assert estimate['core_output_mask'] == row['core_output_mask']
            c = Fraction(row['head']['first_chi']['correlation'])*Fraction(row['last_chi']['correlation'])
            prefix = Fraction(row['exact_prefix_correlation'] if row['rounds']==4 else row['prefix_datf'])
            assert prefix*c*c == Fraction(row['splice_score'])
            if row['rounds'] == 4:
                delta = row['input_difference']['bits']
                col = delta[0][1]; pattern = sum(1 << i for i,j in delta)
                value = int(subprocess.check_output([ROOT/'build/hsilu_dl','prefix-eval',str(col),str(pattern),
                    ';'.join(row['head']['join_mask']['words']),'69'], text=True))
                assert Fraction(value,32) == prefix
    baseline = {}
    for row in csv.DictReader((BASE/'hnw_zero.csv').open()):
        r, rho = int(row['rounds']), float(row['rho'])
        record = baseline.setdefault(r, {'tested':0,'model_nonzero':0,'max_abs':0})
        record['tested'] += 1; record['model_nonzero'] += rho != 0
        record['max_abs'] = max(record['max_abs'],abs(rho))
    assert set(baseline) == set(range(2,7))
    assert all(v['tested']==61504 for v in baseline.values())
    assert all(v['joint_significant'] for v in confirm.values())
    assert all(v['joint_significant'] for v in reference.values())
    assert not any(v['joint_significant'] for v in samples5.values())
    sparse = read(RESULTS/'sparse_diff/candidates.json')
    assert sparse['source_trails']==18 and sparse['distinct_prefix_orbits']==10
    assert sparse['middle_nonzero_matches']==sparse['unique_endpoints']==0
    provenance = {}
    paths = [p for folder in ('include','src','tools','tests') for p in (ROOT/folder).glob('*') if p.is_file()]
    paths += [ROOT/'Makefile', ROOT.parents[1]/'hsilu_ref/hsilu.cpp', ROOT.parents[1]/'hsilu_ref/hsilu.h']
    paths += [ROOT.parents[1]/'hsilu_ref/我的已发布论文：Multiple Rows Mixers and Hsilu.pdf']
    paths += list({Path(s) for p in sparse['prefixes'] for s in p['sources']})
    paths += [p for p in BASE.iterdir() if p.is_file()]
    paths += [RESULTS/'structure_audit.json', RESULTS/'hnw/summary.json', RESULTS/'sparse_diff/candidates.json',
              ROOT.parent/'gaston_dl/tools/datf.py', ROOT.parent/'gaston_dl/tools/datf_k06.py']
    for path in paths:
        provenance[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {'status':'PASS','definition':{'round':'RC -> chi -> L',
               'core':'first RC absorbed by uniform translation; last L omitted',
               'full_mask':'L^{-T}(core mask)','constant_indices':'last R of the 12 reference constants',
               'rho':'E_X (-1)^(lambda dot (F(X) xor F(X xor delta)))','bias':'rho/2'},
               'structural':audit['structural_families'],'coverage':coverage,'baseline_HNW_zero':baseline,
               'four_candidates':list(four.values()),'four_DATF':list(datf4.values()),
               'four_discovery':list(discovery.values()),'four_confirmation':list(confirm.values()),
               'full_reference_confirmation':list(reference.values()),
               'five_candidates':list(five.values()),'five_DATF':list(datf5.values()),
               'five_sampling':list(samples5.values()),'HNW':hnws,
               'sparse_two_round_prefixes':{'source_files':18,'orbits':10,'matches':0,'tested':10*64*384},
               'sampling_note':'Pseudorandom simulations. Wilson and Hoeffding interpretation assumes IID within an endpoint. Endpoints may share states.',
               'provenance_sha256':provenance}
    (RESULTS/'research_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    def fmtci(row):return f'[{row["ci95"][0]:.7f}, {row["ci95"][1]:.7f}]'
    confirm_rows = ''
    for ident in sorted(confirm, key=lambda x:-abs(confirm[x]['rho'])):
        row,sample=four[ident],confirm[ident]
        confirm_rows += f'<tr><td><a href="#{ident}">{ident}</a></td><td>{bitlabel(row["input_difference"])}</td><td>{bitlabel(row["core_output_mask"])}</td><td>{sample["rho"]:.9f}</td><td>{sample["log2_abs"]:.5f}</td><td>{fmtci(sample)}</td><td>{frac(datf4[ident]["correlation"])}</td></tr>'
    details = ''
    for ident in sorted(confirm):
        row = four[ident]
        full_words='\n'.join(f'A{i}: {w}' for i,w in enumerate(row['full_output_mask']['words']))
        details += f'<details id="{ident}"><summary>{ident}：端点与线性尾段</summary><p>Δ={bitlabel(row["input_difference"])}；核心 Λ={bitlabel(row["core_output_mask"])}。精确两 χ 前缀相关度 {row["exact_prefix_correlation"]}，单尾段局部相关度 {row["head"]["first_chi"]["correlation"]} × {row["last_chi"]["correlation"]}，拼接评分 {row["splice_score"]}。</p><p>第三 χ：{bitlabel(row["head"]["first_chi"]["input"])} → {bitlabel(row["head"]["first_chi"]["output"])}；第四 χ：{bitlabel(row["last_chi"]["input"])} → {bitlabel(row["last_chi"]["output"])}。二者之间包含 L。</p><p>完整置换输出掩码 L⁻ᵀΛ，重量 {row["full_output_mask"]["bit_weight"]}：</p><pre>{full_words}</pre></details>'
    comparisons = ''
    for ident in ['R4_pool_000000','R4_pool_000027','R5_pool_000000']:
        row=(four|five)[ident]; d=(datf4|datf5)[ident]
        comparisons += f'<tr><td>{ident}</td><td>{esc(row["splice_score"]) if row["rounds"]==4 else "log₂|评分|≈−58.95959"}</td><td>{frac(d["correlation"]) if row["rounds"]==4 else "log₂|ρ|≈−58.95959"}</td><td>{hnw[ident,"zero"]}</td><td>{hnw[ident,"full"]}</td><td>{confirm[ident]["rho"] if ident in confirm else "未确认"}</td></tr>'
    all_rows = ''
    for ident,row in (four|five).items():
        d=(datf4|datf5)[ident]; sample=(discovery|samples5)[ident]
        log = d['correlation']['log2_abs']
        all_rows += f'<tr><td>{ident}</td><td>{row["rounds"]}</td><td>{bitlabel(row["input_difference"])}</td><td>{bitlabel(row["core_output_mask"])}</td><td>{"模型零" if log is None else f"{log:.6f}"}</td><td>{sample["rho"]:+.9f}</td><td>{"2²⁶" if row["rounds"]==4 else "2²⁴"}</td></tr>'
    refs=''.join(f'<tr><td>{ident}</td><td>{s["rho"]:.9f}</td><td>{fmtci(s)}</td><td>是</td></tr>' for ident,s in reference.items())
    structural=audit['structural_families'][1]
    smask='\n'.join(f'A{i}: {w}' for i,w in enumerate(structural['masks'][0]['full']['words']))
    bridge_rows=''.join(f'<tr><td>{bitlabel(b["alpha"])}</td><td>{bitlabel(b["beta"])}</td><td>{b["alpha"]["active_columns"]} → {b["beta"]["active_columns"]}</td></tr>' for b in audit['linear_bridges']['bridges'])
    base_rows=''.join(f'<tr><td>{r}</td><td>{v["tested"]:,}</td><td>{v["model_nonzero"]:,}</td><td>{v["max_abs"]:g}</td></tr>' for r,v in baseline.items())
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Hsilu 差分线性研究：结构证明、四轮候选与五轮搜索</title>
<style>body{{margin:0;background:#f2f5f7;color:#1e3343;font:16px/1.8 system-ui,sans-serif}}main{{max-width:1200px;margin:auto;padding:32px 24px}}h1{{line-height:1.35}}h2{{font-size:23px}}section{{background:white;padding:25px;margin:22px 0;border:1px solid #d5e0e7;border-radius:12px}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:10px;border-bottom:1px solid #dce5eb;text-align:left;vertical-align:top}}th{{background:#edf3f7}}.scroll{{overflow-x:auto}}pre{{background:#eef3f6;padding:16px;overflow:auto;line-height:1.6}}a{{color:#14618d}}details{{border-top:1px solid #dce5eb;padding:15px 0}}summary{{cursor:pointer;font-weight:600}}.note{{border-left:5px solid #2f7f98;padding-left:16px}}.cards{{display:flex;flex-wrap:wrap;gap:16px}}.card{{background:#eaf2f5;padding:17px;flex:1;min-width:210px}}.number{{font-size:27px;font-weight:700}}nav{{display:flex;gap:15px;flex-wrap:wrap}}footer{{font-size:14px;color:#566b7b}}</style></head><body><main>
<h1>Hsilu 差分线性研究</h1><p>2026-10-18 · 320 比特参考置换 · 明确区分结构证明、近似模型与采样结果</p>
<nav><a href="#definition">边界定义</a><a href="#structure">三轮证明</a><a href="#four">四轮结果</a><a href="#five">五轮搜索</a><a href="#reproduce">复现与数据</a></nav>
<section><div class="cards"><div class="card"><div class="number">3 轮：ρ = 1</div>输入位依赖分析给出的严格结果</div><div class="card"><div class="number">4 轮：ρ̂ ≈ 2⁻⁹</div>独立采样及完整参考置换复核</div><div class="card"><div class="number">5 轮：未确认</div>候选普通 DATF 约 2<sup>−58.96</sup>；采样不足以测量此量级</div></div>
<p class="note">当前可报告三轮结构性区分器与四轮实验非零端点。未证明四轮最优，未证明五轮不存在良好相关度，也没有从减轮置换推出任何完整模式的复杂度结论。</p></section>
<section id="definition"><h2>1. 置换、轮数和输出掩码</h2>
<p>依据本地论文 Table 9 与参考实现，Hsilu 一轮为 RC → χ → L，其中 L 是 MRM 后接行移位。状态是 A₀…A₄ 五个 64 位字；坐标 (i,j) 中 j=0 表示最低位。参考实现使用 12 个常数，减轮取末 R 个；本报告不把 12 轮解释为论文唯一规定的轮数。</p>
<p>搜索核心 Fᵣ 从首个 χ 输入开始，到最后一个 χ 输出结束：χ → (L → RC → χ)^(R−1)。这里省略的是最后一个 L，与 Gaston 项目的入口约定不同。首个 RC 通过完整均匀输入的平移吸收，其他常数保留。</p>
<pre>Pᵣ(X) = L(Fᵣ(X ⊕ c_first))
Δ_full = Δ_core
Λ_full = L⁻ᵀ(Λ_core)
ρ = E[ (−1)^(Λ · (F(X) ⊕ F(X ⊕ Δ))) ]
偏差 ε = ρ / 2</pre>
<p>因此核心端点与转换后的完整置换端点具有同一均匀状态相关度。报告表格默认展示核心掩码；展开候选可见完整输出掩码，不能把稀疏核心掩码直接套到 Pᵣ 的输出。</p>
<p>四轮核心保留常数 69、5a、4b；五轮核心保留 78、69、5a、4b。五轮的三 χ 前缀对应常数索引 7、8、9，首个吸收后余下 78、69；没有误用末三轮的常数序列。</p>
<p>论文报告差分分支数 12、线性分支数 5；这里的活跃数按列/S 盒计。论文线性路径的平方相关度权重与本程序 wL=−log₂|cL| 相差因子 2，不能混用。<a href="https://eprint.iacr.org/2025/1960">原论文及修订说明</a>。</p></section>
<section id="structure"><h2>2. 三轮相关度 1 的结构证明</h2>
<p>取 Δ={{(3,0)}}，第三个 χ 输出处的 (0,29) 与 (1,29) 均不依赖这个输入比特。因此选择这两个输出位中的任意非空子集作为掩码，差分奇偶值对所有输入状态都恒为 0，相关度严格等于 1。</p>
<p>证明使用保守依赖传播：对线性层的 XOR 门，将输入依赖集合取并集；对 χ 的每个输出，取其三个输入的依赖集合并集。这个过程只可能多计依赖，不会漏计依赖。经过各 χ 后的可达位数是 3 → 66 → 318，最后仅上述两位不可达。符号差分 DAG 也直接化简为零，未调用 DATF 的独立性近似来证明。</p>
<p>所有 64 个共同位移都已检查，得到 Δ={{(3,k)}}、核心掩码位于 {{(0,k+29),(1,k+29)}} 的 192 个端点对，位号模 64。依赖证明不依赖轮常数取值。两轮另有 Δ={{(0,0)}}、核心 Λ={{(0,1)}} 的同类严格 ρ=1 结果。</p>
<p>以三轮 Λ={{(0,29)}} 为例，完整输出掩码重量为 134：</p><pre>{smask}</pre>
<p>采样 2²² 对也全部得到偶奇偶值，但严格性来自依赖证明，而不是“没有观察到反例”。这与纯差分或纯线性路径重量较大不矛盾：此处描述的是某个输出奇偶函数对选定输入比特的不敏感性。</p>
<p><a href="structure_audit.json">全部依赖集合、完整掩码与证明检查记录</a></p></section>
<section><h2>3. 两轮线性尾段：最小分支连接的完整枚举</h2>
<p>设尾段为 a₁ → χ → b₁ → L → a₂ → χ → b₂，内部连接满足 b₁=Lᵀ(a₂)。枚举关系两端分别只有一列或两列的情形，共 121,148 项；当两端总活跃数不超过 5 时，至少一端不超过两列，所以覆盖该范围的全部可能。共同移位去重后只有下列三个类别，均为 3→2，没有 1→4、2→3 或 4→1 类别。</p>
<div class="scroll"><table><tr><th>b₁</th><th>a₂</th><th>活跃列数</th></tr>{bridge_rows}</table></div>
<p>补全两侧 χ 的全部非零 LAT 转移后，得到 38,400 条两 χ 线性路径代表、384 个不同连接掩码 γ=Lᵀ(a₁)。全部局部 LAT 值及补全集合均重新核对。</p>
<p>共同移位是 L 的严格对称性；没有据此宣称完整固定常数置换的相关度也具有该对称性。本次四轮与五轮匹配固定这三个尾段方向，遍历输入的全部 64 个位置；仍未穷尽所有绝对尾段移位及更高活跃尾段。</p></section>
<section id="four"><h2>4. 四轮：精确两 χ 前缀 + 两 χ 线性尾段</h2>
<p>首差分限制为一个活跃列，允许该列的全部 31 个非零模式。对每个输入及 γ，条件化首列的 5 个基值，精确求前两个 χ 的差分线性相关度；其他首轮 χ 输出列独立均匀，第二个 χ 的掩码导数为仿射函数，故可通过线性消元/32 项求和得到精确值。</p>
<p>共计算 1,984×384=761,856 个精确前缀，得到 11,805 个非零连接和 102,700 个候选外端点。用 ρ_prefix·cL² 排序只是一种拼接评分；精确前缀不使整条四轮评分自动精确。筛选 16 个多样化候选，先各采样 2²⁶ 对，再预选 4 个端点，以新种子各采样 2²⁸ 对。</p>
<div class="scroll"><table><tr><th>候选</th><th>Δ</th><th>核心 Λ</th><th>独立复核 ρ̂</th><th>log₂|ρ̂|</th><th>逐项 95% Wilson 区间</th><th>完整普通 DATF</th></tr>{confirm_rows}</table></div>
<p>四项独立复核的联合 Hoeffding 阈值为 {next(iter(confirm.values()))['joint_hoeffding_threshold']:.8f}；四项均超过阈值。区间本身是逐项区间；没有因多次筛选而把最大观测值宣称为全局最优。</p>
{details}
<h3>直接调用完整参考置换的交叉复核</h3>
<p>另取两个端点，使用 hsilu_ref/hsilu.cpp、完整输出掩码和不同生成器 std::mt19937_64，各采样 2²⁴ 对。这里实际执行全部四轮最后的 L。两项联合阈值为 {next(iter(reference.values()))['joint_hoeffding_threshold']:.8f}，结果与原核心实验相容：</p>
<div class="scroll"><table><tr><th>候选</th><th>完整参考置换 ρ̂</th><th>逐项 95% Wilson 区间</th><th>超过联合阈值</th></tr>{refs}</table></div>
<p>采样生成器用于可复现仿真；统计阈值采用每个端点内样本独立的模型，不提供信息论随机性证明。不同端点之间允许共享基状态。</p></section>
<section><h2>5. DATF、HNW26 和采样不能混为一个结论</h2>
<div class="scroll"><table><tr><th>端点</th><th>前缀×尾段评分</th><th>完整普通 DATF</th><th>HNW26 零值掩码模型</th><th>HNW26 完整局部模型</th><th>独立实验</th></tr>{comparisons}</table></div>
<p>R4_pool_000027 在两个近似方法中均为零，而采样稳定非零。HNW26 使用精确二进有理算术，只说明模型递推没有浮点舍入误差，不代表整个固定置换的真实相关度被精确求和。</p>
<p>这批数据直接说明：不能把 DATF/HNW26 零值当作删去所有候选的严格判据。高质量尾段拼接可用于产生候选，随后仍须完整端点评估、独立实验或结构性证明。</p>
<details><summary>最初的单活跃列端点扫描，2–6 轮</summary><p>每轮 61,504 个模型代表：输入列固定为 0，枚举 31 个输入模式与 64×31 个输出掩码。此处约简仅用于该模型，固定常数置换不据此认定等价。</p><table><tr><th>轮数</th><th>模型端点</th><th>模型非零</th><th>最大模型 |ρ|</th></tr>{base_rows}</table></details></section>
<section id="five"><h2>6. 五轮搜索及其实际限制</h2>
<p>第一条路线是三 χ 前缀加两 χ 尾段。普通 DATF 完成 761,856 个前缀计算，模型非零 6,005 项，无资源超限。每个输入仅保留评分最高的四个连接，并只补最小权重末 χ，得到 1,000 个不同端点，再选 16 个作完整普通 DATF 与各 2²⁴ 对采样。</p>
<p>最强保留端点为 Δ={{(0,57),(2,57),(3,57)}}、核心 Λ={{(0,10),(2,23)}}；完整普通 DATF 的 log₂|ρ|≈−58.959588，符号为负。HNW26 两种模型均为零。</p>
<p class="note">2²⁴ 对样本无法验证约 2⁻⁵⁹ 的相关度。此次五轮采样只用于排查模型遗漏的大相关度；16 项最大 |ρ̂| 为 {max(abs(s['rho']) for s in samples5.values()):.9f}，均未超过联合阈值 {next(iter(samples5.values()))['joint_hoeffding_threshold']:.9f}。不能据此说模型值正确，也不能说真实值为零。</p>
<p>第二条路线导入已有 SAT 差分路径：18 个文件去重为 10 个两 χ 前缀类别（8 类 6→6，2 类 13→3）。独立 DDT 与正向 L 检查通过；通过正向基向量构造的秩证书确认其两轮局部概率乘积成立。对 10×64×384=245,760 个连接，未找到中间 χ 导数恒定的非零匹配。</p>
<p>这里仅排除了这些前缀与当前尾段池之间的“中间 χ 对均匀输入非零/导数恒定”连接。条件分布上的非恒定导数、更高活跃线性尾段、不同多列输入和其他差分前缀均未被排除。<a href="sparse_diff/candidates.json">差分前缀与匹配记录</a>。</p>
<p>后续优先扩展线性桥的总活跃数到 6 或 7，再考虑跨多列输入；同时保留精确前缀核对与部分模型零值的实验抽查。若引入 MILP/SAT，它适合筛选活动模式与局部权重，不能直接把一般整体相关度或符号抵消当作简单可加的目标。</p></section>
<section><h2>7. 已选 32 个四/五轮端点</h2><p>四轮此表为首次 2²⁶ 对筛选；最终判断使用前面的独立复核。五轮为 2²⁴ 对有限采样。模型零表示数值 0，而不是资源不足或负无穷下溢标记。</p>
<div class="scroll"><table><tr><th>编号</th><th>轮数</th><th>Δ</th><th>核心 Λ</th><th>DATF log₂|ρ|</th><th>首次实验 ρ̂</th><th>样本对数</th></tr>{all_rows}</table></div></section>
<section id="reproduce"><h2>8. 复现、数据与检查范围</h2>
<p><a href="research_summary.json">完整机器可读结果与 SHA-256</a> · <a href="structure_audit.json">结构证书</a> · <a href="baseline/confirmation.csv">四轮独立复核</a> · <a href="baseline/reference_samples.csv">完整参考实现复核</a> · <a href="baseline/five_samples.csv">五轮样本计数</a> · <a href="hnw/summary.json">HNW26 汇总</a></p>
<p>代码和命令说明：<a href="../tools/audit_structure.py">结构枚举</a>、<a href="../tools/pool_search.py">线性尾段与匹配排序</a>、<a href="../tools/search_five.py">五轮前缀搜索</a>、<a href="../tools/evaluate_candidates.py">普通 DATF</a>、<a href="../src/reference_verify.cpp">完整置换独立采样</a>。</p>
<pre>cd /home/yxb/分析/差分线性/hsilu_dl
make test
make build/hnw_endpoint build/reference_verify
./build/hnw_endpoint --selftest
python3 tools/audit_structure.py --output results/structure_audit.json
python3 tools/report.py

# 用新的输出文件复核已有端点；程序拒绝覆盖采样文件
./build/hsilu_dl verify results/baseline/confirmation_jobs.txt results/new_confirmation.csv 268435456 2026101803 4
./build/reference_verify results/baseline/reference_jobs.txt results/new_reference.csv 16777216 2026101821</pre>
<details><summary>从新目录重新枚举四轮/五轮候选</summary><p>下列命令不覆盖现有结果。五轮搜索默认使用四个工作进程；每个端点每次 DATF 计算的预算在脚本中记录。不要把资源超限记录解释为零。</p><pre>python3 tools/audit_structure.py --export-bridges results/reproduce/linear_bridges.json
python3 tools/pool_search.py build --output results/reproduce
./build/hsilu_dl prefix-scan results/reproduce/heads.txt results/reproduce/prefix_matches.csv
python3 tools/pool_search.py rank --output results/reproduce
python3 tools/evaluate_candidates.py --source results/reproduce/pool_selection.json --output results/reproduce/pool_datf.json
./build/hsilu_dl verify results/reproduce/pool_jobs.txt results/reproduce/pool_samples.csv 67108864 2026101802 4
python3 tools/search_five.py --output results/reproduce --workers 4
python3 tools/evaluate_candidates.py --source results/reproduce/five_selection.json --output results/reproduce/five_datf.json
./build/hsilu_dl verify results/reproduce/five_jobs.txt results/reproduce/five_samples.csv 16777216 2026101804 4</pre></details>
<p>核心检查覆盖参考实现 1–12 轮、L 与 Lᵀ、320 基逆映射、2/3/4 χ 符号 DAG 与实际值、五轮前缀常数位置、独立仿射求和的两 χ 精确前缀。HNW26 检查覆盖局部二元 LAT、条件 DDT、线性传播和轮常数。报告生成器重新核对尾段补全集合、候选端点与任务文件、样本计数、阈值及输入输出掩码。</p>
<p>结构枚举未依赖 MILP。已有 SAT 差分路径仅作为输入材料，其来源与逐条复核记录在 JSON 中。普通 DATF 复用 Gaston 项目的通用引擎，Hsilu 的线性层与轮边界单独实现，未修改 Gaston 引擎。</p>
<p>来源：<a href="https://eprint.iacr.org/2025/1960">Multiple Rows Mixers and Hsilu</a>（修订版指出第 5 节证明已更正，本研究没有依赖旧版该节证明）；本地 Table 9 及 hsilu_ref/hsilu.cpp。HNW26 方法依据 <a href="https://eprint.iacr.org/2026/358">论文 2026/358</a>，采用本工作区已核对的局部递推实现。</p></section>
<footer>研究结论的范围：减轮 Hsilu 置换、指定端点与指定搜索域。未给出全局最优性、五轮不可能性或完整 XOF/AEAD 模式的复杂度结论。</footer>
</main></body></html>'''
    (RESULTS/'hsilu_dl_report.html').write_text(page)
    print('PASS report evidence: pool completeness, endpoint boundaries, 4/5-round counts, models, independent samples')
    print(RESULTS/'hsilu_dl_report.html')


if __name__ == '__main__':
    main()

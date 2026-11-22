#!/usr/bin/env python3
"""Audit existing Gaston DL endpoints at the actual XOF/AEAD rate boundary.

All rank/support statements are finite GF(2) computations. They do not imply
that a uniform-state correlation transfers to a message/key distribution.
"""
import argparse
import hashlib
import html
import itertools
import json
import math
from pathlib import Path

from linear_seed_scan import linear

ROOT = Path(__file__).resolve().parents[1]
MASK = (1 << 64) - 1


def pack(state):
    return sum(v << (64 * i) for i, v in enumerate(state))


def unpack(value):
    return tuple((value >> (64 * i)) & MASK for i in range(5))


def words(state):
    return [f'{v:016x}' for v in state]


def bits(state):
    return [[i, j] for i, v in enumerate(state) for j in range(64) if v >> j & 1]


def inverse_basis():
    pivots = {}
    for i in range(320):
        y, x = pack(linear(unpack(1 << i))), 1 << i
        while y:
            k = y.bit_length() - 1
            if k not in pivots:
                pivots[k] = y, x
                break
            a, b = pivots[k]
            y, x = y ^ a, x ^ b
    assert len(pivots) == 320
    inverse = []
    for i in range(320):
        y, x = 1 << i, 0
        while y:
            a, b = pivots[y.bit_length() - 1]
            y, x = y ^ a, x ^ b
        assert pack(linear(unpack(x))) == 1 << i
        inverse.append(x)
    return inverse


def apply_inverse(value, basis):
    result = 0
    while value:
        bit = value & -value
        result ^= basis[bit.bit_length() - 1]
        value ^= bit
    return result


def support_scan(basis):
    result = {rate: {'rate_bits': rate, 'tested': 0, 'feasible': 0,
                     'minimum_capacity_bit_weight': 320} for rate in (64, 128)}
    for weight in (1, 2):
        for chosen in itertools.combinations(range(320), weight):
            delta = sum(1 << i for i in chosen)
            external = apply_inverse(delta, basis)
            for rate, row in result.items():
                weight_outside = (external >> rate).bit_count()
                row['tested'] += 1
                row['feasible'] += weight_outside == 0
                row['minimum_capacity_bit_weight'] = min(
                    row['minimum_capacity_bit_weight'], weight_outside)
    assert all(row['tested'] == 51360 for row in result.values())
    return list(result.values())


def old_endpoints(basis):
    source = ROOT / 'results/hnw26_comparison/confirmed_nonzero_endpoints.json'
    rows = []
    for row in json.loads(source.read_text())['rows']:
        if row['R'] == 4:
            rows.append((row['id'], row['selected'], str(source.relative_to(ROOT))))
    source2 = ROOT / 'results/datf_pools/linear_two_match/hnw26/results.json'
    data = json.loads(source2.read_text())
    rows.append(('L2_2to2_best_DATF', data['candidate'], str(source2.relative_to(ROOT))))
    records = []
    for ident, row, origin in rows:
        delta = tuple(int(x, 16) for x in row['input_difference'])
        mask = tuple(int(x, 16) for x in row['output_mask'])
        external = unpack(apply_inverse(pack(delta), basis))
        assert linear(external) == delta
        records.append({
            'id': ident, 'R': row['R'], 'source': origin,
            'core_input_bits': bits(delta), 'core_input_words': words(delta),
            'external_input_words': words(external),
            'external_input_weight_by_lane': [v.bit_count() for v in external],
            'output_bits': bits(mask), 'output_words': words(mask),
            'xof_input_feasible': not any(external[1:]),
            'aead_input_feasible': not any(external[2:]),
            'xof_output_observable': not any(mask[1:]),
            'aead_output_observable': not any(mask[2:]),
        })
    return records


def check_baseline():
    directory = ROOT / 'results/xof_interface_baseline'
    summary = json.loads((directory / 'summary.json').read_text())
    for group, base in (('source_sha256', ROOT), ('data_sha256', directory)):
        for name, digest in summary[group].items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == digest, name
    sampled = json.loads((directory / 'sampling.json').read_text())
    n = sampled['samples_per_endpoint']
    assert n == 1 << 24
    assert all(x['n0'] + x['n1'] == n for x in sampled['rows'])
    assert len({(x['rounds'], x['input_bit'], x['output_bit'])
                for x in sampled['rows']}) == 7168
    threshold = math.sqrt(2 * math.log(2 * 7168 / .05) / n)
    assert threshold == summary['multiple_testing']['absolute_rho_threshold']
    for result in summary['by_round']:
        values = [(x['n0'] - x['n1']) / n for x in sampled['rows']
                  if x['rounds'] == result['rounds']]
        assert len(values) == 3584
        assert max(map(abs, values)) == result['experimental_max_abs_rho']
        assert sum(abs(v) > threshold for v in values) == 0
    confirmation = json.loads((directory / 'confirmation.json').read_text())
    lookup = {(x['rounds'], x['input_bit'], x['output_bit']): x
              for x in confirmation['rows']}
    for row in summary['independent_confirmation']['rows']:
        raw = lookup[row['rounds'], row['input_bit'], row['output_bit']]
        assert raw['n0'] + raw['n1'] == confirmation['samples_per_endpoint']
        assert row['rho'] == (raw['n0'] - raw['n1']) / confirmation['samples_per_endpoint']
    return summary


def report(output):
    basis = inverse_basis()
    scan = support_scan(basis)
    endpoints = old_endpoints(basis)
    baseline = check_baseline()
    data = {
        'status': 'PASS',
        'definition': 'Gaston full rounds L -> RC -> chi; rate A0 or A0,A1',
        'input_relation': 'Delta_core = L(Delta_external)',
        'output_relation': 'Lambda_external = Lambda_core (no omitted last L)',
        'linear_rank': 320, 'inverse_basis_forward_checks': 320,
        'low_weight_core_input_scan': scan,
        'selected_existing_4_round_endpoints': endpoints,
        'existing_xof_baseline': baseline,
        'limits': [
            'Input impossibility is for identical states before rate injection; no capacity steering or state matching is assumed.',
            'The 51360 inputs cover all core bit weights 1 and 2, not all possible differential prefixes.',
            'Uniform-state DL correlations need not equal message-distribution or fixed-key mode correlations.',
            'No AEAD transcript experiment or key-recovery/forgery complexity is established.',
        ],
        'source_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ['tools/mode_interface_audit.py', 'tools/linear_seed_scan.py']},
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / 'audit.json').write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    esc = html.escape
    scan_html = ''.join(f'<tr><td>{r["rate_bits"]}</td><td>{r["tested"]:,}</td>'
                        f'<td>{r["feasible"]}</td><td>{r["minimum_capacity_bit_weight"]}</td></tr>' for r in scan)
    endpoint_html = ''.join(
        f'<tr><td>{esc(r["id"])}</td><td>{esc(str(r["core_input_bits"]))}</td>'
        f'<td>{esc(str(r["external_input_weight_by_lane"]))}</td><td>{esc(str(r["output_bits"]))}</td>'
        f'<td>{"可" if r["xof_input_feasible"] else "不可"} / {"可" if r["xof_output_observable"] else "不可"}</td>'
        f'<td>{"可" if r["aead_input_feasible"] else "不可"} / {"可" if r["aead_output_observable"] else "不可"}</td></tr>'
        for r in endpoints)
    experiment_html = ''.join(
        f'<tr><td>{r["rounds"]}</td><td>{r["endpoints"]:,}</td><td>{r["DATF_nonzero"]}</td>'
        f'<td>{r["experimental_max_abs_rho"]:.10f}</td><td>{r["above_joint_Hoeffding_threshold"]}</td></tr>'
        for r in baseline['by_round'])
    best = endpoints[-1]
    threshold = baseline['multiple_testing']['absolute_rho_threshold']
    page = f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston：四轮结果到 XOF / AEAD 的接口检查</title>
<style>body{{margin:0;background:#f3f6f5;color:#193c37;font:16px/1.8 system-ui,sans-serif}}main{{max-width:1160px;margin:auto;padding:30px 22px}}section{{background:white;padding:24px;margin:22px 0;border:1px solid #d5e3dd;border-radius:12px}}h1{{line-height:1.35}}h2{{font-size:22px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:10px;text-align:left;border-bottom:1px solid #dce8e1;vertical-align:top}}th{{background:#eaf3ef}}code,pre{{background:#edf4ef}}pre{{padding:16px;overflow:auto}}a{{color:#087161}}.callout{{border-left:5px solid #20806b;padding-left:16px}}.formula{{font-family:monospace;overflow:auto}}footer{{font-size:14px;color:#506b64}}</style></head><body><main>
<h1>Gaston：四轮结果到 XOF / AEAD 的接口检查</h1>
<p>2026-10-18 · 自定义 Gaston 模式 · XOF 末次吸收后的调用允许减为 4 轮</p>
<section><h2>结论</h2><p class="callout">允许末次吸收调用减轮，使四轮输出能直接出现在第一块 XOF 输出中。但已有差分池的入口处于首个 χ 之前。将首个线性层恢复后，所有重量为 1 或 2 的该处差分，都需要改变容量区。因此这批路径不能由同一前置状态下的消息差分直接注入。</p>
<p>已完成的 15 字节消息基线中，三轮与四轮各检查 3,584 个单比特端点，每项 2²⁴ 对。未观察到超过全体联合检验阈值的相关度，预选端点的独立复核也未确认非零。这个结果不排除多比特端点或更小相关度。</p></section>
<section><h2>1. 本次调用位置与完整轮边界</h2>
<p>XOF 速率为 A₀（64 比特），AEAD 速率为 A₀、A₁（128 比特）。初始化与 AEAD 结束调用保留 12 轮。XOF 实验只将最后一个带填充消息块之后的调用改为 4 轮；其他吸收和后续挤出调用保留 12 轮。</p>
<p class="formula">初始化 G¹² → 吸收公共前缀 → G¹² → 吸收末块与填充 → G⁴ → 输出 A₀</p>
<p>Gaston 一轮是 L → RC → χ。旧搜索核心省略第一个 L，并通过均匀状态平移吸收首个 RC；真实模式采样必须执行全部 L 与 RC。</p>
<p class="formula">Δχ = L(Δmsg)，亦即 Δmsg = L⁻¹(Δχ)。输出端没有省略最后一个 L，故 Λout = Λcore。</p>
<p>恢复 L 不增加 χ 层数，也不需要再平方相关度。对完整均匀状态，边界变换是严格等价；对消息产生的受限状态分布，相关度必须重新计算或测量。</p></section>
<section><h2>2. 全部低重量输入的有限枚举</h2>
<p>由正向 L 的 320 个基向量建立矩阵、完成 GF(2) 求逆，再对逆矩阵全部 320 个列向量作正向回代。枚举首个 χ 处重量 1、2 的全部非零差分：320 + C(320,2) = 51,360。</p>
<div class="scroll"><table><tr><th>速率位数</th><th>输入候选数</th><th>L⁻¹ 后仅位于速率区</th><th>容量区非零位数最小值</th></tr>{scan_html}</table></div>
<p>这是针对这些 51,360 个差分与给定行布局的严格有限枚举结论，覆盖原单轮差分池的全部输入端点。即使允许速率区任意多个活跃比特，也不能在保持容量不变时实现它们。结论不涉及其他差分、不同状态匹配、不同速率布局或其他消息接口。</p></section>
<section><h2>3. 已有四轮代表的具体障碍</h2>
<p>坐标为 [行, 位]，位从最低位 0 起。最后两列分别列“输入可注入 / 输出可直接观察”。仅检查这些代表，不将其称为全部四轮端点。</p>
<div class="scroll"><table><tr><th>候选</th><th>旧 Δχ</th><th>L⁻¹ 后各行重量</th><th>输出掩码位</th><th>XOF</th><th>AEAD</th></tr>{endpoint_html}</table></div>
<p>例如已有 L2_2to2 候选 Δχ={{(0,35)}}、Λ={{(1,52),(3,42)}}，置换核心实验曾测得 ρ̂≈0.004011。其实际注入差分共有 149 位，各行重量为 29、30、30、30、30：</p>
<pre>{esc(chr(10).join('A'+str(i)+' = '+w for i,w in enumerate(best['external_input_words'])))}</pre>
<p>输入容量区无法由一次消息吸收改变；输出掩码又包含 A₃，不能从 XOF 第一块或 AEAD 下一密文块直接恢复所需奇偶值。共同循环移位不会改变行支撑，因此不能消除这些接口障碍。</p></section>
<section><h2>4. 可实现消息端点的基线结果</h2>
<p>已保存的实验采用 15 字节消息：8 字节公共前缀与 7 字节末块。每对消息的公共前缀相同，末块仅一位不同；不同样本对重新生成前缀与末块。填充不变，所以差分位置为 0–55，输出掩码位置为 0–63，共 56×64 项。IV 取 0000080000cc0003，并用 Gaston 重新计算初始化。</p>
<p>下表 DATF 是完整均匀 320 比特状态上的普通 DATF；消息实验执行真实吸收流程。两者的分布不同，不能把 DATF 零值解释为消息函数严格零相关。</p>
<div class="scroll"><table><tr><th>末次调用轮数</th><th>端点数</th><th>DATF 非零项</th><th>实验最大 |ρ̂|</th><th>超过联合阈值</th></tr>{experiment_html}</table></div>
<p>每端点 2²⁴ 对样本，7,168 项的 95% Hoeffding 联合阈值为 {threshold:.10f}。统计解释采用端点内独立采样模型；仿真使用可复现的 SplitMix64。对筛选中各轮最大的两项另用新种子复核，也均未超过四项联合阈值。没有显著值不是严格零相关证明。</p>
<p><a href="../xof_interface_baseline/xof_interface_baseline.html">完整消息基线、独立复核与原始计数</a></p></section>
<section><h2>5. AEAD 应单独处理的条件</h2>
<p>若目标是消息处理中一个 G⁴ 调用，其后下一块已知明文与密文可以暴露该次调用输出的速率部分。结束阶段的 G¹² 不会撤回已经输出的密文。因此不能只因结束仍为 12 轮，就排除所有中间调用的统计观察。</p>
<p>但这仍要求输入差分能从速率区注入、输出掩码仅在速率区，并且必须明确如何在合法查询中获得配对的前置状态。跨两次执行复用同一密钥、nonce 和公共前缀可产生相同前置状态，却属于 nonce 复用条件；不能直接当作 nonce 不重复模型下的结果。不同 nonce 经 12 轮初始化后也不能假设容量差分为零。</p>
<p>本报告未执行 AEAD 消息实验，也未建立密钥恢复或标签伪造复杂度。XOF 的消息配对实验不受 AEAD 的 nonce 条件限制，应先完成该方向。</p></section>
<section><h2>6. 下一步搜索目标</h2>
<ol><li>将输入比特上限约束在 Δmsg，而不是 Δχ；计算 Δχ=L(Δmsg)。XOF 输出掩码限制于 A₀，AEAD 限制于 A₀、A₁。</li>
<li>先扩展当前 7 字节末块为至多两个活跃差分位，并允许第一输出字内至多两个掩码位。输入数为 C(56,1)+C(56,2)=1,596，输出数为 C(64,1)+C(64,2)=2,080，共 3,319,680 对端点；这是待搜索规模，并非已完成。</li>
<li>以模型筛选与消息实验共同排序，保留一部分模型零值作对照；候选确定后使用独立样本复核，并校正多重比较。</li>
<li>如需扩大到全部 64 位消息差分，必须重新安排消息长度和目标调用位置。当前 15 字节方案的高 8 位含填充，不能直接当作可控消息位。</li></ol>
<p>如后续确认相关度 ρ，简单独立样本统计量的样本量尺度通常为 O(ρ⁻²)，还需计算消息调用成本、显著性与功效；当前结果不足以给出模式区分复杂度。更不能将其直接当作碰撞、原像、密钥或标签的复杂度。</p></section>
<section><h2>数据与复现</h2><p><a href="audit.json">接口枚举、旧端点与消息基线汇总 JSON</a> · <a href="../../tools/mode_interface_audit.py">接口检查脚本</a> · <a href="../../tests/test_xof_interface.py">消息流程独立检查</a></p>
<pre>cd /home/yxb/分析/差分线性/gaston_dl
python3 tests/test_xof_interface.py
python3 tools/mode_interface_audit.py --output results/mode_interface_audit</pre>
<p>已核对原始数据及源文件 SHA-256、样本计数、阈值、独立复核值。消息实现另对照 Gaston 原始参考实现、填充边界与跨块输出，并由独立 Python 实现检查计数。</p>
<p>接口依据：<a href="https://ascon.isec.tugraz.at/specification.html">Ascon 官方模式说明</a>、<a href="https://github.com/ascon/ascon-c/blob/main/src/hash.c">官方 XOF 吸收与输出实现</a>。本文只借用接口组织方式，研究对象是用户自定义 Gaston 模式。</p></section>
<footer>严格结论：线性边界与有限支持枚举。模型结论：普通 DATF。实验观察：固定消息配置下的相关度计数。待验证：更大端点范围与完整模式的统计区分。</footer>
</main></body></html>'''
    (output / 'gaston_mode_interface.html').write_text(page)
    print(json.dumps({'status': 'PASS', 'scan': scan,
                      'report': str(output / 'gaston_mode_interface.html')}, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    report(parser.parse_args().output)

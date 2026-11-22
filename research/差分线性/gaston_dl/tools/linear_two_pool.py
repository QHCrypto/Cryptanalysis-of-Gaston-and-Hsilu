#!/usr/bin/env python3
"""Exhaust all chi-L-chi trails whose internal L masks have 2 active columns.

The 60,543-pattern scan includes arbitrary nonzero 5-bit column masks; all
surviving bridges have exactly two bits at each end. No cap is imposed on
the outer chi masks. Constants are recorded separately from bare LAT signs.
"""
import argparse
import gzip
import hashlib
import html
import itertools
import json
from collections import Counter
from pathlib import Path

from fixed_tail_search import RC, checked_round, dot, theta_transpose
from linear_seed_scan import (EAST, WEST, TWIN, active, column, embed, encode,
                              linear, rol, rotate, scan_two_column_seeds,
                              transpose, xor)

POOL_KEY = 'linear_two_2to2'
FAMILY = 'L2_2to2'


def words(a):
    return [f'{w:016x}' for w in a]


def bit_weight(a):
    return sum(w.bit_count() for w in a)


def completions(fixed, bound_output, lat):
    positions = [j for j in range(64) if column(fixed, j)]
    choices = [[v for v in range(1, 32)
                if (lat[v][column(fixed, j)] if bound_output else lat[column(fixed, j)][v])]
               for j in positions]
    for values in itertools.product(*choices):
        yield tuple(sum(((v >> i) & 1) << j for j, v in zip(positions, values))
                    for i in range(5))


def build_two_pool():
    from fixed_tail_search import S, lat as lat_entry
    lat = [[lat_entry(a, b) for b in range(32)] for a in range(32)]
    assert sorted(S) == list(range(32))  # zero/nonzero column support is preserved
    seeds, hits = scan_two_column_seeds()
    records, bridges = [], []
    for index, (b1, a2) in enumerate(seeds, 1):
        seed_id = f'T{index:02d}'
        d = tuple(rol(v, -s) for v, s in zip(a2, WEST))
        assert xor(d) == 0 and theta_transpose(d) == d
        assert tuple(rol(v, -s) for v, s in zip(d, EAST)) == b1
        assert bit_weight(b1) == bit_weight(a2) == active(b1) == active(a2) == 2
        # Independent forward-basis dot-product check, not just L^T formula.
        for row, bit in itertools.product(range(5), range(64)):
            assert dot(linear(embed(bit, 1 << row)), a2) == ((b1[row] >> bit) & 1)
        before, after = list(completions(b1, True, lat)), list(completions(a2, False, lat))
        bridges.append({'id': seed_id, 'first_output': encode(b1), 'second_input': encode(a2),
                        'rho_west_transpose': encode(d), 'parity': 0,
                        'first_chi_completions': len(before), 'second_chi_completions': len(after),
                        'trails': len(before)*len(after), 'orbit_size': 64})
        firsts = [(a1, checked_round(a1, b1)) for a1 in before]
        lasts = [(b2, checked_round(a2, b2)) for b2 in after]
        for (a1, r1), (b2, r2) in itertools.product(firsts, lasts):
            w2 = r1['w_squared'] + r2['w_squared']
            sign = r1['sign'] * r2['sign']
            # Tail includes leading L and each reference constant, with chi
            # rounds at indices s, s+1. Bitwise phases make rotation explicit.
            signs = [sign * (-1 if ((a1[0] & RC[s]).bit_count()
                                     + (a2[0] & RC[s+1]).bit_count()) & 1 else 1)
                     for s in range(11)]
            records.append({'rounds': 2, 'bridge_id': seed_id,
                            'input': words(a1), 'first_output': words(b1),
                            'second_input': words(a2), 'output': words(b2),
                            'join_mask': words(transpose(a1)), 'trail': [r1, r2],
                            'weight': w2//2, 'weight_squared': w2,
                            'sign': sign, 'sign_convention': 'bare chi LAT product; constants excluded',
                            'constant_phase_masks': [f'{a1[0]:016x}', f'{a2[0]:016x}'],
                            'sign_by_start_round': signs, 'rotation': 0, 'orbit_size': 64,
                            'input_bit_weight': bit_weight(a1), 'output_bit_weight': bit_weight(b2)})
    records.sort(key=lambda r: (r['weight'], r['bridge_id'], r['input'], r['output']))
    for i, r in enumerate(records):
        r['id'] = f'{FAMILY}_{i:05d}'
    audit = {'status': 'PASS', 'scope': 'All nonzero chi-L-chi trails with 2 active columns on both sides of the internal L, modulo common rotation.',
             'relation': 'b1 = L^T(a2); join_mask = L^T(a1)',
             'round_boundary': 'A tail attached after the middle is L-RC-chi-L-RC-chi; only the first L of the complete evaluated interval is omitted.',
             'bit_index': '(row,bit), zero-based, LSB=0', 'EAST': EAST, 'WEST': WEST, 'TWIN': TWIN,
             'scanned_column_zero_patterns': 63*31*31, 'hits': hits,
             'bridge_rotation_classes': len(seeds), 'expanded_bridges': len(seeds)*64,
             'representatives': len(records), 'expanded_records': len(records)*64,
             'weight_histogram': dict(sorted(Counter(r['weight'] for r in records).items())),
             'weight_squared_histogram': dict(sorted(Counter(r['weight_squared'] for r in records).items())),
             'bare_sign_histogram': dict(sorted(Counter(r['sign'] for r in records).items())),
             'input_bit_weights': sorted({r['input_bit_weight'] for r in records}),
             'output_bit_weights': sorted({r['output_bit_weight'] for r in records}),
             'forward_basis_checks': 320*len(seeds), 'bridges': bridges,
             'limits': ['2->2 counts active chi columns; the internal masks happen also to have bit weight two.',
                        'Outer masks are not restricted to two bits. The cap for single-round pools does not apply here.',
                        'Weights describe individual local LAT products, not a full linear hull or full DATF correlation.',
                        'No global linear branch-number lower bound is proved by this restricted enumeration.',
                        'Only the differential pool needs rotation expansion for relative matching; fixed-constant full correlations require separate validation.']}
    return records, audit


def persist_pool(directory, records, audit):
    path = directory / (POOL_KEY + '.jsonl.gz')
    # Stable gzip header and atomic replacement allow reproducible reruns.
    temp = path.with_suffix('.tmp')
    with temp.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as f:
        for record in records:
            f.write((json.dumps(record, separators=(',', ':'))+'\n').encode())
    temp.replace(path)
    audit_path = directory / (POOL_KEY + '_audit.json')
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2)+'\n')
    return {'file': path.name, 'representatives': len(records),
            'expanded_records': sum(r['orbit_size'] for r in records),
            'weight_histogram': audit['weight_histogram'],
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'rounds': 2, 'family': FAMILY, 'audit': audit_path.name,
            'scope': audit['scope'], 'outer_bit_cap': None,
            'canonicalization': 'Lexicographically smallest (b1,a2) internal bridge under common 64-bit rotation; all outer completions retained.'}


def render_report(directory, audit):
    def bits(x):
        return html.escape(', '.join(f'({r},{b})' for r,b in x['bits']))
    bridge_rows = ''.join(f'<tr><td>{r["id"]}</td><td>{bits(r["first_output"])}</td>'
                          f'<td>{bits(r["second_input"])}</td><td>{bits(r["rho_west_transpose"])}</td>'
                          f'<td>{r["trails"]:,}</td></tr>' for r in audit['bridges'])
    hist_rows = ''.join(f'<tr><td>{w}</td><td>2<sup>−{w}</sup></td><td>{2*int(w)}</td>'
                       f'<td>{n:,}</td><td>{64*n:,}</td></tr>' for w,n in audit['weight_histogram'].items())
    page = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston 两轮 2→2 线性路径池</title><style>
body{margin:0;background:#f1f5f4;color:#19352e;font:16px/1.75 system-ui,sans-serif}main{max-width:1100px;margin:auto;padding:36px 24px 72px}section{background:white;border:1px solid #d6e3de;border-radius:12px;padding:24px;margin-top:24px}h1{line-height:1.3}h2{font-size:23px}a{color:#007665}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:10px;border-bottom:1px solid #d6e3de;text-align:left}th{background:#eaf3ef}.scroll{overflow:auto}.formula{white-space:pre-wrap;padding:18px;background:#eaf3ef;border-radius:8px;font:15px/1.85 monospace}.note{color:#526b62}strong{color:#007665}</style></head><body><main>
<h1>Gaston 两轮 2→2 线性路径池</h1><p>已核对两侧 ρ 移位和 θ 转置的抵消关系，并枚举前后 χ 的全部非零 LAT 转移。</p>
<section><h2>结论与范围</h2><p><strong>10 个线性桥接类别 → 16,000 条两轮路径代表 → 1,024,000 条共同移位展开路径。</strong>每轮恰有两个活跃 χ 列；内部线性层两端也恰有两个活跃比特。外部输入、输出掩码可以包含更多比特。</p><p>路径记作 a₁ → χ → b₁ → L → a₂ → χ → b₂，满足 b₁ = Lᵀa₂。接在中间段之后时，还要保留尾段入口的 L，连接掩码为 γ = Lᵀa₁。仍然只省略完整区间最开始的一个 L。</p></section>
<section><h2>为什么两比特能够抵消</h2><p>L = ρ<sub>west</sub> ∘ θ ∘ ρ<sub>east</sub>。对掩码反向传播，先令 dᵢ = ROTL(a₂,ᵢ, −WESTᵢ)，p = ⊕ᵢdᵢ。</p>
<div class="formula">(θᵀd)ᵢ = dᵢ ⊕ ROTL(p,−23) ⊕ ROTL(p,−24)
              ⊕ ROTL(p,−23−TWINᵢ) ⊕ ROTL(p,−41−TWINᵢ)</div>
<p>当 d 的两个比特位于不同行的同一列时，p = 0，因此 θᵀd = d。最后按 −EASTᵢ 移位即可得到 b₁。这是转置层的抵消；不需要额外满足差分 θ-effect kernel 中的 Q = 0。</p>
<div class="formula">例：a₂ = {(0,0),(1,56)}
ρ_westᵀ：d = {(0,0),(1,0)} → p = 0 → θᵀd = d
ρ_eastᵀ：b₁ = {(0,0),(1,4)}</div><p>一般式：任选行 i &lt; j 及列 k，在 a₂ 的 (i,k+WESTᵢ)、(j,k+WESTⱼ) 放置比特；对应 b₁ 的位置为 (i,k−EASTᵢ)、(j,k−EASTⱼ)，列号均模 64。共有 C(5,2) = 10 个共同移位类别。</p></section>
<section><h2>全部 10 类内部桥接</h2><p>坐标均为零起始 (行,位)，位 0 为 LSB；每行是一个规范化代表。</p><div class="scroll"><table><thead><tr><th>类别</th><th>b₁：L 输入侧掩码</th><th>a₂：L 输出侧掩码</th><th>ρ<sub>west</sub>ᵀ(a₂)</th><th>完整路径代表</th></tr></thead><tbody>'''+bridge_rows+'''</tbody></table></div></section>
<section><h2>补齐 χ 后的路径数量与权重</h2><p>每个桥接类别：第一轮 χ 有 4² = 16 个输入补全，第二轮 χ 有 10² = 100 个输出补全，共 1,600 条。保留全部非零 LAT 转移，没有只保留最优局部权重，也没有对两轮池套用单轮池的外部比特上限。</p><div class="scroll"><table><thead><tr><th>w = −log₂|c|</th><th>|c|</th><th>−log₂c²</th><th>代表数</th><th>移位展开数</th></tr></thead><tbody>'''+hist_rows+'''</tbody></table></div><p class="note">c 是一条路径的局部 LAT 乘积。该表不是完整线性包或完整差分线性相关度。每条记录分别保存 LAT 符号、轮常数相位掩码及各合法起始轮号下的符号。</p></section>
<section><h2>完整性与使用</h2><p>扫描 63 × 31² = 60,543 个“第 0 列活跃”的双列输出掩码，包含任意非零五比特列取值，得到 20 个命中、10 个移位类别。这不是仅枚举所猜测的两比特模式。另以正向 L 的全部 320 个基向量逐一核对每个桥接的转置关系，共 3,200 项。</p><p>相对移位匹配仍只展开差分侧；线性池保留规范代表。固定轮常数下的完整端点相关度需要另行复评。本次完成的是路径池扩充，尚未对新增池执行大规模 DATF 匹配。</p><p><a href="pools.json">池清单</a> · <a href="linear_two_2to2.jsonl.gz">全部两轮路径 JSONL.gz</a> · <a href="linear_two_2to2_audit.json">枚举核验 JSON</a> · <a href="linear_two_2to2_validation.json">独立验证记录</a></p></section></main></body></html>'''
    validation_path = directory / (POOL_KEY + '_validation.json')
    if validation_path.exists():
        validation = json.loads(validation_path.read_text())
        page = page.replace('<h2>完整性与使用</h2>', '<h2>完整性与使用</h2>'
            f'<p><strong>独立验证：{html.escape(validation["status"])}</strong>。从正向 L 的矩阵构造转置，'
            f'穷举全部 {validation["independent_forward_matrix_two_column_patterns"]:,} 个双活跃列掩码，'
            f'恰好得到 {validation["expanded_bridges"]} 个桥接；'
            f'逐条检查 {validation["complete_trails_checked"]:,} 条路径和 '
            f'{validation["constant_start_signs_checked"]:,} 个起始轮号符号。</p>')
    else:
        page = page.replace('<a href="linear_two_2to2_validation.json">独立验证记录</a>',
                            '独立验证记录：运行 make test-linear-two-pool 后生成')
    (directory / (POOL_KEY + '.html')).write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pools', type=Path, required=True, help='Existing pool directory to extend')
    args = parser.parse_args()
    manifest_path = args.pools / 'pools.json'
    meta = json.loads(manifest_path.read_text())
    records, audit = build_two_pool()
    meta['pools'][POOL_KEY] = persist_pool(args.pools, records, audit)
    meta['linear_two_scope'] = audit['scope']
    temp = manifest_path.with_suffix('.tmp')
    temp.write_text(json.dumps(meta, ensure_ascii=False, indent=2)+'\n')
    temp.replace(manifest_path)
    render_report(args.pools, audit)
    print(json.dumps(meta['pools'][POOL_KEY], ensure_ascii=False))


if __name__ == '__main__':
    main()

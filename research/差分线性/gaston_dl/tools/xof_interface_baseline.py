#!/usr/bin/env python3
"""Uniform-state DATF screening and actual-message XOF sampling report.

The two distributions are deliberately recorded separately. Uniform input
translation permits the DATF adapter to absorb the first constant, whereas the
mode sampler executes every linear layer and round constant explicitly.
"""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from datf import BudgetExceeded
from datf_k06 import construct, value_record, verify_dag
from linear_seed_scan import linear

ROOT = Path(__file__).resolve().parents[1]


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def model(directory):
    target = directory / 'datf.json'
    if target.exists():
        raise FileExistsError(target)
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    began = time.monotonic()
    for rounds in (3, 4):
        for bit in range(56):
            delta = linear((1 << bit, 0, 0, 0, 0))
            engine, base, diff, stages = construct(delta, rounds, 12-rounds, 100000)
            if bit == 0:
                verify_dag(engine, base, diff, delta, rounds, 12-rounds, samples=8)
            for out in range(64):
                row = {'rounds': rounds, 'input_bit': bit, 'output_bit': out}
                # Count a separate allowance for each endpoint; completed
                # cached marginals from the same exact DAG may be reused.
                engine.max_bias_nodes = engine.stats['bias_nodes'] + 100000
                try:
                    rho = engine.estimate(diff[0][out])
                    row.update(status='completed', correlation=value_record(rho))
                except BudgetExceeded as exc:
                    row.update(status='resource_limit', correlation=None, reason=str(exc))
                records.append(row)
            if bit % 8 == 7:
                print('DATF', rounds, 'inputs', bit+1, '/ 56', flush=True)
    payload = {
        'method': 'ordinary DATF, unchanged engine; independence approximation',
        'distribution': 'uniform complete 320-bit state, NOT the 15-byte XOF message distribution',
        'mode_first_linear_included': True,
        'adapter': 'Delta_chi=L(Delta_external); first RC absorbed ONLY under uniform-state translation',
        'input_bits': list(range(56)), 'output_bits': list(range(64)),
        'why_56': '7-byte suffix within a 15-byte message; padding at bit 56 is unchanged',
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/datf.py', 'tools/datf_k06.py', 'tools/linear_seed_scan.py')},
        'zero_is_not_proof': True, 'elapsed_seconds': time.monotonic()-began, 'rows': records,
    }
    save(target, payload)
    print('DONE', len(records), 'DATF endpoints', flush=True)


def statistics(row, n):
    rho = (row['n0']-row['n1'])/n
    z = 1.959963984540054
    p = row['n0']/n
    den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    radius = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return {**row, 'rho': rho, 'log2_abs': math.log2(abs(rho)) if rho else None,
            'ci95': [max(-1, 2*(center-radius)-1), min(1, 2*(center+radius)-1)]}


def key(row):
    return row['rounds'], row['input_bit'], row['output_bit']


def report(directory):
    datf = json.loads((directory/'datf.json').read_text())
    sampled = json.loads((directory/'sampling.json').read_text())
    n = sampled['samples_per_endpoint']
    rows = [statistics(row, n) for row in sampled['rows']]
    assert len(rows) == len(datf['rows']) == 7168
    assert {key(row) for row in rows} == {key(row) for row in datf['rows']}
    top = [row for r in (3, 4) for row in sorted(
        (x for x in rows if x['rounds'] == r), key=lambda x: (-abs(x['rho']), key(x)))[:2]]
    # Union-bound Hoeffding threshold; no assumption of independence BETWEEN
    # endpoint tests. IID samples within a test are the statistical model.
    threshold = math.sqrt(2*math.log(2*len(rows)/0.05)/n)
    summary = {
        'profile': {k: v for k, v in sampled.items() if k != 'rows'},
        'interpretation': 'custom Gaston-based XOF experiment; not standardized Ascon',
        'uniform_DATF_not_mode_correlation': True,
        'sample_generation': 'counter-indexed SplitMix64; reproducible simulation, not an information-theoretic randomness proof',
        'multiple_testing': {'tests': len(rows), 'family_alpha': 0.05,
                             'method': 'two-sided Hoeffding union bound under IID within-endpoint sampling',
                             'absolute_rho_threshold': threshold,
                             'cross_endpoint_independence_required': False},
        'by_round': [], 'selected_for_independent_confirmation': top,
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('src/xof_interface_baseline.cpp', 'include/gaston.hpp',
                                       'tools/xof_interface_baseline.py', 'tests/test_xof_interface.py')},
        'data_sha256': {name: hashlib.sha256((directory/name).read_bytes()).hexdigest()
                        for name in ('datf.json', 'sampling.json')},
    }
    for r in (3, 4):
        modeled = [x for x in datf['rows'] if x['rounds'] == r]
        actual = [x for x in rows if x['rounds'] == r]
        summary['by_round'].append({
            'rounds': r, 'endpoints': len(actual),
            'DATF_completed': sum(x['status'] == 'completed' for x in modeled),
            'DATF_nonzero': sum(x['status'] == 'completed' and int(x['correlation']['numerator']) != 0 for x in modeled),
            'DATF_resource_limited': sum(x['status'] != 'completed' for x in modeled),
            'experimental_max_abs_rho': max(abs(x['rho']) for x in actual),
            'above_joint_Hoeffding_threshold': sum(abs(x['rho']) > threshold for x in actual),
        })
    confirmation_path = directory/'confirmation.json'
    if confirmation_path.exists():
        confirm = json.loads(confirmation_path.read_text())
        summary['data_sha256']['confirmation.json'] = hashlib.sha256(confirmation_path.read_bytes()).hexdigest()
        assert confirm['seed'] != sampled['seed']
        lookup = {key(x): x for x in confirm['rows']}
        count = confirm['samples_per_endpoint']
        bound = math.sqrt(2*math.log(2*len(top)/0.05)/count)
        summary['independent_confirmation'] = {
            'seed': confirm['seed'], 'samples_per_endpoint': count,
            'selected_tests': len(top), 'absolute_rho_threshold': bound,
            'rows': [{**statistics(lookup[key(x)], count),
                      'above_threshold': abs(statistics(lookup[key(x)], count)['rho']) > bound}
                     for x in top],
            'scope': 'only these four preselected endpoints; additional output counters are not selected after inspection',
        }
    save(directory/'summary.json', summary)
    save(directory/'confirmation_selection.json', {'selected': top,
         'input_bits_csv': ','.join(str(v) for v in sorted({x['input_bit'] for x in top}))})
    row_html = ''.join('<tr><td>{rounds}</td><td>{endpoints:,}</td><td>{DATF_nonzero}</td>'
                       '<td>{DATF_resource_limited}</td><td>{experimental_max_abs_rho:.8f}</td>'
                       '<td>{above_joint_Hoeffding_threshold}</td></tr>'.format(**x)
                       for x in summary['by_round'])
    candidate_html = ''.join(f'<tr><td>{x["rounds"]}</td><td>(0,{x["input_bit"]})</td>'
                            f'<td>(0,{x["output_bit"]})</td><td>{x["rho"]:+.9f}</td>'
                            f'<td>[{x["ci95"][0]:+.9f}, {x["ci95"][1]:+.9f}]</td></tr>' for x in top)
    confirmation_html = '<p>独立复核尚未运行；当前最大值来自筛选，不能直接作为最终偏差结论。</p>'
    if 'independent_confirmation' in summary:
        c = summary['independent_confirmation']
        confirmation_html = f'<p>独立种子 {c["seed"]}；每端点 {c["samples_per_endpoint"]:,} 对。仅复核预先选定的 4 个端点，联合阈值 {c["absolute_rho_threshold"]:.8f}。</p>'
        confirmation_html += '<table><tr><th>末次调用轮数</th><th>输入位</th><th>输出位</th><th>复核相关度</th><th>超过联合阈值</th></tr>'
        confirmation_html += ''.join(f'<tr><td>{x["rounds"]}</td><td>{x["input_bit"]}</td><td>{x["output_bit"]}</td><td>{x["rho"]:+.9f}</td><td>{"是" if x["above_threshold"] else "否"}</td></tr>' for x in c['rows']) + '</table>'
    page = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gaston-XOF：末次吸收调用减轮基线</title><style>body{{font:16px/1.8 system-ui,sans-serif;color:#173b35;background:#f1f5f4;margin:0}}main{{max-width:1100px;margin:auto;padding:32px 24px}}section{{background:white;border:1px solid #d2e0dc;border-radius:12px;padding:24px;margin:22px 0}}h1{{line-height:1.3}}h2{{font-size:22px}}table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{padding:9px;border-bottom:1px solid #dce7e3;text-align:left}}th{{background:#eaf3ef}}.scroll{{overflow-x:auto}}code{{background:#eaf3ef;padding:2px 5px}}a{{color:#087064}}</style><main>
<h1>Gaston-XOF：末次吸收调用减轮基线</h1><p>自定义研究构造；3 轮作对照，4 轮为目标。轮数指最后一次吸收调用中的完整 Gaston 轮数。</p>
<section><h2>精确定义</h2><p>状态为 5×64 比特，速率为 A₀。IV 沿用 Ascon-XOF128 的数值 <code>0000080000cc0003</code>，作为本次实验配置；从该 IV 使用 Gaston 重新计算初始化，未使用 Ascon 的预计算状态。</p>
<p>初始化 G¹² → 吸收 8 字节公共前缀 → G¹² → 吸收 7 字节后缀及填充 → G³ 或 G⁴ → 输出 A₀。后续挤出调用保留 G¹²，本实验仅观察第一块输出。</p>
<p>消息长度固定 15 字节；同一对消息的前 8 字节相同，后 7 字节恰有一个比特不同。填充字为 <code>0100000000000000</code>。每对的前缀与后缀重新生成。首个线性层和全部轮常数均保留。</p>
<p>差分位为 A₀ 的 0–55，输出掩码为 A₀ 的 0–63；每个轮数 56×64=3,584 个端点。未假定固定轮常数及填充下存在 64 循环移位等价。</p></section>
<section><h2>模型与消息实验</h2><p>DATF 对完整均匀状态使用 Δχ=L(Δext) 的入口变换；这是置换分布下的模型筛选。消息实验直接执行上述 XOF 流程，其受限状态分布不同，不能把两列当作同一相关度的两种精确计算。</p>
<p>每端点样本对数：{n:,}。对全部 7,168 项联合检验，95% Hoeffding 阈值为 |ρ̂| &gt; {threshold:.8f}。该统计阈值使用每个端点内样本独立的模型；端点之间可以共享消息。SplitMix64 用于可复现仿真，不提供信息论随机性证明。</p>
<div class="scroll"><table><tr><th>末次调用轮数</th><th>端点数</th><th>DATF 非零</th><th>DATF 未完成</th><th>实验最大 |ρ̂|</th><th>超过联合阈值</th></tr>{row_html}</table></div>
<p>DATF 零值不是严格零相关证明；没有超过实验阈值也不能排除更小的偏差。本次未搜索多比特端点、其他消息长度或其他前缀分布。</p></section>
<section><h2>筛选最大值与独立复核</h2><p>各轮数取 |ρ̂| 最大的两个端点。下列 95% 区间是逐项 Wilson 区间，未校正筛选偏差。</p><div class="scroll"><table><tr><th>末次调用轮数</th><th>输入差分位</th><th>输出掩码位</th><th>筛选相关度</th><th>逐项 95% 区间</th></tr>{candidate_html}</table></div>{confirmation_html}</section>
<section><h2>复现与数据</h2><p><a href="summary.json">结论、配置与文件校验值</a> · <a href="datf.json">普通 DATF 全部端点</a> · <a href="sampling.json">消息采样原始计数</a> · <a href="confirmation_selection.json">预选复核端点</a></p>
<p>边界验证包括参考 Gaston 对照、0/1/7/8/15/16/23 字节消息、跨块输出、优化采样与完整消息流程对照、首个线性层与常数转换。</p><p>接口依据：<a href="https://ascon.isec.tugraz.at/specification.html">Ascon 官方说明</a>、<a href="https://github.com/ascon/ascon-c/blob/main/src/hash.c">官方模式实现</a>。这不是标准 Ascon 的分析结果。</p></section></main></html>'''
    (directory/'xof_interface_baseline.html').write_text(page)
    print(json.dumps(summary['by_round'], ensure_ascii=False, indent=2))
    print('CONFIRM_INPUT_BITS', ','.join(str(v) for v in sorted({x['input_bit'] for x in top})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('model', 'report'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (model if args.action == 'model' else report)(args.output)

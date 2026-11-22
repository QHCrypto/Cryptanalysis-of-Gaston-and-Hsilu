#!/usr/bin/env python3
"""Fixed linear tails of at least three rounds; backward affine extension.

The one-chi connector is exact. DDT products and a single tail's squared LAT
product rank characteristics, not the full permutation correlation.
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import heapq
import html
import itertools
import json
import math
import random
import subprocess
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

from datf_pool_build import state, words
from datf_paper_reassess import evaluate_endpoint
from fixed_tail_search import S, RC, ExactPrefix, checked_round, dot, state_chi
from linear_seed_scan import EAST, TWIN, active, column, embed, encode, linear, rol, rotate, transpose
from datf_k06 import construct, decode, parity_polynomial, value_record
from datf_quadratic import QuadraticOrdinaryDATF
from datf_paper_audit import audit_components
from datf_fixed_tail import uniform_balance_certificate
from datf import ZERO as POLY_ZERO, add, indices, isolated, support
from six_six_prefix_search import ddt_table

DDT = ddt_table()
ZERO = (0,) * 5
# Trusted exact local arithmetic can produce fractions with >4300 digits.
# This changes serialization limits only, not any DATF transition or estimator.
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def atomic(path, data):
    temp = path.with_suffix(path.suffix + '.partial')
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    temp.replace(path)


def initialize(root):
    """Create a new reproducible search, or reuse verified immutable inputs."""
    root.mkdir(parents=True, exist_ok=True)
    pool = json.loads(Path('results/datf_pools/pools/pools.json').read_text())
    grouped = {}
    for tail in pool['paper_tails']:
        if tail['rounds'] < 3:
            continue
        key = tuple(tail['join_mask'])
        if key not in grouped:
            grouped[key] = {'id': f'G{len(grouped)+1:02d}', 'gamma': list(key), 'tails': []}
        grouped[key]['tails'].append(tail)
    groups = list(grouped.values())
    if (root/'groups.json').exists():
        assert json.loads((root/'groups.json').read_text()) == groups, 'linear pool changed'
    else:
        atomic(root/'groups.json', groups)
    pq = json.loads(Path('results/pq_kernel_six_all/patterns.json').read_text())
    pq_text = '\n'.join(' '.join(r['first_output']['words']) for r in pq['orbits'])+'\n'
    if (root/'pq_states.txt').exists():
        assert (root/'pq_states.txt').read_text() == pq_text, 'PQ seeds changed'
    else:
        (root/'pq_states.txt').write_text(pq_text)
    subprocess.run(['make', 'build/fixed_tail_extend'], check=True)
    pending = []
    for group in groups:
        gamma_path = root/(group['id']+'.txt')
        gamma_path.write_text(' '.join(group['gamma'])+'\n')
        result = root/(group['id']+'_search.json')
        if result.exists():
            old = json.loads(result.read_text())
            assert old['gamma'] == group['gamma']
            assert old['random_information_set_repetitions'] == 128
            assert old['keep_per_bucket'] == 12
        else:
            pending.append([str(Path('build/fixed_tail_extend').resolve()), str(gamma_path),
                            str(root/'pq_states.txt'), str(result), '128', '12'])
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool_:
        for result in pool_.map(lambda cmd: subprocess.run(cmd, check=True), pending):
            assert result.returncode == 0


def xor(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def bit_weight(a):
    return sum(w.bit_count() for w in a)


def coefficient(delta, gamma):
    return tuple((gamma[(k-1) % 5] & delta[(k+1) % 5])
                 ^ (gamma[(k-2) % 5] & delta[(k-1) % 5]) for k in range(5))


def connector(delta, gamma):
    if coefficient(delta, gamma) != ZERO:
        return 0
    return (-1) ** dot(gamma, xor(state_chi(delta), state_chi(ZERO)))


def local_weight(a, b):
    values = [DDT[column(a, j)][column(b, j)] for j in range(64)]
    assert all(values)
    return sum(5 - (v.bit_length()-1) for v in values)


def predecessor_variants(b, limit=64, seed=2026101702):
    """K-best local products plus deterministic diverse valid predecessors.

    No cap on the external input Hamming weight is imposed. Enumeration of
    predecessors is bounded and is not called complete.
    """
    positions = [j for j in range(64) if column(b, j)]
    choices = [
        sorted((5-(DDT[a][column(b, j)].bit_length()-1), a)
               for a in range(1, 32) if DDT[a][column(b, j)])
        for j in positions
    ]
    initial = (0,) * len(positions)

    def score(indices):
        return sum(row[k][0] for row, k in zip(choices, indices))

    queue = [(score(initial), initial)]
    seen = {initial}
    found = {}

    def keep(indices):
        a = ZERO
        for j, row, k in zip(positions, choices, indices):
            a = xor(a, embed(j, row[k][1]))
        w = score(indices)
        assert w == local_weight(a, b)
        found[a] = w

    while queue and len(found) < limit:
        _, indices = heapq.heappop(queue)
        keep(indices)
        for k, row in enumerate(choices):
            if indices[k]+1 >= len(row):
                continue
            nxt = indices[:k] + (indices[k]+1,) + indices[k+1:]
            if nxt not in seen:
                seen.add(nxt)
                heapq.heappush(queue, (score(nxt), nxt))
    rng = random.Random(seed ^ int(hashlib.sha256(str(b).encode()).hexdigest()[:16], 16))
    for _ in range(32):
        # Limit additional local weight, not input bit count.
        ids = tuple(rng.choice([k for k, v in enumerate(row) if v[0] <= row[0][0]+1])
                    for row in choices)
        keep(ids)
    return sorted(found.items(), key=lambda x: (x[1], x[0]))


def record(rd, delta, first, second, gamma, source, pq_inputs):
    if rd == 0:
        middle = delta
        weight = 0
    elif rd == 1:
        middle = linear(first)
        weight = local_weight(delta, first)
    else:
        a2 = linear(first)
        middle = linear(second)
        weight = local_weight(delta, first) + local_weight(a2, second)
    sign = connector(middle, gamma)
    assert sign in (-1, 1)
    return {
        'differential_rounds': rd, 'middle_rounds': 1,
        'input_difference': encode(delta), 'first_output': encode(first),
        'second_input': encode(linear(first)) if rd == 2 else None,
        'second_output': encode(second) if rd == 2 else None,
        'middle_input_difference': encode(middle),
        'join_mask': encode(gamma), 'differential_weight': weight,
        'middle_correlation': sign, 'source': source,
        'input_bit_weight': bit_weight(delta),
        'input_active_columns': active(delta),
        'outside_previous_differential_pool': (
            rd == 0 or (bit_weight(delta) > 2 if rd == 1 else delta not in pq_inputs)
        ),
    }


def fingerprint(record_):
    return hashlib.sha256(json.dumps(record_['input_difference']['words']).encode()).hexdigest()


def select_records(records, cap):
    by_delta = {}
    for r in records:
        key = tuple(r['input_difference']['words'])
        if key not in by_delta or r['differential_weight'] < by_delta[key]['differential_weight']:
            by_delta[key] = r
    order = sorted(by_delta.values(), key=lambda r: (r['differential_weight'], fingerprint(r)))
    buckets = {}
    for r in order:
        n = r['input_bit_weight']
        key = ('new' if r['outside_previous_differential_pool'] else 'old',
               0 if n <= 2 else 1 if n <= 4 else 2 if n <= 8 else 3 if n <= 16 else 4)
        buckets.setdefault(key, []).append(r)
    # Each support range gets some coverage; the rest follow path weight.
    selected, keys = [], set()
    for key in sorted(buckets):
        for r in buckets[key][:max(2, cap//max(1, 2*len(buckets)))]:
            selected.append(r)
            keys.add(tuple(r['input_difference']['words']))
    for r in order:
        if len(selected) >= cap:
            break
        key = tuple(r['input_difference']['words'])
        if key not in keys:
            selected.append(r)
            keys.add(key)
    selected.sort(key=lambda r: (r['differential_weight'], fingerprint(r)))
    return selected[:cap], len(by_delta), len(records)


def group_records(group, search, pq_inputs, cap):
    gamma = state(group['gamma'])
    direct = []
    for j in range(64):
        for a in range(1, 32):
            delta = embed(j, a)
            if connector(delta, gamma):
                direct.append(record(0, delta, ZERO, ZERO, gamma, 'all_single_column_connectors', pq_inputs))
    rng = random.Random(2026101703)
    singles = list(direct)
    for _ in range(1024):
        a, b = rng.sample(singles, 2)
        delta = xor(state(a['input_difference']['words']), state(b['input_difference']['words']))
        if delta != ZERO:
            direct.append(record(0, delta, ZERO, ZERO, gamma, 'sampled_two_column_connectors', pq_inputs))
    result, stats = {}, {}
    result[0], unique, count = select_records(direct, min(cap, 32))
    stats[0] = {'generated_characteristics': count, 'unique_input_differences': unique,
                'selected_input_differences': len(result[0])}
    for rd, key in ((1, 'one_round_seeds'), (2, 'two_round_seeds')):
        records = []
        for seed in search[key]:
            first, second = state(seed['first_output']), state(seed['second_output'])
            for delta, _ in predecessor_variants(first):
                records.append(record(rd, delta, first, second, gamma, seed['source'], pq_inputs))
        result[rd], unique, count = select_records(records, cap)
        stats[rd] = {'generated_characteristics': count, 'unique_input_differences': unique,
                     'selected_input_differences': len(result[rd]),
                     'selected_outside_previous_pool': sum(r['outside_previous_differential_pool'] for r in result[rd])}
    return result, stats


def job_id(delta, mask, rounds, start):
    return hashlib.sha256(json.dumps([delta['words'], mask['words'], rounds, start]).encode()).hexdigest()


def prepare(root, cap):
    groups = json.loads((root/'groups.json').read_text())
    pq = json.loads(Path('results/pq_kernel_six_all/patterns.json').read_text())
    pq_inputs = {rotate(state(p['first_input']['words']), s) for p in pq['orbits'] for s in range(64)}
    jobs, rows, audits = {}, [], []
    for group in groups:
        search = json.loads((root/(group['id']+'_search.json')).read_text())
        records, stats = group_records(group, search, pq_inputs, cap)
        audits.append({'group_id': group['id'], 'generation': stats,
                       **{k: v for k, v in search.items() if not k.endswith('_seeds')}})
        for tail in group['tails']:
            assert tail['rounds'] >= 3
            assert transpose(state(tail['input'])) == state(group['gamma'])
            assert sum(checked_round(tuple(r['input']), tuple(r['output']))['w_squared']
                       for r in tail['trail']) == tail['weight_squared']
            for rd in range(3):
                total = rd+1+tail['rounds']
                start = 12-total
                row = {'id': f'D{rd}_M1_{tail["id"]}', 'group_id': group['id'],
                       'tail_id': tail['id'], 'r_d': rd, 'r_m': 1, 'r_l': tail['rounds'],
                       'rounds': total, 'start_round_index': start,
                       'tail_squared_weight': tail['weight_squared'],
                       'output_mask': encode(state(tail['output'])), 'candidates': []}
                for rank, r in enumerate(records[rd], 1):
                    r = {**r, 'rank': rank,
                         'single_characteristic_weight': r['differential_weight']+tail['weight_squared']}
                    eid = job_id(r['input_difference'], row['output_mask'], total, start)
                    r['endpoint_id'] = eid
                    row['candidates'].append(r)
                    jobs[eid] = {'id': eid, 'rounds': total, 'start_round_index': start,
                                 'input_difference': r['input_difference'], 'output_mask': row['output_mask'],
                                 'kind': 'full_endpoint', 'first_L_omitted': True}
                rows.append(row)
        print('PREPARED', group['id'], stats, flush=True)
    atomic(root/'candidates.json', {'rows': rows, 'group_audits': audits, 'max_per_row': cap,
                                  'scope': 'Fixed unrotated tails; absolute differential positions searched.',
                                  'search_is_not_globally_exhaustive': True})
    with gzip.open(root/'endpoint_jobs.jsonl.gz', 'wt') as f:
        for job in jobs.values():
            f.write(json.dumps(job)+'\n')
    return rows, list(jobs.values())


def evaluate_batch(root, jobs, workers, nodes, label):
    folder = root/label
    folder.mkdir(exist_ok=True)
    backend_paths = [Path('tools')/name for name in (
        'datf.py', 'datf_quadratic.py', 'datf_k06.py',
        'datf_paper_reassess.py', 'datf_fixed_tail.py', 'datf_audit.py')]
    backend = {'backend_sha256': {str(f): hashlib.sha256(f.read_bytes()).hexdigest()
                                 for f in backend_paths},
               'dag_samples': 4, 'max_bias_nodes': nodes}
    if (folder/'configuration.json').exists():
        assert json.loads((folder/'configuration.json').read_text()) == backend, 'DATF backend changed'
    else:
        atomic(folder/'configuration.json', backend)
    known, todo = {}, []
    for job in jobs:
        path = folder/(job['id']+'.summary.json')
        if path.exists():
            known[job['id']] = json.loads(path.read_text())
        else:
            todo.append(job)
    began = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(evaluate_endpoint, job, str(folder), nodes, 4): job['id'] for job in todo}
        for index, f in enumerate(concurrent.futures.as_completed(futures), 1):
            result = f.result()
            known[result['id']] = result
            if index % 16 == 0 or index == len(todo):
                doc = {'stage': label, 'completed': len(known), 'total': len(jobs),
                       'elapsed_seconds': time.monotonic()-began,
                       'model_nonzero': sum(r['status']=='completed' and int(r['correlation']['numerator']) != 0
                                            for r in known.values()),
                       'resource_limited': sum(r['status']=='resource_limit' for r in known.values())}
                atomic(root/'progress.json', doc)
                print('EVALUATE', json.dumps(doc), flush=True)
    return known


def make_prefix_jobs(rows, known):
    jobs = {}
    for row in rows:
        chosen = list(row['candidates'][:2])
        chosen += [r for r in row['candidates']
                   if known[r['endpoint_id']]['status']=='completed'
                   and int(known[r['endpoint_id']]['correlation']['numerator'])]
        # Also check the first new endpoint, if the row starts with an old control.
        chosen += next(([r] for r in row['candidates'] if r['outside_previous_differential_pool']), [])
        for r in chosen:
            rounds = row['r_d']+1
            eid = job_id(r['input_difference'], r['join_mask'], rounds, row['start_round_index'])
            r['prefix_job_id'] = eid
            jobs[eid] = {'id': eid, 'rounds': rounds, 'start_round_index': row['start_round_index'],
                         'input_difference': r['input_difference'], 'output_mask': r['join_mask'],
                         'kind': 'prefix_endpoint', 'first_L_omitted': True}
    return list(jobs.values())


def exact_prefix_audit(root, rows, branch_cap=1048576):
    """Exact two-chi conditional sum; distinct from ordinary DATF."""
    cache, results, seen = {}, [], set()
    for row in rows:
        if row['r_d'] != 1:
            continue
        for c in row['candidates']:
            if 'prefix_job_id' not in c or c['prefix_job_id'] in seen:
                continue
            seen.add(c['prefix_job_id'])
            delta = state(c['input_difference']['words'])
            gamma = state(c['join_mask']['words'])
            count = math.prod(sum(v > 0 for v in DDT[column(delta, j)]) for j in range(64))
            result = {'id': c['prefix_job_id'], 'branch_count': count}
            if count > branch_cap:
                result['status'] = 'not_enumerated'
            else:
                if gamma not in cache:
                    cache[gamma] = ExactPrefix(gamma)
                data = cache[gamma].evaluate(delta, RC[row['start_round_index']+1], details=True)
                exact = Fraction(data['numerator'], data['denominator'])
                parts = data['branches']
                result.update(status='exact_completed', correlation=value_record(exact),
                              positive_contributions=sum(x['exact_signed_count'] > 0 for x in parts),
                              negative_contributions=sum(x['exact_signed_count'] < 0 for x in parts),
                              nonzero_conditional_branches=len(parts), branches=parts)
            results.append(result)
    atomic(root/'exact_prefix_audit.json', {
        'method': 'Exact conditional two-chi sum including base-value dependence',
        'branch_cap': branch_cap, 'results': results})
    return {r['id']: r for r in results}


def theta_effect(b):
    east = tuple(rol(w, s) for w, s in zip(b, EAST))
    p, q = 0, 0
    for i in range(5):
        p ^= east[i]
        q ^= rol(east[i], TWIN[i])
    effect = rol(p ^ rol(p, 1) ^ q ^ rol(q, 18), 23)
    return {'P': f'{p:016x}', 'Q': f'{q:016x}', 'E': f'{effect:016x}'}


def audit_limited_factors(root, known):
    path = root/'supplemental_factor_audits.json'
    existing = json.loads(path.read_text()) if path.exists() else []
    ids = {r['id'] for r in existing}
    for summary in known.values():
        if summary.get('factor_audit_status') == 'PASS' or summary['id'] in ids:
            continue
        if summary['status'] != 'completed':
            continue
        with gzip.open(summary['artifact'], 'rt') as f:
            saved = json.load(f)['result']
        engine, _, diff, _ = construct(
            decode(saved['input_difference']), saved['rounds'], saved['start_round_index'],
            summary['max_bias_nodes'], engine_class=QuadraticOrdinaryDATF)
        target = parity_polynomial(diff, decode(saved['output_mask']))
        trace = []
        rho = engine.estimate(target, trace=trace)
        assert value_record(rho) == saved['correlation']
        audit = None
        # A certified zero factor suffices to check a zero model result,
        # without re-enumerating all earlier nonzero factors.
        polynomial = target
        if not rho:
            for item in trace:
                lone = isolated(polynomial)
                rest = polynomial-frozenset(1 << v for v in indices(lone))
                if int(item['factor_numerator']) == 0:
                    biased = {v for v in indices(support(rest)) if engine.marginal(v)}
                    certificate = uniform_balance_certificate(rest, biased)
                    if certificate:
                        audit = {
                            'status': 'PASS', 'correlation': value_record(Fraction(0)),
                            'method': 'exact balancing-shift certificate for a zero DATF factor',
                            'step': item['step'], 'certificate': certificate,
                            'scope': 'Checks the zero product in the independent DATF model; does not re-enumerate earlier nonzero factors or prove real-permutation zero.'}
                        break
                if any(engine.definitions[v] is None for v in indices(lone)):
                    break
                polynomial = add(*(engine.definitions[v] for v in indices(lone))) if lone else POLY_ZERO
        if audit is None:
            audit = audit_components(engine, target, trace)
        assert audit['correlation'] == saved['correlation']
        existing.append({'id': summary['id'], 'audit': audit})
        atomic(path, existing)
    return {r['id']: r['audit'] for r in existing}


def identify_previous_endpoints(rows):
    """Avoid presenting a different decomposition as a new external endpoint."""
    def key(rounds, start, delta, mask):
        return (rounds, start, state(delta['words']), state(mask['words']))
    selected = {}
    for row in rows:
        for c in row['candidates']:
            if c['full_DATF']['status']=='completed' and fraction(c['full_DATF']['correlation']):
                selected.setdefault(key(row['rounds'], row['start_round_index'],
                                        c['input_difference'], row['output_mask']), []).append(c)
    for file in (Path('results/datf_pools/paper_full24/endpoint_jobs.jsonl.gz'),
                 Path('results/datf_pools/linear_two_match/endpoint_jobs.jsonl.gz')):
        if not file.exists():
            continue
        with gzip.open(file, 'rt') as f:
            for line in f:
                job = json.loads(line)
                k = key(job['rounds'], job['start_round_index'], job['input_difference'], job['output_mask'])
                for c in selected.get(k, []):
                    c.setdefault('previous_endpoint_occurrences', []).append(
                        {'manifest': str(file), 'previous_id': job['id']})


def fraction(record_):
    return Fraction(int(record_['numerator']), int(record_['denominator']))


def fmt_value(record_):
    if record_ is None:
        return '未完成'
    v = fraction(record_)
    if not v:
        return '0（模型值）'
    if len(record_['numerator'])+len(record_['denominator']) > 70:
        sign = '+' if v > 0 else '−'
        return f'{sign}2^({record_["log2_abs"]:.6f})（精确分数见 JSON）'
    return f'{v.numerator}/{v.denominator}；log₂|ρ|={record_["log2_abs"]:.4f}'


def render(root, doc):
    esc = html.escape
    def bits(x):
        return ', '.join(f'({r},{j})' for r, j in x['bits'])
    group_table = ''
    for g in doc['group_audits']:
        counts = g['one_round_kernel_counts_by_bit_weight'][1:]
        group_table += '<tr>'+''.join(f'<td>{x}</td>' for x in (
            g['group_id'], g['connector_constraint_rank'], g['connector_preimage_dimension'],
            ' / '.join(f'{v:,}' for v in counts), g['best_seen_one_round_weight'],
            f'{g["two_round_feasible_first_outputs"]:,} / {g["two_round_first_outputs_checked"]:,}',
            g['two_round_best_seen_weight']))+'</tr>'
    rows_html, details, exact_disagreements = '', '', ''
    for row in doc['rows']:
        best = row['best_nonzero_candidate']
        value = fmt_value(best['full_DATF']['correlation']) if best else (
            '已评估端点均为模型零' if row['completed']==len(row['candidates']) else '尚有未完成项')
        rows_html += '<tr>'+''.join(f'<td>{x}</td>' for x in (
            f'<a href="#{row["id"]}">{esc(row["tail_id"])}</a>',
            f'{row["r_d"]}+1+{row["r_l"]}', row['rounds'],
            row['tail_squared_weight'], row['best_single_characteristic_weight'],
            f'{row["outside_previous_pool"]}/{len(row["candidates"])}',
            f'{row["model_nonzero"]}/{row["completed"]}', esc(value)))+'</tr>'
        candidates = ''
        for c in row['candidates']:
            prefix = c.get('prefix_DATF')
            pv = esc(fmt_value(prefix['correlation'])) if prefix else '未单独计算'
            exact = c.get('exact_prefix')
            ev = (esc(fmt_value(exact['correlation'])).replace('（模型值）', '（精确值）')
                  if exact and exact['status']=='exact_completed' else '未进行精确求和')
            fv = esc(fmt_value(c['full_DATF']['correlation']))
            if exact and exact.get('agrees_with_DATF') is False:
                exact_disagreements += '<tr>'+''.join(f'<td>{x}</td>' for x in (
                    f'<a href="#{row["id"]}">{row["tail_id"]} #{c["rank"]}</a>',
                    row['rounds'], c['input_bit_weight'], ev, pv, fv,
                    f'2<sup>−{c["single_characteristic_weight"]}</sup>'))+'</tr>'
            artifact = Path(c['full_DATF']['artifact']).relative_to(root.resolve()) if c['full_DATF'].get('artifact') else None
            link = f'<a href="{artifact}">DAG / 复评记录</a>' if artifact else ''
            sign = '+' if c['middle_correlation'] > 0 else '−'
            history = ('该完整端点已在旧报告出现，本次提供新的分段结构。'
                       if c.get('previous_endpoint_occurrences') else '')
            candidates += f'''<details class="candidate"><summary>#{c['rank']} · 输入 {c['input_bit_weight']} 比特 / {c['input_active_columns']} 列 · 单路径 {sign}2<sup>−{c['single_characteristic_weight']}</sup> · 完整 DATF：{fv}</summary>
<p>旧差分池之外：{'是' if c['outside_previous_differential_pool'] else '否（对照）'}；候选来源：{esc(c['source'])}。</p>
<p>{history}</p>
<p>输入差分：<code>{esc(bits(c['input_difference']))}</code></p>
<p>输出掩码：<code>{esc(bits(row['output_mask']))}</code></p>
<p>前缀 DATF：{pv}。前缀精确求和：{ev}。</p>
<p>中间单 χ 的精确相关度：{c['middle_correlation']}；差分路径权重：{c['differential_weight']}。两轮差分段的权重表示局部 DDT 乘积，不宣称已验证该多轮差分事件的真实概率。</p>
<p>第一 χ 输出差分：<code>{esc(bits(c['first_output']))}</code></p>
<p>连接 χ 输入差分：<code>{esc(bits(c['middle_input_difference']))}</code></p>
<p>{link}</p></details>'''
        details += f'''<section id="{row['id']}"><h2>{esc(row['tail_id'])}：{row['r_d']}+1+{row['r_l']}</h2>
<p>总轮数 {row['rounds']}，起始轮号 {row['start_round_index']}（零起始），采用参考实现的末 {row['rounds']} 轮。尾段入口的 L 保留。</p>{candidates}</section>'''
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Gaston 固定三轮以上线性尾段：向前扩展</title><style>
body{{margin:0;background:#f0f4f5;color:#17343c;font:16px/1.75 system-ui,sans-serif}}main{{max-width:1280px;margin:auto;padding:36px 24px 80px}}
h1{{line-height:1.3;font-size:32px}}h2{{font-size:23px}}section{{background:white;border:1px solid #cedde0;border-radius:12px;padding:24px;margin:22px 0}}
.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:10px;border-bottom:1px solid #dce5e7;text-align:left;vertical-align:top}}
th{{background:#e8f0f2}}a{{color:#006d83}}code{{overflow-wrap:anywhere}}.candidate{{border:1px solid #d4e0e3;border-radius:8px;padding:12px;margin:10px 0}}
summary{{cursor:pointer;font-weight:600}}.note{{background:#fff5df;padding:16px;border-left:4px solid #bc882c}}.metrics{{display:flex;gap:20px;flex-wrap:wrap}}.metrics div{{background:#eaf3f5;padding:16px;border-radius:8px}}
</style></head><body><main><h1>Gaston：固定三轮及以上线性尾段，向前扩展</h1>
<p>7 条固定尾段，3 种连接掩码；扩展为 0+1+L、1+1+L、2+1+L。只省略完整区间的第一个线性层；其余线性层和参考轮常数保留。</p>
<section><h2>本次结果</h2><div class="metrics"><div>完整端点复评<br><strong>{doc['unique_full_jobs']}</strong></div><div>DATF 模型非零端点<br><strong>{doc['unique_model_nonzero']}</strong></div><div>资源限制未完成<br><strong>{doc['resource_limited']}</strong></div></div>
<p>两个模型非零结果均为四轮的直接连接结构。其中 −9/2⁴⁴ 的端点已在旧报告出现，本次是新的分段表示；另一端点的模型相关度约为 +2<sup>−91.9993</sup>。本次没有得到完整五轮及以上的 DATF 模型非零端点。</p>
<p class="note">局部路径合法、拼接贡献非零、完整 DATF 模型非零、真实置换相关度非零是不同层次。本报告不把局部乘积当作整体相关度，也不把 DATF 的零当作真实相关度为零的证明。没有将完整 DATF 估计标为精确值。</p></section>
<section><h2>关键发现：新前缀精确非零，普通 DATF 却为零</h2>
<p>对选定的两轮前缀进行包含状态条件依赖的精确分支求和。D3_chi1_to_3/4 的新八比特候选出现精确值与普通 DATF 的分歧。因此它们仍是值得保留的五轮、六轮候选，完整相关度尚未确认。</p>
<div class="scroll"><table><thead><tr><th>固定尾段与候选</th><th>总轮数</th><th>输入比特数</th><th>两轮前缀精确 ρ</th><th>前缀 DATF</th><th>完整 DATF</th><th>局部路径乘积绝对值</th></tr></thead><tbody>{exact_disagreements}</tbody></table></div>
<p>排名第一的新前缀首轮差分输入等于输出，有八个单比特活跃列，概率 2<sup>−16</sup>。随后经过 L 仍为八比特，θ 的总作用 E=0，但 P、Q 都非零；它不属于此前限制的 P=Q=0 六比特族。连接 χ 的相关度为 +1；其余首轮差分分支的条件贡献为零，所以完整两轮前缀的精确相关度为 2<sup>−16</sup>。</p>
<p><a href="exact_prefix_audit.json">精确前缀分支记录</a>。精确结论只覆盖两轮前缀；乘上线性路径平方仍是完整区间的启发式评分。</p></section>
<section><h2>搜索方法与覆盖范围</h2>
<p>令 γ=Lᵀu 为固定线性尾段在连接边界的掩码。χ 是二次映射，其导数分量 γ·D<sub>α</sub>χ(x) 是仿射函数；非零相关度等价于全部 x 系数为零，此时相关度严格为 +1 或 −1。</p>
<p>一轮差分前缀：对第一 χ 输出差分 b 搜索 B<sub>γ</sub>Lb=0。完整枚举所有重量 1–4 的 b，用 meet-in-the-middle 消除重复；另外运行 128 次固定种子的随机信息集消元和局部改进。差分输入通过 DDT 反向补全，不限制为旧池的 1–2 比特。</p>
<p>两轮差分前缀：给定第一 χ 输出 b₁，令 a₂=Lb₁；精确求解 χ 在差分 a₂ 下的输出仿射空间与 ker(B<sub>γ</sub>L) 的交。第一输出候选覆盖全部重量 1–2 的状态、全部单活跃列非零模式、11 类 PQ 六比特模式的 64 个移位，以及保留的零空间种子和 64 个稠密对照。每个可行 b₁ 保留一个 b₂ 见证，不枚举全部第二输出。</p>
<p>零空间及交空间的可行性判定精确；全局最小路径权重没有证明。保留种子和反向 DDT 补全有明确数量界限，按路径权重及输入重量分层选取端点复评。普通 DATF 使用原有算法，仅以精确二次 Bias 后端加速，没有采用增强 DATF。</p>
<div class="scroll"><table><thead><tr><th>入口组</th><th>约束秩</th><th>零空间维数</th><th>b 重量 1 / 2 / 3 / 4 的解数</th><th>一轮前缀最佳已见 wD</th><th>两轮可行 b₁ / 检查数</th><th>两轮最佳已见 wD</th></tr></thead><tbody>{group_table}</tbody></table></div>
<p>G01 = D2_repeat_3…6；G02 = D3_chi1_to_3/4；G03 = D3_chi2_to_4。表中最佳已见值属于生成阶段；下表的候选经过端点多样性筛选。</p></section>
<section><h2>固定尾段汇总</h2><p>路径权重 W=wD+2wL，因为连接单 χ 的绝对相关度为 1。0+1+L 表示直接连接对照，前面没有独立差分段。“旧池外”按差分段的外部输入差分判定。</p>
<div class="scroll"><table><thead><tr><th>固定尾段</th><th>分配 D+M+L</th><th>总轮数</th><th>尾段平方权重</th><th>最佳保留路径 W</th><th>旧池外 / 候选</th><th>DATF 非零 / 完成</th><th>最好完整 DATF</th></tr></thead><tbody>{rows_html}</tbody></table></div></section>
{details}<section><h2>数据与复现</h2>
<p><a href="results.json">结果 JSON</a> · <a href="candidates.json">全部保留候选与边界</a> · <a href="validation.json">独立核验</a> · <a href="configuration.json">配置与代码哈希</a> · <a href="G01_search.json">G01 搜索</a> · <a href="G02_search.json">G02 搜索</a> · <a href="G03_search.json">G03 搜索</a> · <a href="exact_prefix_audit.json">精确前缀核对</a></p>
<p>首次独立因子算术审计状态：{esc(str(doc['factor_audit_counts']))}。其中 {doc['supplemental_factor_audits']} 项已完成<a href="supplemental_factor_audits.json">补充核对</a>：采用变量分量的精确加权求和，或给出使某个因子严格为零的平衡移位证书。后一种方法验证整个模型乘积为零，不要求重新穷举之前的非零因子。原始受限记录保留；这些核对不验证独立性假设。</p>
<p>0+1+L 与 1+1+L、2+1+L 是不同扩展结构。只固定池中的尾段规范代表，没有把共同循环移位视为固定常数下真实置换的对称性。未进行对无法分辨微小相关度的无效抽样。</p></section></main></body></html>'''
    (root/'fixed_tail_extension.html').write_text(page)


def finalize(root, rows, known, prefixes):
    exacts = exact_prefix_audit(root, rows)
    supplements = audit_limited_factors(root, known)
    for row in rows:
        for r in row['candidates']:
            r['full_DATF'] = known[r['endpoint_id']]
            if r['endpoint_id'] in supplements:
                r['supplemental_factor_audit_status'] = supplements[r['endpoint_id']]['status']
            if 'prefix_job_id' in r:
                r['prefix_DATF'] = prefixes[r['prefix_job_id']]
                exact = exacts.get(r['prefix_job_id'])
                if exact:
                    exact = {k: v for k, v in exact.items() if k != 'branches'}
                    if exact['status']=='exact_completed':
                        exact['agrees_with_DATF'] = (fraction(exact['correlation']) ==
                                                    fraction(r['prefix_DATF']['correlation']))
                    r['exact_prefix'] = exact
            if row['r_d']:
                r['theta_effect_after_first_chi'] = theta_effect(state(r['first_output']['words']))
        done = [r for r in row['candidates'] if r['full_DATF']['status']=='completed']
        nonzero = [r for r in done if fraction(r['full_DATF']['correlation'])]
        row.update(completed=len(done), model_nonzero=len(nonzero),
                   outside_previous_pool=sum(r['outside_previous_differential_pool'] for r in row['candidates']),
                   best_single_characteristic_weight=min(r['single_characteristic_weight'] for r in row['candidates']),
                   best_nonzero_candidate=max(nonzero, key=lambda r: abs(fraction(r['full_DATF']['correlation'])), default=None))
    identify_previous_endpoints(rows)
    candidate_doc = json.loads((root/'candidates.json').read_text())
    doc = {'rows': rows, 'group_audits': candidate_doc['group_audits'],
           'unique_full_jobs': len(known), 'unique_prefix_jobs': len(prefixes),
           'unique_model_nonzero': sum(r['status']=='completed' and fraction(r['correlation']) != 0 for r in known.values()),
           'resource_limited': sum(r['status']!='completed' for r in known.values()),
           'factor_audit_counts': dict(Counter(r.get('factor_audit_status') for r in known.values())),
           'supplemental_factor_audits': len(supplements),
           'first_L_omitted': True, 'enhanced_DATF': False,
           'correlation_scope': 'ordinary DATF model; local DDT/DLCT arithmetic independently checked'}
    atomic(root/'results.json', doc)
    render(root, doc)
    print('REPORT', root/'fixed_tail_extension.html', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=Path('results/fixed_tail_extension'))
    p.add_argument('--max-per-row', type=int, default=48)
    p.add_argument('--workers', type=int, default=12)
    p.add_argument('--max-bias-nodes', type=int, default=1000000)
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--reuse-candidates', action='store_true',
                   help='reuse the saved, independently checked endpoint selection')
    args = p.parse_args()
    root = args.output.resolve()
    initialize(root)
    if args.reuse_candidates:
        saved = json.loads((root/'candidates.json').read_text())
        assert saved['max_per_row'] == args.max_per_row
        rows = saved['rows']
        with gzip.open(root/'endpoint_jobs.jsonl.gz', 'rt') as f:
            jobs = [json.loads(line) for line in f]
    else:
        rows, jobs = prepare(root, args.max_per_row)
    sources = [Path(__file__), Path('src/fixed_tail_extend.cpp'), Path('tools/datf.py'),
               Path('tools/datf_quadratic.py'), Path('tools/datf_k06.py'),
               Path('tools/fixed_tail_search.py'), Path('tools/linear_seed_scan.py'),
               Path('tools/datf_paper_reassess.py'), Path('tools/datf_fixed_tail.py'),
               Path('tools/datf_paper_audit.py'),
               Path('tools/datf_audit.py'), Path('include/gaston.hpp'),
               Path('tests/test_fixed_tail_extend.py'),
               Path('results/pq_kernel_six_all/patterns.json'),
               Path('results/datf_pools/pools/pools.json')]
    atomic(root/'configuration.json', {
        'workers': args.workers, 'max_per_row': args.max_per_row,
        'max_bias_nodes': args.max_bias_nodes, 'first_L_omitted': True,
        'linear_tail_min_rounds': 3, 'differential_rounds': [0, 1, 2], 'middle_rounds': 1,
        'ordinary_DATF': True, 'MILP_used': False, 'uncomputed_is_not_zero': True,
        'integer_decimal_serialization_limit': 'disabled for trusted exact fractions; no arithmetic rule changed',
        'compiler': subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
        'python': sys.version,
        'code_and_input_sha256': {str(f.resolve()): hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
    })
    if args.prepare_only:
        print('READY', len(jobs), 'unique full jobs')
        return
    known = evaluate_batch(root, jobs, args.workers, args.max_bias_nodes, 'full_jobs')
    prefix_jobs = make_prefix_jobs(rows, known)
    prefixes = evaluate_batch(root, prefix_jobs, args.workers, args.max_bias_nodes, 'prefix_jobs')
    finalize(root, rows, known, prefixes)


if __name__ == '__main__':
    main()

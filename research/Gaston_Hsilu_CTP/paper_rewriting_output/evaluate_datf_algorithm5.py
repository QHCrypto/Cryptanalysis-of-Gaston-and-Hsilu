#!/usr/bin/env python3
"""Apply the paper's DATF Algorithm 5 to the two complete four-layer cores."""
import argparse
import hashlib
import importlib.util
import json
import random
import resource
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2] / '差分线性'
BACKEND = ROOT / 'gaston_dl/results/g4_middle_algorithm5_20261027'
sys.path.insert(0, str(BACKEND))
from algorithm5 import Algorithm5
from algorithm3 import ResourceLimit
from datf import evaluate
from datf_k06 import construct, parity_polynomial, value_record, verify_dag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', choices=['G4a', 'H4a'], required=True)
    args = ap.parse_args()
    output = Path(__file__).parent / 'dl_datf_algorithm5_20261111'
    output.mkdir(exist_ok=True)
    source = ROOT / '路径分段汇总_20261027/results.json'
    row = next(r for r in json.loads(source.read_text())['rows'] if r['code'] == args.target)
    delta = tuple(int(w, 16) for w in row['delta_core'])
    mask = tuple(int(w, 16) for w in row['lambda_core'])
    rounds, start = row['rounds'], row['start_round_index']
    resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
    def timeout(signum, frame):
        raise ResourceLimit('hard wall-time limit')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(180)
    started = time.perf_counter()
    result = {
        'id': args.target, 'algorithm': row['algorithm'], 'rounds': rounds,
        'start_round_index': start, 'delta_core': row['delta_core'], 'lambda_core': row['lambda_core'],
        'method': 'DATF Algorithm 5: Trade-Off Estimation, Journal of Cryptology 39:7 (2026), page 19',
        'scope': 'Direct evaluation of complete four-layer core; no differential-prefix or linear-tail multiplication',
        'quantity': 'signed correlation rho = 2 * paper bias epsilon',
        'transitional_boundary': 'before each chi after the first, following the internal linear map and round constant',
        'algorithm3_attempt': row['datf_algorithm3'],
        'limits': {'seconds': 150, 'max_nodes': 2000000, 'max_terms': 5000000, 'address_space_bytes': 4 << 30},
        'status': 'running', 'correlation': None,
    }
    alg = None
    try:
        if args.target == 'G4a':
            spec = importlib.util.spec_from_file_location('existing_algorithm5_validation', BACKEND / 'run.py')
            validation = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(validation)
            result['backend_validation'] = validation.validate_weighted_backend()
            engine, base, diff, _ = construct(delta, rounds, start)
            verify_dag(engine, base, diff, delta, rounds, start, samples=16)
        else:
            sys.path.insert(0, str(ROOT / 'hsilu_dl/tools'))
            import hsilu
            engine, base, diff = hsilu.construct(delta, rounds, start=start)
            rng = random.Random(20261111)
            for _ in range(16):
                state = tuple(rng.getrandbits(64) for _ in range(5))
                values = engine.evaluate_dag([(w >> j) & 1 for w in state for j in range(64)])
                out0 = hsilu.core(state, rounds, start)
                out1 = hsilu.core(tuple(x ^ d for x, d in zip(state, delta)), rounds, start)
                assert tuple(sum(evaluate(f, values) << j for j, f in enumerate(r)) for r in base) == out0
                assert tuple(sum(evaluate(f, values) << j for j, f in enumerate(r)) for r in diff) == tuple(x ^ y for x, y in zip(out0, out1))
        result['full_state_dag_checks'] = 16
        target = parity_polynomial(diff, mask)
        dag_path = output / (args.target + '_dag.json')
        dag_path.write_text(json.dumps({'target_terms': [hex(t) for t in sorted(target)], 'dag': engine.serialized_dag()}) + '\n')
        result['dag_file'] = str(dag_path)
        result['dag_sha256'] = hashlib.sha256(dag_path.read_bytes()).hexdigest()
        alg = Algorithm5(engine, seconds=150, max_nodes=2000000, max_terms=5000000)
        correlation = alg.estimate(target)
        result.update(status='completed', correlation=value_record(correlation), bias=value_record(correlation / 2))
        result['trace'] = alg.trace
        result['marginal_correlations'] = {str(k): value_record(v) for k, v in alg.distributions.items()}
    except (ResourceLimit, MemoryError, RecursionError) as exc:
        result.update(status='resource_limit', reason=repr(exc), correlation=None)
    finally:
        signal.alarm(0)
        result['elapsed_seconds'] = time.perf_counter() - started
        if alg:
            result['statistics'] = alg.statistics()
            result['phase'] = alg.phase
        files = [Path(__file__), source, BACKEND / 'algorithm5.py',
                 ROOT / 'DATF算法3复核_20261021/algorithm3.py',
                 ROOT / 'gaston_dl/tools/datf.py', ROOT / 'gaston_dl/tools/datf_k06.py',
                 ROOT / 'hsilu_dl/tools/hsilu.py']
        result['sources'] = [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
        (output / (args.target + '.json')).write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(args.target, result['status'], result['correlation'], round(result['elapsed_seconds'], 3))


if __name__ == '__main__':
    main()

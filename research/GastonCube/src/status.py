"""Read saved certificates and the current pilot process; never infer a coefficient."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    config = json.loads((ROOT/'cube_config.json').read_text())
    validation = json.loads((ROOT/'results/validation_summary.json').read_text())
    result = json.loads((ROOT/'results/r6_even_k0_o0_count.json').read_text())
    job_file = ROOT/'results/active_job.json'
    job = json.loads(job_file.read_text()) if job_file.exists() else {}
    pid = job.get('pid')
    process_file = Path(f'/proc/{pid}/cmdline') if pid else None
    alive = False
    if process_file and process_file.exists():
        alive = process_file.read_bytes().decode().rstrip('\0').split('\0') == job.get('argv')
    status = dict(nonce_indices=config['nonce_indices'],
                  selected_cube='even original nonce coordinates 0,2,...,124',
                  target='[k0 * product(n_i, i in I)] S0[0] after six rounds',
                  process_alive=alive, process_id=pid if alive else None,
                  pilot_key_degree_split=job.get('model_key_degree_split',False),
                  completed_batches=result['batches'], trails_counted=result['trails_counted'],
                  pilot_complete=result['complete'], pilot_coefficient=result['coefficient'],
                  certified_six_round_one_terms=([dict(key=0,output=0)]
                        if result['complete'] and result['coefficient']==1 else []),
                  full_128_by_128_matrix_complete=False,
                  validation_passed=validation['all_passed'],
                  note='Only complete results are coefficients. Initial cube selection used the legacy circuit and four candidates; it is not proven fastest. The new per-round key split is validated separately.')
    (ROOT/'results/research_status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps(status,indent=2))


if __name__ == '__main__':
    main()

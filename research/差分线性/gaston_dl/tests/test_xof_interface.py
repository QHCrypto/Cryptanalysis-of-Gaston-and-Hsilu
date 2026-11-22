#!/usr/bin/env python3
"""Independent direct-bit counts check the sampler's histogram accumulation."""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from fixed_tail_search import RC, state_chi
from linear_seed_scan import linear

MASK = (1 << 64)-1


def mix(z):
    z = ((z ^ (z >> 30))*0xbf58476d1ce4e5b9) & MASK
    z = ((z ^ (z >> 27))*0x94d049bb133111eb) & MASK
    return z ^ (z >> 31)


def permute(s, rounds):
    for rc in RC[12-rounds:]:
        s = linear(s)
        s = state_chi((s[0] ^ rc, *s[1:]))
    return s


def main():
    binary = ROOT/'build/xof_interface_baseline'
    subprocess.run([binary, '--selftest'], check=True)
    with tempfile.TemporaryDirectory() as temp:
        payloads = []
        for threads in (1, 3):
            dest = Path(temp)/f'{threads}.json'
            subprocess.run([binary, dest, '31', '1234', str(threads), '0,55'], check=True)
            payloads.append(json.loads(dest.read_text()))
        assert payloads[0]['rows'] == payloads[1]['rows']
        actual = {(r['rounds'], r['input_bit'], r['output_bit']): r['n1'] for r in payloads[0]['rows']}
        expected = dict.fromkeys(actual, 0)
        init = permute((0x0000080000cc0003, 0, 0, 0, 0), 12)
        for k in range(31):
            words = [mix((1234+(5*k+i+1)*0x9e3779b97f4a7c15) & MASK) for i in range(2)]
            s = permute((init[0] ^ words[0], *init[1:]), 12)
            s = (s[0] ^ (words[1] & ((1 << 56)-1)) ^ (1 << 56), *s[1:])
            for rounds in (3, 4):
                base = permute(s, rounds)[0]
                for bit in (0, 55):
                    diff = base ^ permute((s[0] ^ (1 << bit), *s[1:]), rounds)[0]
                    for out in range(64):
                        expected[rounds, bit, out] += (diff >> out) & 1
        assert expected == actual
    print('PASS independent Python message-pair counts; thread-independent sampling')


if __name__ == '__main__':
    main()

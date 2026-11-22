"""Check the Section 2 parameter transcription against existing implementations.

Tables are parsed from the numerical inputs extracted from the manuscript. Equations are transcribed below and
compared with independently maintained Python and native implementations.
Finite mode cross-checks are implementation checks, not security proofs.
"""
from pathlib import Path
import ctypes
import hashlib
import json
import random
import re
import sys

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parent
SOURCE = PAPER / 'verification_inputs/paper_parameters.tex'
tex = SOURCE.read_text()
sys.path.insert(0, str(ROOT / 'HsiluCube/src'))
sys.path.insert(0, str(ROOT / 'nist_compatibility_20261109'))
import permutations as reference
import nist_modes as modes


def table(label):
    pos = tex.index(r'\label{' + label + '}')
    return tex[pos:tex.index(r'\end{table}', pos)]


def row(block, name):
    match = re.search(re.escape(name) + r'\s*&([^\n]+)', block)
    return list(map(int, re.findall(r'\d+', match[1])))


gtable = table('tab:gaston-parameters')
GE, GW, GT = [row(gtable, '$' + x + '_i$') for x in ('e', 'w', 't')]
assert (GE, GW, GT) == (reference.GE, reference.GW, reference.GT)
steps = {}
for line in table('tab:hsilu-parameters').splitlines():
    if re.match(r'^\d+&', line):
        values = list(map(int, re.findall(r'\d+', line)))
        for offset in range(0, len(values), 4):
            k, a, b, m = values[offset:offset+4]
            assert a < k + 5 and b < k + 5
            steps[k] = (a, b, m)
assert len(steps) == 13
assert [steps[k] for k in range(13)] == list(zip(reference.HN1, reference.HN2, reference.HM))
hmatch = re.search(r'\(h_0,h_1,h_2,h_3,h_4\)&=\(([^)]+)\)', tex)
HR = list(map(int, hmatch[1].split(',')))
assert HR == reference.HR
# The submission defines chi algebraically; the earlier S-box table was removed.
# Independently tabulated chi5, with row zero as the least-significant bit.
SBOX = [0,9,18,11,5,12,22,15,10,3,24,1,13,4,30,7,
        20,21,6,23,17,16,2,19,26,27,8,25,29,28,14,31]
assert len(SBOX) == len(set(SBOX)) == 32
for x, y in enumerate(SBOX):
    bits = [(x >> j) & 1 for j in range(5)]
    assert y == sum(b << j for j, b in enumerate(reference.chi(bits)))
degrees = []
for bit in range(5):
    anf = [(y >> bit) & 1 for y in SBOX]
    for j in range(5):
        for x in range(32):
            if (x >> j) & 1:
                anf[x] ^= anf[x ^ (1 << j)]
    degrees.append(max(x.bit_count() for x, coefficient in enumerate(anf) if coefficient))
assert degrees == [2] * 5
array = tex[tex.index(r'\begin{array}{c|rrrrrr}'):tex.index(r'\end{array}')]
RC = [int(v, 16) for v in re.findall(r'\\mathtt\{([0-9a-f]{2})\}', array)]
assert RC == reference.RC
MASK = (1 << 64) - 1


def rol(x, n):
    return ((x << n) | (x >> ((64 - n) % 64))) & MASK


def linear(state, name):
    if name == 'gaston':
        a = [rol(x, e) for x, e in zip(state, GE)]
        p = q = 0
        for x, t in zip(a, GT):
            p ^= x
            q ^= rol(x, t)
        effect = p ^ rol(p, 1) ^ q ^ rol(q, 18)
        return [rol(x ^ rol(effect, 23), w) for x, w in zip(a, GW)]
    x = list(state)
    for k in range(13):
        a, b, m = steps[k]
        x.append(x[a] ^ rol(x[b], m))
    return [rol(x[13+i], HR[i]) for i in range(5)]


def chi(state):
    output = [0] * 5
    for j in range(64):
        y = SBOX[sum(((state[i] >> j) & 1) << i for i in range(5))]
        for i in range(5):
            output[i] |= ((y >> i) & 1) << j
    return output


def permute(state, name, rounds):
    x = list(state)
    for c in RC[12-rounds:]:
        if name == 'gaston':
            x = linear(x, name)
        x[0] ^= c
        x = chi(x)
        if name == 'hsilu':
            x = linear(x, name)
    return x


def binary_rank(vectors):
    pivots = {}
    for v in vectors:
        while v:
            lead = v.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = v
                break
            v ^= pivots[lead]
    return len(pivots)


linear_ranks = {}
for name in ('gaston', 'hsilu'):
    images = []
    for j in range(320):
        state = [0] * 5
        state[j // 64] = 1 << (j % 64)
        output = linear(state, name)
        assert output == getattr(reference, name + '_linear')(state)
        images.append(sum(x << (64*i) for i, x in enumerate(output)))
    linear_ranks[name] = binary_rank(images)
assert set(linear_ranks.values()) == {320}

lib = ctypes.CDLL(str(ROOT / 'nist_compatibility_20261109/native.so'))
U = ctypes.c_uint64
for name in ('hperm', 'gperm'):
    getattr(lib, name).argtypes = [ctypes.POINTER(U), ctypes.c_int]
rng = random.Random(20261111)
permutation_checks = 0
for name in ('gaston', 'hsilu'):
    for rounds in range(1, 13):
        for state in [[0]*5, [MASK]*5] + [[rng.getrandbits(64) for _ in range(5)] for _ in range(8)]:
            output = permute(state, name, rounds)
            native = (U*5)(*state)
            getattr(lib, 'gperm' if name == 'gaston' else 'hperm')(native, rounds)
            assert output == list(native) == reference.permute(state, name, rounds)
            permutation_checks += 1
for state, output in reference.GASTON_KAT:
    assert permute(state, 'gaston', 12) == output

IV_X = int(re.search(r'IV_X=\\hex\{([0-9a-f]+)\}', tex)[1], 16)
IV_A = int(re.search(r'IV_A=\\hex\{([0-9a-f]+)\}', tex)[1], 16)
assert (IV_X, IV_A) == (modes.XOF_IV, modes.AEAD_IV)


def xof(message, length, name, rounds, output_length):
    bits = [(message >> i) & 1 for i in range(length)] + [1]
    bits += [0] * ((-len(bits)) % 64)
    s = permute([IV_X, 0, 0, 0, 0], name, rounds)
    for offset in range(0, len(bits), 64):
        s[0] ^= sum(bits[offset+j] << j for j in range(64))
        s = permute(s, name, rounds)
    output = []
    while len(output) < output_length:
        output.extend((s[0] >> j) & 1 for j in range(64))
        if len(output) < output_length:
            s = permute(s, name, rounds)
    return sum(b << i for i, b in enumerate(output[:output_length]))


def aead(key, nonce, ad, plaintext, name, rounds):
    k = [int.from_bytes(key[i:i+8], 'little') for i in (0, 8)]
    n = [int.from_bytes(nonce[i:i+8], 'little') for i in (0, 8)]
    s = permute([IV_A] + k + n, name, rounds)
    s[3] ^= k[0]
    s[4] ^= k[1]
    def absorb(block):
        s[0] ^= int.from_bytes(block[:8], 'little')
        s[1] ^= int.from_bytes(block[8:], 'little')
    if ad:
        padded = ad + b'\x01' + bytes((-(len(ad)+1)) % 16)
        for offset in range(0, len(padded), 16):
            absorb(padded[offset:offset+16])
            s = permute(s, name, 8)
    s[4] ^= 1 << 63
    ciphertext = bytearray()
    for offset in range(0, len(plaintext)//16*16, 16):
        absorb(plaintext[offset:offset+16])
        ciphertext.extend(b''.join(w.to_bytes(8, 'little') for w in s[:2]))
        s = permute(s, name, 8)
    tail = plaintext[len(plaintext)//16*16:]
    absorb(tail + b'\x01' + bytes(15-len(tail)))
    ciphertext.extend(b''.join(w.to_bytes(8, 'little') for w in s[:2])[:len(tail)])
    s[2] ^= k[0]
    s[3] ^= k[1]
    s = permute(s, name, 12)
    tag = b''.join((s[3+i] ^ k[i]).to_bytes(8, 'little') for i in range(2))
    return bytes(ciphertext), tag


xof_checks = aead_checks = 0
lengths = (0, 1, 7, 8, 62, 63, 64, 65, 127, 128, 190, 191, 192)
for name in ('gaston', 'hsilu'):
    for rounds in (2, 3, 6, 7, 12):
        for length in lengths:
            message = rng.getrandbits(length)
            assert xof(message, length, name, rounds, 193) == modes.xof_bits(message, length, name, rounds, 193)
            xof_checks += 1
    for rounds in (6, 7, 12):
        for ad_length in (0, 1, 15, 16, 17, 32):
            for plaintext_length in (0, 1, 15, 16, 17, 32):
                args = (rng.randbytes(16), rng.randbytes(16), rng.randbytes(ad_length), rng.randbytes(plaintext_length), name, rounds)
                assert aead(*args) == modes.aead(*args)
                aead_checks += 1

report = {
    'date': '2026-11-11', 'source': str(SOURCE),
    'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'parameters_parsed_from_latex': True, 'chi_lookup_entries': len(SBOX),
    'chi_coordinate_degrees': degrees, 'gaston_rotation_triples': 5,
    'hsilu_circuit_steps': len(steps), 'round_constants': len(RC),
    'linear_basis_comparisons_per_permutation': 320, 'linear_ranks': linear_ranks,
    'permutation_comparisons_with_python_and_native': permutation_checks,
    'round_counts_checked': list(range(1, 13)), 'gaston_original_vectors': 3,
    'bitstring_xof_comparisons': xof_checks, 'xof_message_lengths_bits': lengths,
    'aead_comparisons': aead_checks, 'aead_ad_and_plaintext_lengths_bytes': [0, 1, 15, 16, 17, 32],
    'all_passed': True,
    'scope': 'Exact finite parameter, lookup and linear-map checks; sampled nonlinear permutation/mode agreement. No new security claim.',
}
target = Path(__file__).with_name('section2_validation.json')
target.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
print(json.dumps(report, indent=2, ensure_ascii=False))

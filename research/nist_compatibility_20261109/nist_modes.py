"""SP 800-232 outer modes with an explicit permutation argument.

Integer bit i stores the i-th bit of the input string (LSB-first). Gaston
and Hsilu keep their original word-level permutations and last-round constants.
The IV stays at its full-mode value when permutation calls are reduced.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'HsiluCube/src'))
from permutations import permute as research_permute

MASK = (1 << 64) - 1
XOF_IV = 0x0000080000CC0003
AEAD_IV = 0x00001000808C0001


def rotr(x, n):
    return ((x >> n) | (x << (64-n))) & MASK


def ascon_permute(state, rounds):
    a = list(state)
    for rc in [0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b][-rounds:]:
        a[2] ^= rc
        a[0] ^= a[4]; a[4] ^= a[3]; a[2] ^= a[1]
        t = [(~a[i] & a[(i+1) % 5]) & MASK for i in range(5)]
        a = [a[i] ^ t[(i+1) % 5] for i in range(5)]
        a[1] ^= a[0]; a[0] ^= a[4]; a[3] ^= a[2]; a[2] ^= MASK
        a = [v ^ rotr(v, x) ^ rotr(v, y) for v, (x,y) in zip(a, [(19,28),(61,39),(1,6),(10,17),(7,41)])]
    return a


def permute(state, name, rounds):
    if not 1 <= rounds <= 12:
        raise ValueError('round count outside the audited range')
    if name == 'ascon':
        return ascon_permute(state, rounds)
    return research_permute(state, name, rounds)


def prefix_state(words, name, rounds):
    s = permute([XOF_IV, 0, 0, 0, 0], name, rounds)
    for w in words:
        s[0] ^= w
        s = permute(s, name, rounds)
    return s


def xof_bits(message, bit_length, name, rounds=12, output_bits=128):
    if bit_length < 0 or message < 0 or message >> bit_length or output_bits <= 0:
        raise ValueError('invalid bit string')
    full, tail = divmod(bit_length, 64)
    s = prefix_state([(message >> (64*i)) & MASK for i in range(full)], name, rounds)
    s[0] ^= (message >> (64*full)) | (1 << tail)
    s = permute(s, name, rounds)
    out = 0
    for i in range((output_bits + 63)//64):
        if i:
            s = permute(s, name, rounds)
        out |= s[0] << (64*i)
    return out & ((1 << output_bits)-1)


def xof190(prefix, tail, name, rounds):
    assert len(prefix) == 2 and 0 <= tail < 1 << 62
    message = prefix[0] | (prefix[1] << 64) | (tail << 128)
    out = xof_bits(message, 190, name, rounds)
    return out & MASK, out >> 64


def aead(key, nonce, ad, plaintext, name, initialization_rounds=12,
         processing_rounds=8, finalization_rounds=12):
    """NIST Algorithm 3 for arbitrary byte lengths; returns ciphertext, tag."""
    assert len(key) == len(nonce) == 16
    k0, k1 = [int.from_bytes(key[i:i+8], 'little') for i in (0,8)]
    n0, n1 = [int.from_bytes(nonce[i:i+8], 'little') for i in (0,8)]
    s = permute([AEAD_IV,k0,k1,n0,n1], name, initialization_rounds)
    s[3] ^= k0; s[4] ^= k1

    def absorb(block):
        v = int.from_bytes(block, 'little')
        s[0] ^= v & MASK; s[1] ^= v >> 64

    if ad:
        for i in range(0, len(ad)//16*16, 16):
            absorb(ad[i:i+16]); s = permute(s, name, processing_rounds)
        tail = ad[len(ad)//16*16:]
        absorb(tail + b'\x01' + bytes(15-len(tail)))
        s = permute(s, name, processing_rounds)
    s[4] ^= 1 << 63
    ciphertext = bytearray()
    for i in range(0, len(plaintext)//16*16, 16):
        absorb(plaintext[i:i+16])
        ciphertext += s[0].to_bytes(8,'little') + s[1].to_bytes(8,'little')
        s = permute(s, name, processing_rounds)
    tail = plaintext[len(plaintext)//16*16:]
    absorb(tail + b'\x01' + bytes(15-len(tail)))
    ciphertext += (s[0].to_bytes(8,'little') + s[1].to_bytes(8,'little'))[:len(tail)]
    # NIST final key injection: 0^128 || K || 0^64.
    s[2] ^= k0; s[3] ^= k1
    s = permute(s, name, finalization_rounds)
    tag = (s[3]^k0).to_bytes(8,'little') + (s[4]^k1).to_bytes(8,'little')
    return bytes(ciphertext), tag

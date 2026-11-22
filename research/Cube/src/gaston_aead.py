"""Word-oriented Gaston-AEAD128-style outer mode.

This is a research model, not a standards submission.  It follows the NIST
Ascon-AEAD128 outer schedule and replaces every Ascon permutation call by
Gaston.  The implementation is deliberately explicit about the little-endian
64-bit words used by the NIST specification.  The helper currently targets
full 128-bit message blocks; that is sufficient for the first-ciphertext cube
interface used in the reduced-round experiments.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from permutations import permute

IV = 0x00001000808C0001
MASK64 = (1 << 64) - 1


def _words128(value):
    """Split a 128-bit integer into the state words S0,S1."""
    return value & MASK64, (value >> 64) & MASK64


def _block_words(block):
    if len(block) != 16:
        raise ValueError("full-block helper requires exactly 16 bytes")
    return (int.from_bytes(block[:8], "little"),
            int.from_bytes(block[8:], "little"))


def _words_block(a, b):
    return a.to_bytes(8, "little") + b.to_bytes(8, "little")


def encrypt_full_blocks(key, nonce, plaintext, *, initialization_rounds=12,
                        processing_rounds=8, finalization_rounds=12):
    """Encrypt a multiple-of-16-byte message with empty associated data.

    ``key`` and ``nonce`` are 16-byte strings.  The returned pair is
    ``(ciphertext, tag)``.  For a multiple-of-rate message, the empty final
    fragment still contributes the standard padding bit before finalization.
    """
    if len(key) != 16 or len(nonce) != 16:
        raise ValueError("key and nonce must be 16 bytes")
    if len(plaintext) % 16:
        raise ValueError("this research helper accepts full 128-bit blocks only")

    k0, k1 = _block_words(key)
    n0, n1 = _block_words(nonce)
    state = permute([IV, k0, k1, n0, n1], "gaston", initialization_rounds)

    # Initialization key injection is into the 192-bit capacity suffix.
    state[3] ^= k0
    state[4] ^= k1

    # Empty AD: domain separation is the high bit of S4 in the integer-word
    # representation specified by NIST SP 800-232.
    state[4] ^= 0x8000000000000000

    ciphertext = bytearray()
    for offset in range(0, len(plaintext), 16):
        p0, p1 = _block_words(plaintext[offset:offset + 16])
        c0, c1 = state[0] ^ p0, state[1] ^ p1
        ciphertext += _words_block(c0, c1)
        state[0] ^= p0
        state[1] ^= p1
        state = permute(state, "gaston", processing_rounds)

    # parse() returns an empty final fragment for an exact multiple.  Its
    # 128-bit padding is represented by the low bit of S0.
    state[0] ^= 1
    # NIST SP 800-232 Algorithm 3: 0^128 || K || 0^64.
    state[2] ^= k0
    state[3] ^= k1
    state = permute(state, "gaston", finalization_rounds)
    tag = _words_block(state[3] ^ k0, state[4] ^ k1)
    return bytes(ciphertext), tag


def first_ciphertext_rate_block(key_words, nonce_words, *, rounds=12,
                                plaintext_words=(0, 0)):
    """Return the first 128-bit rate block before the following p[8] call.

    This is the exact observable used by the reduced-round cube experiments
    when associated data is empty and the first plaintext block is complete.
    """
    state = permute([IV, *key_words, *nonce_words], "gaston", rounds)
    return state[0] ^ plaintext_words[0], state[1] ^ plaintext_words[1]


def self_check():
    key = bytes.fromhex("00112233445566778899aabbccddeeff")
    nonce = bytes.fromhex("ffeeddccbbaa99887766554433221100")
    plaintext = bytes(16)
    ciphertext, tag = encrypt_full_blocks(key, nonce, plaintext)
    k = _block_words(key)
    n = _block_words(nonce)
    direct = first_ciphertext_rate_block(k, n)
    assert ciphertext[:16] == _words_block(*direct)
    assert len(tag) == 16
    return {"first_block_matches_initialization_rate": True,
            "ciphertext_bytes": len(ciphertext), "tag_bytes": len(tag)}


if __name__ == "__main__":
    import json
    print(json.dumps(self_check(), indent=2))

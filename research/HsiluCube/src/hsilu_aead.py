"""Word-oriented Hsilu-AEAD128 research model.

The outer schedule follows the NIST Ascon-AEAD128 layout with empty
associated data and full 128-bit message blocks.  Every Ascon permutation
call is replaced by the Hsilu permutation from ``hsilu_ref``.  This module is
an executable consistency model for the cube experiments; it is not a claim
about a standardised Hsilu-based construction.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from permutations import permute

IV = 0x00001000808C0001
MASK64 = (1 << 64) - 1


def _block_words(block):
    if len(block) != 16:
        raise ValueError("full-block helper requires exactly 16 bytes")
    return (int.from_bytes(block[:8], "little"),
            int.from_bytes(block[8:], "little"))


def _words_block(a, b):
    return a.to_bytes(8, "little") + b.to_bytes(8, "little")


def encrypt_full_blocks(key, nonce, plaintext, *, initialization_rounds=12,
                        processing_rounds=8, finalization_rounds=12):
    """Encrypt a multiple-of-16-byte message with empty associated data."""
    if len(key) != 16 or len(nonce) != 16:
        raise ValueError("key and nonce must be 16 bytes")
    if len(plaintext) % 16:
        raise ValueError("this research helper accepts full 128-bit blocks only")

    k0, k1 = _block_words(key)
    n0, n1 = _block_words(nonce)
    state = permute([IV, k0, k1, n0, n1], "hsilu", initialization_rounds)
    state[3] ^= k0
    state[4] ^= k1
    state[4] ^= 0x8000000000000000  # empty-AD domain separation

    ciphertext = bytearray()
    for offset in range(0, len(plaintext), 16):
        p0, p1 = _block_words(plaintext[offset:offset + 16])
        ciphertext += _words_block(state[0] ^ p0, state[1] ^ p1)
        state[0] ^= p0
        state[1] ^= p1
        state = permute(state, "hsilu", processing_rounds)

    # Exact-multiple messages still absorb the final padding bit.
    state[0] ^= 1
    # NIST SP 800-232 Algorithm 3: 0^128 || K || 0^64.
    state[2] ^= k0
    state[3] ^= k1
    state = permute(state, "hsilu", finalization_rounds)
    tag = _words_block(state[3] ^ k0, state[4] ^ k1)
    return bytes(ciphertext), tag


def first_ciphertext_rate_block(key_words, nonce_words, *, rounds=12,
                                plaintext_words=(0, 0)):
    """Return the first 128-bit rate block observed by the cube interface."""
    state = permute([IV, *key_words, *nonce_words], "hsilu", rounds)
    return state[0] ^ plaintext_words[0], state[1] ^ plaintext_words[1]


def self_check():
    key = bytes.fromhex("00112233445566778899aabbccddeeff")
    nonce = bytes.fromhex("ffeeddccbbaa99887766554433221100")
    ciphertext, tag = encrypt_full_blocks(key, nonce, bytes(16))
    direct = first_ciphertext_rate_block(_block_words(key), _block_words(nonce))
    assert ciphertext[:16] == _words_block(*direct)
    assert len(tag) == 16
    return {"first_block_matches_initialization_rate": True,
            "ciphertext_bytes": len(ciphertext), "tag_bytes": len(tag)}


if __name__ == "__main__":
    import json
    print(json.dumps(self_check(), indent=2))

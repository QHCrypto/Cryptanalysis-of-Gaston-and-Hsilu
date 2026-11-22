"""Independent full Boolean Gaston reference; bit 0 is the word LSB.

Actual Gaston has 64-bit words. Smaller widths are validation models only.
"""
IV = 0x00001000808C0001
E = (0, 60, 22, 27, 4)
W = (0, 56, 31, 46, 43)
T = (25, 32, 52, 60, 63)
RC = (0xF0, 0xE1, 0xD2, 0xC3, 0xB4, 0xA5,
      0x96, 0x87, 0x78, 0x69, 0x5A, 0x4B)


def rol(x, shift, width=64):
    shift %= width
    return ((x << shift) | (x >> ((width - shift) % width))) & ((1 << width) - 1)


def linear(s, width=64):
    a = [rol(x, e, width) for x, e in zip(s, E)]
    p = a[0] ^ a[1] ^ a[2] ^ a[3] ^ a[4]
    p ^= rol(p, 1, width)
    q = 0
    for x, t in zip(a, T):
        q ^= rol(x, t, width)
    q ^= rol(q, 18, width)
    p = rol(p ^ q, 23, width)
    return [rol(x ^ p, w, width) for x, w in zip(a, W)]


def permute(s, rounds=6, width=64, start=None):
    mask = (1 << width) - 1
    a = [x & mask for x in s]
    start = 12 - rounds if start is None else start
    for rc in RC[start:start + rounds]:
        b = linear(a, width)
        b[0] ^= rc & mask
        a = [(b[i] ^ ((~b[(i+1) % 5]) & b[(i+2) % 5])) & mask
             for i in range(5)]
    return a


def coefficient(rounds, key, nonce, output, width=64, iv=IV, start=None):
    """Exact ANF coefficient by exhaustive Boolean evaluation (small tests).

    Other key/nonce inputs are zero. Includes both values of the key variable.
    """
    value = 0
    for assignment in range(1 << (len(nonce) + 1)):
        s = [iv, 0, 0, 0, 0]
        s[1 + key // width] = (assignment & 1) << (key % width)
        for j, n in enumerate(nonce):
            s[3 + n // width] |= ((assignment >> (j+1)) & 1) << (n % width)
        value ^= (permute(s, rounds, width, start)[output // width]
                  >> (output % width)) & 1
    return value

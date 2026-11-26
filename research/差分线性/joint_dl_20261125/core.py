"""Small standalone Gaston/Hsilu core, native LSB convention."""
from functools import reduce
from operator import xor

MASK = (1 << 64) - 1
ZERO = (0,) * 5
RC = (0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b)
E = (0, 60, 22, 27, 4)
W = (0, 56, 31, 46, 43)
T = (25, 32, 52, 60, 63)
H = (0, 57, 13, 31, 20)
N1 = (0, 3, 1, 0, 5, 5, 7, 11, 12, 12, 13, 14, 15)
N2 = (4, 4, 2, 1, 2, 3, 6, 8, 9, 6, 10, 11, 9)
SH = (47, 51, 54, 30, 4, 55, 0, 14, 59, 4, 10, 45, 53)


def rot(x, s):
    s %= 64
    return ((x << s) | (x >> ((-s) % 64))) & MASK


def add(*states):
    return tuple(reduce(xor, lane) for lane in zip(*states))


def dot(x, y):
    return reduce(xor, (a & b for a, b in zip(x, y))).bit_count() & 1


def pack(x):
    return sum(v << (64 * i) for i, v in enumerate(x))


def embed(j, p):
    return tuple(((p >> i) & 1) << j for i in range(5))


def column(x, j):
    return sum(((x[i] >> j) & 1) << i for i in range(5))


def chi(x):
    return tuple(x[i] ^ ((~x[(i + 1) % 5] & MASK) & x[(i + 2) % 5]) for i in range(5))


def sbox(x):
    return column(chi(embed(0, x)), 0)


def polar(x, d):
    return tuple((x[(i + 1) % 5] & d[(i + 2) % 5]) ^
                 (d[(i + 1) % 5] & x[(i + 2) % 5]) for i in range(5))


def polar_transpose(w, d):
    return tuple((w[(i - 1) % 5] & d[(i + 1) % 5]) ^
                 (w[(i - 2) % 5] & d[(i - 1) % 5]) for i in range(5))


def lg(x):
    a = tuple(rot(v, e) for v, e in zip(x, E))
    p = reduce(xor, a)
    q = reduce(xor, (rot(v, t) for v, t in zip(a, T)))
    z = rot(p ^ rot(p, 1) ^ q ^ rot(q, 18), 23)
    return tuple(rot(v ^ z, w) for v, w in zip(a, W))


def lgt(beta):
    d = tuple(rot(v, -s) for v, s in zip(beta, W))
    p = reduce(xor, d)
    return tuple(rot(d[i] ^ rot(p, -23) ^ rot(p, -24) ^
                     rot(p, -23-T[i]) ^ rot(p, -41-T[i]), -E[i]) for i in range(5))


def lh(x):
    v = list(x)
    for a, b, s in zip(N1, N2, SH):
        v.append(v[a] ^ rot(v[b], s))
    return tuple(rot(v[13+i], H[i]) for i in range(5))


def lht(beta):
    g = [0] * 18
    for i in range(5):
        g[13+i] = rot(beta[i], -H[i])
    for k in range(12, -1, -1):
        g[N1[k]] ^= g[5+k]
        g[N2[k]] ^= rot(g[5+k], -SH[k])
    return tuple(g[:5])


def permutation(x, name, rounds=4):
    L = lg if name == 'Gaston' else lh
    x = chi(x)
    for i in range(1, rounds):
        x = chi(add(L(x), (RC[12-rounds+i], 0, 0, 0, 0)))
    return x


ENDPOINTS = {
    'Gaston': (embed(1, 1), add(embed(0, 2), embed(54, 8))),
    'Hsilu': (embed(18, 8), add(embed(9, 17), embed(2, 4))),
}

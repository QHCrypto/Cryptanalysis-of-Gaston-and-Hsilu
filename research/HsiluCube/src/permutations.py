"""Independent word-level transcriptions of the original specifications.

Bit j is the coefficient of 2**j; all rotations are left rotations.
No material from previous Gaston analysis directories is used.
"""
MASK = (1 << 64) - 1
RC = [0xF0, 0xE1, 0xD2, 0xC3, 0xB4, 0xA5,
      0x96, 0x87, 0x78, 0x69, 0x5A, 0x4B]
GE = [0, 60, 22, 27, 4]
GW = [0, 56, 31, 46, 43]
GT = [25, 32, 52, 60, 63]
HN1 = [0, 3, 1, 0, 5, 5, 7, 11, 12, 12, 13, 14, 15]
HN2 = [4, 4, 2, 1, 2, 3, 6, 8, 9, 6, 10, 11, 9]
HM = [47, 51, 54, 30, 4, 55, 0, 14, 59, 4, 10, 45, 53]
HR = [0, 57, 13, 31, 20]


def rol(x, n):
    n %= 64
    return ((x << n) | (x >> ((64 - n) % 64))) & MASK


def chi(s):
    return [s[i] ^ ((~s[(i + 1) % 5]) & s[(i + 2) % 5])
            for i in range(5)]


def gaston_effect(s):
    a = [rol(x, e) for x, e in zip(s, GE)]
    p = q = 0
    for x, t in zip(a, GT):
        p ^= x
        q ^= rol(x, t)
    return rol(p ^ rol(p, 1) ^ q ^ rol(q, 18), 23)


def gaston_linear(s):
    effect = gaston_effect(s)
    return [rol(rol(x, e) ^ effect, w) for x, e, w in zip(s, GE, GW)]


def hsilu_linear(s):
    x = list(s)
    for a, b, m in zip(HN1, HN2, HM):
        x.append(x[a] ^ rol(x[b], m))
    return [rol(x[13 + j], HR[j]) for j in range(5)]


def permute_window(s, name, start, rounds):
    if (name not in ("gaston", "hsilu") or not 0 <= start <= 12
            or not 0 <= rounds <= 12 or start + rounds > 12):
        raise ValueError((name, start, rounds))
    a = list(s)
    for rc in RC[start:start + rounds]:
        if name == "gaston":
            a = gaston_linear(a)
        a[0] ^= rc
        a = chi(a)
        if name == "hsilu":
            a = hsilu_linear(a)
    return a


def permute(s, name, rounds):
    if not 0 <= rounds <= 12:
        raise ValueError((name, rounds))
    return permute_window(s, name, 12 - rounds, rounds)


GASTON_KAT = [
    ([0] * 5,
     [0x88B326096BEBC635, 0x6CA8FB64BC5CE6CA, 0xF1CE3840D8190713,
      0x54D70067438689B5, 0xF17FE863F958F32B]),
    ([0x1F4AD9906DA6A254, 0x4B84D7F83F2BDDFA, 0x468A0853578A00E3,
      0x6C05A0506DF7F66E, 0x4EFB22112453C964],
     [0x1BA89B5B5C4583B6, 0x22135709AE53417D, 0x9847B975E9EC9F3D,
      0xCE042DF2A402591D, 0x563EC68FC30307EA]),
    ([MASK, 0x0123456789ABCDEF, 0xFEDCBA9876543210,
      0xAAAAAAAAAAAAAAAA, 0x0101010101010101],
     [0x3117D51B14937067, 0x338F17F773C13F79, 0xDFB86E0868D252AB,
      0x0D461D35EB863DE7, 0x08BCE3E354C7231A]),
]


def self_check():
    for source, target in GASTON_KAT:
        assert permute(source, "gaston", 12) == target
    return {"gaston_original_paper_vectors": len(GASTON_KAT)}


if __name__ == "__main__":
    import json
    print(json.dumps(self_check(), indent=2))

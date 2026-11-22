#pragma once

#include <array>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace gdl {
using State = std::array<uint64_t, 5>;
inline constexpr std::array<uint64_t, 12> RC = {
    0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b};
inline constexpr int EAST[5] = {0, 60, 22, 27, 4};
inline constexpr int WEST[5] = {0, 56, 31, 46, 43};
inline constexpr int THETA[5] = {25, 32, 52, 60, 63};

inline uint64_t rotl(uint64_t x, unsigned r) {
    r &= 63;
    return (x << r) | (x >> ((64 - r) & 63));
}
inline int parity(unsigned x) { return __builtin_parity(x); }
inline int dot(const State &a, const State &b) {
    uint64_t w = 0;
    for (int i = 0; i < 5; ++i) w ^= a[i] & b[i];
    return __builtin_parityll(w);
}
inline State operator^(State a, const State &b) {
    for (int i = 0; i < 5; ++i) a[i] ^= b[i];
    return a;
}
inline State embed(int col, int pattern) {
    State a{};
    for (int i = 0; i < 5; ++i)
        if ((pattern >> i) & 1) a[i] = uint64_t(1) << col;
    return a;
}
inline unsigned column(const State &a, int col) {
    unsigned v = 0;
    for (int i = 0; i < 5; ++i) v |= unsigned((a[i] >> col) & 1) << i;
    return v;
}
inline State rotate(State a, int r) {
    for (auto &x : a) x = rotl(x, r);
    return a;
}
inline State linear(State a) {
    for (int i = 0; i < 5; ++i) a[i] = rotl(a[i], EAST[i]);
    uint64_t p = 0, q = 0;
    for (int i = 0; i < 5; ++i) { p ^= a[i]; q ^= rotl(a[i], THETA[i]); }
    p ^= rotl(p, 1);
    q ^= rotl(q, 18);
    p = rotl(p ^ q, 23);
    for (int i = 0; i < 5; ++i) a[i] = rotl(a[i] ^ p, WEST[i]);
    return a;
}
inline State chi(const State &a) {
    State b{};
    for (int i = 0; i < 5; ++i) b[i] = a[i] ^ (~a[(i + 1) % 5] & a[(i + 2) % 5]);
    return b;
}
inline unsigned chi5(unsigned a) {
    unsigned b = 0;
    for (int i = 0; i < 5; ++i)
        b |= (((a >> i) & 1) ^ ((1 ^ ((a >> ((i + 1) % 5)) & 1)) &
                                  ((a >> ((i + 2) % 5)) & 1))) << i;
    return b;
}
inline State round(State a, uint64_t rc) {
    a = linear(a);
    a[0] ^= rc;
    return chi(a);
}

// Input is immediately before the first chi. Its preceding linear layer is
// omitted; its constant is absorbed by a translation of the uniform input.
// All subsequent reference round constants remain in the evaluated function.
inline State core(State a, int rounds, int start, bool constants = true) {
    a = chi(a);
    for (int r = 1; r < rounds; ++r) a = round(a, constants ? RC[start + r] : 0);
    return a;
}
inline uint64_t splitmix_final(uint64_t z) {
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
// Counter-indexed generation makes a run reproducible across thread counts.
inline State sample_state(uint64_t seed, uint64_t index) {
    State s{};
    for (uint64_t i = 0; i < 5; ++i)
        s[i] = splitmix_final(seed + (5 * index + i + 1) * UINT64_C(0x9e3779b97f4a7c15));
    return s;
}
inline std::string hex_state(const State &a) {
    std::ostringstream s;
    s << std::hex << std::setfill('0');
    for (int i = 0; i < 5; ++i) {
        if (i) s << ';';
        s << std::setw(16) << a[i];
    }
    return s.str();
}
} // namespace gdl

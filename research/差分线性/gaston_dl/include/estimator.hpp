#pragma once

#include "gaston.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace gdl {
using Real = long double;
using Local = std::array<Real, 32>;
using Correlations = std::array<Local, 64>;

struct Tables {
    int ddt[32][32]{};       // [input difference][output difference], counts / 32
    int lat[32][32]{};       // [output mask][input mask], signed counts / 32
    Real dlct[32][32]{};     // [input difference][output mask], normalized correlation
    std::array<std::vector<std::pair<int, Real>>, 32> kernel;
    Tables() {
        for (int a = 0; a < 32; ++a)
            for (int x = 0; x < 32; ++x) ++ddt[a][chi5(x) ^ chi5(x ^ a)];
        for (int v = 0; v < 32; ++v)
            for (int u = 0; u < 32; ++u) {
                for (int x = 0; x < 32; ++x)
                    lat[v][u] += parity((v & chi5(x)) ^ (u & x)) ? -1 : 1;
                if (lat[v][u]) kernel[v].push_back({u, Real(lat[v][u] * lat[v][u]) / 1024});
            }
        for (int a = 0; a < 32; ++a)
            for (int v = 0; v < 32; ++v)
                for (int b = 0; b < 32; ++b)
                    dlct[a][v] += Real(ddt[a][b]) / 32 * (parity(v & b) ? -1 : 1);
    }
};

struct LinearMasks {
    std::array<State, 320> rows{};
    // terms[col * 32 + v]: nonzero 5-bit blocks of L^T embed(col,v).
    std::array<std::vector<std::pair<int, int>>, 2048> terms;
    LinearMasks() {
        for (int in = 0; in < 320; ++in) {
            const auto out = linear(embed(in % 64, 1 << (in / 64)));
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 64; ++j)
                    if ((out[i] >> j) & 1) rows[64 * i + j][in / 64] |= uint64_t(1) << (in % 64);
        }
        for (int j = 0; j < 64; ++j)
            for (int v = 1; v < 32; ++v) {
                const auto m = transpose(embed(j, v));
                for (int k = 0; k < 64; ++k) {
                    int u = column(m, k);
                    if (u) terms[j * 32 + v].push_back({k, u});
                }
            }
    }
    State transpose(const State &a) const {
        State b{};
        for (int i = 0; i < 5; ++i) {
            uint64_t x = a[i];
            while (x) {
                const int j = __builtin_ctzll(x);
                b = b ^ rows[64 * i + j];
                x &= x - 1;
            }
        }
        return b;
    }
};

struct Estimate {
    std::vector<Correlations> rounds; // index is the number of chi layers
    bool underflow = false;
};

class Estimator {
public:
    Tables tables;
    LinearMasks masks;

    Correlations apply_chi(const Correlations &x) const {
        Correlations y{};
        for (int j = 0; j < 64; ++j)
            for (int v = 0; v < 32; ++v)
                for (auto [u, weight] : tables.kernel[v]) y[j][v] += weight * x[j][u];
        return y;
    }
    Correlations apply_linear(const Correlations &x, bool &underflow) const {
        Correlations y{};
        for (int j = 0; j < 64; ++j) {
            y[j][0] = 1;
            for (int v = 1; v < 32; ++v) {
                Real value = 1;
                for (auto [k, u] : masks.terms[j * 32 + v]) {
                    const Real next = value * x[k][u];
                    if (next == 0 && value != 0 && x[k][u] != 0) underflow = true;
                    value = next;
                    if (value == 0) break;
                }
                y[j][v] = value;
            }
        }
        return y;
    }
    // HNW zero-value-mask optimization, rd=1, rl=0. Exact first-DDT
    // probabilities; later conditional value distributions and cross-column
    // dependencies are approximated. This is not an exact fixed-permutation hull.
    Estimate estimate(int input_col, int input_pattern, int max_rounds) const {
        std::feclearexcept(FE_UNDERFLOW);
        Estimate result;
        result.rounds.resize(max_rounds + 1);
        for (int j = 0; j < 64; ++j)
            for (int v = 0; v < 32; ++v)
                result.rounds[1][j][v] = j == input_col ? tables.dlct[input_pattern][v] : 1;
        if (max_rounds == 1) return result;
        for (int b = 0; b < 32; ++b) {
            if (!tables.ddt[input_pattern][b]) continue;
            const Real weight = Real(tables.ddt[input_pattern][b]) / 32;
            // L on a fixed difference is exact. Do not average DDT branches
            // before the subsequent nonlinear product approximation.
            const State d = linear(embed(input_col, b));
            Correlations x{};
            for (int j = 0; j < 64; ++j)
                for (int v = 0; v < 32; ++v) x[j][v] = parity(v & column(d, j)) ? -1 : 1;
            for (int r = 2; r <= max_rounds; ++r) {
                x = apply_chi(x);
                for (int j = 0; j < 64; ++j)
                    for (int v = 0; v < 32; ++v) result.rounds[r][j][v] += weight * x[j][v];
                if (r < max_rounds) x = apply_linear(x, result.underflow);
            }
        }
        result.underflow |= std::fetestexcept(FE_UNDERFLOW) != 0;
        return result;
    }
};
} // namespace gdl

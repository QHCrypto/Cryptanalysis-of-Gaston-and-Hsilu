#pragma once

// Exact evaluation of the independence factors in ordinary DATF Algorithm 1
// for one or two chi layers. This is not enhanced DATF. For two layers every
// first-cut beta is an affine function of the original uniform input, so its
// marginal is either uniform or deterministic. All remaining Bias calls are
// uniform quadratic sums, evaluated by elimination instead of enumeration.
#include "gaston.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

namespace gdl {
using Term = uint64_t;
using Poly = std::set<Term>;
inline Term monomial(int a, int b = -1) {
    if (a < 0) return b < 0 ? 0 : Term(b + 1);
    if (b < 0 || a == b) return Term(a + 1);
    if (a > b) std::swap(a, b);
    return (Term(a + 1) << 32) | Term(b + 1);
}
inline void toggle(Poly &p, Term t) {
    const auto it = p.find(t);
    if (it == p.end()) p.insert(t); else p.erase(it);
}
inline void add_poly(Poly &p, const Poly &q) { for (auto t : q) toggle(p, t); }
inline void multiply_affine(Poly &p, const Poly &a, const Poly &b) {
    for (auto x : a) for (auto y : b) {
        if ((x >> 32) || (y >> 32)) throw std::runtime_error("non-affine multiplication");
        toggle(p, monomial(int(x) - 1, int(y) - 1));
    }
}
struct CorrelationPower {
    int sign = 0, weight = 0;
    uint64_t rc_mask = 0;
};

inline CorrelationPower quadratic_sum(const Poly &p) {
    std::set<int> vars;
    for (auto t : p) if (t) {
        vars.insert(int(uint32_t(t)) - 1);
        if (t >> 32) vars.insert(int(t >> 32) - 1);
    }
    if (vars.empty()) return {p.count(0) ? -1 : 1, 0, 0};
    const int n = int(vars.size()), nw = (n + 63) / 64;
    std::map<int, int> ids;
    int idx = 0; for (int v : vars) ids[v] = idx++;
    std::vector<std::vector<uint64_t>> rows(n, std::vector<uint64_t>(nw));
    std::vector<unsigned char> lin(n), alive(n, 1);
    int constant = int(p.count(0)), w = 0;
    for (auto t : p) if (t) {
        int b = ids.at(int(uint32_t(t)) - 1);
        if (!(t >> 32)) lin[b] ^= 1;
        else {
            int a = ids.at(int(t >> 32) - 1);
            rows[a][b / 64] ^= uint64_t(1) << (b % 64);
            rows[b][a / 64] ^= uint64_t(1) << (a % 64);
        }
    }
    while (true) {
        int u = -1, v = -1;
        for (int i = 0; i < n && u < 0; ++i) if (alive[i]) {
            bool connected = false;
            for (int k = 0; k < nw; ++k) if (rows[i][k]) {
                u = i; v = 64 * k + __builtin_ctzll(rows[i][k]); connected = true; break;
            }
            if (!connected && lin[i]) return {};
        }
        if (u < 0) {
            for (int i = 0; i < n; ++i) if (alive[i] && lin[i]) return {};
            return {constant ? -1 : 1, w, 0};
        }
        auto a = rows[u], b = rows[v];
        a[v / 64] &= ~(uint64_t(1) << (v % 64));
        b[u / 64] &= ~(uint64_t(1) << (u % 64));
        const int la = lin[u], lb = lin[v];
        constant ^= la & lb;
        alive[u] = alive[v] = 0;
        for (int j = 0; j < n; ++j) if (alive[j]) {
            const int aj = int((a[j / 64] >> (j % 64)) & 1);
            const int bj = int((b[j / 64] >> (j % 64)) & 1);
            lin[j] ^= (la & bj) ^ (lb & aj) ^ (aj & bj);
            rows[j][u / 64] &= ~(uint64_t(1) << (u % 64));
            rows[j][v / 64] &= ~(uint64_t(1) << (v % 64));
            if (aj) for (int k = 0; k < nw; ++k) rows[j][k] ^= b[k];
            if (bj) for (int k = 0; k < nw; ++k) rows[j][k] ^= a[k];
            rows[j][j / 64] &= ~(uint64_t(1) << (j % 64));
        }
        std::fill(rows[u].begin(), rows[u].end(), 0);
        std::fill(rows[v].begin(), rows[v].end(), 0);
        ++w;
    }
}

inline CorrelationPower middle_one(const State &delta, const State &gamma) {
    for (int i = 0; i < 5; ++i)
        if ((gamma[(i + 4) % 5] & delta[(i + 1) % 5]) ^
            (gamma[(i + 3) % 5] & delta[(i + 4) % 5])) return {};
    return {dot(gamma, chi(delta)) ? -1 : 1, 0, 0};
}

struct MiddleTables {
    std::array<std::vector<int>, 320> forward;
    std::array<Poly, 320> alpha_def;
    MiddleTables() {
        for (int i = 0; i < 320; ++i) {
            const auto v = linear(embed(i % 64, 1 << (i / 64)));
            for (int row = 0; row < 5; ++row) {
                uint64_t bits = v[row];
                while (bits) {
                    int col = __builtin_ctzll(bits); bits &= bits - 1;
                    forward[i].push_back(64 * row + col);
                }
            }
            const int row = i / 64, col = i % 64;
            Poly p{monomial(i), monomial(((row + 2) % 5) * 64 + col),
                   monomial(((row + 1) % 5) * 64 + col, ((row + 2) % 5) * 64 + col)};
            for (int out : forward[i]) add_poly(alpha_def[out], p);
        }
    }
};

class MiddleTwo {
    const MiddleTables &tables;
    std::array<Poly, 320> difference;
    std::vector<Poly> betas;
public:
    MiddleTwo(const State &delta, const MiddleTables &t) : tables(t) {
        std::array<Poly, 320> affine;
        for (int row = 0; row < 5; ++row) for (int col = 0; col < 64; ++col) {
            const int i = 64 * row + col, k = (row + 1) % 5, l = (row + 2) % 5;
            const int di = int((delta[row] >> col) & 1), dk = int((delta[k] >> col) & 1), dl = int((delta[l] >> col) & 1);
            Poly p;
            if (di ^ dl ^ (dk & dl)) p.insert(0);
            if (dl) p.insert(monomial(k * 64 + col));
            if (dk) p.insert(monomial(l * 64 + col));
            for (int out : tables.forward[i]) add_poly(affine[out], p);
        }
        std::map<Poly, int> lookup;
        for (int i = 0; i < 320; ++i) {
            Poly nonconstant = affine[i];
            const bool c = nonconstant.erase(0);
            if (nonconstant.size() <= 1) difference[i] = affine[i];
            else {
                auto it = lookup.find(nonconstant);
                if (it == lookup.end()) {
                    int id = 640 + int(betas.size());
                    lookup[nonconstant] = id; betas.push_back(nonconstant);
                    difference[i].insert(monomial(id));
                } else difference[i].insert(monomial(it->second));
                if (c) difference[i].insert(0);
            }
        }
    }

    CorrelationPower evaluate(const State &gamma) const {
        Poly f;
        for (int row = 0; row < 5; ++row) {
            uint64_t bits = gamma[row];
            while (bits) {
                const int col = __builtin_ctzll(bits); bits &= bits - 1;
                const int i = row * 64 + col, k = ((row + 1) % 5) * 64 + col, l = ((row + 2) % 5) * 64 + col;
                add_poly(f, difference[i]); add_poly(f, difference[l]);
                multiply_affine(f, Poly{monomial(320 + k)}, difference[l]);
                multiply_affine(f, difference[k], Poly{monomial(320 + l)});
                multiply_affine(f, difference[k], difference[l]);
            }
        }
        std::set<int> linear_vars, nonlinear_vars;
        for (auto term : f) if (term) {
            if (term >> 32) {
                nonlinear_vars.insert(int(term >> 32) - 1);
                nonlinear_vars.insert(int(uint32_t(term)) - 1);
            } else linear_vars.insert(int(term) - 1);
        }
        Poly remainder;
        uint64_t rc_mask = 0;
        for (int v : linear_vars) if (!nonlinear_vars.count(v)) {
            f.erase(monomial(v));
            if (v < 320) return {}; // Original Algorithm 1's isolated uniform input.
            if (v < 640) {
                add_poly(remainder, tables.alpha_def[v - 320]);
                if (v < 384) rc_mask ^= uint64_t(1) << (v - 320);
            } else add_poly(remainder, betas.at(v - 640));
        }
        auto a = quadratic_sum(f);
        if (!a.sign) return {};
        auto b = quadratic_sum(remainder);
        if (!b.sign) return {};
        return {a.sign * b.sign, a.weight + b.weight, rc_mask};
    }
};
} // namespace gdl

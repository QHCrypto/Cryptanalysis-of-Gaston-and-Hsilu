#include "experiment.hpp"
#include <iostream>
#include <set>
extern "C" {
#include "gaston.h"
void gaston_round(gaston_state_t *, uint64_t);
}
using namespace gdl;

static void require(bool ok, const char *message) {
    if (!ok) throw std::runtime_error(message);
}
static int rank(std::vector<State> rows) {
    int r = 0;
    for (int bit = 0; bit < 320; ++bit) {
        size_t p = r;
        while (p < rows.size() && !((rows[p][bit / 64] >> (bit % 64)) & 1)) ++p;
        if (p == rows.size()) continue;
        std::swap(rows[p], rows[r]);
        for (size_t j = r + 1; j < rows.size(); ++j)
            if ((rows[j][bit / 64] >> (bit % 64)) & 1) rows[j] = rows[j] ^ rows[r];
        ++r;
    }
    return r;
}
static void test_reference() {
    for (int k = 0; k < 128; ++k) {
        const State a = sample_state(101, k);
        for (int r = 0; r < 12; ++r) {
            gaston_state_t b{};
            for (int i = 0; i < 5; ++i) b.x[i] = a[i];
            gaston_round(&b, RC[r]);
            const State out = round(a, RC[r]);
            for (int i = 0; i < 5; ++i) require(out[i] == b.x[i], "round differs from reference");
        }
        for (int rounds = 1; rounds <= 12; ++rounds) {
            const int start = 12 - rounds;
            gaston_state_t b{};
            for (int i = 0; i < 5; ++i) b.x[i] = a[i];
            for (int r = start; r < 12; ++r) gaston_round(&b, RC[r]);
            State first_input = linear(a);
            first_input[0] ^= RC[start];
            const State out = core(first_input, rounds, start);
            for (int i = 0; i < 5; ++i) require(out[i] == b.x[i], "first-chi boundary differs from reference");
        }
    }
    std::cout << "PASS reference round and omitted-first-linear boundary (all round counts)\n";
}
static void test_tables(const Tables &t) {
    std::set<unsigned> outputs;
    for (unsigned x = 0; x < 32; ++x) {
        outputs.insert(chi5(x));
        for (int j = 0; j < 64; ++j) require(column(chi(embed(j, x)), j) == chi5(x), "chi5 bit order");
    }
    require(outputs.size() == 32, "chi is not a permutation");
    for (int a = 0; a < 32; ++a) {
        int total = 0;
        for (int b = 0; b < 32; ++b) total += t.ddt[a][b];
        require(total == 32, "DDT normalization");
        for (int v = 0; v < 32; ++v) {
            int exact = 0;
            for (int x = 0; x < 32; ++x) exact += parity(v & (chi5(x) ^ chi5(x ^ a))) ? -1 : 1;
            Real transformed = 0;
            for (auto [u, weight] : t.kernel[v]) transformed += weight * (parity(u & a) ? -1 : 1);
            require(transformed == Real(exact) / 32 && transformed == t.dlct[a][v], "LAT-squared/DDT/direct-derivative mismatch");
        }
    }
    std::cout << "PASS chi/DDT/LAT/direct derivative identities (all local endpoints)\n";
}
static void test_linear(const LinearMasks &m) {
    require(rank({m.rows.begin(), m.rows.end()}) == 320, "linear layer is not invertible");
    for (int k = 0; k < 256; ++k) {
        const State x = sample_state(201, 2 * k), v = sample_state(201, 2 * k + 1);
        require(dot(linear(x), v) == dot(x, m.transpose(v)), "transpose identity");
        require(linear(rotate(x, k % 64)) == rotate(linear(x), k % 64), "linear rotation equivariance");
        require(chi(rotate(x, k % 64)) == rotate(chi(x), k % 64), "chi rotation equivariance");
    }
    // The second-chi target column is uniform conditioned on the entire first
    // active column if the map from all OTHER columns has rank 5. This is what
    // makes the two-round single-output-box estimate exact for this L.
    for (int dc = 0; dc < 64; ++dc)
        for (int mc = 0; mc < 64; ++mc) {
            std::vector<State> rows;
            for (int i = 0; i < 5; ++i) {
                State row = m.rows[64 * i + mc];
                for (auto &word : row) word &= ~(uint64_t(1) << dc);
                rows.push_back(row);
            }
            require(rank(rows) == 5, "two-round uniformity rank condition failed");
        }
    std::cout << "PASS rank(L)=320, transpose, rotations, all 4096 two-round conditional rank checks\n";
}
static void test_estimator(const Estimator &e) {
    for (int a = 1; a < 32; ++a) {
        const auto r = e.estimate(0, a, 3);
        require(!r.underflow, "unexpected short-round underflow");
        for (int j = 0; j < 64; ++j)
            for (int v = 1; v < 32; ++v) {
                // Independent enumeration of first input and second local base
                // value, justified by the rank condition checked above.
                int exact = 0;
                for (int x = 0; x < 32; ++x) {
                    const int d = column(linear(embed(0, chi5(x) ^ chi5(x ^ a))), j);
                    for (int y = 0; y < 32; ++y) exact += parity(v & (chi5(y) ^ chi5(y ^ d))) ? -1 : 1;
                }
                require(r.rounds[2][j][v] == Real(exact) / 1024, "two-round exhaustive endpoint mismatch");
            }
    }
    for (int a : {1, 3, 17, 31}) {
        const auto base = e.estimate(0, a, 3);
        for (int shift : {1, 7, 63}) {
            const auto shifted = e.estimate(shift, a, 3);
            for (int r = 1; r <= 3; ++r)
                for (int j = 0; j < 64; ++j)
                    for (int v = 0; v < 32; ++v)
                        require(std::abs(base.rounds[r][j][v] - shifted.rounds[r][(j + shift) % 64][v]) < 1e-18L,
                                "model symmetry reduction mismatch");
        }
    }
    std::cout << "PASS two-round exhaustive local enumeration and rotated three-round estimates\n";
}
static void test_sampler() {
    const auto a = measure(0, 1, 0, 1, 3, 9, true, 4097, 123, 1);
    const auto b = measure(0, 1, 0, 1, 3, 9, true, 4097, 123, 4);
    require(a.n0 == b.n0 && a.n1 == b.n1 && a.n0 + a.n1 == 4097, "thread-count dependent sampling");
    const auto positive = measure(0, 1, 0, 1, 1, 11, true, 1024, 1, 3);
    require(positive.correlation == -1, "one-round known negative derivative");
    require(positive.low == -1 && positive.high > -1, "boundary confidence interval");
    std::cout << "PASS reproducible sample partitioning, signed correlation, boundary Wilson interval\n";
}
static void test_model_zero_mechanism(const Estimator &e) {
    uint64_t zero_products = 0, products = 0;
    for (int b = 1; b < 32; ++b) {
        const State delta = linear(embed(0, b));
        Correlations x{};
        for (int j = 0; j < 64; ++j)
            for (int v = 0; v < 32; ++v) x[j][v] = parity(column(delta, j) & v) ? -1 : 1;
        x = e.apply_chi(x); // second chi
        bool underflow = false;
        x = e.apply_linear(x, underflow);
        x = e.apply_chi(x); // third chi
        for (const auto &local : x)
            for (Real value : local) require(16 * value == std::round(16 * value), "unexpected denominator before fourth-round collapse");
        for (int j = 0; j < 64; ++j)
            for (int v = 1; v < 32; ++v) {
                bool zero_factor = false;
                for (auto [k, u] : e.masks.terms[j * 32 + v]) zero_factor |= x[k][u] == 0;
                zero_products += zero_factor;
                ++products;
            }
        require(!underflow, "model zero diagnosis underflow");
    }
    require(zero_products == products, "not all fourth-round linear products have exact zero factors");
    std::cout << "PASS fourth-round model collapse: all " << products
              << " nontrivial linear products have exact zero factors (not underflow)\n";
}
int main() {
    try {
        test_reference();
        Estimator e;
        test_tables(e.tables);
        test_linear(e.masks);
        test_estimator(e);
        test_sampler();
        test_model_zero_mechanism(e);
        std::cout << "All tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}

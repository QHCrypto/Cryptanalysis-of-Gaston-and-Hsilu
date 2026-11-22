#pragma once

#include "estimator.hpp"
#include <thread>

namespace gdl {
struct Measurement {
    uint64_t n0 = 0, n1 = 0;
    Real correlation = 0, low = 0, high = 0;
};
inline Measurement measure_states(const State &d, const State &mask, int rounds, int start,
                                  bool constants, uint64_t samples, uint64_t seed, int threads) {
    if (samples == 0 || threads < 1 || rounds < 1 || start < 0 || start + rounds > 12)
        throw std::runtime_error("invalid sampling configuration");
    std::vector<uint64_t> counts(threads);
    std::vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t] {
            const uint64_t q = samples / threads, rem = samples % threads;
            const uint64_t begin = q * t + std::min<uint64_t>(t, rem);
            const uint64_t end = begin + q + (uint64_t(t) < rem);
            uint64_t n0 = 0;
            for (uint64_t k = begin; k < end; ++k) {
                const State a = sample_state(seed, k);
                const State y = core(a, rounds, start, constants) ^ core(a ^ d, rounds, start, constants);
                n0 += !dot(y, mask);
            }
            counts[t] = n0;
        });
    }
    for (auto &w : workers) w.join();
    Measurement m;
    for (auto n : counts) m.n0 += n;
    m.n1 = samples - m.n0;
    // Subtract in floating point to avoid unsigned integer wraparound.
    m.correlation = (Real(m.n0) - Real(m.n1)) / Real(samples);
    // Two-sided 95% Wilson interval for p=Pr[parity=0], mapped by rho=2p-1.
    // Unlike a plug-in standard error, this is non-degenerate at rho=+/-1.
    const Real z = 1.959963984540054L, n = Real(samples), p = Real(m.n0) / n;
    const Real den = 1 + z * z / n;
    const Real center = (p + z * z / (2 * n)) / den;
    const Real radius = z * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den;
    m.low = std::max(Real(-1), 2 * (center - radius) - 1);
    m.high = std::min(Real(1), 2 * (center + radius) - 1);
    return m;
}
inline Measurement measure(int dc, int dp, int mc, int mp, int rounds, int start,
                           bool constants, uint64_t samples, uint64_t seed, int threads) {
    return measure_states(embed(dc, dp), embed(mc, mp), rounds, start,
                          constants, samples, seed, threads);
}
} // namespace gdl

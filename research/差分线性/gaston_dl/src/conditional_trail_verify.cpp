#include "experiment.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>

using namespace gdl;

struct Counts {
    uint64_t positive = 0, negative = 0, bad_first = 0, bad_middle = 0;
};

int main(int argc, char **argv) {
    try {
        if (argc != 6) {
            std::fprintf(stderr, "Usage: %s jobs.txt output.csv conditional_samples seed threads\n", argv[0]);
            return 1;
        }
        const uint64_t samples = std::stoull(argv[3]), seed = std::stoull(argv[4]);
        const int threads = std::stoi(argv[5]);
        if (!samples || threads < 1) throw std::runtime_error("invalid sampling parameters");
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("output already exists");
        std::ifstream jobs(argv[1]);
        if (!jobs) throw std::runtime_error("cannot open jobs file");
        std::ofstream csv(argv[2]);
        if (!csv) throw std::runtime_error("cannot create output");
        csv << "id,start,conditional_samples,seed,threads,p_first_exact,p_second_theory,second_event_hits,positive_hits,negative_hits,p_second_estimate,ci95_low,ci95_high,theory_inside_ci95,bad_first,bad_middle,input_words,first_output_words,second_output_words,join_words\n";
        csv << std::setprecision(18);
        Tables tables;
        std::string line;
        unsigned index = 0;
        while (std::getline(jobs, line)) {
            if (line.empty()) continue;
            std::istringstream in(line);
            std::string id, word;
            int start, expected_sign;
            State d0{}, b1{}, b2{}, gamma{};
            if (!(in >> id >> start >> expected_sign) || start < 0 || start + 3 > 12 || std::abs(expected_sign) != 1)
                throw std::runtime_error("invalid job header");
            for (auto *state : {&d0, &b1, &b2, &gamma}) for (auto &v : *state) {
                if (!(in >> word)) throw std::runtime_error("missing word");
                size_t n = 0;
                v = std::stoull(word, &n, 16);
                if (n != word.size()) throw std::runtime_error("invalid word");
            }
            if (in >> word) throw std::runtime_error("extra job field");
            std::vector<int> positions;
            std::vector<std::vector<unsigned>> allowed;
            Real p1 = 1, p2 = 1, middle = 1;
            const State a2 = linear(b1), a3 = linear(b2);
            for (int j = 0; j < 64; ++j) {
                const unsigned a = column(d0, j), b = column(b1, j);
                p1 *= Real(tables.ddt[a][b]) / 32;
                p2 *= Real(tables.ddt[column(a2, j)][column(b2, j)]) / 32;
                middle *= tables.dlct[column(a3, j)][column(gamma, j)];
                if (a) {
                    std::vector<unsigned> values;
                    for (unsigned x = 0; x < 32; ++x) if ((chi5(x) ^ chi5(x ^ a)) == b) values.push_back(x);
                    if (values.empty() || (values.size() & (values.size() - 1)))
                        throw std::runtime_error("conditional local sample count must be a power of two");
                    positions.push_back(j);
                    allowed.push_back(values);
                }
            }
            if (!p1 || !p2 || middle != expected_sign) throw std::runtime_error("invalid DDT or middle connection");
            const uint64_t run_seed = seed + index++;
            std::vector<Counts> counts(threads);
            std::vector<std::thread> workers;
            for (int t = 0; t < threads; ++t) workers.emplace_back([&, t] {
                const uint64_t q = samples / threads, rem = samples % threads;
                const uint64_t begin = q * t + std::min<uint64_t>(t, rem);
                const uint64_t end = begin + q + (uint64_t(t) < rem);
                Counts result;
                for (uint64_t k = begin; k < end; ++k) {
                    const auto draw = [&](uint64_t offset) {
                        return splitmix_final(run_seed + ((5 + positions.size()) * k + offset + 1) * UINT64_C(0x9e3779b97f4a7c15));
                    };
                    State x{};
                    for (int i = 0; i < 5; ++i) x[i] = draw(i);
                    for (size_t c = 0; c < positions.size(); ++c) {
                        const unsigned v = allowed[c][draw(5 + c) & (allowed[c].size() - 1)];
                        const uint64_t bit = UINT64_C(1) << positions[c];
                        for (int i = 0; i < 5; ++i) x[i] = (x[i] & ~bit) | (((v >> i) & 1) ? bit : 0);
                    }
                    State y = chi(x), other = chi(x ^ d0);
                    if ((y ^ other) != b1) { ++result.bad_first; continue; }
                    y = round(y, RC[start + 1]);
                    other = round(other, RC[start + 1]);
                    if ((y ^ other) != b2) continue;
                    y = round(y, RC[start + 2]);
                    other = round(other, RC[start + 2]);
                    const int sign = dot(y ^ other, gamma) ? -1 : 1;
                    if (sign == 1) ++result.positive; else ++result.negative;
                    if (sign != expected_sign) ++result.bad_middle;
                }
                counts[t] = result;
            });
            for (auto &worker : workers) worker.join();
            Counts total;
            for (const auto &c : counts) {
                total.positive += c.positive; total.negative += c.negative;
                total.bad_first += c.bad_first; total.bad_middle += c.bad_middle;
            }
            const uint64_t hits = total.positive + total.negative;
            const Real p = Real(hits) / samples, z = 1.959963984540054L, n = Real(samples);
            const Real den = 1 + z * z / n, center = (p + z * z / (2 * n)) / den;
            const Real radius = z * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den;
            const Real low = std::max(Real(0), center - radius), high = std::min(Real(1), center + radius);
            csv << id << ',' << start << ',' << samples << ',' << run_seed << ',' << threads << ','
                << p1 << ',' << p2 << ',' << hits << ',' << total.positive << ',' << total.negative << ','
                << p << ',' << low << ',' << high << ',' << (low <= p2 && p2 <= high) << ','
                << total.bad_first << ',' << total.bad_middle << ',' << hex_state(d0) << ',' << hex_state(b1) << ','
                << hex_state(b2) << ',' << hex_state(gamma) << '\n';
            csv.flush();
            std::printf("%s conditional hits=%llu/%llu p2=%.9Lf theory=%.9Lf signs(+,-)=(%llu,%llu) bad=(%llu,%llu)\n",
                        id.c_str(), (unsigned long long)hits, (unsigned long long)samples, p, p2,
                        (unsigned long long)total.positive, (unsigned long long)total.negative,
                        (unsigned long long)total.bad_first, (unsigned long long)total.bad_middle);
            std::fflush(stdout);
            if (total.bad_first || total.bad_middle) throw std::runtime_error("conditional trail verification failed");
        }
        if (!index || !csv) throw std::runtime_error("no jobs or failed CSV write");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
}

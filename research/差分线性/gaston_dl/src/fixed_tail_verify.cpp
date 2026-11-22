#include "experiment.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>

// Jobs: id rounds start exact_numerator exact_denominator delta[5] mask[5].
// The mask is a complete 320-bit mask, not a single active output column.
int main(int argc, char **argv) {
    try {
        if (argc != 6) {
            std::fprintf(stderr, "Usage: %s jobs.txt output.csv samples seed threads\n", argv[0]);
            return 1;
        }
        const uint64_t samples = std::stoull(argv[3]), seed = std::stoull(argv[4]);
        const int threads = std::stoi(argv[5]);
        if (!samples || threads < 1) throw std::runtime_error("invalid sampling parameters");
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("output already exists");
        std::ifstream jobs(argv[1]);
        if (!jobs) throw std::runtime_error("cannot open jobs file");
        std::ofstream csv(argv[2]);
        if (!csv) throw std::runtime_error("cannot create output file");
        csv << "id,rounds,start,exact_numerator,exact_denominator,samples,seed,threads,n0,n1,rho,ci95_low,ci95_high,exact_inside_ci95,input_words,mask_words\n";
        csv << std::setprecision(18);
        std::string line;
        unsigned index = 0;
        while (std::getline(jobs, line)) {
            if (line.empty()) continue;
            std::istringstream in(line);
            std::string id, word;
            int rounds, start;
            int64_t numerator;
            uint64_t denominator;
            gdl::State delta{}, mask{};
            if (!(in >> id >> rounds >> start >> numerator >> denominator) || !denominator)
                throw std::runtime_error("invalid job header");
            for (auto *state : {&delta, &mask}) for (auto &v : *state) {
                if (!(in >> word)) throw std::runtime_error("missing state word");
                size_t n = 0;
                v = std::stoull(word, &n, 16);
                if (n != word.size()) throw std::runtime_error("invalid state word");
            }
            if (in >> word) throw std::runtime_error("extra job field");
            const uint64_t run_seed = seed + index++;
            const auto m = gdl::measure_states(delta, mask, rounds, start, true, samples, run_seed, threads);
            const auto expected = gdl::Real(numerator) / denominator;
            const bool covered = m.low <= expected && expected <= m.high;
            csv << id << ',' << rounds << ',' << start << ',' << numerator << ',' << denominator << ','
                << samples << ',' << run_seed << ',' << threads << ',' << m.n0 << ',' << m.n1 << ','
                << m.correlation << ',' << m.low << ',' << m.high << ',' << covered << ','
                << gdl::hex_state(delta) << ',' << gdl::hex_state(mask) << '\n';
            csv.flush();
            std::printf("%s prefix rho=%.9Lf exact=%.9Lf CI=[%.9Lf,%.9Lf] covered=%d\n",
                        id.c_str(), m.correlation, expected, m.low, m.high, covered);
            std::fflush(stdout);
        }
        if (!index) throw std::runtime_error("empty jobs file");
        if (!csv) throw std::runtime_error("failed to write output");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
}

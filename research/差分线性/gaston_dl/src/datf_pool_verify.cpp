#include "experiment.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

// Jobs: id rounds start model_numerator model_denominator delta[5] mask[5].
// Complete endpoint sampling, using the same checked core as other experiments.
int main(int argc, char **argv) {
    try {
        if (argc != 6) throw std::runtime_error("usage: datf_pool_verify jobs output samples seed threads");
        const uint64_t samples = std::stoull(argv[3]), seed = std::stoull(argv[4]);
        const int threads = std::stoi(argv[5]);
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("output exists");
        std::ifstream in(argv[1]);
        std::ofstream out(argv[2]);
        if (!in || !out) throw std::runtime_error("cannot open jobs/output");
        out << "id,rounds,start,model_numerator,model_denominator,samples,seed,threads,n0,n1,rho,ci95_low,ci95_high,model_inside_ci95,input_words,mask_words\n";
        out << std::setprecision(18);
        std::string line;
        uint64_t index = 0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream f(line);
            std::string id, extra;
            int rounds, start;
            int64_t numerator;
            uint64_t denominator;
            gdl::State delta{}, mask{};
            if (!(f >> id >> rounds >> start >> numerator >> denominator) || !denominator)
                throw std::runtime_error("invalid job");
            for (auto *state : {&delta, &mask}) for (auto &w : *state)
                if (!(f >> std::hex >> w)) throw std::runtime_error("invalid state");
            if (f >> extra) throw std::runtime_error("extra job fields");
            const uint64_t run_seed = seed + index++;
            const auto m = gdl::measure_states(delta, mask, rounds, start, true, samples, run_seed, threads);
            const auto model = gdl::Real(numerator)/denominator;
            const bool covered = m.low <= model && model <= m.high;
            out << id << ',' << rounds << ',' << start << ',' << numerator << ',' << denominator << ','
                << samples << ',' << run_seed << ',' << threads << ',' << m.n0 << ',' << m.n1 << ','
                << m.correlation << ',' << m.low << ',' << m.high << ',' << covered << ','
                << gdl::hex_state(delta) << ',' << gdl::hex_state(mask) << '\n';
            out.flush();
            std::cout << id << " rho=" << m.correlation << " CI95=[" << m.low << ',' << m.high
                      << "] model_inside=" << covered << std::endl;
        }
        if (!index || !out) throw std::runtime_error("empty jobs or output failure");
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

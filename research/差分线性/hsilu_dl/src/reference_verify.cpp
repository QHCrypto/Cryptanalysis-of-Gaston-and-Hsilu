// Independent full-permutation sampling: reference hsilu_round, dense output
// masks, and std::mt19937_64 instead of the core sampler's counter generator.
#include "hsilu.h"
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
void hsilu_round(hsilu_state_t*, uint64_t);
using State = std::array<uint64_t, 5>;
constexpr uint64_t constants[12] = {0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b};
State parse(const std::string& words) {
    std::stringstream stream(words); std::string word; State a{};
    for (auto& x : a) {
        if (!std::getline(stream, word, ';')) throw std::runtime_error("five words required");
        x = std::stoull(word, nullptr, 16);
    }
    return a;
}
State reference(State a, int rounds) {
    hsilu_state_t s{};
    for (int i=0; i<5; ++i) s.x[i] = a[i];
    for (int r=12-rounds; r<12; ++r) hsilu_round(&s, constants[r]);
    for (int i=0; i<5; ++i) a[i] = s.x[i];
    return a;
}
int main(int argc, char** argv) try {
    if (argc != 5) throw std::runtime_error("reference_verify JOBS OUTPUT SAMPLES SEED");
    if (std::filesystem::exists(argv[2])) throw std::runtime_error("output exists");
    const uint64_t n = std::stoull(argv[3]), seed = std::stoull(argv[4]);
    if (!n) throw std::runtime_error("empty sample");
    std::ifstream jobs(argv[1]); std::ofstream output(argv[2]);
    if (!jobs || !output) throw std::runtime_error("file error");
    output << "id,rounds,samples,seed,n0,n1,rho\n" << std::setprecision(20);
    std::string id, d, m; int rounds; unsigned count = 0;
    while (jobs >> id >> rounds >> d >> m) {
        if (rounds < 1 || rounds > 12) throw std::runtime_error("invalid round count");
        const auto delta = parse(d), mask = parse(m);
        std::mt19937_64 rng(seed); uint64_t n0 = 0;
        for (uint64_t k=0; k<n; ++k) {
            State a{}, b{};
            for (int i=0; i<5; ++i) {a[i]=rng(); b[i]=a[i]^delta[i];}
            a=reference(a,rounds); b=reference(b,rounds);
            uint64_t parity_word = 0;
            for (int i=0; i<5; ++i) parity_word ^= (a[i]^b[i])&mask[i];
            n0 += !__builtin_parityll(parity_word);
        }
        output << id << ',' << rounds << ',' << n << ',' << seed << ',' << n0 << ',' << n-n0
               << ',' << (2*static_cast<long double>(n0)-n)/n << '\n';
        output.flush(); ++count;
        std::cerr << "REFERENCE VERIFIED " << id << '\n';
    }
    if (!jobs.eof() || !count || !output) throw std::runtime_error("incomplete jobs/output");
} catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}

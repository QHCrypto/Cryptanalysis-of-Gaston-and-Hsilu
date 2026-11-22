// Exact interval sums for a rotation-invariant Gaston affine nonce cube.
// Links the original Gaston C round implementation. An incomplete interval
// is only a partial XOR sum and is never reported as a rank certificate.
#include <array>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>

extern "C" {
#include "gaston.h"
void gaston_round(gaston_state_t*, uint64_t);
}

constexpr uint64_t IV = UINT64_C(0x00001000808c0001);
constexpr uint64_t RC[12] = {0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b};

uint64_t number(const char* s) {
    if (!s[0] || s[0]=='-') throw std::runtime_error("expected an unsigned integer");
    size_t used = 0;
    auto v = std::stoull(s, &used, 0);
    if (used != std::string(s).size()) throw std::runtime_error("invalid integer");
    return v;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::fprintf(stderr, "usage: %s rounds begin count output.json\n", argv[0]);
        return 2;
    }
    try {
        auto rounds = number(argv[1]);
        if (rounds < 1 || rounds > 6) throw std::runtime_error("rounds must be 1..6");
        unsigned period = 1u << rounds;
        uint64_t total = UINT64_C(1) << (period - 1);
        uint64_t begin = number(argv[2]), count = number(argv[3]);
        if (begin >= total || count == 0 || count > total - begin)
            throw std::runtime_error("interval outside the nonce cube");
        std::array<std::array<uint64_t,2>,3> sums{};
        auto start = std::chrono::steady_clock::now();
        for (uint64_t index=begin; index < begin+count; ++index) {
            uint64_t word = (index << 1) | uint64_t(__builtin_parityll(index));
            uint64_t nonce = 0;
            for (unsigned offset=0; offset<64; offset+=period) nonce |= word << offset;
            for (unsigned k=0; k<3; ++k) {
                gaston_state_t state{{IV, uint64_t(k==1), uint64_t(k==2), nonce, 0}};
                for (unsigned r=12-rounds; r<12; ++r) gaston_round(&state, RC[r]);
                for (unsigned w=0; w<2; ++w) sums[k][w] ^= state.x[w];
            }
        }
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        FILE* out=std::fopen(argv[4],"w");
        if (!out) throw std::runtime_error("cannot open output file");
        std::fprintf(out,"{\n  \"implementation\": \"original_gaston_c\",\n  \"word_bits\": 64,\n  \"rounds\": %" PRIu64 ",\n  \"period\": %u,\n  \"cube_dimension\": %u,\n",rounds,period,period-1);
        std::fprintf(out,"  \"iv_hex\": \"%016" PRIx64 "\",\n  \"begin\": %" PRIu64 ",\n  \"count\": %" PRIu64 ",\n  \"total_points\": %" PRIu64 ",\n",IV,begin,count,total);
        std::fprintf(out,"  \"interval_is_entire_cube\": %s,\n  \"seconds\": %.9f,\n  \"sums_hex\": [\n",begin==0&&count==total?"true":"false",seconds);
        for (unsigned k=0;k<3;++k)
            std::fprintf(out,"    [\"%016" PRIx64 "\", \"%016" PRIx64 "\"]%s\n",sums[k][0],sums[k][1],k==2?"":",");
        std::fprintf(out,"  ],\n  \"scope\": \"Interval XOR sums for keys 0,e_(0,0),e_(1,0); no full-cube coefficient follows from an incomplete interval.\"\n}\n");
        std::fclose(out);
        std::printf("rounds=%" PRIu64 " dimension=%u points=%" PRIu64 "/%" PRIu64 " seconds=%.6f complete=%s\n",rounds,period-1,count,total,seconds,begin==0&&count==total?"yes":"no");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"Error: %s\n",e.what());
        return 1;
    }
}

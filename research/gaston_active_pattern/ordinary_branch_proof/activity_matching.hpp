#pragma once
#include <array>
#include <cstdint>
struct Matcher {
    std::array<uint64_t,64> edges{};
    std::array<int,64> mate;
    uint64_t seen=0;
    bool augment(int x) {
        auto choices=edges[x]&~seen;
        while(choices) {
            const int y=__builtin_ctzll(choices);choices&=choices-1;
            if(seen>>y&1)continue;
            seen|=1ULL<<y;
            if(mate[y]<0 || augment(mate[y])){mate[y]=x;return true;}
        }
        return false;
    }
    unsigned lower(uint64_t e) {
        edges.fill(0);mate.fill(-1);
        constexpr int east[]={0,60,22,27,4},west[]={0,56,31,46,43};
        while(e) {
            int j=__builtin_ctzll(e);e&=e-1;
            for(int r=0;r<5;++r)
                edges[(j+23-east[r]+64)%64]|=1ULL<<((j+23+west[r])%64);
        }
        unsigned n=0;
        for(int i=0;i<64;++i)if(edges[i]) {
            seen=0;if(augment(i) && ++n>15)return n;
        }
        return n;
    }
};

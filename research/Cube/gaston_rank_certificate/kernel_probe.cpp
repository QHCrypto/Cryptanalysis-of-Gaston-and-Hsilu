// Exact cube sums for user-specified key directions on a coordinate facet.
#include <array>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
extern "C" {
#include "gaston.h"
void gaston_round(gaston_state_t*, uint64_t);
}
constexpr uint64_t RC[12]={0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b};
uint64_t number(const char* s) {
    if(!s[0]||s[0]=='-')throw std::runtime_error("unsigned integer required");
    size_t used=0;auto v=std::stoull(s,&used,0);
    if(used!=std::string(s).size())throw std::runtime_error("invalid integer");
    return v;
}
int main(int argc,char** argv) {
    if(argc!=6){std::fprintf(stderr,"usage: %s rounds begin count keys.txt output.json\n",argv[0]);return 2;}
    try {
        auto rounds=number(argv[1]);
        if(rounds<1||rounds>6)throw std::runtime_error("rounds must be 1..6");
        unsigned period=1u<<rounds;
        uint64_t total=UINT64_C(1)<<(period-1),begin=number(argv[2]),count=number(argv[3]);
        if(begin>=total||!count||count>total-begin)throw std::runtime_error("invalid interval");
        std::ifstream input(argv[4]);
        std::vector<std::array<uint64_t,2>> keys;
        std::array<uint64_t,2> key;
        while(input>>std::hex>>key[0]>>key[1])keys.push_back(key);
        if(keys.empty()||!input.eof())throw std::runtime_error("invalid key file");
        std::vector<std::array<uint64_t,2>> sums(keys.size());
        auto start=std::chrono::steady_clock::now();
        for(uint64_t i=begin;i<begin+count;++i) {
            uint64_t word=i<<1,nonce=0;
            for(unsigned offset=0;offset<64;offset+=period)nonce|=word<<offset;
            for(size_t k=0;k<keys.size();++k) {
                gaston_state_t state{{UINT64_C(0x00001000808c0001),keys[k][0],keys[k][1],nonce,0}};
                for(unsigned r=12-rounds;r<12;++r)gaston_round(&state,RC[r]);
                for(unsigned w=0;w<2;++w)sums[k][w]^=state.x[w];
            }
        }
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        FILE* out=std::fopen(argv[5],"w");if(!out)throw std::runtime_error("cannot open output");
        std::fprintf(out,"{\n  \"implementation\": \"original_gaston_c\",\n  \"word_bits\": 64,\n  \"rounds\": %" PRIu64 ",\n  \"period\": %u,\n  \"cube_dimension\": %u,\n  \"cube\": \"N1=0; N0 repeats a period-bit word with bit zero fixed to zero\",\n",rounds,period,period-1);
        std::fprintf(out,"  \"begin\": %" PRIu64 ",\n  \"count\": %" PRIu64 ",\n  \"total_points\": %" PRIu64 ",\n  \"coverage_complete\": %s,\n  \"seconds\": %.9f,\n  \"keys_hex\": [\n",begin,count,total,begin==0&&count==total?"true":"false",seconds);
        for(size_t k=0;k<keys.size();++k)std::fprintf(out,"    [\"%016" PRIx64 "\",\"%016" PRIx64 "\"]%s\n",keys[k][0],keys[k][1],k+1==keys.size()?"":",");
        std::fprintf(out,"  ],\n  \"sums_hex\": [\n");
        for(size_t k=0;k<keys.size();++k)std::fprintf(out,"    [\"%016" PRIx64 "\",\"%016" PRIx64 "\"]%s\n",sums[k][0],sums[k][1],k+1==keys.size()?"":",");
        std::fprintf(out,"  ]\n}\n");std::fclose(out);
        std::printf("rounds=%" PRIu64 " points=%" PRIu64 " keys=%zu seconds=%.6f complete=%s\n",rounds,count,keys.size(),seconds,begin==0&&count==total?"yes":"no");
    }catch(const std::exception& e){std::fprintf(stderr,"Error: %s\n",e.what());return 1;}
}

// Enumerate E rotation orbits of weights 2,4,6. A matching in the forced
// input/output activity graph lower-bounds wc(X)+wc(Y).
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <string>
#include "activity_matching.hpp"

uint64_t rol(uint64_t v,int d) {return d?(v<<d)|(v>>(64-d)):v;}
int main(int argc,char**argv) {
    if(argc<2||argc>3){std::fprintf(stderr,"Usage: %s OUTPUT.tsv [PRUNED_MATCHINGS.bin]\n",argv[0]);return 1;}
    std::ofstream file(argv[1]);file<<"effect_weight\teffect_hex\tmatching_lower_bound\n";
    std::ofstream certificates;
    if(argc==3){certificates.open(argv[2],std::ios::binary);certificates.write("GSTNMAT1",8);}
    for(int k:{2,4,6}) {
        uint64_t total=0,retained=0;std::array<uint64_t,31> histogram{};
        Matcher matcher;
        auto visit=[&](uint64_t mask) {
            uint64_t bits=mask;
            while(bits) {
                int j=__builtin_ctzll(bits);bits&=bits-1;
                if(rol(mask,(64-j)%64)<mask)return;
            }
            ++total;
            unsigned bound=matcher.lower(mask);++histogram[bound];
            if(bound<=15) {
                ++retained;
                file<<k<<"\t"<<std::hex<<mask<<std::dec<<"\t"<<bound<<"\n";
            } else if(certificates.is_open()) {
                // Little-endian mask followed by 16 little-endian 12-bit edges:
                // low six bits are the input column, next six the output column.
                for(int byte=0;byte<8;++byte)certificates.put(char(mask>>(8*byte)));
                unsigned edges=0;
                for(unsigned b=0;b<64;++b)if(matcher.mate[b]>=0) {
                    unsigned pair=unsigned(matcher.mate[b])|(b<<6);
                    certificates.put(char(pair));certificates.put(char(pair>>8));++edges;
                }
                if(edges!=16)std::abort();
            }
        };
        auto enumerate=[&](auto&& self,int pos,int left,uint64_t mask)->void {
            if(!left){visit(mask);return;}
            for(int j=pos;j<=64-left;++j)self(self,j+1,left-1,mask|(1ULL<<j));
        };
        enumerate(enumerate,1,k-1,1);
        std::printf("weight=%d orbits=%llu retained=%llu\n",k,
                    (unsigned long long)total,(unsigned long long)retained);
        for(unsigned n=0;n<histogram.size();++n)if(histogram[n])
            std::printf("  matching_bound=%u count=%llu\n",n,(unsigned long long)histogram[n]);
        std::fflush(stdout);
    }
    return !file||(certificates.is_open()&&!certificates);
}

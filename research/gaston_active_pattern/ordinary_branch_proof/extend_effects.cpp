// Every >=8-bit E contains a 6-bit subset. If that subset was not already
// excluded by a size-16 matching, rotate it to one of the three retained
// weight-six representatives, then select any two further bits.
#include "activity_matching.hpp"
#include <fstream>
#include <sstream>
#include <cstdio>
#include <stdexcept>
#include <string>
int main() {
    std::ifstream input("effects_total15.tsv");
    std::ofstream output("extension_matchings.bin",std::ios::binary);
    output.write("GSTNEXT1",8);
    std::string line;std::getline(input,line);
    unsigned bases=0,cases=0;
    auto write64=[&](uint64_t value){for(int i=0;i<8;++i)output.put(char(value>>(8*i)));};
    while(std::getline(input,line)) {
        std::istringstream stream(line);
        unsigned weight,bound;std::string hex;
        stream>>weight>>hex>>bound;
        if(weight!=6)continue;
        ++bases;const uint64_t base=std::stoull(hex,nullptr,16);
        for(unsigned a=0;a<64;++a)if(!(base>>a&1))
        for(unsigned b=a+1;b<64;++b)if(!(base>>b&1)) {
            const uint64_t e=base|(1ULL<<a)|(1ULL<<b);
            Matcher matcher;
            if(matcher.lower(e)!=16)throw std::runtime_error("unexcluded extension");
            write64(base);write64(e);unsigned n=0;
            for(unsigned y=0;y<64;++y)if(matcher.mate[y]>=0) {
                unsigned edge=unsigned(matcher.mate[y])|(y<<6);
                output.put(char(edge));output.put(char(edge>>8));++n;
            }
            if(n!=16)throw std::runtime_error("wrong matching size");
            ++cases;
        }
    }
    if(bases!=3||cases!=4959||!output)throw std::runtime_error("incomplete extension enumeration");
    std::printf("All %u two-bit extensions of %u weight-six representatives have a size-16 matching.\n",cases,bases);
}

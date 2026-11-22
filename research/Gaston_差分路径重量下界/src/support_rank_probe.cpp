// Exact GF(2) necessary-condition probe for a common middle activity mask.
// No MILP model or solver parameters are modified.
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using Word=uint64_t;
using State=std::array<Word,5>;
using Columns=std::array<std::array<Word,5>,64>;
constexpr std::array<int,5> EAST{0,60,22,27,4},WEST{0,56,31,46,43},TWIN{25,32,52,60,63};
Word rol(Word x,int s) {s=(s%64+64)%64;return s?(x<<s)|(x>>(64-s)):x;}
Word effect(const State& a) {
    Word p=0,q=0;
    for(int i=0;i<5;++i){p^=a[i];q^=rol(a[i],TWIN[i]);}
    return p^rol(p,1)^q^rol(q,18);
}
struct Basis {
    std::array<Word,64> b{};
    int rank=0;
    int add(Word x) {
        while(x){int p=63-__builtin_clzll(x);if(b[p])x^=b[p];else{b[p]=x;++rank;return p;}}
        return -1;
    }
};
Columns left_columns{},right_columns{};
Word canonical(Word mask) {
    Word best=mask,bits=mask;
    while(bits){int z=__builtin_ctzll(bits);best=std::min(best,rol(mask,-z));bits&=bits-1;}
    return best;
}
int restricted_rank(const Columns& columns,Word mask) {
    Basis b;
    while(mask){int z=__builtin_ctzll(mask);for(Word v:columns[z])b.add(v);mask&=mask-1;}
    return b.rank;
}
void init() {
    for(int z=0;z<64;++z)for(int i=0;i<5;++i){
        State a{};a[i]=rol(Word(1)<<z,-WEST[i]);left_columns[z][i]=effect(a);
        a[i]=rol(Word(1)<<z,EAST[i]);right_columns[z][i]=effect(a);
    }
    // Every vortex provides valid left/right kernel examples in the correct coordinates.
    for(int i=0;i<5;++i)for(int j=i+1;j<5;++j)for(int k=j+1;k<5;++k){
        State a{};a[i]=(Word(1)<<TWIN[j])^(Word(1)<<TWIN[k]);
        a[j]=(Word(1)<<TWIN[i])^(Word(1)<<TWIN[k]);a[k]=(Word(1)<<TWIN[i])^(Word(1)<<TWIN[j]);
        if(effect(a))throw std::runtime_error("bad vortex");
        Word left=0,right=0;
        for(int r=0;r<5;++r){left|=rol(a[r],WEST[r]);right|=rol(a[r],-EAST[r]);}
        if(restricted_rank(left_columns,left)>=5*__builtin_popcountll(left) ||
           restricted_rank(right_columns,right)>=5*__builtin_popcountll(right))throw std::runtime_error("coordinate check failed");
    }
}
int main(int argc,char** argv) {
    try {
        unsigned long long samples=20000;
        bool full=false;int active=6;
        if(argc==2 && std::string(argv[1])=="--full-six")full=true;
        else if(argc==3 && std::string(argv[1])=="--full") {full=true;active=std::stoi(argv[2]);if(active<6||active>8)throw std::runtime_error("full active count must be 6..8");}
        else if(argc==3 && std::string(argv[1])=="--samples")samples=std::stoull(argv[2]);
        else if(argc!=1)throw std::runtime_error("Use --samples N or --full-six");
        init();auto started=std::chrono::steady_clock::now();
        uint64_t checked=0,left_dependent=0,both_dependent=0;
        std::set<Word> left_orbits,both_orbits;
        auto visit=[&](Word mask,int rank_left){
            ++checked;
            if(rank_left==5*active)return;
            ++left_dependent;left_orbits.insert(canonical(mask));
            if(restricted_rank(right_columns,mask)<5*active){++both_dependent;both_orbits.insert(canonical(mask));}
        };
        if(full){
            Basis b;for(auto v:left_columns[0])b.add(v);
            auto search=[&](auto&& self,int start,int remain,Word mask)->void {
                if(!remain){visit(mask,b.rank);return;}
                for(int z=start;z<=64-remain;++z){
                    std::array<int,5> pivots;
                    for(int i=0;i<5;++i)pivots[i]=b.add(left_columns[z][i]);
                    self(self,z+1,remain-1,mask|(Word(1)<<z));
                    for(int p:pivots)if(p>=0){b.b[p]=0;--b.rank;}
                }
            };
            search(search,1,active-1,1);
            uint64_t expected=1;for(int j=1;j<active;++j)expected=expected*(64-j)/j;
            if(checked!=expected)throw std::runtime_error("coverage count mismatch");
        } else {
            std::mt19937_64 rng(20261024);
            for(uint64_t n=0;n<samples;++n){Word mask=1;while(__builtin_popcountll(mask)<6)mask|=Word(1)<<(1+rng()%63);visit(mask,restricted_rank(left_columns,mask));}
        }
        double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::printf("{\n  \"scope\":\"necessary rank conditions only; rank deficiency does not imply full activity or low trail weight\",\n  \"complete_anchored_enumeration\":%s,\n  \"middle_active_columns\":%d,\n  \"checked\":%llu,\n  \"left_dependent\":%llu,\n  \"both_dependent\":%llu,\n  \"left_orbits\":%zu,\n  \"both_orbits\":%zu,\n  \"seconds\":%.6f,\n  \"both_orbit_masks\":[",full?"true":"false",active,(unsigned long long)checked,(unsigned long long)left_dependent,(unsigned long long)both_dependent,left_orbits.size(),both_orbits.size(),elapsed);
        bool first=true;for(Word mask:both_orbits){std::printf("%s\"%016llx\"",first?"":",",(unsigned long long)mask);first=false;}
        std::printf("],\n  \"left_orbit_masks\":[");
        first=true;for(Word mask:left_orbits){std::printf("%s\"%016llx\"",first?"":",",(unsigned long long)mask);first=false;}
        std::printf("]\n}\n");return 0;
    }catch(const std::exception& e){std::fprintf(stderr,"Error: %s\n",e.what());return 1;}
}

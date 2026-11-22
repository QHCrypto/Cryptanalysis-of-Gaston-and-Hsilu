// Exhaustive rank/consistency audit of q(b,w)=0 after the X2 conditions.
// Input is produced by verify_prefix.py; no permutation search is performed.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc,char**argv){
    if(argc!=2){std::cerr<<"Usage: enumerate_branches branch_system.txt\n";return 1;}
    std::ifstream input(argv[1]);
    unsigned guesses,free,eqs;input>>guesses>>free>>eqs;
    if(!input||guesses>30||free>32||eqs>32)throw std::runtime_error("Invalid dimensions");
    std::vector<uint64_t> row(eqs);
    for(auto &v:row)input>>v;
    std::vector<std::vector<uint64_t>> delta(eqs,std::vector<uint64_t>(guesses));
    auto neighbor=delta;
    for(unsigned r=0;r<eqs;++r){
        for(auto&v:delta[r])input>>v;
        for(auto&v:neighbor[r])input>>v;
    }
    if(!input)throw std::runtime_error("Truncated system");
    const uint64_t mask=(1ULL<<free)-1,total=1ULL<<guesses;
    std::array<uint64_t,33> all{},consistent{},inconsistent{};
    uint64_t branch=0,backgrounds=0,hash=14695981039346656037ULL;
    uint64_t spread_hash=14695981039346656037ULL;
    for(uint64_t seq=0;seq<total;++seq){
        if(seq){
            const unsigned bit=__builtin_ctzll(seq);
            for(unsigned r=0;r<eqs;++r)
                row[r]^=delta[r][bit]^(uint64_t(__builtin_parityll(branch&neighbor[r][bit]))<<free);
            branch^=1ULL<<bit;
        }
        if(seq<1024)for(auto v:row){hash^=v;hash*=1099511628211ULL;}
        if(seq%8191==0||seq+1==total)
            for(auto v:row){spread_hash^=v;spread_hash*=1099511628211ULL;}
        std::array<uint64_t,33> pivots{};
        unsigned rank=0;bool good=true;
        for(auto value:row){
            uint64_t coefficients;
            while((coefficients=value&mask)){
                const unsigned bit=63-__builtin_clzll(coefficients);
                if(pivots[bit])value^=pivots[bit];
                else{pivots[bit]=value;++rank;break;}
            }
            if(!(value&mask)&&(value>>free))good=false;
        }
        ++all[rank];
        if(good){++consistent[rank];backgrounds+=1ULL<<(free-rank);}
        else ++inconsistent[rank];
    }
    auto print=[&](const char*name,const auto&a){
        std::cout<<"  \""<<name<<"\": {";
        bool first=true;
        for(unsigned r=0;r<=free;++r)if(a[r]){
            if(!first) std::cout<<", ";
            first=false;
            std::cout<<'"'<<r<<"\": "<<a[r];
        }
        std::cout<<"},\n";
    };
    std::cout<<"{\n  \"independent_guesses\": "<<guesses<<",\n  \"free_variables\": "<<free
             <<",\n  \"branches_enumerated\": "<<total<<",\n  \"backgrounds\": "<<backgrounds<<",\n";
    print("coefficient_rank_histogram",all);print("consistent_rank_histogram",consistent);
    print("inconsistent_rank_histogram",inconsistent);
    std::cout<<"  \"first_1024_augmented_rows_fnv64\": \""<<hash<<"\",\n"
             <<"  \"spread_augmented_rows_fnv64\": \""<<spread_hash<<"\"\n}\n";
}

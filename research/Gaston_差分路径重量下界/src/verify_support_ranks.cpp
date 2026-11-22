// Independent verifier: matrices from official reference C, fresh elimination
// on every lexicographic support (no recursive incremental basis or rollback).
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <vector>
using Word=uint64_t;
using Table=std::array<std::array<Word,5>,64>;
int rank_support(const Table& table,const std::vector<int>& s) {
    std::array<Word,64> basis{};
    int rank=0;
    for(int z:s)for(Word x:table[z]){
        while(x){int p=__builtin_ctzll(x);if(basis[p])x^=basis[p];else{basis[p]=x;++rank;break;}}
    }
    return rank;
}
Word rol(Word x,int n){n&=63;return n?(x<<n)|(x>>(64-n)):x;}
Word canonical(const std::vector<int>& s){Word mask=0;for(int z:s)mask|=Word(1)<<z;Word best=mask;for(int z:s)best=std::min(best,rol(mask,64-z));return best;}
int main(int argc,char** argv){try{
    if(argc<2 || argc>3)throw std::runtime_error("matrix file and optional active count required");
    int active=argc==3?std::stoi(argv[2]):6;if(active<6||active>8)throw std::runtime_error("active count must be 6..8");
    std::ifstream input(argv[1]);Table left{},right{};
    for(auto* table:{&left,&right})for(auto& column:*table)for(auto& value:column){std::string word;if(!(input>>word))throw std::runtime_error("truncated matrix");value=std::stoull(word,nullptr,16);}
    uint64_t checked=0,left_count=0,both_count=0;std::set<Word> left_orbits,both_orbits;
    std::vector<int> s(active);for(int j=0;j<active;++j)s[j]=j;
    while(true){
        ++checked;
        if(rank_support(left,s)<5*active){++left_count;left_orbits.insert(canonical(s));if(rank_support(right,s)<5*active){++both_count;both_orbits.insert(canonical(s));}}
        int j=active-1;while(j>=1 && s[j]==64-active+j)--j;if(j<1)break;
        ++s[j];for(int k=j+1;k<active;++k)s[k]=s[k-1]+1;
    }
    uint64_t expected=1;for(int j=1;j<active;++j)expected=expected*(64-j)/j;
    if(checked!=expected)throw std::runtime_error("coverage mismatch");
    std::printf("{\"checked\":%llu,\"left_dependent\":%llu,\"both_dependent\":%llu,\"left_orbit_masks\":[",(unsigned long long)checked,(unsigned long long)left_count,(unsigned long long)both_count);
    bool first=true;for(Word mask:left_orbits){std::printf("%s\"%016llx\"",first?"":",",(unsigned long long)mask);first=false;}
    std::printf("],\"both_orbit_masks\":[");first=true;for(Word mask:both_orbits){std::printf("%s\"%016llx\"",first?"":",",(unsigned long long)mask);first=false;}
    std::printf("]}\n");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}

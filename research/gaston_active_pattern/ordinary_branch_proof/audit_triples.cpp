#include "../direct_enumeration/gaston_common.hpp"
#include <numeric>

unsigned column(const State& s, unsigned j) {
    unsigned v=0;
    for(unsigned r=0;r<5;++r)v|=((s[r]>>j)&1)<<r;
    return v;
}
bool inv15(const std::array<unsigned,15>& columns) {
    std::array<unsigned,15> basis{};
    for(auto col:columns) {
        while(col) {
            unsigned p=31-__builtin_clz(col);
            if(basis[p])col^=basis[p];
            else {basis[p]=col;break;}
        }
        if(!col)return false;
    }
    return true;
}
int main(int argc,char**argv) {
    Options o;o.input_fixed=argc>1 && std::string(argv[1])=="input";
    const auto table=make_table(o);
    unsigned histogram[22]={},total=0;
    const unsigned limit=argc>2?std::stoul(argv[2]):1953;
    for(unsigned b=1;b<63;++b)for(unsigned c=b+1;c<64;++c) {
        if(total>=limit)break;
        std::array<unsigned,15> cols{};
        unsigned projections[64][15]{};
        for(unsigned k=0;k<15;++k)for(unsigned j=0;j<64;++j)
            projections[j][k]=column(table[k<5?0:(k<10?b:c)][1u<<(k%5)],j);
        uint64_t used=0;unsigned found=0;
        for(unsigned i=0;i<64 && found<10;++i)if(!(used>>i&1))
        for(unsigned j=i+1;j<64 && !(used>>i&1);++j)if(!(used>>j&1))
        for(unsigned k=j+1;k<64 && !(used>>i&1);++k)if(!(used>>k&1)) {
            for(unsigned v=0;v<15;++v)cols[v]=projections[i][v]|(projections[j][v]<<5)|(projections[k][v]<<10);
            if(inv15(cols)) {used|=(1ULL<<i)|(1ULL<<j)|(1ULL<<k);++found;break;}
        }
        ++histogram[found];++total;
        if(total%100==0) {std::printf("configs=%u; last=%u,%u count=%u\n",total,b,c,found);std::fflush(stdout);}
    }
    for(unsigned i=0;i<22;++i)if(histogram[i])std::printf("disjoint triples %u: %u configs\n",i,histogram[i]);
}

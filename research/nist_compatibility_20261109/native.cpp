#include <cstdint>
#include "../hsilu_ref/hsilu.h"
#include "../gaston分析/three_round_algebra.hpp"
void hsilu_round(hsilu_state_t*, uint64_t);
static const uint64_t RC[12]={0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b};
static void hp(hsilu_state_t& a,int q){for(int i=12-q;i<12;++i)hsilu_round(&a,RC[i]);}
static uint64_t random_word(uint64_t& s){uint64_t z=(s+=0x9e3779b97f4a7c15ULL);z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
static gaston3::State g_y2(uint64_t m,gaston3::State s){
    s[0]^=m;s=gaston3::linear(s);s[0]^=0x5a;s=gaston3::chi(s);
    s=gaston3::linear(s);s[0]^=0x4b;return s;
}
static uint64_t target_bits(uint64_t w){
    const int j[17]={0,1,6,17,18,19,24,29,30,41,42,47,52,53,59,60,63};
    uint64_t out=0;for(int i=0;i<17;++i)out|=uint64_t(gaston3::bit(w,j[i]))<<i;return out;
}
extern "C" {
uint64_t prefix_scan(uint64_t iv,int q,uint64_t first,uint64_t start,uint64_t trials,
                     uint64_t m1,uint64_t m4,uint64_t* found){
    hsilu_state_t a={{iv,0,0,0,0}};hp(a,q);a.x[0]^=first;hp(a,q);
    uint64_t hits=0;
    for(uint64_t j=0;j<trials;++j){auto b=a;b.x[0]^=start+j;hp(b,q);
        if(!(b.x[1]&m1)&&(b.x[4]&m4)==m4){
            if(!hits){found[0]=first;found[1]=start+j;for(int i=0;i<5;++i)found[2+i]=b.x[i];}++hits;
        }
    }return hits;
}
void hperm(uint64_t* a,int q){hsilu_state_t s;for(int i=0;i<5;++i)s.x[i]=a[i];hp(s,q);for(int i=0;i<5;++i)a[i]=s.x[i];}
void gperm(uint64_t* a,int q){gaston3::State s;for(int i=0;i<5;++i)s[i]=a[i];for(int i=12-q;i<12;++i){s=gaston3::linear(s);s[0]^=RC[i];s=gaston3::chi(s);}for(int i=0;i<5;++i)a[i]=s[i];}
void gaston_conditioned_sample(const uint64_t* input,uint64_t target,uint64_t seed,
                              uint64_t trials,uint64_t* stats,uint64_t* hist){
    using namespace gaston3;
    const int freej[17]={0,4,5,10,11,22,23,27,33,34,40,41,45,51,56,57,63};
    bool free[64]={};for(int j:freej)free[j]=true;
    std::vector<Equation> fixed;for(int j=0;j<64;++j)if(!free[j])fixed.push_back({urow(j),false});
    auto space=eliminate(fixed);State s;for(int i=0;i<5;++i)s[i]=input[i];
    for(uint64_t t=0;t<trials;++t){
        uint64_t m=random_word(seed),constant=target_bits(g_y2(m,s)[0]);
        std::vector<Equation> rows(17);
        for(int i=0;i<17;++i)rows[i].rhs=(target^constant)>>i&1;
        for(int k=0;k<17;++k){uint64_t col=target_bits(g_y2(m^space.basis[k],s)[0])^constant;
            for(int i=0;i<17;++i)rows[i].row|=((col>>i)&1)<<k;
        }
        auto sol=eliminate(rows,17);++hist[sol.rank];++stats[0];if(!sol.consistent)continue;++stats[1];
        for(uint64_t z=0;z<(uint64_t(1)<<sol.basis.size());++z){
            uint64_t assignment=point(sol,z),v=m;
            for(int k=0;k<17;++k)if(assignment>>k&1)v^=space.basis[k];
            if(!v){++stats[4];continue;} // the sole invalid final padded word
            auto y=g_y2(v,s);if(target_bits(y[0])!=target){++stats[5];continue;}
            ++stats[2];if(target_bits(chi(y)[0])==target)++stats[3];
        }
    }
}
}

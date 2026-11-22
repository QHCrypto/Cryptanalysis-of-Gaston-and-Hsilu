#include "datf_middle.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace gdl;
int value(const Poly &p, const std::vector<int> &x) {
    int r=0;
    for(auto t:p) {
        if(!t)r^=1;
        else if(t>>32)r^=x.at(int(t>>32)-1)&x.at(int(uint32_t(t))-1);
        else r^=x.at(int(t)-1);
    }
    return r;
}
int main() {
    std::mt19937_64 rng(20261016);
    for(int trial=0;trial<240;++trial) {
        int n=1+trial%10;Poly p;
        if(rng()&1)p.insert(0);
        for(int i=0;i<n;++i) {
            if(rng()&1)p.insert(monomial(i));
            for(int j=i+1;j<n;++j)if((rng()%4)==0)p.insert(monomial(i,j));
        }
        int sum=0;
        for(int mask=0;mask<(1<<n);++mask) {
            std::vector<int>x(n);for(int i=0;i<n;++i)x[i]=(mask>>i)&1;
            sum+=value(p,x)?-1:1;
        }
        auto c=quadratic_sum(p);
        assert(c.sign==0?sum==0:sum==c.sign*(1<<(n-c.weight)));
    }
    for(int a=0;a<32;++a)for(int b=0;b<32;++b) {
        int sum=0;for(int x=0;x<32;++x)sum+=parity(b&(chi5(x)^chi5(x^a)))?-1:1;
        auto c=middle_one(embed(0,a),embed(0,b));assert(sum==32*c.sign);
    }
    MiddleTables tables;
    for(int k=0;k<64;++k) {
        State s{};for(auto &w:s)w=rng();auto expected=linear(chi(s));
        std::vector<int>x(320);for(int i=0;i<320;++i)x[i]=int((s[i/64]>>(i%64))&1);
        for(int i=0;i<320;++i)assert(value(tables.alpha_def[i],x)==int((expected[i/64]>>(i%64))&1));
    }
    std::cout<<"PASS: 240 exhaustive quadratic sums; 1024 exact one-chi correlations; 64 full linear-chi symbolic states\n";
}

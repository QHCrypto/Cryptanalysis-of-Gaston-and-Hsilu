#pragma once
#include "hsilu.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
namespace hdl {
inline std::string hex_state(const State& s){return words(s);}
struct Tables {
    int ddt[32][32]{},lat[32][32]{};
    Tables(){
        for(int a=0;a<32;++a)for(int x=0;x<32;++x)++ddt[a][chi5(x)^chi5(x^a)];
        for(int v=0;v<32;++v)for(int u=0;u<32;++u)for(int x=0;x<32;++x)
            lat[v][u]+=parity((v&chi5(x))^(u&x))?-1:1;
    }
};
struct LinearMasks {
    std::array<std::vector<std::pair<int,int>>,2048> terms;
    LinearMasks(){for(int j=0;j<64;++j)for(int v=1;v<32;++v){
        const auto a=transpose(embed(j,v));
        for(int k=0;k<64;++k)if(column(a,k))terms[j*32+v].push_back({k,column(a,k)});
    }}
};
}

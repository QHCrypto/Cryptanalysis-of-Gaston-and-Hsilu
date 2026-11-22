#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using U = uint64_t;
using State = std::array<U, 5>;
constexpr int RC[12] = {0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b};
constexpr int GE[5] = {0,60,22,27,4}, GW[5] = {0,56,31,46,43}, GT[5] = {25,32,52,60,63};
constexpr int A[13] = {0,3,1,0,5,5,7,11,12,12,13,14,15};
constexpr int B[13] = {4,4,2,1,2,3,6,8,9,6,10,11,9};
constexpr int M[13] = {47,51,54,30,4,55,0,14,59,4,10,45,53};
constexpr int HR[5] = {0,57,13,31,20};
static inline U rol(U x, unsigned n) { return (x << n) | (x >> ((64-n)&63)); }

State permute_window(State s, bool gaston, int start, int rounds) {
    if (start < 0 || rounds < 0 || start + rounds > 12) return s;
    for (int q=start; q<start+rounds; ++q) {
        if (gaston) {
            U p=0, v=0;
            for (int j=0;j<5;++j) { s[j]=rol(s[j],GE[j]); p^=s[j]; v^=rol(s[j],GT[j]); }
            U e=rol(p^rol(p,1)^v^rol(v,18),23);
            for (int j=0;j<5;++j) s[j]=rol(s[j]^e,GW[j]);
        }
        s[0]^=RC[q];
        State t=s;
        for (int j=0;j<5;++j) s[j]=t[j]^((~t[(j+1)%5])&t[(j+2)%5]);
        if (!gaston) {
            U x[18]={s[0],s[1],s[2],s[3],s[4]};
            for(int j=0;j<13;++j) x[j+5]=x[A[j]]^rol(x[B[j]],M[j]);
            for(int j=0;j<5;++j) s[j]=rol(x[j+13],HR[j]);
        }
    }
    return s;
}

State permute(State s, bool gaston, int rounds) {
    return permute_window(s, gaston, 12-rounds, rounds);
}

State linear_cube_window(bool gaston, int start, int rounds, U iv, U k0, U k1,
                         U base0, U base1,
                         const std::vector<std::array<U,2>>& basis);

State cube(bool gaston, int rounds, U iv, U k0, U k1,
           U base0, U base1, const std::vector<int>& vars) {
    State sum{};
    U n0=base0, n1=base1;
    for (U a=0; a<(U(1)<<vars.size()); ++a) {
        if (a) {
            int p=vars[__builtin_ctzll(a)];
            if(p<64) n0^=U(1)<<p; else n1^=U(1)<<(p-64);
        }
        State s=permute({iv,k0,k1,n0,n1},gaston,rounds);
        for(int j=0;j<5;++j) sum[j]^=s[j];
    }
    return sum;
}

State linear_cube(bool gaston, int rounds, U iv, U k0, U k1,
                  U base0, U base1,
                  const std::vector<std::array<U,2>>& basis) {
    return linear_cube_window(gaston, 12-rounds, rounds, iv, k0, k1,
                              base0, base1, basis);
}

State linear_cube_window(bool gaston, int start, int rounds, U iv, U k0, U k1,
                         U base0, U base1,
                         const std::vector<std::array<U,2>>& basis) {
    State sum{};
    U n0=base0, n1=base1;
    for (U a=0; a<(U(1)<<basis.size()); ++a) {
        if (a) {
            const auto &v=basis[__builtin_ctzll(a)];
            n0^=v[0];
            n1^=v[1];
        }
        State s=permute_window({iv,k0,k1,n0,n1},gaston,start,rounds);
        for(int j=0;j<5;++j) sum[j]^=s[j];
    }
    return sum;
}

void print(const State &s) {
    for(U x:s) std::cout<<std::hex<<std::setw(16)<<std::setfill('0')<<x<<' ';
    std::cout<<'\n';
}

// Persistent local protocol, hex state words, decimal round/dimension/indices.
// p name rounds x0 x1 x2 x3 x4
// c name rounds iv k0 k1 n0 n1 dimension indices...
// t name rounds iv k0 k1 n0 n1 dimension indices... start_key_bit count
// t returns the base sum followed by one sum for each toggled key bit.
// l name rounds iv k0 k1 n0 n1 dimension basis0_lo basis0_hi ...
// q name rounds start iv k0 k1 n0 n1 dimension basis0_lo basis0_hi ...
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string command, name;
    while(std::cin>>command>>name) {
        int rounds;
        std::cin>>std::dec>>rounds;
        if((name!="gaston" && name!="hsilu") || rounds<0 || rounds>12) return 2;
        if(command=="p") {
            State s;
            for(U &x:s) std::cin>>std::hex>>x;
            print(permute(s,name=="gaston",rounds));
        } else if(command=="w") {
            int start;
            State s;
            std::cin>>std::dec>>start;
            for(U &x:s) std::cin>>std::hex>>x;
            if(start<0 || start+rounds>12) return 3;
            print(permute_window(s,name=="gaston",start,rounds));
        } else if(command=="c" || command=="t") {
            U iv,k0,k1,n0,n1;
            int d;
            std::cin>>std::hex>>iv>>k0>>k1>>n0>>n1>>std::dec>>d;
            if(d<0 || d>24) return 3;
            std::vector<int> vars(d);
            std::array<bool,128> seen{};
            for(int &p:vars) {
                std::cin>>p;
                if(p<0 || p>=128 || seen[p]) return 4;
                seen[p]=true;
            }
            print(cube(name=="gaston",rounds,iv,k0,k1,n0,n1,vars));
            if(command=="t") {
                int start,count;
                std::cin>>start>>count;
                if(start<0 || count<0 || start+count>128) return 5;
                for(int j=start;j<start+count;++j)
                    print(cube(name=="gaston",rounds,iv,
                               k0^(j<64 ? U(1)<<j : 0),
                               k1^(j>=64 ? U(1)<<(j-64) : 0),n0,n1,vars));
            }
        } else if(command=="l") {
            U iv2,k02,k12,n02,n12;
            int d;
            std::cin>>std::hex>>iv2>>k02>>k12>>n02>>n12>>std::dec>>d;
            if(d<0 || d>24) return 3;
            std::vector<std::array<U,2>> basis(d);
            for(auto &v:basis) std::cin>>std::hex>>v[0]>>v[1];
            print(linear_cube(name=="gaston",rounds,iv2,k02,k12,n02,n12,basis));
        } else if(command=="q") {
            int start;
            std::cin>>std::dec>>start;
            U iv2,k02,k12,n02,n12;
            int d;
            std::cin>>std::hex>>iv2>>k02>>k12>>n02>>n12>>std::dec>>d;
            if(start<0 || start+rounds>12 || d<0 || d>24) return 3;
            std::vector<std::array<U,2>> basis(d);
            for(auto &v:basis) std::cin>>std::hex>>v[0]>>v[1];
            print(linear_cube_window(name=="gaston",start,rounds,iv2,k02,k12,n02,n12,basis));
        } else return 6;
        std::cout.flush();
    }
    return std::cin.eof() ? 0 : 7;
}

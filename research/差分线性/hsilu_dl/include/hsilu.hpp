#pragma once
#include <array>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
namespace hdl {
using State=std::array<uint64_t,5>;
inline constexpr int N1[13]={0,3,1,0,5,5,7,11,12,12,13,14,15};
inline constexpr int N2[13]={4,4,2,1,2,3,6,8,9,6,10,11,9};
inline constexpr int SH[13]={47,51,54,30,4,55,0,14,59,4,10,45,53};
inline constexpr int ROT[5]={0,57,13,31,20};
inline constexpr uint64_t RC[12]={0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b};
inline uint64_t rotl(uint64_t x,unsigned n){n&=63;return (x<<n)|(x>>((-n)&63));}
inline State operator^(State a,const State& b){for(int i=0;i<5;++i)a[i]^=b[i];return a;}
inline State embed(int col,int pattern){State s{};for(int i=0;i<5;++i)if(pattern&(1<<i))s[i]=UINT64_C(1)<<col;return s;}
inline int column(const State& s,int j){int v=0;for(int i=0;i<5;++i)v|=((s[i]>>j)&1)<<i;return v;}
inline int parity(unsigned x){return __builtin_parity(x);}
inline int dot(const State& a,const State& b){uint64_t t=0;for(int i=0;i<5;++i)t^=a[i]&b[i];return __builtin_parityll(t);}
inline int active(const State& s){uint64_t v=0;for(auto w:s)v|=w;return __builtin_popcountll(v);}
inline State linear(const State& a){
    uint64_t x[18]{};for(int i=0;i<5;++i)x[i]=a[i];
    for(int j=0;j<13;++j)x[5+j]=x[N1[j]]^rotl(x[N2[j]],SH[j]);
    State out{};for(int i=0;i<5;++i)out[i]=rotl(x[13+i],ROT[i]);return out;
}
inline State transpose(const State& a){
    uint64_t x[18]{};for(int i=0;i<5;++i)x[13+i]=rotl(a[i],-ROT[i]);
    for(int j=12;j>=0;--j){x[N1[j]]^=x[5+j];x[N2[j]]^=rotl(x[5+j],-SH[j]);}
    return State{x[0],x[1],x[2],x[3],x[4]};
}
inline State chi(const State& s){State t{};for(int i=0;i<5;++i)t[i]=s[i]^(~s[(i+1)%5]&s[(i+2)%5]);return t;}
inline int chi5(int a){int b=0;for(int i=0;i<5;++i)b|=(((a>>i)&1)^((1^((a>>((i+1)%5))&1))&((a>>((i+2)%5))&1)))<<i;return b;}
inline State round(State a,uint64_t rc){a[0]^=rc;return linear(chi(a));}
inline State full(State a,int R){for(int r=12-R;r<12;++r)a=round(a,RC[r]);return a;}
// Core begins after the first iota and ends before the last linear layer.
inline State core(State a,int R){a=chi(a);for(int r=1;r<R;++r){a=linear(a);a[0]^=RC[12-R+r];a=chi(a);}return a;}
inline uint64_t mix(uint64_t z){z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);return z^(z>>31);}
inline State sample(uint64_t seed,uint64_t k){State a{};for(int i=0;i<5;++i)a[i]=mix(seed+(5*k+i+1)*UINT64_C(0x9e3779b97f4a7c15));return a;}
inline std::string words(const State& a){std::ostringstream s;s<<std::hex<<std::setfill('0');for(int i=0;i<5;++i){if(i)s<<';';s<<std::setw(16)<<a[i];}return s.str();}
}

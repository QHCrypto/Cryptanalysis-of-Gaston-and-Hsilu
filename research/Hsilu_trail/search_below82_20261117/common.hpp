#pragma once
#include <array>
#include <vector>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
using U=uint64_t;
using S=std::array<U,5>;

inline U rol(U x,int n){n&=63;return (x<<n)|(x>>((-n)&63));}
inline S operator^(const S&a,const S&b){S c;for(int i=0;i<5;i++)c[i]=a[i]^b[i];return c;}
inline S rotate(S a,int n){for(auto &x:a)x=rol(x,n);return a;}
inline int hw(const S&a){int w=0;for(U x:a)w+=__builtin_popcountll(x);return w;}
inline int active(const S&a){U v=0;for(U x:a)v|=x;return __builtin_popcountll(v);}
inline int col(const S&a,int z){int c=0;for(int i=0;i<5;i++)c|=((a[i]>>z)&1)<<i;return c;}
inline S linear(S a){
 constexpr int n1[13]={0,3,1,0,5,5,7,11,12,12,13,14,15};
 constexpr int n2[13]={4,4,2,1,2,3,6,8,9,6,10,11,9};
 constexpr int m[13]={47,51,54,30,4,55,0,14,59,4,10,45,53};
 constexpr int r[5]={0,57,13,31,20};
 U x[18]{};for(int i=0;i<5;i++)x[i]=a[i];
 for(int i=0;i<13;i++)x[5+i]=x[n1[i]]^rol(x[n2[i]],m[i]);
 for(int i=0;i<5;i++)a[i]=rol(x[13+i],r[i]);
 return a;
}
inline S canonical(S a){S best=a;for(int k=1;k<64;k++)best=std::min(best,rotate(a,k));return best;}
inline int chi(int a){int b=0;for(int i=0;i<5;i++)b|=(((a>>i)&1)^((1^((a>>((i+1)%5))&1))&((a>>((i+2)%5))&1)))<<i;return b;}
inline int weight(const S&a){
    U any=0,two=0,three=0;
    for(int i=0;i<5;i++){two|=any&a[i];any|=a[i];three|=a[i]&a[(i+1)%5]&a[(i+2)%5];}
    return 2*__builtin_popcountll(any)+__builtin_popcountll(two)+__builtin_popcountll(three);
}
struct Tables{
    int ddt[32][32]{},w[32]{},rev[32]{};
    std::vector<int> basis[32],outputs[32];int offset[32]{};
    Tables(){
        for(int a=0;a<32;a++)for(int x=0;x<32;x++)ddt[a][chi(x)^chi(x^a)]++;
        for(int a=0;a<32;a++){
            for(int b=0;b<32;b++)if(ddt[a][b])outputs[a].push_back(b);
            offset[a]=outputs[a][0];w[a]=__builtin_ctz((unsigned)outputs[a].size());
            int piv[5]{};
            for(int b:outputs[a]){int v=b^offset[a];for(int i=4;i>=0;i--)if((v>>i)&1){if(piv[i])v^=piv[i];else{piv[i]=v;basis[a].push_back(v);break;}}}
            assert((int)basis[a].size()==w[a]);
            S s{};for(int i=0;i<5;i++)s[i]=(a>>i)&1;assert(weight(s)==w[a]);
        }
        for(int b=0;b<32;b++){rev[b]=99;for(int a=0;a<32;a++)if(ddt[a][b])rev[b]=std::min(rev[b],w[a]);}
    }
    int reverse_weight(const S&s)const{int w=0;for(int z=0;z<64;z++)w+=rev[col(s,z)];return w;}
    std::pair<S,std::vector<S>> affine(const S&a)const{
        S off{};std::vector<S> bs;
        for(int z=0;z<64;z++){
            int v=col(a,z);for(int i=0;i<5;i++)off[i]|=U((offset[v]>>i)&1)<<z;
            for(int b:basis[v]){S s{};for(int i=0;i<5;i++)s[i]=U((b>>i)&1)<<z;bs.push_back(s);}
        }return {off,bs};
    }
};
inline std::string json_state(const S&s){std::ostringstream o;o<<'[';for(int i=0;i<5;i++){if(i)o<<',';o<<'"'<<std::hex<<std::setw(16)<<std::setfill('0')<<s[i]<<'"';}return o.str()+']';}

inline int reverse_weight_fast(const S&a){U any=0,extra=0;for(int i=0;i<5;i++){any|=a[i];extra|=a[i]&a[(i+1)%5]&a[(i+3)%5];}return 2*__builtin_popcountll(any)+__builtin_popcountll(extra);}
inline const S INVERSE[5]={{0x3bdd40f15eda019aULL,0x7b680668ef6503c5ULL,0xa019a3bd940f15adULL,0xb477b281e2b5b403ULL,0xf6503c56b680748eULL},{0x6a9019932bb45d3dULL,0xaed134f5aa40664cULL,0x44d3d6b9019912bbULL,0x7ad720332257689aULL,0xe406644aed124f5aULL},{0xb246ba90e08cf012ULL,0x8233c04ac91aca43ULL,0x4f012b246b290e08ULL,0x25608d6521c10960ULL,0x11aca438212c04a8ULL},{0x500e51c461084084ULL,0x8421021140394611ULL,0xc4084500e5184610ULL,0x8a01ca308c21881ULL,0x394611843102114ULL},{0x568ac7d21c918e8ULL,0x872463a015a2b1b4ULL,0x908e80568ac6d21cULL,0xd00ad158da439210ULL,0x522b1b4872421a01ULL}};
inline S inverse(const S&a){S b{};for(int i=0;i<5;i++){U w=a[i];while(w){int z=__builtin_ctzll(w);w&=w-1;b=b^rotate(INVERSE[i],z);}}return b;}

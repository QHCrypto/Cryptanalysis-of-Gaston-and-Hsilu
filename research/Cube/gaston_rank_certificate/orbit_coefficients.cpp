// Exact rotation-orbit summation for low coefficients of the parity-cube
// circulant matrix. Incomplete coverage never certifies a matrix coefficient.
#include <array>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
extern "C" {
#include "gaston.h"
}

using Word = uint64_t;
using State = std::array<Word, 5>;
constexpr unsigned E[5] = {0,60,22,27,4};
constexpr unsigned W[5] = {0,56,31,46,43};
constexpr unsigned T[5] = {25,32,52,60,63};

State skeleton(State a, unsigned rounds) {
    for (unsigned r=0;r<rounds;++r) {
        Word p=0,q=0;
        for(unsigned i=0;i<5;++i) {
            a[i]=rotate_left64(a[i],E[i]);p^=a[i];
            q^=rotate_left64(a[i],T[i]);
        }
        Word effect=rotate_left64(p^rotate_left64(p,1)^q^rotate_left64(q,18),23);
        for(unsigned i=0;i<5;++i)a[i]=rotate_left64(a[i]^effect,W[i]);
        State b;
        for(unsigned i=0;i<5;++i)b[i]=a[(i+1)%5]&a[(i+2)%5];
        a=b;
    }
    return a;
}

std::array<State,2> skeleton_derivatives(Word nonce,Word key,unsigned rounds) {
    // Ordinary formal derivative before Boolean reduction of the key.
    // Repeated key occurrences cannot contribute to degree 64 with 63
    // distinct nonce variables. Thus this derivative gives the exact
    // full-cube matrix while making every partial contribution key-linear.
    std::array<State,3> a{{State{0,0,0,nonce,0},State{0,key,0,0,0},State{0,0,key,0,0}}};
    for(unsigned r=0;r<rounds;++r) {
        for(auto& s:a) {
            Word p=0,q=0;
            for(unsigned i=0;i<5;++i) {
                s[i]=rotate_left64(s[i],E[i]);p^=s[i];q^=rotate_left64(s[i],T[i]);
            }
            Word effect=rotate_left64(p^rotate_left64(p,1)^q^rotate_left64(q,18),23);
            for(unsigned i=0;i<5;++i)s[i]=rotate_left64(s[i]^effect,W[i]);
        }
        std::array<State,3> b;
        for(unsigned i=0;i<5;++i) {
            unsigned j=(i+1)%5,k=(i+2)%5;
            b[0][i]=a[0][j]&a[0][k];
            for(unsigned d=1;d<3;++d)b[d][i]=(a[d][j]&a[0][k])^(a[0][j]&a[d][k]);
        }
        a=b;
    }
    return {a[1],a[2]};
}

Word number(const char* s) {
    if(!s[0]||s[0]=='-')throw std::runtime_error("unsigned integer required");
    size_t n=0;Word v=std::stoull(s,&n,0);
    if(n!=std::string(s).size())throw std::runtime_error("invalid integer");
    return v;
}

int main(int argc,char** argv) {
    if(argc!=8) {
        std::fprintf(stderr,"usage: %s rounds nonce_period m even_filter max_orbits output.json mode\n"
                     "mode: sum, dual, or census; max_orbits=0 means complete traversal\n",argv[0]);
        return 2;
    }
    try {
        unsigned rounds=number(argv[1]),period=number(argv[2]),m=number(argv[3]);
        bool even=number(argv[4])!=0;
        Word limit=number(argv[5]);
        std::string mode=argv[7];
        if(rounds<1||rounds>6||period<2||period>64||(period&(period-1))||
           m==0||m>8||(m&(m-1))||m>period||(mode!="sum"&&mode!="dual"&&mode!="census"))
            throw std::runtime_error("invalid parameters");
        unsigned blocks=period/m,alphabet=1u<<m;
        Word lowmask=(Word(1)<<m)-1,key=0;
        for(unsigned j=0;j<64;j+=m)key|=Word(1)<<j;
        std::vector<unsigned> a(blocks+1);
        std::array<std::array<Word,2>,3> sums{};
        std::array<Word,65> histogram{};
        Word orbits=0,points=0;
        bool stopped=false;
        auto start=std::chrono::steady_clock::now();
        std::function<void(unsigned,unsigned)> visit=[&](unsigned pos,unsigned p) {
            if(stopped)return;
            if(pos>blocks) {
                if(blocks%p)return;
                Word base=0;
                for(unsigned j=0;j<blocks;++j)base|=Word(a[j+1])<<(m*j);
                if(even&&__builtin_parityll(base))return;
                // The limit is checked before the next selected orbit, so a
                // run ending naturally at exactly the limit is still complete.
                if(limit&&orbits==limit){stopped=true;return;}
                ++orbits;points+=p;++histogram[p];
                if(mode=="census")return;
                Word nonce=0;
                for(unsigned j=0;j<64;j+=period)nonce|=base<<j;
                std::array<State,2> derivatives{};
                if(mode=="dual")derivatives=skeleton_derivatives(nonce,key,rounds);
                for(unsigned k=(mode=="dual"?1:0);k<3;++k) {
                    State s=mode=="dual"?derivatives[k-1]:
                        skeleton({0,k==1?key:0,k==2?key:0,nonce,0},rounds);
                    // Sum over R_m^j, 0<=j<p. Inputs and output have period mp.
                    // One m-bit block of the orbit norm is enough.
                    for(unsigned w=0;w<2;++w) {
                        Word folded=0;
                        for(unsigned j=0;j<p;++j)folded^=(s[w]>>(m*j))&lowmask;
                        sums[k][w]^=folded;
                    }
                }
                return;
            }
            a[pos]=a[pos-p];visit(pos+1,p);
            for(unsigned j=a[pos-p]+1;j<alphabet&&!stopped;++j) {
                a[pos]=j;visit(pos+1,pos);
            }
        };
        visit(1,1);
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        bool target_cube=period==(1u<<rounds)&&even;
        bool certificate=target_cube&&!stopped&&mode!="census";
        FILE* out=std::fopen(argv[6],"w");
        if(!out)throw std::runtime_error("cannot open output");
        std::fprintf(out,"{\n  \"rounds\": %u,\n  \"word_bits\": 64,\n  \"nonce_period\": %u,\n  \"m\": %u,\n"
                     "  \"even_filter\": %s,\n  \"mode\": \"%s\",\n  \"orbits_processed\": %" PRIu64 ",\n"
                     "  \"points_covered\": %" PRIu64 ",\n  \"enumeration_complete\": %s,\n"
                     "  \"target_cube\": %s,\n  \"coefficient_certificate\": %s,\n  \"seconds\": %.9f,\n",
                     rounds,period,m,even?"true":"false",mode.c_str(),orbits,points,
                     stopped?"false":"true",target_cube?"true":"false",certificate?"true":"false",seconds);
        std::fprintf(out,"  \"orbit_length_histogram\": {");
        bool first=true;
        for(unsigned j=1;j<=blocks;++j)if(histogram[j]) {
            std::fprintf(out,"%s\"%u\": %" PRIu64,first?"":", ",j,histogram[j]);first=false;
        }
        std::fprintf(out,"},\n  \"folded_sums_hex\": [");
        for(unsigned k=0;k<3;++k)std::fprintf(out,"%s[\"%" PRIx64 "\",\"%" PRIx64 "\"]",
                                            k?",":"",sums[k][0],sums[k][1]);
        std::fprintf(out,"],\n  \"low_coefficient_matrices\": ");
        if(!certificate)std::fprintf(out,"null");
        else {
            std::fprintf(out,"[");
            for(unsigned i=0;i<m;++i) {
                if(i)std::fprintf(out,",");
                std::fprintf(out,"[");
                for(unsigned w=0;w<2;++w) {
                    if(w)std::fprintf(out,",");
                    std::fprintf(out,"[");
                    for(unsigned k=1;k<=2;++k) {
                        Word v=sums[k][w]^sums[0][w];unsigned bit=0;
                        for(unsigned j=0;j<m;++j)if((j&i)==i)bit^=(v>>j)&1;
                        std::fprintf(out,"%s%u",k==2?",":"",bit);
                    }
                    std::fprintf(out,"]");
                }
                std::fprintf(out,"]");
            }
            std::fprintf(out,"]");
        }
        std::fprintf(out,"\n}\n");std::fclose(out);
        std::printf("r=%u period=%u m=%u orbits=%" PRIu64 " points=%" PRIu64
                    " seconds=%.6f complete=%s coefficient_certificate=%s\n",
                    rounds,period,m,orbits,points,seconds,stopped?"no":"yes",certificate?"yes":"no");
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Error: %s\n",e.what());return 1;
    }
}

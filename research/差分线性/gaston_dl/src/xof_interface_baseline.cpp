// Research profile: Ascon-XOF byte interface with Gaston permutations.
// Only the permutation after the final padded block is reduced.
#include "gaston.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>
extern "C" {
#include "gaston.h"
void gaston_round(gaston_state_t *, uint64_t);
}
using namespace gdl;
constexpr uint64_t IV = UINT64_C(0x0000080000cc0003);
constexpr uint64_t TAIL_MASK = (UINT64_C(1) << 56) - 1;

State permute(State s, int n) {
    for (int r=12-n; r<12; ++r) s=round(s,RC[r]);
    return s;
}
State reference_permute(State s, int n) {
    gaston_state_t t{};
    for(int i=0;i<5;++i)t.x[i]=s[i];
    for(int r=12-n;r<12;++r)gaston_round(&t,RC[r]);
    for(int i=0;i<5;++i)s[i]=t.x[i];
    return s;
}
using Permutation = State (*)(State,int);
std::vector<uint8_t> xof(const std::vector<uint8_t>& message, size_t outlen,
                         int final_rounds, Permutation f=permute) {
    State s=f(State{IV,0,0,0,0},12);
    size_t offset=0;
    while(message.size()-offset>=8){
        for(int j=0;j<8;++j)s[0]^=uint64_t(message[offset+j])<<(8*j);
        s=f(s,12);offset+=8;
    }
    const size_t rest=message.size()-offset;
    for(size_t j=0;j<rest;++j)s[0]^=uint64_t(message[offset+j])<<(8*j);
    s[0]^=UINT64_C(1)<<(8*rest);
    s=f(s,final_rounds);
    std::vector<uint8_t> out;
    while(out.size()<outlen){
        for(int j=0;j<8 && out.size()<outlen;++j)out.push_back(uint8_t(s[0]>>(8*j)));
        if(out.size()<outlen)s=f(s,12);
    }
    return out;
}
State before_last(const State& initialized,uint64_t seed,uint64_t k){
    auto random=sample_state(seed,k);
    State s=initialized;s[0]^=random[0];s=permute(s,12);
    s[0]^=(random[1]&TAIL_MASK)^(UINT64_C(1)<<56);
    return s;
}
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
void selftest(){
    for(size_t length: {size_t(0),size_t(1),size_t(7),size_t(8),size_t(15),size_t(16),size_t(23)}){
        std::vector<uint8_t> m(length);
        for(size_t j=0;j<length;++j)m[j]=uint8_t(37*j+11);
        for(int R:{3,4,12})for(size_t n:{size_t(0),size_t(8),size_t(17)})
            require(xof(m,n,R)==xof(m,n,R,reference_permute),"reference XOF mismatch");
    }
    const State init=permute(State{IV,0,0,0,0},12);
    for(uint64_t k=0;k<32;++k){
        auto random=sample_state(71,k);std::vector<uint8_t> m(15);
        for(int j=0;j<8;++j)m[j]=uint8_t(random[0]>>(8*j));
        for(int j=0;j<7;++j)m[8+j]=uint8_t(random[1]>>(8*j));
        const auto s=before_last(init,71,k);
        for(int R:{3,4}){
            auto encoded=xof(m,8,R,reference_permute);uint64_t expected=0;
            for(int j=0;j<8;++j)expected|=uint64_t(encoded[j])<<(8*j);
            require(permute(s,R)[0]==expected,"optimized 15-byte message mismatch");
            auto first=linear(s);first[0]^=RC[12-R];
            require(core(first,R,12-R)==permute(s,R),"first linear/constant boundary mismatch");
            for(int bit:{0,7,31,55}){
                auto changed=m;changed[8+bit/8]^=uint8_t(1<<(bit%8));
                auto pair=s;pair[0]^=UINT64_C(1)<<bit;
                auto y=xof(changed,8,R,reference_permute);uint64_t v=0;
                for(int j=0;j<8;++j)v|=uint64_t(y[j])<<(8*j);
                require(permute(pair,R)[0]==v,"legal message difference mismatch");
            }
        }
    }
    std::cout<<"PASS reference permutation, padding boundaries, squeeze continuation, optimized message pairs, first L and RC\n";
}
struct Counts { int rounds, bit; std::array<uint64_t,64> ones{}; };
int main(int argc,char**argv)try{
    if(argc==2 && std::string(argv[1])=="--selftest"){selftest();return 0;}
    if(argc!=6)throw std::runtime_error("usage: xof_interface_baseline OUTPUT.json SAMPLES SEED THREADS INPUT_BITS_CSV|all");
    const std::string output=argv[1];
    require(!std::ifstream(output).good(),"output already exists");
    const uint64_t n=std::stoull(argv[2]),seed=std::stoull(argv[3]);
    const int threads=std::stoi(argv[4]);
    require(n>0 && n<=(UINT64_C(1)<<26) && threads>0 && threads<=32,"invalid resources");
    std::vector<int> bits;
    if(std::string(argv[5])=="all")for(int b=0;b<56;++b)bits.push_back(b);
    else{std::stringstream ss(argv[5]);std::string token;while(std::getline(ss,token,','))bits.push_back(std::stoi(token));}
    std::sort(bits.begin(),bits.end());bits.erase(std::unique(bits.begin(),bits.end()),bits.end());
    require(!bits.empty() && bits.front()>=0 && bits.back()<56,"difference touches padding or is empty");
    const auto began=std::chrono::steady_clock::now();
    const auto initialized=permute(State{IV,0,0,0,0},12);
    std::vector<Counts> results;
    for(int R:{3,4})for(int bit:bits)results.push_back(Counts{R,bit,{}});
    // Fixed-size chunks bound memory independently of the requested sample count.
    constexpr uint64_t chunk=32768;
    for(uint64_t offset=0;offset<n;offset+=chunk){
        const uint64_t size=std::min(chunk,n-offset);
        std::vector<State> states(size);
        std::array<std::vector<uint64_t>,2> base;
        for(auto& v:base)v.resize(size);
        for(uint64_t k=0;k<size;++k){
            states[k]=before_last(initialized,seed,offset+k);
            base[0][k]=permute(states[k],3)[0];base[1][k]=permute(states[k],4)[0];
        }
        std::atomic<size_t> next{0};std::vector<std::thread> workers;
        for(int t=0;t<threads;++t)workers.emplace_back([&]{
            while(true){
                const size_t index=next.fetch_add(1);if(index>=results.size())break;
                auto& result=results[index];
                std::array<std::array<uint64_t,256>,8> histogram{};
                for(uint64_t k=0;k<size;++k){
                    auto s=states[k];s[0]^=UINT64_C(1)<<result.bit;
                    uint64_t d=permute(s,result.rounds)[0]^base[result.rounds-3][k];
                    for(int byte=0;byte<8;++byte){++histogram[byte][d&255];d>>=8;}
                }
                for(int byte=0;byte<8;++byte)for(int value=0;value<256;++value)
                    for(int bit=0;bit<8;++bit)if(value&(1<<bit))result.ones[8*byte+bit]+=histogram[byte][value];
            }
        });
        for(auto& thread:workers)thread.join();
        if((offset+size)%(UINT64_C(1)<<20)==0 || offset+size==n)
            std::cerr<<"sampled "<<offset+size<<" / "<<n<<" base messages\n";
    }
    std::ofstream out(output);require(out.good(),"cannot open output");
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    out<<"{\"profile\":\"Gaston-XOF final-absorption-only reduction\",\"iv\":\"0000080000cc0003\","
       <<"\"message_bytes\":15,\"random_prefix_bytes\":8,\"random_suffix_bytes\":7,"
       <<"\"padding_word\":\"0100000000000000\",\"initialization_rounds\":12,\"prefix_rounds\":12,"
       <<"\"squeeze_continuation_rounds\":12,\"first_linear_included\":true,"
       <<"\"samples_per_endpoint\":"<<n<<",\"seed\":"<<seed<<",\"threads\":"<<threads
       <<",\"elapsed_seconds\":"<<seconds<<",\"rows\":[";
    bool comma=false;
    for(const auto& r:results)for(int b=0;b<64;++b){
        if(comma)out<<',';
        comma=true;
        out<<"{\"rounds\":"<<r.rounds<<",\"input_bit\":"<<r.bit<<",\"output_bit\":"<<b
           <<",\"n0\":"<<n-r.ones[b]<<",\"n1\":"<<r.ones[b]<<'}';
    }
    out<<"]}\n";require(out.good(),"output write failed");
    std::cout<<"DONE "<<results.size()*64<<" endpoints, "<<seconds<<" seconds\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

#include "hsilu.hpp"
#include "hsilu.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>
void hsilu_round(hsilu_state_t*,uint64_t);
using namespace hdl;
using Real=long double;
using Local=std::array<Real,32>;
using Corr=std::array<Local,64>;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
struct Tables{
    int ddt[32][32]{},lat[32][32]{};
    std::array<std::vector<std::pair<int,Real>>,32> kernel;
    std::array<std::vector<std::pair<int,int>>,2048> masks;
    Tables(){
        for(int a=0;a<32;++a)for(int x=0;x<32;++x)++ddt[a][chi5(x)^chi5(x^a)];
        for(int v=0;v<32;++v)for(int u=0;u<32;++u){
            for(int x=0;x<32;++x)lat[v][u]+=parity((v&chi5(x))^(u&x))?-1:1;
            if(lat[v][u])kernel[v].push_back({u,Real(lat[v][u]*lat[v][u])/1024});
        }
        for(int j=0;j<64;++j)for(int v=1;v<32;++v){
            auto a=transpose(embed(j,v));
            for(int k=0;k<64;++k)if(column(a,k))masks[j*32+v].push_back({k,column(a,k)});
        }
    }
    Corr step_chi(const Corr& x)const{
        Corr y{};
        for(int j=0;j<64;++j)for(int v=0;v<32;++v)
            for(auto [u,w]:kernel[v])y[j][v]+=w*x[j][u];
        return y;
    }
    Corr step_linear(const Corr& x)const{
        Corr y{};for(int j=0;j<64;++j){y[j][0]=1;for(int v=1;v<32;++v){
            Real c=1;for(auto [k,u]:masks[j*32+v]){c*=x[k][u];if(c==0)break;}y[j][v]=c;
        }}return y;
    }
    std::vector<Corr> estimate(int pattern,int R)const{
        std::vector<Corr> result(R+1);
        for(int b=1;b<32;++b)if(ddt[pattern][b]){
            const Real p=Real(ddt[pattern][b])/32;
            const auto d=linear(embed(0,b));Corr x{};
            for(int j=0;j<64;++j)for(int v=0;v<32;++v)x[j][v]=parity(v&column(d,j))?-1:1;
            for(int r=2;r<=R;++r){
                x=step_chi(x);
                for(int j=0;j<64;++j)for(int v=1;v<32;++v)result[r][j][v]+=p*x[j][v];
                if(r<R)x=step_linear(x);
            }
        }return result;
    }
};
void selftest(){
    for(int k=0;k<64;++k){auto a=sample(19,k),b=sample(32,k);
        require(dot(linear(a),b)==dot(a,transpose(b)),"transpose identity");
        for(int R=1;R<=12;++R){
            hsilu_state_t s{};for(int i=0;i<5;++i)s.x[i]=a[i];
            for(int r=12-R;r<12;++r)hsilu_round(&s,RC[r]);
            auto out=full(a,R);for(int i=0;i<5;++i)require(out[i]==s.x[i],"reference mismatch");
            auto shifted=a;shifted[0]^=RC[12-R];
            require(linear(core(shifted,R))==out,"core boundary mismatch");
        }
    }
    Tables t;for(int a=0;a<32;++a){int sum=0;for(auto x:t.ddt[a])sum+=x;require(sum==32,"DDT sum");}
    std::cout<<"PASS reference rounds 1..12, core/output boundary, transpose, DDT\n";
}
void search(const std::string& path,int maxR){
    require(maxR>=2&&maxR<=12,"invalid max rounds");
    require(!std::filesystem::exists(path),"output exists");
    std::ofstream out(path);require(out.good(),"cannot create output");Tables tables;
    out<<"rounds,input_pattern,output_col,output_pattern,rho\n"<<std::setprecision(20);
    for(int a=1;a<32;++a){auto estimates=tables.estimate(a,maxR);
        for(int r=2;r<=maxR;++r)for(int j=0;j<64;++j)for(int v=1;v<32;++v)
            out<<r<<','<<a<<','<<j<<','<<v<<','<<estimates[r][j][v]<<'\n';
        std::cerr<<"HNW zero-mask input "<<a<<"/31\n";
    }
}
struct Job{std::string id;int rounds;State delta,mask;};
State parse(const std::string& s){State a{};std::stringstream in(s);std::string token;for(int i=0;i<5;++i){require(bool(std::getline(in,token,';')),"five words required");a[i]=std::stoull(token,nullptr,16);}return a;}
// Exact two-chi correlation for a difference in one input column.
// Condition on its five input bits. All other first-chi output bits are
// independent uniform; the second-chi derivative is affine in L(first output).
int exact_prefix2(int col,int pattern,const State& gamma,uint64_t second_rc){
    int total=0;
    for(int x=0;x<32;++x){
        const int y=chi5(x),b=y^chi5(x^pattern);
        const auto delta=linear(embed(col,b));State a{},constant{};
        for(int i=0;i<5;++i){
            a[i]=(gamma[(i+4)%5]&delta[(i+1)%5])^(gamma[(i+3)%5]&delta[(i+4)%5]);
            constant[i]=delta[i]^delta[(i+2)%5]^(delta[(i+1)%5]&delta[(i+2)%5]);
        }
        const auto pulled=transpose(a);bool external=false;
        for(auto w:pulled)external|=bool(w&~(UINT64_C(1)<<col));
        if(external)continue;
        const int sign=dot(gamma,constant)^__builtin_parityll(a[0]&second_rc)^parity(column(pulled,col)&y);
        total+=sign?-1:1;
    }
    return total;
}
void prefix_scan(const std::string& source,const std::string& output){
    require(!std::filesystem::exists(output),"output exists");
    std::ifstream in(source);std::ofstream out(output);
    require(in.good()&&out.good(),"prefix scan file error");
    out<<"head_id,input_col,input_pattern,prefix_numerator,prefix_denominator,head_weight\n";
    std::string id,mask;int weight,count=0;
    while(in>>id>>weight>>mask){
        const auto gamma=parse(mask);
        for(int c=0;c<64;++c)for(int p=1;p<32;++p){
            const auto num=exact_prefix2(c,p,gamma,RC[9]); // R=4, chi indices 8,9,10,11.
            if(num)out<<id<<','<<c<<','<<p<<','<<num<<",32,"<<weight<<'\n';
        }
        if(++count%32==0)std::cerr<<"PREFIX HEADS "<<count<<'\n';
    }
}
void verify(const std::string& jobs,const std::string& output,uint64_t n,uint64_t seed,int threads){
    require(n>0&&threads>0&&threads<=32,"invalid sampling resources");
    require(!std::filesystem::exists(output),"output exists");
    std::ifstream in(jobs);require(in.good(),"cannot read jobs");std::vector<Job> work;
    std::string id,d,m;int R;while(in>>id>>R>>d>>m){require(R>=1&&R<=12,"invalid rounds");work.push_back({id,R,parse(d),parse(m)});}
    require(in.eof()&&!work.empty(),"invalid/empty jobs");
    std::ofstream out(output);require(out.good(),"cannot create output");
    out<<"id,rounds,samples,seed,n0,n1,rho\n"<<std::setprecision(20);
    for(const auto& job:work){std::vector<uint64_t> counts(threads);std::vector<std::thread> workers;
        for(int t=0;t<threads;++t)workers.emplace_back([&,t]{uint64_t n0=0;
            const auto q=n/threads,rem=n%threads,begin=q*t+std::min<uint64_t>(t,rem),end=begin+q+(uint64_t(t)<rem);
            for(uint64_t k=begin;k<end;++k){auto a=sample(seed,k);const auto b=core(a,job.rounds)^core(a^job.delta,job.rounds);n0+=!dot(b,job.mask);}counts[t]=n0;
        });for(auto& t:workers)t.join();uint64_t n0=0;for(auto c:counts)n0+=c;
        out<<job.id<<','<<job.rounds<<','<<n<<','<<seed<<','<<n0<<','<<n-n0<<','<<(2*Real(n0)-Real(n))/Real(n)<<'\n';out.flush();
        std::cerr<<"VERIFIED "<<job.id<<'\n';
    }
}
int main(int argc,char**argv)try{
    require(argc>=2,"commands: selftest | search OUT MAX_ROUNDS | verify JOBS OUT N SEED THREADS");
    std::string mode=argv[1];
    if(mode=="selftest"&&argc==2)selftest();
    else if(mode=="search"&&argc==4)search(argv[2],std::stoi(argv[3]));
    else if(mode=="prefix-scan"&&argc==4)prefix_scan(argv[2],argv[3]);
    else if(mode=="prefix-eval"&&argc==6)std::cout<<exact_prefix2(std::stoi(argv[2]),std::stoi(argv[3]),parse(argv[4]),std::stoull(argv[5],nullptr,16))<<'\n';
    else if(mode=="verify"&&argc==7)verify(argv[2],argv[3],std::stoull(argv[4]),std::stoull(argv[5]),std::stoi(argv[6]));
    else throw std::runtime_error("invalid command arguments");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

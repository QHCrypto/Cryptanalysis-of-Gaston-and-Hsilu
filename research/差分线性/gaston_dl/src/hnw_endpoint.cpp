#include "estimator.hpp"
#include <gmpxx.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <iostream>
#include <limits>

// Exact dyadic arithmetic: the independence approximation is in the model,
// never in floating point or in a threshold for small correlations.
using namespace gdl;
using Int = mpz_class;
struct Dyadic {
    Int n=0; int e=0;
    Dyadic()=default;
    Dyadic(Int a,int b=0):n(std::move(a)),e(b){norm();}
    void norm(){if(n==0){e=0;return;}Int a=n;if(a<0)a=-a;int k=std::min(e,int(mpz_scan1(a.get_mpz_t(),0)));n>>=k;e-=k;}
    Dyadic operator*(const Dyadic& b)const{return Dyadic(n*b.n,e+b.e);}
    Dyadic& operator+=(const Dyadic& b){int k=std::max(e,b.e);n=(n<<(k-e))+(b.n<<(k-b.e));e=k;norm();return *this;}
    bool operator==(const Dyadic& b)const{return n==b.n&&e==b.e;}
};
void json(std::ostream& out,const Dyadic& q){
    out<<"{\"numerator\":\""<<q.n<<"\",\"denominator\":\""<<(Int(1)<<q.e)<<"\",\"denominator_exponent\":"<<q.e<<",\"log2_abs\":";
    if(q.n==0)out<<"null";
    else{Int a=q.n;if(a<0)a=-a;int shift=std::max(0,int(mpz_sizeinbase(a.get_mpz_t(),2))-1-62);out<<std::setprecision(18)<<(std::log2(Int(a>>shift).get_d())+shift-q.e);}
    out<<'}';
}
struct Column {
    std::vector<Int> n; int e=0;
    explicit Column(int size=32):n(size){}
    void norm(){int k=e;for(const auto& v:n)if(v!=0){Int a=v;if(a<0)a=-a;k=std::min(k,int(mpz_scan1(a.get_mpz_t(),0)));}if(k){for(auto&v:n)v>>=k;e-=k;}}
    Dyadic get(int i)const{return Dyadic(n[i],e);}
};
using Vector = std::vector<Column>;
void fwt(std::vector<Int>& x){for(size_t half=1;half<x.size();half*=2)for(size_t base=0;base<x.size();base+=2*half)for(size_t j=0;j<half;++j){Int a=x[base+j],b=x[base+j+half];x[base+j]=a+b;x[base+j+half]=a-b;}}

struct HNW {
    Tables tables; LinearMasks masks;
    bool full; int size;
    explicit HNW(bool f):full(f),size(f?1024:32){}
    Column first(int a,int b)const{
        Column q(size);
        if(!full){for(int v=0;v<32;++v)q.n[v]=parity(v&b)?-1:1;return q;}
        const int count=tables.ddt[a][b];assert(count && !(count&(count-1)));q.e=__builtin_ctz(unsigned(count));
        for(int x=0;x<32;++x)if((chi5(x)^chi5(x^a))==unsigned(b))q.n[chi5(x)*32+b]+=1;
        fwt(q.n);q.norm();return q;
    }
    Column nonlinear(const Column& x)const{
        Column y(size);y.e=x.e+10;
        if(full){
            auto p=x.n;fwt(p);
            for(int a=0;a<32;++a)for(int d=0;d<32;++d)y.n[chi5(a)*32+(chi5(a)^chi5(a^d))]=p[a*32+d];
            fwt(y.n);
        }else{
            for(int v=0;v<32;++v)for(int u=0;u<32;++u)if(tables.lat[v][u])y.n[v]+=tables.lat[v][u]*tables.lat[v][u]*x.n[u];
        }
        y.norm();return y;
    }
    Vector linear_layer(const Vector& x)const{
        Vector y(64,Column(size));
        for(int col=0;col<64;++col){
            std::vector<Dyadic> entries(size);entries[0]=Dyadic(1);int maxe=0;
            for(int v0=0;v0<(full?32:1);++v0)for(int v1=0;v1<32;++v1){
                if(!v0&&!v1)continue;
                const auto& a=masks.terms[col*32+v0];const auto& b=masks.terms[col*32+v1];
                size_t ia=0,ib=0;Int value=1;int exponent=0;
                while(ia<a.size()||ib<b.size()){
                    int k=std::min(ia<a.size()?a[ia].first:64,ib<b.size()?b[ib].first:64),u0=0,u1=0;
                    if(ia<a.size()&&a[ia].first==k)u0=a[ia++].second;
                    if(ib<b.size()&&b[ib].first==k)u1=b[ib++].second;
                    const auto& f=x[k].n[u0*32+u1];
                    if(f==0){value=0;break;}value*=f;exponent+=x[k].e;
                }
                entries[v0*32+v1]=Dyadic(value,exponent);maxe=std::max(maxe,entries[v0*32+v1].e);
            }
            y[col].e=maxe;
            for(int v=0;v<size;++v)y[col].n[v]=entries[v].n<<(maxe-entries[v].e);
            y[col].norm();
        }
        return y;
    }
    void constant(Vector& x,uint64_t rc)const{
        if(!full)return;
        for(int j=0;j<64;++j)if((rc>>j)&1)for(int a=1;a<32;a+=2)for(int b=0;b<32;++b)x[j].n[a*32+b]*=-1;
    }
    bool absorbed(const Vector& x)const{
        for(const auto& col:x){assert(col.get(0)==Dyadic(1));for(int i=1;i<size;++i)if(col.n[i]!=0)return false;}
        return true;
    }
    Dyadic output(const Vector& x,const State& mask)const{
        Dyadic q(1);for(int j=0;j<64;++j){q=q*x[j].get(column(mask,j));if(q.n==0)break;}return q;
    }
    // A branch-independent certificate, not another independence assumption.
    // For a fixed first DDT branch the difference is deterministic, so the
    // support in the base-value mask is independent of the difference mask.
    int uniform_base_certificate(const State& delta)const{
        assert(full);
        bool possible[64][32]{};
        for(int j=0;j<64;++j){int a=column(delta,j);
            for(int b=0;b<32;++b)if(tables.ddt[a][b]){
                auto q=first(a,b);
                for(int u=0;u<32;++u)possible[j][u]|=q.n[u*32]!=0;
            }
        }
        int killed=0;
        for(int j=0;j<64;++j)for(int a=1;a<32;++a){
            for(auto [k,u]:masks.terms[j*32+a])if(!possible[k][u]){++killed;break;}
        }
        return killed;
    }
};

void selftest(){
    HNW full(true),zero(false);std::mt19937 rng(20261016);
    // FWHT implementation against the independent two-wise LAT formula.
    for(int trial=0;trial<4;++trial){
        Column q(1024);q.e=6;for(int n=0;n<64;++n)q.n[rng()%1024]+=1;fwt(q.n);q.norm();auto y=full.nonlinear(q);
        for(int a=0;a<32;++a)for(int b=0;b<32;++b){Int sum=0;for(int u0=0;u0<32;++u0)for(int u1=0;u1<32;++u1){int c=full.tables.lat[a^b][u0^u1]*full.tables.lat[b][u1];if(c)sum+=c*q.n[u0*32+u1];}assert(y.get(a*32+b)==Dyadic(sum,q.e+10));}
    }
    // Enumerate every (input difference, output difference) branch directly.
    int branches=0;
    for(int a=0;a<32;++a)for(int b=0;b<32;++b)if(full.tables.ddt[a][b]){
        auto x=full.first(a,b);auto z=zero.first(a,b);for(int v=0;v<32;++v)assert(x.get(v)==z.get(v));++branches;
    }
    // The optimized zero-mask nonlinear update agrees with the a=0 slice of
    // the full update when all nonzero value-mask input coordinates are zero.
    for(int d=0;d<32;++d){Column q(1024);auto z=zero.first(1,d);for(int v=0;v<32;++v)q.n[v]=z.n[v];auto fy=full.nonlinear(q),zy=zero.nonlinear(z);for(int v=0;v<32;++v)assert(fy.get(v)==zy.get(v));}
    // Independent deterministic-state check of every two-wise linear mask.
    State value{},difference{};for(int i=0;i<5;++i){value[i]=(uint64_t(rng())<<32)|rng();difference[i]=(uint64_t(rng())<<32)|rng();}
    Vector deterministic(64,Column(1024));
    for(int j=0;j<64;++j)for(int a=0;a<32;++a)for(int b=0;b<32;++b)deterministic[j].n[a*32+b]=parity((a&column(value,j))^(b&column(difference,j)))?-1:1;
    auto propagated=full.linear_layer(deterministic);auto lv=linear(value),ld=linear(difference);
    for(int j=0;j<64;++j)for(int a=0;a<32;++a)for(int b=0;b<32;++b)assert(propagated[j].get(a*32+b)==Dyadic(parity((a&column(lv,j))^(b&column(ld,j)))?-1:1));
    full.constant(propagated,RC[8]);lv[0]^=RC[8];
    for(int j=0;j<64;++j)for(int a=0;a<32;++a)for(int b=0;b<32;++b)assert(propagated[j].get(a*32+b)==Dyadic(parity((a&column(lv,j))^(b&column(ld,j)))?-1:1));
    // Single-active output agrees with the previous Table-13 implementation.
    Estimator old;
    for(int d:{1,4,7,31}){int col=3;auto reference=old.estimate(col,d,3);Dyadic total;
        for(int b=0;b<32;++b)if(zero.tables.ddt[d][b]){Vector q(64,Column(32));for(int j=0;j<64;++j)q[j]=zero.first(j==col?d:0,j==col?b:0);q=zero.linear_layer(q);for(auto&v:q)v=zero.nonlinear(v);q=zero.linear_layer(q);for(auto&v:q)v=zero.nonlinear(v);total+=Dyadic(zero.tables.ddt[d][b],5)*zero.output(q,embed(19,8));}
        long double exact=std::ldexp(total.n.get_d(),-total.e);assert(exact==reference.rounds[3][19][8]);
    }
    std::cout<<"{\"status\":\"PASS\",\"two_wise_LAT_coordinates_checked\":4096,\"conditional_DDT_branches_checked\":"<<branches<<",\"zero_slice_coordinates_checked\":1024,\"two_wise_linear_coordinates_checked\":65536,\"round_constant_coordinates_checked\":65536,\"prior_estimator_checks\":4}\n";
}

int main(int argc,char**argv){
 try{
    if(argc==2&&std::string(argv[1])=="--selftest"){selftest();return 0;}
    if(argc!=15&&argc!=17&&argc!=22)throw std::runtime_error("usage: hnw_endpoint full|zero R START DELTA[5] MASK[5] OUTPUT.json [BEGIN END [PREFIX2_MASK[5]]]");
    std::string mode=argv[1];if(mode!="full"&&mode!="zero")throw std::runtime_error("invalid mode");
    int R=std::stoi(argv[2]),start=std::stoi(argv[3]);if(R<1||start<0||start+R>12)throw std::runtime_error("invalid interval");
    State delta{},mask{};for(int i=0;i<5;++i){delta[i]=std::stoull(argv[4+i],nullptr,16);mask[i]=std::stoull(argv[9+i],nullptr,16);}
    if(delta==State{}||mask==State{})throw std::runtime_error("nonzero endpoints required");
    HNW hnw(mode=="full");std::vector<int> cols;uint64_t expected=1;
    for(int j=0;j<64;++j)if(column(delta,j)){cols.push_back(j);int n=0;for(int b=0;b<32;++b)n+=hnw.tables.ddt[column(delta,j)][b]!=0;if(expected>std::numeric_limits<uint64_t>::max()/n)throw std::runtime_error("branch count overflow");expected*=n;}
    uint64_t begin=argc>=17?std::stoull(argv[15]):0,end=argc>=17?std::stoull(argv[16]):expected;
    if(begin>=end||end>expected)throw std::runtime_error("invalid branch range");
    const bool prefix_readout=argc==22;
    State prefix_mask{};if(prefix_readout){for(int i=0;i<5;++i)prefix_mask[i]=std::stoull(argv[17+i],nullptr,16);if(R<2||prefix_mask==State{})throw std::runtime_error("invalid prefix readout");}
    // The legacy and range-only modes deliberately keep the unoptimized path
    // so the certified shortcut can be checked against an independent run.
    int killed=prefix_readout&&hnw.full?hnw.uniform_base_certificate(delta):-1;
    bool fast_second=prefix_readout&&(!hnw.full||killed==64*31);
    std::vector<Column> second_lookup;
    if(fast_second)for(int d=0;d<32;++d){Column q(hnw.size);for(int v=0;v<32;++v)q.n[v]=parity(v&d)?-1:1;second_lookup.push_back(hnw.nonlinear(q));}
    std::ofstream out(argv[14]);if(!out)throw std::runtime_error("cannot write output");
    auto began=std::chrono::steady_clock::now();Dyadic total,mass,prefix_total;uint64_t done=0,nonzero=0,prefix_nonzero=0,ordinal=0;std::map<int,uint64_t> extinction;
    out<<"{\"mode\":\""<<mode<<"\",\"rounds\":"<<R<<",\"start_round_index\":"<<start<<",\"rd\":1,\"rl\":0,\"input_words\":\""<<hex_state(delta)<<"\",\"output_words\":\""<<hex_state(mask)<<"\",\"arithmetic\":\"exact dyadic GMP integers\",\"expected_first_DDT_branches\":"<<expected<<",\"branch_begin\":"<<begin<<",\"branch_end\":"<<end<<",\"first_linear_killed_value_masks\":"<<killed<<",\"certified_fast_second_round\":"<<(fast_second?"true":"false");
    if(prefix_readout)out<<",\"prefix_rounds\":2,\"prefix_mask_words\":\""<<hex_state(prefix_mask)<<'"';
    out<<",\"branches\":[\n";
    State diff{};
    std::function<void(size_t,Dyadic)> visit=[&](size_t depth,Dyadic probability){
        if(depth<cols.size()){
            int col=cols[depth],a=column(delta,col);
            for(int b=0;b<32;++b)if(hnw.tables.ddt[a][b]){diff=diff^embed(col,column(diff,col))^embed(col,b);visit(depth+1,probability*Dyadic(hnw.tables.ddt[a][b],5));}
            return;
        }
        uint64_t index=ordinal++;if(index<begin||index>=end)return;
        Vector q(64,Column(hnw.size));Dyadic prefix_rho;
        if(fast_second){auto d=linear(diff);for(int j=0;j<64;++j)q[j]=second_lookup[column(d,j)];prefix_rho=hnw.output(q,prefix_mask);}
        else for(int j=0;j<64;++j)q[j]=hnw.first(column(delta,j),column(diff,j));
        int absorbed_round=0;
        for(int r=fast_second?3:2;r<=R;++r){q=hnw.linear_layer(q);hnw.constant(q,RC[start+r-1]);
            if(hnw.absorbed(q)){absorbed_round=r;break;}
            for(auto&v:q)v=hnw.nonlinear(v);
            if(r==2&&prefix_readout)prefix_rho=hnw.output(q,prefix_mask);
            if(hnw.absorbed(q)){absorbed_round=r;break;}
        }
        Dyadic rho=hnw.output(q,mask),contribution=probability*rho;total+=contribution;mass+=probability;nonzero+=rho.n!=0;if(absorbed_round)++extinction[absorbed_round];
        Dyadic prefix_contribution=probability*prefix_rho;if(prefix_readout){prefix_total+=prefix_contribution;prefix_nonzero+=prefix_rho.n!=0;}
        if(done++)out<<",\n";
        out<<"{\"branch_index\":"<<index<<",\"first_output_difference\":\""<<hex_state(diff)<<"\",\"probability\":";json(out,probability);out<<",\"conditional_model_correlation\":";json(out,rho);out<<",\"contribution\":";json(out,contribution);out<<",\"absorbed_by_round\":"<<absorbed_round;
        if(prefix_readout){out<<",\"prefix_conditional_model_correlation\":";json(out,prefix_rho);out<<",\"prefix_contribution\":";json(out,prefix_contribution);}out<<'}';
        if(done%256==0)std::cerr<<mode<<" R"<<R<<" branches "<<done<<'/'<<(end-begin)<<'\n';
    };
    visit(0,Dyadic(1));assert(ordinal==expected&&done==end-begin);if(begin==0&&end==expected)assert(mass==Dyadic(1));
    out<<"\n],\"probability_mass\":";json(out,mass);
    if(prefix_readout){out<<",\"prefix_correlation\":";json(out,prefix_total);out<<",\"nonzero_prefix_branch_correlations\":"<<prefix_nonzero;}
    out<<",\"correlation\":";json(out,total);out<<",\"bias\":";json(out,total*Dyadic(1,1));out<<",\"completed_branches\":"<<done<<",\"nonzero_branch_correlations\":"<<nonzero<<",\"absorbed_histogram\":{";bool comma=false;for(auto [r,count]:extinction){if(comma)out<<',';out<<'"'<<r<<"\":"<<count;comma=true;}out<<"},\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"status\":\"completed\"}\n";
    std::cout<<mode<<" R"<<R<<" ";json(std::cout,total);std::cout<<" branches="<<done<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}

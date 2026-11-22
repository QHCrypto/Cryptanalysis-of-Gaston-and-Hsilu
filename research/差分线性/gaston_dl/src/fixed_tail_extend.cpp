#include "gaston.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <vector>

// Fixed-tail backward search. All feasibility decisions are exact over GF(2).
// Minimum-weight search is exhaustive only in the explicitly counted domains.
using V = gdl::State;
V xx(V a, const V &b) { for (int i=0;i<5;++i) a[i]^=b[i]; return a; }
bool nz(const V &a) { return a[0]|a[1]|a[2]|a[3]|a[4]; }
V unit(int k) { V a{}; a[k/64]=uint64_t(1)<<(k%64); return a; }
int wt(const V &a) { int n=0; for(auto w:a)n+=__builtin_popcountll(w); return n; }
int ac(const V &a) { return __builtin_popcountll(a[0]|a[1]|a[2]|a[3]|a[4]); }
int pivot(const V &a) {
    for(int r=4;r>=0;--r) if(a[r]) return 64*r+63-__builtin_clzll(a[r]);
    return -1;
}
uint64_t fingerprint(const V &v) {
    uint64_t z=0; for(auto w:v) z=gdl::splitmix_final(z^w^0x9e3779b97f4a7c15ULL);
    return z;
}
V deriv_linear(const V &a, const V &x) {
    V y{};
    for(int i=0;i<5;++i)
        y[i]=(a[(i+2)%5]&x[(i+1)%5])^(a[(i+1)%5]&x[(i+2)%5]);
    return y;
}
V derivative_coeff(const V &delta, const V &mask) {
    V h{};
    for(int k=0;k<5;++k)
        h[k]=(mask[(k+4)%5]&delta[(k+1)%5])^
             (mask[(k+3)%5]&delta[(k+4)%5]);
    return h;
}
struct Basis {
    std::array<V,320> rows{}, witnesses{};
    int rank=0;
    bool add(V a,V witness) {
        while(nz(a)) {
            int k=pivot(a);
            if(nz(rows[k])) { a=xx(a,rows[k]); witness=xx(witness,witnesses[k]); }
            else { rows[k]=a; witnesses[k]=witness; ++rank; return true; }
        }
        return false;
    }
    bool solve(V a,V &answer) const {
        answer={};
        while(nz(a)) {
            int k=pivot(a);
            if(!nz(rows[k])) return false;
            a=xx(a,rows[k]); answer=xx(answer,witnesses[k]);
        }
        return true;
    }
};
int ddt[32][32], min_weight[32], diff_weight[32];
void tables() {
    for(int a=0;a<32;++a) for(int x=0;x<32;++x)
        ++ddt[a][gdl::chi5(x)^gdl::chi5(x^a)];
    for(int b=0;b<32;++b) {
        int best=0; for(int a=0;a<32;++a)best=std::max(best,ddt[a][b]);
        min_weight[b]=5-__builtin_ctz(unsigned(best));
    }
    for(int a=0;a<32;++a) {
        int count=0;
        for(int b=0;b<32;++b) if(ddt[a][b]) {
            if(count && count!=ddt[a][b])throw std::runtime_error("nonuniform derivative image");
            count=ddt[a][b];
        }
        diff_weight[a]=5-__builtin_ctz(unsigned(count));
        std::set<unsigned> span{0};
        const unsigned base=gdl::chi5(a)^gdl::chi5(0);
        for(int k=0;k<5;++k) {
            const unsigned v=gdl::chi5(1<<k)^gdl::chi5(a^(1<<k))^base;
            auto old=span; for(auto w:old)span.insert(w^v);
        }
        for(int b=0;b<32;++b)
            if(bool(ddt[a][b])!=bool(span.count(b^base)))
                throw std::runtime_error("affine image validation failed");
    }
    // Independent truth-table check of the wordwise connector equation.
    for(int a=0;a<32;++a)for(int m=0;m<32;++m) {
        V A=gdl::embed(0,a), M=gdl::embed(0,m);
        const unsigned h=gdl::column(derivative_coeff(A,M),0);
        const int c=gdl::parity(m&(gdl::chi5(a)^gdl::chi5(0)));
        for(int x=0;x<32;++x)
            if(gdl::parity(m&(gdl::chi5(x)^gdl::chi5(x^a)))!=(c^gdl::parity(h&x)))
                throw std::runtime_error("connector validation failed");
    }
}
int backward_weight(const V &b) {
    int s=0; for(int j=0;j<64;++j)s+=min_weight[gdl::column(b,j)]; return s;
}
int forward_weight(const V &a) {
    int s=0; for(int j=0;j<64;++j)s+=diff_weight[gdl::column(a,j)]; return s;
}
struct Candidate {
    V first{},second{};
    int weight=0;
    std::string source;
};
bool less_candidate(const Candidate&a,const Candidate&b) {
    if(a.weight!=b.weight)return a.weight<b.weight;
    auto ha=fingerprint(a.first), hb=fingerprint(b.first);
    if(ha!=hb)return ha<hb;
    if(a.first!=b.first)return a.first<b.first;
    return a.second<b.second;
}
int bucket(const V&a) {
    int n=wt(a),c=ac(a);
    int cb=c<=4?c:(c<=8?5:(c<=16?6:(c<=32?7:8)));
    return cb*8+std::min(n,7);
}
struct Keeper {
    size_t per_bucket;
    std::map<int,std::vector<Candidate>> groups;
    explicit Keeper(size_t n):per_bucket(n){}
    void add(Candidate a) {
        auto &v=groups[bucket(a.first)];
        if(v.size()>=per_bucket && !less_candidate(a,v.back()))return;
        auto p=std::lower_bound(v.begin(),v.end(),a,less_candidate);
        if(p!=v.end() && p->first==a.first && p->second==a.second)return;
        v.insert(p,a);
        if(v.size()>per_bucket)v.pop_back();
    }
    std::vector<Candidate> all() const {
        std::vector<Candidate> a;
        for(auto &g:groups)a.insert(a.end(),g.second.begin(),g.second.end());
        std::sort(a.begin(),a.end(),less_candidate);return a;
    }
};
struct Search {
    V gamma;
    std::array<V,320> columns;
    Keeper one,two;
    std::array<uint64_t,5> counts{};
    uint64_t null_trials=0, two_checked=0,two_feasible=0;
    int rank=0, bounded_min=10000, unrestricted_min=10000,two_min=10000;
    explicit Search(V g,size_t keep):gamma(g),one(keep),two(keep) {
        Basis b;
        for(int k=0;k<320;++k) {
            columns[k]=derivative_coeff(gdl::linear(unit(k)),gamma);
            b.add(columns[k],unit(k));
        }
        rank=b.rank;
    }
    V apply(const V&a)const {
        V b{};
        for(int r=0;r<5;++r) {
            uint64_t w=a[r];
            while(w){int j=__builtin_ctzll(w);b=xx(b,columns[64*r+j]);w&=w-1;}
        }
        return b;
    }
    void accept_one(const V&v,const std::string&source,bool bounded) {
        if(!nz(v)||nz(apply(v)))throw std::runtime_error("invalid kernel vector");
        const int w=backward_weight(v);
        if(bounded){++counts[wt(v)];bounded_min=std::min(bounded_min,w);}
        unrestricted_min=std::min(unrestricted_min,w);
        one.add({v,{},w,source});
    }
    void sparse(int bits=320) {
        std::map<V,std::vector<std::pair<int,int>>> pairs;
        for(int a=0;a<bits;++a) {
            if(!nz(columns[a]))accept_one(unit(a),"exhaustive_bit_weight_1_to_4",true);
            for(int b=a+1;b<bits;++b) {
                V syndrome=xx(columns[a],columns[b]);
                pairs[syndrome].emplace_back(a,b);
                if(!nz(syndrome))accept_one(xx(unit(a),unit(b)),"exhaustive_bit_weight_1_to_4",true);
            }
        }
        // Exactly one canonical partition per support: a | (b,c), (a,b) | (c,d).
        for(int a=0;a<bits;++a) {
            auto it=pairs.find(columns[a]);if(it==pairs.end())continue;
            for(auto [b,c]:it->second) if(a<b)
                accept_one(xx(unit(a),xx(unit(b),unit(c))),"exhaustive_bit_weight_1_to_4",true);
        }
        for(auto &item:pairs) {
            const auto &p=item.second;
            for(size_t i=0;i<p.size();++i)for(size_t j=i+1;j<p.size();++j)
                if(p[i].second<p[j].first)
                    accept_one(xx(xx(unit(p[i].first),unit(p[i].second)),
                                  xx(unit(p[j].first),unit(p[j].second))),
                               "exhaustive_bit_weight_1_to_4",true);
        }
    }
    void randomized_kernel(unsigned repetitions,std::mt19937_64 &rng) {
        std::array<int,320> order;std::iota(order.begin(),order.end(),0);
        for(unsigned trial=0;trial<repetitions;++trial) {
            std::shuffle(order.begin(),order.end(),rng);
            Basis basis;std::vector<V> kernel;
            for(int k:order) {
                V v{};
                if(basis.solve(columns[k],v)) {
                    v=xx(v,unit(k));kernel.push_back(v);
                    accept_one(v,"random_information_set_kernel",false);++null_trials;
                }else basis.add(columns[k],unit(k));
            }
            // Pairwise local improvement is heuristic, never a completeness claim.
            for(int q=0;q<64 && !kernel.empty();++q) {
                V v=xx(kernel[rng()%kernel.size()],kernel[rng()%kernel.size()]);
                if(!nz(v))continue;
                int best=backward_weight(v);
                for(auto &u:kernel) {
                    V next=xx(v,u);if(!nz(next))continue;
                    int w=backward_weight(next);
                    if(w<best){v=next;best=w;}
                }
                accept_one(v,"random_kernel_greedy",false);++null_trials;
            }
        }
    }
    void try_two(const V&b1,const std::string&source) {
        ++two_checked;
        const V a2=gdl::linear(b1);
        const V base=xx(gdl::chi(a2),gdl::chi(V{}));
        Basis basis;
        for(int r=0;r<5;++r) for(int j=0;j<64;++j) {
            V direction=deriv_linear(a2,unit(64*r+j));
            if(nz(direction))basis.add(apply(direction),direction);
        }
        V correction{};
        if(!basis.solve(apply(base),correction))return;
        const V b2=xx(base,correction);
        if(nz(apply(b2)))throw std::runtime_error("bad affine intersection");
        for(int j=0;j<64;++j)
            if(!ddt[gdl::column(a2,j)][gdl::column(b2,j)])
                throw std::runtime_error("invalid second DDT");
        ++two_feasible;
        int w=backward_weight(b1)+forward_weight(a2);two_min=std::min(two_min,w);
        two.add({b1,b2,w,source});
    }
};
void state_json(std::ostream&out,const V&a) {
    out<<'[';char buf[32];
    for(int i=0;i<5;++i){if(i)out<<',';std::snprintf(buf,sizeof(buf),"\"%016llx\"",(unsigned long long)a[i]);out<<buf;}
    out<<']';
}
void candidates_json(std::ostream&out,const std::vector<Candidate>&records) {
    out<<'[';bool first=true;
    for(const auto&r:records) {
        if(!first)out<<',';
        first=false;
        out<<"{\"first_output\":";state_json(out,r.first);
        out<<",\"second_output\":";state_json(out,r.second);
        out<<",\"minimum_ddt_weight\":"<<r.weight<<",\"source\":\""<<r.source<<"\"}";
    }
    out<<']';
}
int main(int argc,char**argv) {
    try {
        tables();
        if(argc==2 && std::string(argv[1])=="--self-test") {
            V g{0x123,0x456,0x789,0xabc,0xdef};Search s(g,10000);
            // A synthetic rank-four matrix supplies many nontrivial collisions,
            // including supports of weight three and four.
            for(int k=0;k<12;++k){s.columns[k]={};s.columns[k][0]=(11*k+3)&15;}
            s.sparse(12);
            std::array<uint64_t,5> expected{};
            for(unsigned mask=1;mask<(1u<<12);++mask)if(__builtin_popcount(mask)<=4) {
                V v{};v[0]=mask;if(!nz(s.apply(v)))++expected[__builtin_popcount(mask)];
            }
            if(expected!=s.counts)throw std::runtime_error("MITM completeness check failed");
            if(!expected[3] || !expected[4])throw std::runtime_error("weak MITM test domain");
            Search real(g,10000);
            std::mt19937_64 rng(20261017);
            for(int k=0;k<64;++k) {
                V a{};for(auto&w:a)w=rng();
                V x{};for(auto&w:x)w=rng();
                const V actual=xx(gdl::chi(x),gdl::chi(xx(x,a)));
                if(actual!=xx(xx(gdl::chi(a),gdl::chi(V{})),deriv_linear(a,x)))
                    throw std::runtime_error("full-state affine derivative failed");
                real.try_two(a,"selftest");
            }
            std::printf("PASS: all chi DDT affine images, 32768 connector truth values, MITM vs exhaustive, full-state affine intersections.\n");
            return 0;
        }
        if(argc!=6)throw std::runtime_error("usage: fixed_tail_extend gamma.txt pq_states.txt output.json random_trials keep_per_bucket");
        V gamma{};std::ifstream input(argv[1]);
        for(auto&w:gamma)if(!(input>>std::hex>>w))throw std::runtime_error("invalid gamma");
        std::ifstream pq(argv[2]);if(!pq)throw std::runtime_error("cannot open PQ seeds");
        std::ifstream existing(argv[3]);if(existing)throw std::runtime_error("output already exists");
        const auto began=std::chrono::steady_clock::now();
        std::mt19937_64 rng(2026101701);
        Search search(gamma,std::stoul(argv[5]));
        std::printf("kernel rank=%d dimension=%d; exhaustive output bit weight <=4\n",search.rank,320-search.rank);std::fflush(stdout);
        search.sparse();
        std::printf("bounded counts: %llu %llu %llu %llu; minimum DDT weight=%d\n",
                    (unsigned long long)search.counts[1],(unsigned long long)search.counts[2],
                    (unsigned long long)search.counts[3],(unsigned long long)search.counts[4],search.bounded_min);std::fflush(stdout);
        search.randomized_kernel(std::stoul(argv[4]),rng);
        auto one=search.one.all();
        std::set<V> first_outputs;
        for(int a=0;a<320;++a) {
            first_outputs.insert(unit(a));
            for(int b=a+1;b<320;++b)first_outputs.insert(xx(unit(a),unit(b)));
        }
        for(int j=0;j<64;++j)for(int a=1;a<32;++a)first_outputs.insert(gdl::embed(j,a));
        V seed{};while(pq>>std::hex>>seed[0]) {
            for(int r=1;r<5;++r)if(!(pq>>seed[r]))throw std::runtime_error("invalid PQ state");
            for(int shift=0;shift<64;++shift)first_outputs.insert(gdl::rotate(seed,shift));
        }
        for(const auto&r:one)first_outputs.insert(r.first);
        // Dense controls ensure the affine intersection method is also exercised
        // outside sparse families; these are not advertised as exhaustive.
        for(int k=0;k<64;++k){V a{};for(auto&w:a)w=rng();first_outputs.insert(a);}
        int n=0;
        for(const auto &b:first_outputs) {
            search.try_two(b,"bounded_first_output_or_PQ_or_kernel_or_dense_seed");
            if(++n%10000==0){std::printf("two-round checked=%d feasible=%llu best_weight=%d\n",n,(unsigned long long)search.two_feasible,search.two_min);std::fflush(stdout);}
        }
        std::ofstream out(argv[3]);if(!out)throw std::runtime_error("cannot create output");
        out<<"{\"method\":\"GF2 kernel and affine DDT intersection\",\"seed\":2026101701,\"gamma\":";
        state_json(out,gamma);
        out<<",\"connector_constraint_rank\":"<<search.rank<<",\"connector_preimage_dimension\":"<<320-search.rank;
        out<<",\"one_round_exhaustive_output_bit_cap\":4,\"one_round_kernel_counts_by_bit_weight\":[";
        for(int i=0;i<5;++i){if(i)out<<',';out<<search.counts[i];}
        out<<"],\"bounded_minimum_ddt_weight\":";
        if(search.bounded_min==10000)out<<"null";else out<<search.bounded_min;
        out<<",\"best_seen_one_round_weight\":"<<search.unrestricted_min
           <<",\"random_information_set_repetitions\":"<<argv[4]
           <<",\"heuristic_kernel_vectors_examined\":"<<search.null_trials
           <<",\"two_round_first_outputs_checked\":"<<search.two_checked
           <<",\"two_round_feasible_first_outputs\":"<<search.two_feasible
           <<",\"two_round_best_seen_weight\":"<<search.two_min
           <<",\"keep_per_bucket\":"<<argv[5]<<",\"one_round_seeds\":";
        candidates_json(out,one);out<<",\"two_round_seeds\":";candidates_json(out,search.two.all());
        out<<",\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<"}\n";
        if(!out)throw std::runtime_error("write failed");
        std::printf("DONE one-round seeds=%zu two-round seeds=%zu output=%s\n",one.size(),search.two.all().size(),argv[3]);
    }catch(const std::exception&e){std::fprintf(stderr,"%s\n",e.what());return 1;}
}

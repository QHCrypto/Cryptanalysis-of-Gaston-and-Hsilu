#include "datf_middle.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <tuple>
using namespace gdl;

struct Diff { std::string id, family; int w, orbit; State input, output; };
struct Lin { std::string id, family; int rounds, w2; State gamma, output; };
struct Hit { int weight, sign, shift, mw; size_t di, li; uint64_t phase; };
struct Summary {
    std::string family, tail;
    int middle = 0, bound = 0;
    uint64_t evaluated = 0, nonzero = 0;
    bool exhausted = false, certified = false;
    std::vector<Hit> hits;
};
struct PooledLinear {
    std::vector<size_t> indices, scout;
    uint64_t budget = 0;
    Summary summary;
};
bool better(const Hit &a, const Hit &b) {
    return std::tie(a.weight,a.di,a.shift,a.li) < std::tie(b.weight,b.di,b.shift,b.li);
}
void keep(Summary &s, Hit h) {
    ++s.nonzero;
    if (s.hits.size() < 32 || better(h,s.hits.back())) {
        s.hits.push_back(h); std::sort(s.hits.begin(),s.hits.end(),better);
        if (s.hits.size()>32) s.hits.resize(32);
    }
    s.certified = !s.hits.empty() && s.hits[0].weight == s.bound;
}
State read_state(std::istream &f) { State a{}; for(auto &w:a) f>>std::hex>>w; f>>std::dec; return a; }
int main(int argc,char **argv) {
    if(argc<6 || argc>9) { std::cerr<<"usage: datf_pool_match DIFFERENTIALS LINEARS PREFIX M2_REPS_PER_FAMILY L1_PAIR_BUDGET [L1_M2_REPS=64] [L2_PAIR_BUDGET=200000] [L2_M2_REPS=64]\n"; return 2; }
    const size_t m2_limit=std::stoull(argv[4]);
    const uint64_t l1_budget=std::stoull(argv[5]);
    const size_t l1_m2_limit=argc>=7?std::stoull(argv[6]):64;
    const uint64_t l2_budget=argc>=8?std::stoull(argv[7]):200000;
    const size_t l2_m2_limit=argc>=9?std::stoull(argv[8]):64;
    if(!m2_limit || !l1_budget || !l1_m2_limit || !l2_budget || !l2_m2_limit) return 2;
    std::ifstream df(argv[1]),lf(argv[2]);
    if(!df||!lf) return 2;
    std::vector<Diff> ds; std::vector<Lin> ls;
    while(true) { Diff d; if(!(df>>d.id>>d.family>>d.w>>d.orbit))break; d.input=read_state(df);d.output=read_state(df);ds.push_back(d); }
    while(true) { Lin l; if(!(lf>>l.id>>l.family>>l.rounds>>l.w2))break;l.gamma=read_state(lf);l.output=read_state(lf);ls.push_back(l); }
    const std::string prefix=argv[3];
    std::ofstream hits(prefix+"_hits.tsv"),summary(prefix+"_summary.tsv"),selected(prefix+"_selected.tsv");
    hits<<"family\tmiddle\ttail\tdiff_id\trotation\tlinear_id\tweight\tsign\tmiddle_weight\trc_phase_mask\n";
    summary<<"family\tmiddle\ttail\tevaluated\tnonzero\texhausted\toptimum_certified\tupper_bound_weight\tbest_weight\n";
    selected<<"family\tdiff_id\n";
    MiddleTables tables;
    const auto began=std::chrono::steady_clock::now();
    for(const std::string family: {"D1","D2"}) {
        std::vector<size_t> all,chosen;
        for(size_t i=0;i<ds.size();++i)if(ds[i].family==family)all.push_back(i);
        const int min_dw=family=="D1"?2:24;
        // Input files contain an explicit deterministic scouting priority.
        chosen.assign(all.begin(),all.begin()+std::min(m2_limit,all.size()));
        for(size_t i:chosen)selected<<family<<'\t'<<ds[i].id<<'\n';
        for(int middle:{1,2}) {
            std::vector<Summary> stats;
            std::vector<size_t> fixed;
            std::map<std::string,PooledLinear> pools;
            for(size_t j=0;j<ls.size();++j) {
                if(ls[j].family=="L1" || ls[j].family=="L2_2to2")pools[ls[j].family].indices.push_back(j);
                else { fixed.push_back(j); Summary s;s.family=family;s.middle=middle;s.tail=ls[j].id;s.bound=min_dw+ls[j].w2;stats.push_back(s); }
            }
            for(auto &entry:pools) {
                auto &p=entry.second;
                const size_t limit=entry.first=="L1"?l1_m2_limit:l2_m2_limit;
                p.budget=entry.first=="L1"?l1_budget:l2_budget;
                p.scout.assign(p.indices.begin(),p.indices.begin()+std::min(limit,p.indices.size()));
                auto &s=p.summary;s.family=family;s.middle=middle;s.tail=entry.first+"_pool";
                int min_lw=ls[p.indices.front()].w2;
                for(size_t li:p.indices)min_lw=std::min(min_lw,ls[li].w2);
                s.bound=min_dw+min_lw;
            }
            const auto &indices=middle==1?all:chosen;
            // All fixed tails with the same gamma reuse one DATF computation.
            std::map<State,std::vector<size_t>> mask_groups;
            for(size_t k=0;k<fixed.size();++k)mask_groups[ls[fixed[k]].gamma].push_back(k);
            uint64_t visited=0;
            for(size_t di:indices) {
                const auto &d=ds[di];
                std::vector<PooledLinear*> scouting;
                for(auto &entry:pools) {
                    auto &p=entry.second;
                    if(middle==2 && !p.summary.certified &&
                       p.summary.evaluated+d.orbit*p.scout.size()<=p.budget)scouting.push_back(&p);
                }
                for(int shift=0;shift<d.orbit;++shift) {
                    const State delta=linear(rotate(d.output,shift));
                    std::unique_ptr<MiddleTwo> compiled;
                    if(middle==2)compiled=std::make_unique<MiddleTwo>(delta,tables);
                    for(const auto &entry:mask_groups) {
                        const auto c=middle==1?middle_one(delta,entry.first):compiled->evaluate(entry.first);
                        for(size_t k:entry.second) {
                            auto &s=stats[k];const size_t li=fixed[k];const auto &l=ls[li];++s.evaluated;
                            if(c.sign) {
                                int sign=c.sign;
                                if(middle==2 && __builtin_parityll(c.rc_mask & RC[11-l.rounds]))sign=-sign;
                                keep(s,{d.w+l.w2+c.weight,sign,shift,c.weight,di,li,c.rc_mask});
                            }
                        }
                    }
                    // A hard universal |rho_middle|<=1 bound can certify the
                    // best individual splice score; this is not a hull bound.
                    for(auto *p:scouting) {
                        auto &one=p->summary;
                        if(one.certified)continue;
                        for(size_t li:p->scout) {
                            const auto &l=ls[li];
                            if(!one.hits.empty() && d.w+l.w2>one.hits.front().weight)continue;
                            ++one.evaluated;
                            auto c=middle==1?middle_one(delta,l.gamma):compiled->evaluate(l.gamma);
                            if(c.sign) {
                                int sign=c.sign;
                                if(middle==2&&__builtin_parityll(c.rc_mask & RC[11-l.rounds]))sign=-sign;
                                keep(one,{d.w+l.w2+c.weight,sign,shift,c.weight,di,li,c.rc_mask});
                                if(one.certified)break;
                            }
                        }
                    }
                    ++visited;
                }
                if(visited%262144==0)std::cerr<<family<<" middle="<<middle<<" differences="<<visited<<" elapsed="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<'\n';
            }
            for(auto &entry:pools) {
              auto &one=entry.second.summary;
              const auto &ones=entry.second.indices;
              if(middle==1) {
                // For one chi, every nonzero middle correlation is +/-1.
                // Search weight strata in ascending order. The first hit
                // certifies the best score over the complete two pools.
                std::set<int> strata;
                std::set<int> dw;
                std::map<int,std::vector<size_t>> by_weight;
                for(size_t di:all)dw.insert(ds[di].w);
                for(size_t li:ones)by_weight[ls[li].w2].push_back(li);
                for(int w:dw)for(const auto &entry:by_weight)strata.insert(w+entry.first);
                bool found=false;
                for(int stratum:strata) {
                    for(size_t di:all) {
                        const auto &d=ds[di];
                        const auto it=by_weight.find(stratum-d.w);
                        if(it==by_weight.end())continue;
                        const auto &layer=it->second;
                        for(int shift=0;shift<d.orbit && !found;++shift) {
                            State delta=linear(rotate(d.output,shift));
                            for(size_t li:layer) {
                                ++one.evaluated;auto c=middle_one(delta,ls[li].gamma);
                                if(c.sign) { keep(one,{stratum,c.sign,shift,0,di,li,0});found=true;break; }
                            }
                        }
                        if(found)break;
                    }
                    if(found) { one.bound=stratum;one.certified=true;break; }
                }
              }
            // This count includes skipped candidates only when no bound or
            // budget terminated the linear-pool search. Otherwise report a
            // certificate separately, never imply exhaustive enumeration.
            one.exhausted=indices.size()==all.size() && one.evaluated==visited*ones.size();
              stats.push_back(one);
            }
            for(size_t k=0;k<fixed.size();++k)stats[k].exhausted=indices.size()==all.size();
            for(const auto &s:stats) {
                summary<<s.family<<'\t'<<s.middle<<'\t'<<s.tail<<'\t'<<s.evaluated<<'\t'<<s.nonzero<<'\t'<<s.exhausted<<'\t'<<s.certified<<'\t'<<s.bound<<'\t';
                if(!s.hits.empty())summary<<s.hits[0].weight;
                summary<<'\n';
                for(const auto &h:s.hits)hits<<s.family<<'\t'<<s.middle<<'\t'<<s.tail<<'\t'<<ds[h.di].id<<'\t'<<h.shift<<'\t'<<ls[h.li].id<<'\t'<<h.weight<<'\t'<<h.sign<<'\t'<<h.mw<<'\t'<<std::hex<<h.phase<<std::dec<<'\n';
                std::cerr<<s.family<<" M"<<s.middle<<' '<<s.tail<<" checked="<<s.evaluated<<" nonzero="<<s.nonzero<<" best="<<(s.hits.empty()?-1:s.hits[0].weight)<<" certified="<<s.certified<<'\n';
            }
            summary.flush();hits.flush();
        }
    }
    std::ofstream timing(prefix+"_timing.txt");timing<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<'\n';
}

#include "datf_middle.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

using namespace gdl;
struct Diff {std::string id; int family,weight,orbit; State input,output; bool middle2;};
struct Group {std::string id; State gamma; uint64_t members;};
State read_state(std::istream &f) {State a{};for(auto &w:a)f>>std::hex>>w;f>>std::dec;return a;}
int main(int argc,char **argv) {
 try {
    if(argc!=7)throw std::runtime_error("usage: datf_linear_two_shard DIFFS GROUPS PREFIX SHARD SHARDS M2_REPS_PER_FAMILY (0=all)");
    const std::string prefix=argv[3];const int shard=std::stoi(argv[4]),shards=std::stoi(argv[5]);
    const size_t m2_limit=std::stoull(argv[6]);
    if(shards<1||shard<0||shard>=shards)throw std::runtime_error("invalid shard");
    if(std::filesystem::exists(prefix+"_summary.json"))throw std::runtime_error("completed shard exists");
    std::ifstream df(argv[1]),gf(argv[2]);if(!df||!gf)throw std::runtime_error("cannot open inputs");
    std::vector<Diff> ds;std::vector<Group> groups;size_t order=0,counts[2]{};
    for(;;) {
        Diff d;std::string family;if(!(df>>d.id>>family>>d.weight>>d.orbit))break;
        if(family!="D1"&&family!="D2")throw std::runtime_error("invalid differential family");
        d.family=family=="D1"?0:1;d.input=read_state(df);d.output=read_state(df);
        if(!df)throw std::runtime_error("invalid differential record");
        d.middle2=!m2_limit||counts[d.family]<m2_limit;++counts[d.family];
        if(order++%shards==size_t(shard))ds.push_back(d);
    }
    for(;;) {Group g;if(!(gf>>g.id>>g.members))break;g.gamma=read_state(gf);
        if(!gf)throw std::runtime_error("invalid group record");
        groups.push_back(g);}
    if(groups.empty())throw std::runtime_error("empty linear groups");
    uint64_t members=0;for(auto &g:groups)members+=g.members;
    uint64_t states[2][2]{},nonzero_groups[2][2]{},nonzero_pairs[2][2]{},representatives[2][2]{};
    std::ofstream hits(prefix+"_hits.tsv.partial");
    hits<<"family\tmiddle\tgroup\tdiff_id\trotation\tsign\tmiddle_weight\trc_phase_mask\n";
    MiddleTables tables;auto began=std::chrono::steady_clock::now();
    auto progress=[&](size_t done) {
        std::ofstream f(prefix+"_progress.json.partial");
        f<<"{\"shard\":"<<shard<<",\"representatives_done\":"<<done<<",\"representatives_total\":"<<ds.size()
         <<",\"middle1_states\":"<<states[0][0]+states[1][0]<<",\"middle2_states\":"<<states[0][1]+states[1][1]
         <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<"}\n";
        f.close();std::filesystem::rename(prefix+"_progress.json.partial",prefix+"_progress.json");hits.flush();
    };
    progress(0);
    for(size_t index=0;index<ds.size();++index) {
        const auto &d=ds[index];
        for(int shift=0;shift<d.orbit;++shift) {
            auto delta=linear(rotate(d.output,shift));
            std::unique_ptr<MiddleTwo> two;if(d.middle2)two=std::make_unique<MiddleTwo>(delta,tables);
            for(const auto &g:groups) {
                CorrelationPower values[2]={middle_one(delta,g.gamma),two?two->evaluate(g.gamma):CorrelationPower{}};
                for(int m=0;m<(d.middle2?2:1);++m) {
                    auto c=values[m];if(!c.sign)continue;
                    if(m==1&&__builtin_parityll(c.rc_mask&RC[9]))c.sign=-c.sign; // two-chi tail: 11-r_l
                    ++nonzero_groups[d.family][m];nonzero_pairs[d.family][m]+=g.members;
                    hits<<'D'<<d.family+1<<'\t'<<m+1<<'\t'<<g.id<<'\t'<<d.id<<'\t'<<shift<<'\t'
                        <<c.sign<<'\t'<<c.weight<<'\t'<<std::hex<<c.rc_mask<<std::dec<<'\n';
                }
            }
            ++states[d.family][0];if(d.middle2)++states[d.family][1];
        }
        ++representatives[d.family][0];if(d.middle2)++representatives[d.family][1];
        if((index+1)%16==0)progress(index+1);
    }
    progress(ds.size());hits.close();if(!hits)throw std::runtime_error("hit write failed");
    std::filesystem::rename(prefix+"_hits.tsv.partial",prefix+"_hits.tsv");
    std::ofstream out(prefix+"_summary.json.partial");
    out<<"{\"shard\":"<<shard<<",\"shards\":"<<shards<<",\"seconds\":"
       <<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"rows\":[";
    bool first=true;
    for(int f=0;f<2;++f)for(int m=0;m<2;++m) {
        if(!first)out<<',';
        first=false;
        out<<"{\"family\":\"D"<<f+1<<"\",\"middle\":"<<m+1<<",\"representatives\":"<<representatives[f][m]
           <<",\"states\":"<<states[f][m]<<",\"group_evaluations\":"<<states[f][m]*groups.size()
           <<",\"logical_pairs\":"<<states[f][m]*members<<",\"nonzero_groups\":"<<nonzero_groups[f][m]
           <<",\"nonzero_pairs\":"<<nonzero_pairs[f][m]<<'}';
    }
    out<<"]}\n";out.close();if(!out)throw std::runtime_error("summary write failed");
    std::filesystem::rename(prefix+"_summary.json.partial",prefix+"_summary.json");
    std::cout<<"DONE "<<shard<<std::endl;
 }catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}

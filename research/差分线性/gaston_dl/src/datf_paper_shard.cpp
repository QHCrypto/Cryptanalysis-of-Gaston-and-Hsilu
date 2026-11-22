#include "datf_middle.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace gdl;
struct Diff {std::string id,family;int weight,orbit;State input,output;};
struct Tail {std::string id;int rounds,weight;State gamma,output;};
State read_state(std::istream& f){State a{};for(auto &w:a)f>>std::hex>>w;f>>std::dec;return a;}
int main(int argc,char**argv){
 try {
  if(argc!=6)throw std::runtime_error("usage: datf_paper_shard DIFFERENTIALS LINEARS PREFIX SHARD SHARDS");
  int shard=std::stoi(argv[4]),shards=std::stoi(argv[5]);
  if(shards<1||shard<0||shard>=shards)throw std::runtime_error("invalid shard");
  std::string prefix=argv[3];
  if(std::filesystem::exists(prefix+"_summary.json"))throw std::runtime_error("completed shard exists");
  std::ifstream df(argv[1]),lf(argv[2]);
  if(!df||!lf)throw std::runtime_error("cannot open inputs");
  std::vector<Diff> ds;std::vector<Tail> tails;size_t global=0;
  for(;;){Diff d;if(!(df>>d.id>>d.family>>d.weight>>d.orbit))break;d.input=read_state(df);d.output=read_state(df);
   if(!df)throw std::runtime_error("invalid differential input");
   if(global++%shards==size_t(shard))ds.push_back(d);
  }
  for(;;){Tail t;std::string family;if(!(lf>>t.id>>family>>t.rounds>>t.weight))break;t.gamma=read_state(lf);t.output=read_state(lf);
   if(!lf)throw std::runtime_error("invalid linear input");
   if(family=="paper")tails.push_back(t);
  }
  if(tails.size()!=11)throw std::runtime_error("expected 11 paper tails");
  std::map<State,std::vector<size_t>> groups;
  for(size_t j=0;j<tails.size();++j)groups[tails[j].gamma].push_back(j);
  uint64_t evaluated[2][2][11]{},nonzero[2][2][11]{},states[2]{},representatives[2]{};
  std::ofstream hits(prefix+"_hits.tsv.partial");
  if(!hits)throw std::runtime_error("cannot open hit output");
  hits<<"family\tmiddle\ttail\tdiff_id\trotation\tlinear_id\tweight\tsign\tmiddle_weight\trc_phase_mask\n";
  MiddleTables tables;auto began=std::chrono::steady_clock::now();
  auto progress=[&](size_t done){
   std::ofstream f(prefix+"_progress.json.partial");
   f<<"{\"shard\":"<<shard<<",\"representatives_done\":"<<done<<",\"representatives_total\":"<<ds.size()
    <<",\"D1_states\":"<<states[0]<<",\"D2_states\":"<<states[1]
    <<",\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<"}\n";
   f.close();std::filesystem::rename(prefix+"_progress.json.partial",prefix+"_progress.json");hits.flush();
  };
  progress(0);
  for(size_t index=0;index<ds.size();++index){
   const auto &d=ds[index];int family=d.family=="D1"?0:1;
   if(d.family!="D1"&&d.family!="D2")throw std::runtime_error("unknown family");
   for(int shift=0;shift<d.orbit;++shift){
    auto delta=linear(rotate(d.output,shift));MiddleTwo two(delta,tables);
    for(auto &group:groups){
     CorrelationPower cs[2]={middle_one(delta,group.first),two.evaluate(group.first)};
     for(int m=0;m<2;++m)for(size_t j:group.second){
      ++evaluated[family][m][j];auto c=cs[m];if(!c.sign)continue;++nonzero[family][m][j];
      if(m==1&&__builtin_parityll(c.rc_mask&RC[11-tails[j].rounds]))c.sign=-c.sign;
      hits<<d.family<<'\t'<<m+1<<'\t'<<tails[j].id<<'\t'<<d.id<<'\t'<<shift<<'\t'<<tails[j].id<<'\t'
       <<d.weight+tails[j].weight+c.weight<<'\t'<<c.sign<<'\t'<<c.weight<<'\t'<<std::hex<<c.rc_mask<<std::dec<<'\n';
     }
    }
    ++states[family];
   }
   ++representatives[family];if((index+1)%64==0)progress(index+1);
  }
  progress(ds.size());hits.close();if(!hits)throw std::runtime_error("hit write failed");
  std::filesystem::rename(prefix+"_hits.tsv.partial",prefix+"_hits.tsv");
  std::ofstream summary(prefix+"_summary.json.partial");
  summary<<"{\"shard\":"<<shard<<",\"shards\":"<<shards<<",\"D1_representatives\":"<<representatives[0]
   <<",\"D2_representatives\":"<<representatives[1]<<",\"D1_states\":"<<states[0]<<",\"D2_states\":"<<states[1]
   <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"rows\":[";
  bool first=true;
  for(int family=0;family<2;++family)for(int m=0;m<2;++m)for(size_t j=0;j<tails.size();++j){
   if(!first)summary<<',';
   first=false;
   summary<<"{\"family\":\"D"<<family+1<<"\",\"middle\":"<<m+1<<",\"tail\":\""<<tails[j].id
    <<"\",\"evaluated\":"<<evaluated[family][m][j]<<",\"nonzero\":"<<nonzero[family][m][j]<<'}';
  }
  summary<<"]}\n";summary.close();if(!summary)throw std::runtime_error("summary write failed");
  std::filesystem::rename(prefix+"_summary.json.partial",prefix+"_summary.json");
  std::cout<<"completed shard "<<shard<<" states "<<states[0]+states[1]<<std::endl;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
}

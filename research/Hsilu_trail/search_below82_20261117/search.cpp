#include "common.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include <chrono>
Tables tab;
int best=82;
uint64_t extensions=0,cores=0,scanned=0;
std::ofstream out,catalog;
std::set<S> seen,catseen;
void save(const S&b0,const S&a1,const S&b1,const S&a2,int total){
 out<<"{\"weight\":"<<total<<",\"weights\":["<<reverse_weight_fast(b0)<<','<<weight(a1)<<','<<weight(a2)<<"],\"b0\":"<<json_state(b0)<<",\"a1\":"<<json_state(a1)<<",\"b1\":"<<json_state(b1)<<",\"a2\":"<<json_state(a2)<<"}\n";out.flush();
 std::cerr<<"candidate "<<total<<" ("<<reverse_weight_fast(b0)<<','<<weight(a1)<<','<<weight(a2)<<") cores="<<cores<<" extensions="<<extensions<<std::endl;
}
int extend(const S&b0,const S&a1,int maxdim=24){
 int base=reverse_weight_fast(b0)+weight(a1);
 if(base+2>best)return 999;
 if(!seen.insert(canonical(a1)).second)return 999;
 auto [off,bs]=tab.affine(a1);int d=bs.size();if(d>maxdim)return 999;
 cores++;S a=linear(off),b=off,besta=a,bestb=b;std::vector<S> images;for(auto s:bs)images.push_back(linear(s));
 int local=weight(a);uint64_t end=1ULL<<d;
 for(uint64_t i=1;i<end;i++){int j=__builtin_ctzll(i);a=a^images[j];b=b^bs[j];int w=weight(a);if(w<local){local=w;besta=a;bestb=b;}}
 extensions+=end;
 if(base+local<=best){best=base+local;save(b0,a1,bestb,besta,best);}
 return base+local;
}
void record(const S&b0,const S&a1){
 int cost=reverse_weight_fast(b0)+weight(a1);
 if(cost<=60){S c=canonical(a1);if(catseen.insert(c).second){for(auto v:c)catalog<<std::hex<<v<<' ';catalog<<'\n';catalog.flush();}}
 extend(b0,a1);
}
int main(int argc,char**argv){
 std::string mode=argc>1?argv[1]:"seeds";
 out.open("results/"+mode+".jsonl");catalog.open("results/"+mode+"_catalog.txt");
 for(int a=0;a<32;a++){S s{};for(int i=0;i<5;i++)s[i]=a>>i&1;assert(reverse_weight_fast(s)==tab.reverse_weight(s));}
 for(int i=0;i<320;i++){S s{};s[i/64]=1ULL<<(i%64);assert(linear(inverse(s))==s);}
 if(mode=="seeds" || mode=="cuts"){
  std::ifstream in(mode=="seeds"?"seeds.txt":"cut_seeds.txt");S a{};while(in>>std::hex>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]){std::cerr<<"seed weights "<<reverse_weight_fast(inverse(a))<<','<<weight(a)<<std::endl;record(inverse(a),a);}
 }else if(mode=="sparse3"){
  S pc[64][32]{},ic[64][32]{};int w[32];
  for(int z=0;z<64;z++)for(int a=1;a<32;a++){for(int i=0;i<5;i++)pc[z][a][i]=U(a>>i&1)<<z;ic[z][a]=inverse(pc[z][a]);w[a]=tab.w[a];}
  for(int a=1;a<32;a++)record(ic[0][a],pc[0][a]);
  for(int z=1;z<64;z++)for(int a=1;a<32;a++)for(int b=1;b<32;b++){
   S y=ic[0][a]^ic[z][b],x=pc[0][a]^pc[z][b];if(reverse_weight_fast(y)+w[a]+w[b]+2<=best)record(y,x);
  }
  for(int z=1;z<63;z++){
   for(int k=z+1;k<64;k++)for(int a=1;a<32;a++)for(int b=1;b<32;b++){
    S y=ic[0][a]^ic[z][b],x=pc[0][a]^pc[z][b];int cost=w[a]+w[b];
    for(int c=1;c<32;c++){scanned++;S yy=y^ic[k][c];if(reverse_weight_fast(yy)+cost+w[c]+2<=best)record(yy,x^pc[k][c]);}
   }
   std::cerr<<"z="<<z<<" scanned="<<scanned<<" cores="<<cores<<" best="<<best<<std::endl;
  }
 }
 std::cout<<"{\"mode\":\""<<mode<<"\",\"scanned\":"<<scanned<<",\"extended_cores\":"<<cores<<",\"extensions\":"<<extensions<<",\"catalog\":"<<catseen.size()<<",\"best\":"<<best<<"}"<<std::endl;
}

#pragma once
// Shared exact saturated-degree MP engine for Gaston.
// The degree 2^r homogeneous component is computed in GF(2)[x]/(x_i^2).
// This graded algebra is NOT the full Boolean ANF algebra. It is exact ONLY
// for the maximal degree 2^r component of this quadratic-round circuit.
#include <algorithm>
#include <array>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <gurobi_c++.h>

namespace gaston_mp {
volatile std::sig_atomic_t stopped=0;
GRBModel* live=nullptr;
void interrupt(int) { stopped=1; if(live)live->terminate(); }
struct Guard { explicit Guard(GRBModel& m){live=&m;} ~Guard(){live=nullptr;} };
enum Kind {ZERO,SOURCE,SUM,PRODUCT};
struct Node {Kind kind; int degree; std::vector<int> in;};
using State=std::array<std::vector<int>,5>;
constexpr int E[5]={0,60,22,27,4}, W[5]={0,56,31,46,43}, T[5]={25,32,52,60,63};
struct Circuit {
    std::vector<Node> nodes{{ZERO,0,{}}};
    std::map<std::pair<Kind,std::vector<int>>,int> memo;
    int source(){int id=nodes.size();nodes.push_back({SOURCE,1,{}});return id;}
    int intern(Kind kind,int degree,std::vector<int> in){
        auto key=std::make_pair(kind,in);auto it=memo.find(key);
        if(it!=memo.end())return it->second;
        int id=nodes.size();nodes.push_back({kind,degree,std::move(in)});memo.emplace(key,id);return id;
    }
    int sum(const std::vector<int>& terms){
        std::vector<int> flat,out;
        for(int t:terms)if(t){
            if(nodes[t].kind==SUM)flat.insert(flat.end(),nodes[t].in.begin(),nodes[t].in.end());
            else flat.push_back(t);
        }
        std::sort(flat.begin(),flat.end());
        for(size_t i=0;i<flat.size();){size_t j=i+1;while(j<flat.size()&&flat[j]==flat[i])++j;
            if((j-i)&1)out.push_back(flat[i]);
            i=j;}
        if(out.empty())return 0;
        if(out.size()==1)return out[0];
        int degree=nodes[out[0]].degree;
        for(int x:out)if(nodes[x].degree!=degree)throw std::runtime_error("inhomogeneous sum");
        return intern(SUM,degree,std::move(out));
    }
    int product(int a,int b){
        if(!a||!b||a==b)return 0; // Every positive-degree element squares to zero.
        if(a>b)std::swap(a,b);
        return intern(PRODUCT,nodes[a].degree+nodes[b].degree,{a,b});
    }
    int quadratic(int a,int b){
        if(!a||!b)return 0;
        auto aa=nodes[a].kind==SUM?nodes[a].in:std::vector<int>{a};
        auto bb=nodes[b].kind==SUM?nodes[b].in:std::vector<int>{b};
        std::vector<int> out;for(int x:aa)for(int y:bb)out.push_back(product(x,y));
        return sum(out);
    }
    State linear(const State& s,int width){
        State out;for(auto& row:out)row.resize(width);
        auto back=[&](int b,int shift){return (b+4*width-shift%width)%width;};
        for(int i=0;i<5;++i)for(int b=0;b<width;++b){
            int q=back(b,W[i]);std::vector<int> terms{s[i][back(q,E[i])]};
            for(int j=0;j<5;++j)for(int shift:{23,24,23+T[j],41+T[j]})
                terms.push_back(s[j][back(q,E[j]+shift)]);
            out[i][b]=sum(terms);
        }
        return out;
    }
};
int number(const std::string& s){size_t n=0;int v=std::stoi(s,&n);if(n!=s.size())throw std::runtime_error("invalid integer");return v;}
std::vector<int> indices(const std::string& s){
    std::vector<int> v;std::istringstream in(s);std::string p;
    while(std::getline(in,p,',')){auto d=p.find('-');if(d==std::string::npos)v.push_back(number(p));
        else {int a=number(p.substr(0,d)),b=number(p.substr(d+1));
            if(a<0||b<a||b>127)throw std::runtime_error("invalid range");
            for(int i=a;i<=b;++i)v.push_back(i);}}
    std::sort(v.begin(),v.end());
    if(std::adjacent_find(v.begin(),v.end())!=v.end())throw std::runtime_error("duplicate nonce index");
    return v;
}
struct Options {int rounds=6,width=64,key=0,output=0,pool=256,threads=4;bool resume=false,expand=true,key_split=true,parity_cube=false;
    std::string mode="screen",result="results/result.json",checkpoint;std::vector<int> nonce;};
Options parse(int argc,char** argv,bool parity_cube){Options o;o.parity_cube=parity_cube;
    for(int i=1;i<argc;++i){std::string k=argv[i];if(k=="--resume"){o.resume=true;continue;}
        if(k=="--help"){
            std::printf("Gaston MP: --rounds 6 --width 64 --key 0|all --nonce 0-62 --output-bit 0\n"
                        "  --mode screen|count|build --pool 256 --threads 4 --expand-first 1 --key-split 1\n"
                        "  key-split=1 keeps only key-degree 0/1 at every round; 0 selects the legacy circuit.\n"
                        "  --result results/result.json --checkpoint results/trails.txt [--resume]\n"
                        "  Coordinate program: --nonce gives original N0/N1 bit indices.\n"
                        "  Parity program: --nonce is forbidden; N0 repeats an even-parity 2^rounds-bit block, N1=0.\n"
                        "  No TimeLimit. A screen witness is NOT a coefficient-one certificate.\n");std::exit(0);}
        if(i+1==argc)throw std::runtime_error("missing option value");
        std::string v=argv[++i];
        if(k=="--rounds")o.rounds=number(v);else if(k=="--width")o.width=number(v);
        else if(k=="--key")o.key=v=="all"?-1:number(v);else if(k=="--nonce"){
            if(o.parity_cube)throw std::runtime_error("parity cube uses directions e0+ej, not coordinate --nonce");
            o.nonce=indices(v);}
        else if(k=="--output-bit")o.output=number(v);else if(k=="--pool")o.pool=number(v);
        else if(k=="--threads")o.threads=number(v);else if(k=="--mode")o.mode=v;
        else if(k=="--result")o.result=v;else if(k=="--checkpoint")o.checkpoint=v;
        else if(k=="--expand-first")o.expand=number(v)!=0;
        else if(k=="--key-split")o.key_split=number(v)!=0;else throw std::runtime_error("unknown option: "+k);
    }
    if(o.rounds<1||o.rounds>6)throw std::runtime_error("invalid rounds");
    if(o.parity_cube){
        int period=1<<o.rounds;
        if(o.width<period||o.width%period)throw std::runtime_error("parity period 2^rounds must divide word width");
        for(int j=1;j<period;++j)o.nonce.push_back(j);
    }else if(o.nonce.empty())for(int j=0;j<(1<<o.rounds)-1;++j)o.nonce.push_back(j);
    if(o.rounds<1||o.rounds>6||o.width<2||o.width>64||o.key < -1||o.key>=2*o.width||
       o.output<0||o.output>=5*o.width||o.nonce.size()!=size_t((1<<o.rounds)-1)||
       o.nonce.front()<0||o.nonce.back()>=2*o.width||o.pool<1||o.threads<1||
       (o.mode!="screen"&&o.mode!="count"&&o.mode!="build"))throw std::runtime_error("invalid options");
    if(o.resume&&(o.checkpoint.empty()||o.mode!="count"))throw std::runtime_error("resume requires count and checkpoint");
    return o;
}
State zero_state(int width){State s;for(auto& row:s)row.assign(width,0);return s;}
int build_target(Circuit& c,const Options& o,std::vector<int>& keys,std::vector<int>& nonce){
    State a=zero_state(o.width),b=zero_state(o.width);
    for(int i=0;i<2*o.width;++i)if(o.key<0||i==o.key){keys[i]=c.source();b[1+i/o.width][i%o.width]=keys[i];}
    if(o.parity_cube){
        for(size_t j=0;j<o.nonce.size();++j)nonce.push_back(c.source());
        int pivot=c.sum(nonce),period=1<<o.rounds;
        for(int j=0;j<o.width;++j)a[3][j]=(j%period)==0?pivot:nonce[j%period-1];
    }else for(int i:o.nonce){int id=c.source();nonce.push_back(id);a[3+i/o.width][i%o.width]=id;}
    if(!o.key_split){
        // Legacy saturated-total-degree circuit: key degree is filtered only in MILP.
        for(int i=0;i<5;++i)for(int j=0;j<o.width;++j)a[i][j]=c.sum({a[i][j],b[i][j]});
        for(int t=0;t<o.rounds;++t){auto u=c.linear(a,o.width);
            for(int i=0;i<5;++i)for(int j=0;j<o.width;++j)
                a[i][j]=t==0&&o.expand?c.quadratic(u[(i+1)%5][j],u[(i+2)%5][j]):c.product(u[(i+1)%5][j],u[(i+2)%5][j]);}
        return a[o.output/o.width][o.output%o.width];
    }
    // Bidegrees at round t: a=(key 0, nonce 2^t), b=(key 1, nonce 2^t-1).
    // In the associated graded algebra, repeated sources vanish. Therefore key
    // degree never decreases along a surviving product and all b*b terms vanish
    // after projection onto key degree <=1. This is exact for total degree 2^r.
    for(int t=0;t<o.rounds;++t){
        auto u=c.linear(a,o.width),v=c.linear(b,o.width);
        auto mul=[&](int x,int y){return t==0&&o.expand?c.quadratic(x,y):c.product(x,y);};
        for(int i=0;i<5;++i)for(int j=0;j<o.width;++j){int p=(i+1)%5,q=(i+2)%5;
            a[i][j]=(size_t(1<<(t+1))>o.nonce.size())?0:mul(u[p][j],u[q][j]);
            b[i][j]=c.sum({mul(u[p][j],v[q][j]),mul(v[p][j],u[q][j])});
        }
    }
    return b[o.output/o.width][o.output%o.width];
}
struct Trail {int key=-1;std::vector<int> selected;};
struct Model {
    const Circuit& c;GRBModel& m;int root;
    std::vector<bool> relevant;std::vector<GRBVar> node,choices;
    std::vector<std::pair<int,int>> choice_edge;std::vector<int> keys,nonce;
    Model(const Circuit& cc,GRBModel& mm,int rr,const std::vector<int>& kk,const std::vector<int>& nn):
        c(cc),m(mm),root(rr),relevant(cc.nodes.size(),false),node(cc.nodes.size()),keys(kk),nonce(nn){
        mark(root);for(int k:keys)mark(k);for(int n:nonce)mark(n);
        std::vector<std::vector<GRBVar>> uses(c.nodes.size());
        for(size_t i=1;i<c.nodes.size();++i)if(relevant[i])node[i]=m.addVar(0,1,0,GRB_BINARY,"n"+std::to_string(i));
        for(size_t i=1;i<c.nodes.size();++i)if(relevant[i]){
            if(c.nodes[i].kind==SUM){GRBLinExpr sum=0;
                for(int child:c.nodes[i].in){auto e=m.addVar(0,1,0,GRB_BINARY,"e"+std::to_string(choices.size()));
                    choices.push_back(e);choice_edge.push_back({int(i),child});uses[child].push_back(e);sum+=e;}
                m.addConstr(sum==node[i]); // Select a summand, not a value-XOR.
            }else if(c.nodes[i].kind==PRODUCT)for(int child:c.nodes[i].in)uses[child].push_back(node[i]);
        }
        for(size_t i=1;i<c.nodes.size();++i)if(relevant[i]){GRBLinExpr demand=i==size_t(root)?1:0;
            for(auto e:uses[i])demand+=e;
            m.addConstr(demand==node[i]);}
        GRBLinExpr count=0;for(int k:keys)if(k)count+=node[k];m.addConstr(count==1);
        for(int n:nonce)m.addConstr(node[n]==1);
        m.setObjective(GRBLinExpr(0),GRB_MINIMIZE);
    }
    void mark(int id){if(!id||relevant[id])return;relevant[id]=true;for(int x:c.nodes[id].in)mark(x);}
    Trail solution(){Trail t;for(size_t i=0;i<keys.size();++i)if(keys[i]&&node[keys[i]].get(GRB_DoubleAttr_Xn)>.5)t.key=i;
        for(size_t i=0;i<choices.size();++i)if(choices[i].get(GRB_DoubleAttr_Xn)>.5)t.selected.push_back(i);
        if(t.key<0)throw std::runtime_error("solution missing key");
        return t;}
    void exclude(const Trail& t){GRBLinExpr sum=0;for(int i:t.selected){
        if(i<0||i>=int(choices.size()))throw std::runtime_error("invalid checkpoint edge");
        sum+=choices[i];}
        // A rooted acyclic SUM/PRODUCT trace is uniquely fixed by its selected
        // SUM edges. Source/node activity then follows from exact use equations.
        m.addConstr(sum<=int(t.selected.size())-1);
    }
    uint64_t hash()const{uint64_t h=1469598103934665603ULL;auto add=[&](uint64_t x){h^=x;h*=1099511628211ULL;};
        add(root);for(const auto& n:c.nodes){add(n.kind);add(n.degree);add(n.in.size());for(int i:n.in)add(i);}
        for(int k:keys)add(k);
        for(int n:nonce)add(n);
        return h;}
};
struct Result {bool complete=false;std::string classification="not_started";uint64_t trails=0,batches=0;
    int status=0,solcount=0,variables=0,constraints=0;size_t circuit=0;double seconds=0;
    std::vector<int> parity;Trail witness;};
void report(const Options& o,const Result& r){
    std::string tmp=o.result+".tmp";FILE* f=std::fopen(tmp.c_str(),"w");if(!f)throw std::runtime_error("cannot write result");
    std::fprintf(f,"{\n  \"schema\": \"GastonCube-v1\", \"rounds\": %d, \"word_bits\": %d,\n"
        "  \"input\": \"IV,K0,K1,N0,N1\",\n"
        "  \"cube_kind\": \"%s\", \"cube_dimension\": %zu, \"cube_period\": %d,\n"
        "  \"output_bit\": %d, \"key_index\": %d,\n",o.rounds,o.width,
        o.parity_cube?"even_parity_periodic":"independent_coordinates",o.nonce.size(),
        o.parity_cube?1<<o.rounds:0,o.output,o.key);
    if(o.parity_cube)std::fprintf(f,"  \"nonce_indices\": null, \"cube_variables\": [");
    else std::fprintf(f,"  \"nonce_indices\": [");
    for(size_t i=0;i<o.nonce.size();++i)std::fprintf(f,"%s%d",i?",":"",o.nonce[i]);
    std::fprintf(f,"],\n  \"key_degree_split\": %s,\n  \"mode\": \"%s\", \"classification\": \"%s\", \"complete\": %s,\n"
        "  \"coefficient\": ",o.key_split?"true":"false",o.mode.c_str(),r.classification.c_str(),r.complete?"true":"false");
    if(r.complete&&o.key>=0)std::fprintf(f,"%d",r.parity[o.key]);else std::fprintf(f,"null");
    std::fprintf(f,",\n  \"coefficients\": ");
    if(r.complete){std::fprintf(f,"{");bool first=true;for(size_t i=0;i<r.parity.size();++i)if(o.key<0||int(i)==o.key){
        std::fprintf(f,"%s\"%zu\":%d",first?"":",",i,r.parity[i]);first=false;}std::fprintf(f,"}");}else std::fprintf(f,"null");
    std::fprintf(f,",\n  \"trails_counted\": %llu, \"batches\": %llu, \"solver_status\": %d,\n"
        "  \"solcount\": %d, \"circuit_nodes\": %zu, \"binary_variables\": %d, \"constraints\": %d,\n"
        "  \"seconds\": %.6f, \"witness_key\": ",static_cast<unsigned long long>(r.trails),
        static_cast<unsigned long long>(r.batches),r.status,r.solcount,r.circuit,r.variables,r.constraints,r.seconds);
    if(r.witness.key>=0)std::fprintf(f,"%d",r.witness.key);else std::fprintf(f,"null");
    std::fprintf(f,",\n  \"scope\": \"A feasible trail is not a nonzero coefficient. Coefficients are certified only when complete=true.\"\n}\n");
    if(std::fclose(f)||std::rename(tmp.c_str(),o.result.c_str()))throw std::runtime_error("result write failed");
}
int run(int argc,char** argv,bool parity_cube=false){
    std::signal(SIGINT,interrupt);std::signal(SIGTERM,interrupt);
    try {Options o=parse(argc,argv,parity_cube);Result r;r.parity.resize(2*o.width);auto start=std::chrono::steady_clock::now();
        Circuit c;std::vector<int> keys(2*o.width),nonce;
        int root=build_target(c,o,keys,nonce);r.circuit=c.nodes.size();
        if(!root){r.complete=true;r.classification="symbolic_zero";report(o,r);std::printf("Exact symbolic coefficient zero.\n");return 0;}
        GRBEnv env;env.set(GRB_IntParam_Threads,o.threads);GRBModel model(env);Guard guard(model);Model mp(c,model,root,keys,nonce);
        model.update(); // Needed here for reporting sizes before optimize/build-only.
        r.variables=model.get(GRB_IntAttr_NumVars);r.constraints=model.get(GRB_IntAttr_NumConstrs);
        std::ostringstream signature;signature<<(o.parity_cube?"GASTONPARITY_V1 ":"GASTONCUBE_V1 ")<<o.rounds<<' '<<o.width<<' '<<o.key<<' '<<o.output<<' '<<mp.hash();
        std::ofstream checkpoint;
        if(o.resume){std::ifstream in(o.checkpoint);std::string line;std::getline(in,line);
            if(line!=signature.str())throw std::runtime_error("checkpoint signature mismatch");
            while(std::getline(in,line)){std::istringstream row(line);Trail t;size_t n;
                if(!(row>>t.key>>n)||t.key<0||t.key>=2*o.width){
                    throw std::runtime_error("truncated checkpoint");
                }
                for(size_t i=0;i<n;++i){int e;if(!(row>>e))throw std::runtime_error("truncated checkpoint");t.selected.push_back(e);}
                std::string extra;if(row>>extra)throw std::runtime_error("bad checkpoint record");mp.exclude(t);r.parity[t.key]^=1;++r.trails;}
        }
        if(!o.checkpoint.empty()){
            if(!o.resume&&std::ifstream(o.checkpoint).good())throw std::runtime_error("checkpoint exists; use --resume");
            checkpoint.open(o.checkpoint,o.resume?std::ios::app:std::ios::out);if(!checkpoint)throw std::runtime_error("cannot open checkpoint");
            if(!o.resume)checkpoint<<signature.str()<<'\n';
            checkpoint.flush();
        }
        std::printf("Target: rounds=%d, degree=%d, key=%d, nonce_degree=%zu, output=%d, key_degree_split=%d\n",o.rounds,1<<o.rounds,o.key,o.nonce.size(),o.output,int(o.key_split));
        std::printf("Circuit=%zu nodes; MILP=%d binary variables, %d constraints.\n",r.circuit,r.variables,r.constraints);std::fflush(stdout);
        r.classification=o.mode=="build"?"model_built":"running";report(o,r);if(o.mode=="build")return 0;
        if(o.mode=="count"){model.set(GRB_IntParam_PoolSearchMode,2);model.set(GRB_IntParam_PoolSolutions,o.pool);}
        while(!stopped){model.optimize();r.status=model.get(GRB_IntAttr_Status);r.solcount=model.get(GRB_IntAttr_SolCount);++r.batches;
            std::vector<Trail> batch;
            if(r.solcount>0){std::printf("Status=%d Objective=%.0f Runtime=%.6f s SolCount=%d\n",r.status,model.get(GRB_DoubleAttr_ObjVal),model.get(GRB_DoubleAttr_Runtime),r.solcount);
                model.set(GRB_IntParam_SolutionNumber,0);r.witness=mp.solution();
                std::printf("Incumbent trail: key=%d, selected_summands=%zu; ANF coefficient remains unknown until complete.\n",r.witness.key,r.witness.selected.size());
                if(o.mode=="count")for(int i=0;i<r.solcount;++i){model.set(GRB_IntParam_SolutionNumber,i);batch.push_back(mp.solution());}
            }
            if(o.mode=="count"){
                for(const auto& t:batch){if(r.trails==std::numeric_limits<uint64_t>::max())throw std::runtime_error("counter overflow");
                    if(checkpoint.is_open()){checkpoint<<t.key<<' '<<t.selected.size();for(int e:t.selected)checkpoint<<' '<<e;checkpoint<<'\n';}
                    r.parity[t.key]^=1;++r.trails;}
                if(checkpoint.is_open()){checkpoint.flush();if(!checkpoint)throw std::runtime_error("checkpoint write failed");}
                if(r.status==GRB_INFEASIBLE||(r.status==GRB_OPTIMAL&&r.solcount<o.pool)){r.complete=true;r.classification="exact_parity";}
                else r.classification=stopped?"interrupted_partial_count":"partial_count";
            }else if(r.status==GRB_INFEASIBLE){r.complete=true;r.classification="infeasible_zero";}
            else r.classification=r.solcount?"feasible_only":"unresolved";
            r.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();report(o,r);
            if(o.mode!="count"||r.complete||stopped||r.status!=GRB_OPTIMAL)break;
            for(const auto& t:batch)mp.exclude(t);
        }
        if(r.complete)std::printf("Certified complete parity: counted=%llu; target coefficient=%s\n",static_cast<unsigned long long>(r.trails),o.key<0?"see coefficients map":r.parity[o.key]?"1":"0");
        return stopped?130:(r.complete||(o.mode=="screen"&&r.solcount>0)?0:3);
    }catch(const GRBException& e){std::fprintf(stderr,"Gurobi error %d: %s\n",e.getErrorCode(),e.getMessage().c_str());return 1;
    }catch(const std::exception& e){std::fprintf(stderr,"Error: %s\n",e.what());return 1;}
}

} // namespace gaston_mp

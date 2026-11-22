// Safe quartic screening by 3BDPwoU-style monomial-trail feasibility.
// INFEASIBLE certifies a zero full-cube coefficient. Feasibility does not
// certify a nonzero coefficient: different trails may cancel modulo two.
// At nonce degree 4*2^(rounds-1), lower-degree branches cannot contribute.
// Reparameterize the first-round coefficient by lambda_j=1+K1_j, so its
// variable-bearing state is (0,lambda*x,0,x,0), before the linear layer.
#include <algorithm>
#include <array>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <gurobi_c++.h>

namespace {
GRBModel* active_model=nullptr;
volatile std::sig_atomic_t interrupted=0;
void handle_interrupt(int) { interrupted=1; if(active_model) active_model->terminate(); }
struct Guard { explicit Guard(GRBModel& m) {active_model=&m;} ~Guard(){active_model=nullptr;} };
using Word=std::vector<int>;
using State=std::array<Word,5>;
constexpr int N1[13]={0,3,1,0,5,5,7,11,12,12,13,14,15};
constexpr int N2[13]={4,4,2,1,2,3,6,8,9,6,10,11,9};
constexpr int M[13]={47,51,54,30,4,55,0,14,59,4,10,45,53};
constexpr int R[5]={0,57,13,31,20};

struct Circuit {
    GRBModel& model;
    int n,zero;
    std::vector<GRBVar> node;
    std::vector<std::vector<GRBVar>> uses;
    Circuit(GRBModel& m,int width):model(m),n(width) {
        zero=fresh(); model.addConstr(node[zero]==0);
    }
    int fresh() {
        int id=static_cast<int>(node.size());
        node.push_back(model.addVar(0,1,0,GRB_BINARY)); uses.emplace_back(); return id;
    }
    int multiply(int a,int b) {
        if(a==zero||b==zero)return zero;
        if(a==b)return a;
        int z=fresh(); uses[a].push_back(node[z]); uses[b].push_back(node[z]);
        return z;
    }
    int add(int a,int b) {
        if(a==zero)return b;
        if(b==zero)return a;
        if(a==b)return zero;
        int z=fresh();
        auto u=model.addVar(0,1,0,GRB_BINARY),v=model.addVar(0,1,0,GRB_BINARY);
        uses[a].push_back(u);uses[b].push_back(v);
        // XOR on MONOMIAL EXPONENTS selects one input; this is not bit parity.
        model.addConstr(node[z]==u+v); return z;
    }
    State linear(const State& state) {
        std::vector<Word> x(state.begin(),state.end());
        for(int t=0;t<13;++t) {
            Word w(n);
            for(int b=0;b<n;++b)w[b]=add(x[N1[t]][b],x[N2[t]][(b-M[t]%n+n)%n]);
            x.push_back(w);
        }
        State out;
        for(int i=0;i<5;++i) {
            out[i].resize(n);
            for(int b=0;b<n;++b) {
                // Distinct boundary copies are essential: an internal linear
                // node can also feed another output of this same layer.
                int source=x[13+i][(b-R[i]%n+n)%n];
                out[i][b]=fresh();uses[source].push_back(node[out[i][b]]);
            }
        }
        return out;
    }
    void weight_cut(const State& state,int weight) {
        GRBLinExpr sum=0;
        for(const auto& w:state)for(int id:w)sum+=node[id];
        model.addConstr(sum==weight);
    }
    void connect() {
        for(size_t i=0;i<node.size();++i) {
            GRBLinExpr sum=0;
            for(auto edge:uses[i]) {model.addConstr(edge<=node[i]);sum+=edge;}
            model.addConstr(sum>=node[i]); // COPY is OR, as in Appendix A.
            // At the saturated nonce degree, a reused positive-degree node
            // would replace f*f by f and lose degree. Thus every useful trail
            // is disjoint; lambda sources have only one outgoing use anyway.
            model.addConstr(sum<=node[i]);
        }
    }
};
}

int main(int argc,char** argv) {
    std::signal(SIGINT,handle_interrupt);
    if(argc<3) {std::fprintf(stderr,"usage: %s terms.txt results.tsv [start [count]]\n",argv[0]);return 2;}
    try {
        std::ifstream input(argv[1]);
        int n,rounds,count;
        if(!(input>>n>>rounds>>count)||n<4||n>64||rounds<1||rounds>5||n!=4*(1<<(rounds-1)))
            throw std::runtime_error("invalid term-file header");
        std::vector<std::array<int,4>> terms(count);
        for(auto& t:terms) {
            for(int& p:t)if(!(input>>p)||p<0||p>=5*n)throw std::runtime_error("invalid term");
            if(!std::is_sorted(t.begin(),t.end())||std::adjacent_find(t.begin(),t.end())!=t.end())
                throw std::runtime_error("quartic must have four distinct sorted indices");
        }
        int start=argc>3?std::atoi(argv[3]):0;
        int end=argc>4?std::min(count,start+std::atoi(argv[4])):count;
        if(start<0||start>=end)throw std::runtime_error("invalid range");
        GRBEnv env;
        env.set(GRB_IntParam_OutputFlag,1);env.set(GRB_IntParam_LogToConsole,1);
        env.set(GRB_IntParam_Threads,4);
        GRBModel model(env);Guard guard(model);Circuit c(model,n);
        State initial;
        for(auto& w:initial)w.assign(n,c.zero);
        std::vector<int> nonce(n),lambda(n);
        for(int b=0;b<n;++b) {
            nonce[b]=c.fresh();lambda[b]=c.fresh();
            model.addConstr(c.node[nonce[b]]==1);
            initial[1][b]=c.multiply(nonce[b],lambda[b]);initial[3][b]=nonce[b];
        }
        State state=c.linear(initial);
        // These valid maximum-degree cuts also prevent loss from shared
        // intermediate factors: to reach degree n all 4*2^(r-1) leaves differ.
        c.weight_cut(state,n);
        for(int r=1;r<rounds;++r) {
            State next;
            for(int i=0;i<5;++i) {
                next[i].resize(n);
                for(int b=0;b<n;++b)next[i][b]=c.multiply(state[(i+1)%5][b],state[(i+2)%5][b]);
            }
            state=c.linear(next);c.weight_cut(state,n>>r);
        }
        std::vector<GRBVar> output(5*n);
        for(int i=0;i<5;++i)for(int b=0;b<n;++b) {
            output[i*n+b]=model.addVar(0,0,0,GRB_BINARY);
            c.uses[state[i][b]].push_back(output[i*n+b]);
        }
        c.connect();model.setObjective(GRBLinExpr(0),GRB_MINIMIZE);
        bool write_header=start==0||!std::ifstream(argv[2]).good();
        FILE* report=std::fopen(argv[2],start==0?"w":"a");
        if(!report)throw std::runtime_error("cannot open result file");
        if(write_header)std::fprintf(report,"id\tp\tq\ts\tt\tstatus\tsolcount\tseconds\tclassification\tlambda_support\n");
        std::array<int,4> previous={-1,-1,-1,-1};
        for(int index=start;index<end&&!interrupted;++index) {
            for(int p:previous)if(p>=0) {output[p].set(GRB_DoubleAttr_LB,0);output[p].set(GRB_DoubleAttr_UB,0);}
            for(int p:terms[index]) {output[p].set(GRB_DoubleAttr_LB,1);output[p].set(GRB_DoubleAttr_UB,1);}
            std::printf("SCREEN %d/%d: %d %d %d %d\n",index+1,count,terms[index][0],terms[index][1],terms[index][2],terms[index][3]);
            std::fflush(stdout);
            model.optimize();
            int status=model.get(GRB_IntAttr_Status),solutions=model.get(GRB_IntAttr_SolCount);
            double seconds=model.get(GRB_DoubleAttr_Runtime);
            const char* label=status==GRB_INFEASIBLE?"excluded":solutions>0?"retained":"unresolved";
            std::string support;
            if(solutions>0) {
                std::printf("Status=%d Objective=%.0f Runtime=%.6f s; retained trail, not a nonzero coefficient.\n",status,model.get(GRB_DoubleAttr_ObjVal),seconds);
                std::printf("lambda exponent support:");
                for(int b=0;b<n;++b)if(c.node[lambda[b]].get(GRB_DoubleAttr_X)>.5) {
                    std::printf(" %d",b);support+=(support.empty()?"":",")+std::to_string(b);
                }
                std::printf("\n");
            }
            std::fprintf(report,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%s\t%s\n",index,terms[index][0],terms[index][1],terms[index][2],terms[index][3],status,solutions,seconds,label,support.c_str());
            std::fflush(report);previous=terms[index];
        }
        std::fclose(report);
        return interrupted?130:0;
    } catch(const GRBException& e) {
        std::fprintf(stderr,"Gurobi error %d: %s\n",e.getErrorCode(),e.getMessage().c_str());return 1;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Error: %s\n",e.what());return 1;
    }
}

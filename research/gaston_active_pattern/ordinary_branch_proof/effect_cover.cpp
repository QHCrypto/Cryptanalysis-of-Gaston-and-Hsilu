// Necessary activity-only conditions; infeasibility implies infeasibility
// of the full bitwise propagation, while feasibility alone proves nothing.
#include <gurobi_c++.h>
#include <array>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <string>
namespace {
GRBModel* active=nullptr;
void stop(int){if(active)active->terminate();}
struct Guard{explicit Guard(GRBModel& m){active=&m;}~Guard(){active=nullptr;}};
}
int main(int argc,char**argv) {
    if(argc!=4){std::fprintf(stderr,"Usage: %s E_WEIGHT_MIN TOTAL_MAX PREFIX\n",argv[0]);return 1;}
    std::signal(SIGINT,stop);
    try{
        const int ew=std::stoi(argv[1]),bound=std::stoi(argv[2]);
        const std::string prefix=argv[3];
        GRBEnv env;GRBModel model(env);Guard guard(model);
        model.set(GRB_IntParam_Threads,2);
        model.set(GRB_DoubleParam_MIPGap,0);
        model.set(GRB_DoubleParam_MIPGapAbs,0);
        std::array<GRBVar,64>a,b,e;
        GRBLinExpr n=0,ne=0;
        for(int j=0;j<64;++j) {
            a[j]=model.addVar(0,1,0,GRB_BINARY,"a_"+std::to_string(j));
            b[j]=model.addVar(0,1,0,GRB_BINARY,"b_"+std::to_string(j));
            e[j]=model.addVar(0,1,0,GRB_BINARY,"e_"+std::to_string(j));
            n+=a[j]+b[j];ne+=e[j];
        }
        constexpr int east[]={0,60,22,27,4},west[]={0,56,31,46,43};
        for(int j=0;j<64;++j)for(int i=0;i<5;++i)
            model.addConstr(e[j]<=a[(j+23-east[i]+64)%64]+b[(j+23+west[i])%64]);
        model.addConstr(ne>=ew);model.addConstr(n<=bound);model.addConstr(e[0]==1);
        model.setObjective(GRBLinExpr(0),GRB_MINIMIZE);
        model.write(prefix+".lp");model.optimize();
        std::ofstream f(prefix+".json");
        f<<"{\"status\":"<<model.get(GRB_IntAttr_Status)<<",\"effect_weight_min\":"<<ew
         <<",\"total_max\":"<<bound<<",\"runtime_seconds\":"<<model.get(GRB_DoubleAttr_Runtime)<<"}\n";
        if(model.get(GRB_IntAttr_SolCount)>0) {
            std::printf("Status %d; objective %.0f; runtime %.6f s\n",model.get(GRB_IntAttr_Status),
                        model.get(GRB_DoubleAttr_ObjVal),model.get(GRB_DoubleAttr_Runtime));
            model.write(prefix+".sol");
        }
    }catch(const GRBException&e){std::fprintf(stderr,"%d: %s\n",e.getErrorCode(),e.getMessage().c_str());return 1;}
    catch(const std::exception&e){std::fprintf(stderr,"%s\n",e.what());return 1;}
}

// Exact total-activity decision instance. XOR hulls are read from base.lp.
#include <gurobi_c++.h>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <string>
#include <stdexcept>
namespace {
GRBModel* active=nullptr;
void stop(int) {if(active)active->terminate();}
struct Guard {explicit Guard(GRBModel& m){active=&m;} ~Guard(){active=nullptr;}};
}
int main(int argc,char**argv) {
    if(argc!=7) {
        std::fprintf(stderr,"Usage: %s BASE.lp TOTAL_MAX OUTPUT_EXACT_OR_0 EFFECT_HEX_OR_any THREADS PREFIX\n",argv[0]);
        return 1;
    }
    std::signal(SIGINT,stop);
    try {
        const int total=std::stoi(argv[2]),out=std::stoi(argv[3]);
        const std::string effect=argv[4],prefix=argv[6];
        GRBEnv env;
        GRBModel model(env,argv[1]);
        Guard guard(model);
        model.set(GRB_IntParam_Threads,std::stoi(argv[5]));
        model.set(GRB_IntParam_MIPFocus,3);
        model.set(GRB_DoubleParam_MIPGap,0);
        model.set(GRB_DoubleParam_MIPGapAbs,0);
        model.set(GRB_DoubleParam_FeasibilityTol,1e-9);
        model.set(GRB_DoubleParam_IntFeasTol,1e-9);
        model.remove(model.getConstrByName("fixed_count"));
        GRBLinExpr ni=0,no=0,ne=0;
        for(int j=0;j<64;++j) {
            ni+=model.getVarByName("active_input_"+std::to_string(j));
            no+=model.getVarByName("active_output_"+std::to_string(j));
            ne+=model.getVarByName("effect_"+std::to_string(j));
        }
        model.addConstr(ni+no<=total,"total_activity_bound");
        if(out)model.addConstr(no==out,"fixed_output");
        model.addConstr(ne<=ni+no,"effect_weight_upper");
        model.addConstr(ne>=2,"nonzero_even_effect");
        if(effect=="any")model.addConstr(model.getVarByName("effect_0")==1,"rotation_representative");
        else {
            const auto e=std::stoull(effect,nullptr,16);
            if(!e)throw std::invalid_argument("E must be nonzero");
            for(int j=0;j<64;++j)model.addConstr(model.getVarByName("effect_"+std::to_string(j))==double(e>>j&1));
        }
        model.setObjective(GRBLinExpr(0),GRB_MINIMIZE);
        model.write(prefix+".lp");
        model.optimize();
        const int status=model.get(GRB_IntAttr_Status),solutions=model.get(GRB_IntAttr_SolCount);
        std::ofstream f(prefix+".status.json");
        f.precision(17);
        f<<"{\"total_max\":"<<total<<",\"output_exact\":"<<out<<",\"effect_search\":\""<<effect
         <<"\",\"status\":"<<status<<",\"solutions\":"<<solutions<<",\"infeasible_proven\":"
         <<(status==GRB_INFEASIBLE?"true":"false")<<",\"runtime_seconds\":"
         <<model.get(GRB_DoubleAttr_Runtime)<<"}\n";
        if(solutions) {
            std::printf("Status: %d; objective: %.0f; runtime: %.6f s\n",status,
                        model.get(GRB_DoubleAttr_ObjVal),model.get(GRB_DoubleAttr_Runtime));
            for(const char* side:{"input","output"}) {
                std::printf("%s:",side);
                for(int i=0;i<5;++i) {
                    unsigned long long w=0;
                    for(int j=0;j<64;++j)
                        if(model.getVarByName(std::string(side)+"_"+std::to_string(i)+"_"+std::to_string(j)).get(GRB_DoubleAttr_X)>.5)
                            w|=1ULL<<j;
                    std::printf(" 0x%016llx",w);
                }
                std::printf("\n");
            }
            model.write(prefix+".sol");
        }
        if(!f)throw std::runtime_error("Cannot write status");
        return status==GRB_INFEASIBLE||solutions?0:2;
    } catch(const GRBException& e) {
        std::fprintf(stderr,"Gurobi error %d: %s\n",e.getErrorCode(),e.getMessage().c_str());return 1;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Error: %s\n",e.what());return 1;
    }
}

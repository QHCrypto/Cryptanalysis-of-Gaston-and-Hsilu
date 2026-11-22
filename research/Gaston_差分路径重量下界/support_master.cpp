// Stream distinct feasible middle supports of the existing lossy three-round model.
#include <gurobi_c++.h>
#include <array>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
GRBModel* active = nullptr;
volatile std::sig_atomic_t stopped = 0;
void interrupt(int) {
    stopped = 1;
    if (active) active->terminate();
}
struct Guard {
    explicit Guard(GRBModel& model) { active = &model; }
    ~Guard() { active = nullptr; }
};
}

int main(int argc, char** argv) {
    // The caller exports the model with its weight cap and without rotation constraints.
    if (argc != 4 && argc != 5) {
        std::fprintf(stderr, "usage: support_master model.mps threads proposals.jsonl [first_effect_weight]\n");
        return 2;
    }
    try {
        std::signal(SIGINT, interrupt);
        std::signal(SIGTERM, interrupt);
        GRBEnv env;
        GRBModel model(env, argv[1]);
        Guard guard(model);
        std::ofstream proposals(argv[3]);
        proposals.exceptions(std::ios::failbit | std::ios::badbit);
        model.set(GRB_IntParam_Threads, std::stoi(argv[2]));
        model.set(GRB_IntParam_MIPFocus, 1);
        model.set(GRB_IntParam_SolutionLimit, 1);
        std::array<GRBVar,64> support;
        const int first_effect = argc==5 ? std::stoi(argv[4]) : 0;
        if (first_effect!=0 && first_effect!=2) throw std::invalid_argument("first effect must be 0 or 2");
        GRBLinExpr effect1 = 0, effect2 = 0;
        for (int z=0; z<64; ++z) {
            auto id = std::to_string(z);
            auto bit = model.getVarByName("L0_e_"+id);
            if (!first_effect) bit.set(GRB_DoubleAttr_UB, 0.0);
            else effect1 += bit;
            effect2 += model.getVarByName("L1_e_"+id);
            support[z] = model.getVarByName("chi1_active_"+id);
        }
        model.addConstr(effect2 <= 4, "second_effect_cap");
        if (first_effect) model.addConstr(effect1 == first_effect, "first_effect_weight");
        while (!stopped) {
            model.optimize();
            int status = model.get(GRB_IntAttr_Status);
            if (!model.get(GRB_IntAttr_SolCount)) {
                proposals << "{\"master_status\":" << status
                          << ",\"exhausted\":" << (status==GRB_INFEASIBLE ? "true" : "false") << "}" << std::endl;
                break;
            }
            uint64_t mask = 0;
            GRBLinExpr different = 0;
            for (int z=0; z<64; ++z) {
                if (support[z].get(GRB_DoubleAttr_X) > 0.5) {
                    mask |= uint64_t(1)<<z;
                    different += 1-support[z];
                } else different += support[z];
            }
            char hex[17];
            std::snprintf(hex, sizeof(hex), "%016llx", static_cast<unsigned long long>(mask));
            double weight = model.get(GRB_DoubleAttr_ObjVal);
            std::printf("status=%d relaxed_weight=%.0f seconds=%.3f support=%s\n", status, weight,
                        model.get(GRB_DoubleAttr_Runtime), hex);
            proposals << "{\"mask\":\"" << hex << "\",\"relaxed_weight\":" << weight << "}" << std::endl;
            std::string reply;
            if (stopped || !std::getline(std::cin,reply)) break;
            if (reply != "next") throw std::runtime_error("invalid support acknowledgement");
            // Enumerate each support once. Unresolved SAT checks are owned by the caller.
            model.addConstr(different >= 1);
        }
    } catch (const GRBException& e) {
        std::fprintf(stderr,"Gurobi error %d: %s\n",e.getErrorCode(),e.getMessage().c_str());
        return 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"%s\n",e.what());
        return 1;
    }
}

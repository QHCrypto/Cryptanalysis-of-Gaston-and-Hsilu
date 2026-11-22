// Exact decision problem built from the existing, exported Gaston MILP.
// Optional fixed E only restricts the counterexample search: its infeasibility
// is NOT a proof for unrestricted nonzero E.
#include <gurobi_c++.h>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <string>

namespace {
GRBModel* active = nullptr;
void interrupt(int) { if (active) active->terminate(); }
struct Guard {
    explicit Guard(GRBModel& m) { active = &m; }
    ~Guard() { active = nullptr; }
};
}

int main(int argc, char** argv) {
    if (argc != 7) {
        std::fprintf(stderr, "Usage: %s BASE.lp OUTPUT INPUT_MAX EFFECT_HEX_OR_any THREADS PREFIX\n", argv[0]);
        return 1;
    }
    std::signal(SIGINT, interrupt);
    try {
        const int s = std::stoi(argv[2]), p = std::stoi(argv[3]);
        const std::string effect = argv[4], prefix = argv[6];
        if (s < 1 || s > 64 || p < 1 || p > 64) return 1;
        GRBEnv env;
        GRBModel model(env, argv[1]);
        Guard guard(model);
        model.set(GRB_IntParam_Threads, std::stoi(argv[5]));
        model.set(GRB_IntParam_MIPFocus, 1);
        model.set(GRB_IntParam_SolutionLimit, 1);
        model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
        model.set(GRB_DoubleParam_IntFeasTol, 1e-9);
        model.getConstrByName("fixed_count").set(GRB_DoubleAttr_RHS, s);
        GRBLinExpr nin = 0, ne = 0;
        for (int j = 0; j < 64; ++j) {
            nin += model.getVarByName("active_input_" + std::to_string(j));
            ne += model.getVarByName("effect_" + std::to_string(j));
        }
        model.addConstr(nin <= p, "counterexample_input_max");
        // For every row, wt(E) <= wt(X_i)+wt(Y_i) <= p+s.
        model.addConstr(ne <= p+s, "effect_weight_valid_upper");
        model.addConstr(ne >= 2, "effect_nonzero_even_weight");
        if (effect == "any") {
            // Every nonzero E admits a simultaneous row rotation with E[0]=1.
            // BASE must not contain an additional input/output anchor.
            model.addConstr(model.getVarByName("effect_0") == 1, "effect_rotation_anchor");
        } else {
            const auto e = std::stoull(effect, nullptr, 16);
            if (!e) return 1;
            for (int j = 0; j < 64; ++j)
                model.addConstr(model.getVarByName("effect_" + std::to_string(j)) ==
                                double((e >> j) & 1), "fixed_effect_" + std::to_string(j));
        }
        model.setObjective(GRBLinExpr(0), GRB_MINIMIZE);
        model.write(prefix + ".lp");
        model.optimize();
        const int status = model.get(GRB_IntAttr_Status);
        std::ofstream meta(prefix + ".status.json");
        meta << "{\"status\":" << status << ",\"solutions\":" << model.get(GRB_IntAttr_SolCount)
             << ",\"output_active\":" << s << ",\"input_max\":" << p
             << ",\"effect_search\":\"" << effect << "\",\"runtime_seconds\":"
             << model.get(GRB_DoubleAttr_Runtime) << "}\n";
        if (model.get(GRB_IntAttr_SolCount) > 0) {
            std::printf("Status: %d; feasibility objective: %.0f; runtime: %.6f s\n",
                        status, model.get(GRB_DoubleAttr_ObjVal), model.get(GRB_DoubleAttr_Runtime));
            for (const char* side : {"input", "output"}) {
                std::printf("%s:", side);
                for (int i = 0; i < 5; ++i) {
                    unsigned long long w = 0;
                    for (int j = 0; j < 64; ++j)
                        if (model.getVarByName(std::string(side) + "_" + std::to_string(i) +
                                              "_" + std::to_string(j)).get(GRB_DoubleAttr_X) > .5)
                            w |= 1ULL << j;
                    std::printf(" 0x%016llx", w);
                }
                std::printf("\n");
            }
            model.write(prefix + ".sol");
        }
    } catch (const GRBException& e) {
        std::fprintf(stderr, "Gurobi error %d: %s\n", e.getErrorCode(), e.getMessage().c_str());
        return 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }
}

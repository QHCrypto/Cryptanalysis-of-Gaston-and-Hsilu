// Fu et al. (2024), Algorithm 2 的原始消息比特猜测口径。
#include "model_common.hpp"

int main(int argc, char **argv) {
    using namespace gaston;
    std::signal(SIGINT, handle_sigint);
    try {
        const auto o = options(argc, argv, "original");
        GRBEnv env;
        GRBModel model(env);
        ActiveModelGuard active_guard(model);
        model.set(GRB_StringAttr_ModelName, "Gaston_original_message_bits");
        const auto guess = variables(model, "m");
        const auto selected = variables(model, "l");
        std::array<std::array<GRBVar, N>, 5> known;
        for (int row = 0; row < 5; ++row)
            for (int pos = 0; pos < N; ++pos)
                known[row][pos] = model.addVar(0, 1, 0, GRB_BINARY,
                                                "k" + std::to_string(row) + "_" +
                                                    name("", pos).substr(1));

        GRBLinExpr objective = 0, budget = 0;
        for (int i = 0; i < N; ++i) {
            objective += selected[i];
            budget += guess[i] + selected[i];
        }
        model.setObjective(objective, GRB_MAXIMIZE);
        model.addConstr(budget <= o.budget, "budget");

        // Y^1_0[s] = m[s] xor u[s] xor constant;
        // Y^1_i[s] = u[s+shift_i] xor constant (i=1..4).
        // k=1 当且仅当对应仿射形式的全部原始消息比特都被猜测。
        for (int row = 0; row < 5; ++row)
            for (int pos = 0; pos < N; ++pos) {
                const auto s = support(row, pos);
                GRBLinExpr sum = 0;
                for (int t : s) {
                    model.addConstr(known[row][pos] <= guess[t]);
                    sum += guess[t];
                }
                model.addConstr(known[row][pos] >= sum - (int(s.size()) - 1));
            }

        // 式(15): 每个 Y^2_0[j] 包含 21 个第一轮 chi 输出。
        // 每个乘积至少固定一个因子，作为逐项仿射化的充分条件。
        for (int j = 0; j < N; ++j)
            for (int row = 0; row < 5; ++row) {
                const auto factors = product_rows(row);
                for (int q = 0; q < D_size[row]; ++q) {
                    const int pos = mod64(j + D[row][q]);
                    model.addConstr(selected[j] <= known[factors[0]][pos] +
                                                       known[factors[1]][pos]);
                }
            }

        finish(model, o, guess, selected, "猜测的 m");
    } catch (const GRBException &e) {
        std::fprintf(stderr, "Gurobi 错误 %d: %s\n", e.getErrorCode(),
                     e.getMessage().c_str());
        return 1;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "错误: %s\n", e.what());
        return 1;
    }
}

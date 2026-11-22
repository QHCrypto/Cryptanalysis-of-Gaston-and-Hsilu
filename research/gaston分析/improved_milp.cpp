// 分享对话提出的猜测共享线性字 u=Pm 的覆盖模型。
#include "model_common.hpp"
#include <set>

int main(int argc, char **argv) {
    using namespace gaston;
    std::signal(SIGINT, handle_sigint);
    try {
        const auto o = options(argc, argv, "improved");
        GRBEnv env;
        GRBModel model(env);
        ActiveModelGuard active_guard(model);
        model.set(GRB_StringAttr_ModelName, "Gaston_improved_shared_u");
        auto guess = variables(model, "u");
        auto selected = variables(model, "l");
        GRBLinExpr objective = 0, budget = 0;
        for (int i = 0; i < N; ++i) {
            objective += selected[i];
            budget += guess[i] + selected[i];
        }
        model.setObjective(objective, GRB_MAXIMIZE);
        model.addConstr(budget <= o.budget, "budget");

        // 对 chi 第 0..2 行，两个乘积因子都只依赖一个 u 坐标；
        // 第 3..4 行的另一个因子是 Y^1_0，仅猜 u 无法将其固定。
        const std::array<std::array<int, 2>, 5> choices = {{{56, 31},
                                                              {31, 46},
                                                              {46, 43},
                                                              {43, -1},
                                                              {56, -1}}};
        for (int j = 0; j < N; ++j)
            for (int row = 0; row < 5; ++row)
                for (int q = 0; q < D_size[row]; ++q) {
                    const int pos = mod64(j + D[row][q]);
                    GRBLinExpr coverage = guess[mod64(pos + choices[row][0])];
                    if (choices[row][1] >= 0)
                        coverage += guess[mod64(pos + choices[row][1])];
                    model.addConstr(selected[j] <= coverage);
                }

        // 47/17 是分享对话的可行性见证，用作求解初值。
        const std::set<int> unguessed = {0, 4, 5, 10, 11, 22, 23, 27, 33,
                                         34, 40, 41, 45, 51, 56, 57, 63};
        const std::set<int> targets = {0, 1, 6, 17, 18, 19, 24, 29, 30,
                                       41, 42, 47, 52, 53, 59, 60, 63};
        if (o.budget >= 64) {
            for (int i = 0; i < N; ++i)
                guess[i].set(GRB_DoubleAttr_Start, unguessed.count(i) ? 0 : 1);
            for (int j = 0; j < N; ++j)
                selected[j].set(GRB_DoubleAttr_Start, targets.count(j) ? 1 : 0);
        }
        finish(model, o, guess, selected, "猜测的 u");
    } catch (const GRBException &e) {
        std::fprintf(stderr, "Gurobi 错误 %d: %s\n", e.getErrorCode(),
                     e.getMessage().c_str());
        return 1;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "错误: %s\n", e.what());
        return 1;
    }
}

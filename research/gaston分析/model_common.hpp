#pragma once
#include <gurobi_c++.h>
#include <array>
#include <csignal>
#include <cmath>
#include <cstdio>
#include <exception>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace gaston {
constexpr int N = 64;
constexpr std::array<std::array<int, 5>, 5> D = {{{0, 23, 24, 48, 2},
                                                  {19, 20, 51, 5, -1},
                                                  {45, 46, 33, 51, -1},
                                                  {50, 51, 46, 0, -1},
                                                  {27, 28, 26, 44, -1}}};
constexpr std::array<int, 5> D_size = {5, 4, 4, 4, 4};
constexpr std::array<int, 5> u_shift = {0, 56, 31, 46, 43};
constexpr std::array<int, 4> p_shift = {23, 24, 48, 2};
inline int mod64(int x) { return (x % N + N) % N; }
inline std::string name(const std::string &prefix, int i) {
    return prefix + "_" + (i < 10 ? "0" : "") + std::to_string(i);
}
inline std::vector<int> support(int row, int pos) {
    std::vector<int> s;
    if (row == 0) s.push_back(pos);
    for (int d : p_shift) s.push_back(mod64(pos + u_shift[row] + d));
    return s;
}
inline std::array<int, 2> product_rows(int row) {
    return {(row + 1) % 5, (row + 2) % 5};
}
struct Options {
    std::optional<double> time_limit;
    int budget = 64;
    std::string prefix;
};
inline Options options(int argc, char **argv, const std::string &default_prefix) {
    if (argc > 4)
        throw std::invalid_argument("用法: PROGRAM [时间限制秒|none] [预算0..64] [输出前缀]");
    Options o;
    o.prefix = default_prefix;
    if (argc >= 2 && std::string(argv[1]) != "none") {
        std::size_t used = 0;
        const std::string raw = argv[1];
        const double value = std::stod(raw, &used);
        if (used != raw.size() || !std::isfinite(value) || value <= 0)
            throw std::invalid_argument("时间限制必须为正数，或使用 none");
        o.time_limit = value;
    }
    if (argc >= 3) {
        std::size_t used = 0;
        const std::string raw = argv[2];
        o.budget = std::stoi(raw, &used);
        if (used != raw.size()) throw std::invalid_argument("预算必须是 0..64 的整数");
    }
    if (argc >= 4) o.prefix = argv[3];
    if (o.budget < 0 || o.budget > 64)
        throw std::invalid_argument("预算必须在 0..64 之间");
    if (o.prefix.empty()) throw std::invalid_argument("输出前缀不能为空");
    return o;
}

inline GRBModel *active_model = nullptr;
inline void handle_sigint(int) {
    if (active_model != nullptr) active_model->terminate();
}
class ActiveModelGuard {
public:
    explicit ActiveModelGuard(GRBModel &model) { active_model = &model; }
    ~ActiveModelGuard() { active_model = nullptr; }
    ActiveModelGuard(const ActiveModelGuard &) = delete;
    ActiveModelGuard &operator=(const ActiveModelGuard &) = delete;
};

inline std::array<GRBVar, N> variables(GRBModel &model, const std::string &prefix) {
    std::array<GRBVar, N> result;
    for (int i = 0; i < N; ++i)
        result[i] = model.addVar(0, 1, 0, GRB_BINARY, name(prefix, i));
    return result;
}
inline void finish(GRBModel &model, const Options &o,
                   const std::array<GRBVar, N> &guess,
                   const std::array<GRBVar, N> &selected,
                   const std::string &guess_label) {
    if (o.time_limit) model.set(GRB_DoubleParam_TimeLimit, *o.time_limit);
    // 导出求解前的完整 LP，因此这里确实需要应用待处理的模型修改。
    model.update();
    model.write(o.prefix + ".lp");
    model.optimize();
    const int status = model.get(GRB_IntAttr_Status);
    const int solutions = model.get(GRB_IntAttr_SolCount);
    std::printf("状态码=%d 求解时间=%.3f 秒 解数=%d\n", status,
                model.get(GRB_DoubleAttr_Runtime), solutions);
    if (solutions == 0) return;
    model.write(o.prefix + ".sol");
    int g_count = 0, l_count = 0;
    std::printf("%s位置:", guess_label.c_str());
    for (int i = 0; i < N; ++i)
        if (guess[i].get(GRB_DoubleAttr_X) > 0.5) {
            std::printf(" %d", i);
            ++g_count;
        }
    std::printf("\n线性化目标位置:");
    for (int j = 0; j < N; ++j)
        if (selected[j].get(GRB_DoubleAttr_X) > 0.5) {
            std::printf(" %d", j);
            ++l_count;
        }
    std::printf("\n猜测数=%d 目标数=%d 预算使用=%d\n", g_count, l_count,
                g_count + l_count);
    std::printf("当前可行目标值=%.10g", model.get(GRB_DoubleAttr_ObjVal));
    if (status == GRB_OPTIMAL) std::printf("（已证明最优）");
    const double bound = model.get(GRB_DoubleAttr_ObjBound);
    const double gap = model.get(GRB_DoubleAttr_MIPGap);
    if (std::isfinite(bound)) std::printf(" 最优界=%.10g", bound);
    if (std::isfinite(gap)) std::printf(" MIPGap=%.6g", gap);
    std::printf("\n模型: %s.lp; 解: %s.sol\n", o.prefix.c_str(), o.prefix.c_str());
}
} // namespace gaston

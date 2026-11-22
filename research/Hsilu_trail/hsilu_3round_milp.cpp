// Hsilu: two chi layers and three linear layers.
//
// The modeled topology is
//
//     state0 -- L0 -- chi0 -- L1 -- chi1 -- L2 -- state5.
//
// The first two chi layers use one of the saved chi5 inner models from
// chi5_fixed_budget/models/best_{13,14,15}.json.  The output of L2 is the
// input to the third chi layer, so the third chi box is not constrained here;
// only its input-row cost is included in the objective.  The objective is an
// estimate: an input row weight 1, 2, or at least 3 is charged 2, 3, or 4,
// respectively.  The reported exact weight is recomputed from complete
// chi5 DDT rows after an incumbent has been found.

#include <gurobi_c++.h>

#include "../hsilu_ref/hsilu.h"

#include <array>
#include <cctype>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int ROWS = 5;
constexpr int BITS = 64;
constexpr int CHI_ROWS = 64;
constexpr int CHI_BITS = 5;

using StateVars = std::array<std::array<GRBVar, BITS>, ROWS>;
using StateValue = std::array<uint64_t, ROWS>;
using ChiRow = std::array<int, 11>;

GRBModel *active_model = nullptr;

void handle_sigint(int) {
    if (active_model != nullptr) active_model->terminate();
}

class ActiveModelGuard {
public:
    explicit ActiveModelGuard(GRBModel &model) { active_model = &model; }
    ~ActiveModelGuard() { active_model = nullptr; }
    ActiveModelGuard(const ActiveModelGuard &) = delete;
    ActiveModelGuard &operator=(const ActiveModelGuard &) = delete;
};

struct Options {
    std::string model_path = "chi5_fixed_budget/models/best_13.json";
    std::string prefix = "hsilu_3round";
};

void usage(const char *program) {
    std::printf(
        "用法: %s [--model PATH] [--prefix PREFIX]\n"
        "默认模型: chi5_fixed_budget/models/best_13.json\n"
        "拓扑: L0 -> chi0 -> L1 -> chi1 -> L2\n",
        program);
}

Options parse_options(int argc, char **argv) {
    Options result;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            std::exit(0);
        }
        if (arg == "--model" || arg == "--prefix") {
            if (i + 1 >= argc) throw std::invalid_argument(arg + " 缺少参数");
            const std::string value = argv[++i];
            if (arg == "--model") result.model_path = value;
            else result.prefix = value;
            continue;
        }
        throw std::invalid_argument("未知参数: " + arg);
    }
    if (result.model_path.empty() || result.prefix.empty())
        throw std::invalid_argument("模型路径和输出前缀不能为空");
    return result;
}

// Parse only the nested integer array following the JSON key "rows".  The
// saved model files are intentionally simple and contain exactly 11 integers
// per row, so no external JSON library is needed by this Gurobi program.
std::vector<ChiRow> load_chi_rows(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("无法打开 χ5 模型文件: " + path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    const std::string text = buffer.str();
    const std::size_t key = text.find("\"rows\"");
    if (key == std::string::npos) throw std::runtime_error("模型文件没有 rows 字段");
    const std::size_t open = text.find('[', key);
    if (open == std::string::npos) throw std::runtime_error("rows 字段不是数组");
    int depth = 0;
    std::size_t close = std::string::npos;
    for (std::size_t i = open; i < text.size(); ++i) {
        if (text[i] == '[') ++depth;
        else if (text[i] == ']') {
            --depth;
            if (depth == 0) {
                close = i;
                break;
            }
        }
    }
    if (close == std::string::npos) throw std::runtime_error("rows 数组未闭合");

    std::vector<int> numbers;
    for (std::size_t i = open + 1; i < close;) {
        if (text[i] == '-' || std::isdigit(static_cast<unsigned char>(text[i]))) {
            const char *start = text.c_str() + i;
            char *end = nullptr;
            const long value = std::strtol(start, &end, 10);
            if (end == start) throw std::runtime_error("rows 中存在无法解析的整数");
            if (value < std::numeric_limits<int>::min() ||
                value > std::numeric_limits<int>::max())
                throw std::runtime_error("rows 中的系数超出 int 范围");
            numbers.push_back(static_cast<int>(value));
            i = static_cast<std::size_t>(end - text.c_str());
        } else {
            ++i;
        }
    }
    if (numbers.empty() || numbers.size() % 11 != 0)
        throw std::runtime_error("rows 中的整数数量不是 11 的整数倍");
    std::vector<ChiRow> rows(numbers.size() / 11);
    for (std::size_t r = 0; r < rows.size(); ++r)
        for (int j = 0; j < 11; ++j) rows[r][j] = numbers[11 * r + j];
    return rows;
}

std::string bit_name(const std::string &prefix, int row, int bit) {
    return prefix + "_r" + std::to_string(row) + "_b" + std::to_string(bit);
}

StateVars add_state(GRBModel &model, const std::string &prefix) {
    StateVars result;
    for (int row = 0; row < ROWS; ++row)
        for (int bit = 0; bit < BITS; ++bit)
            result[row][bit] = model.addVar(0, 1, 0, GRB_BINARY,
                                            bit_name(prefix, row, bit));
    return result;
}

struct LinearAux {
    std::array<std::array<GRBVar, BITS>, NXOR> generated;
};

void add_xor(GRBModel &model, const GRBVar &a, const GRBVar &b,
             const GRBVar &c, const std::string &name) {
    model.addConstr(a + b + c <= 2, name + "_sum");
    model.addConstr(a + b >= c, name + "_a_b_ge_c");
    model.addConstr(a + c >= b, name + "_a_c_ge_b");
    model.addConstr(b + c >= a, name + "_b_c_ge_a");
}

GRBVar word_bit(const StateVars &input, const LinearAux &aux, int index, int bit) {
    if (index < ROWS) return input[index][bit];
    return aux.generated[index - ROWS][bit];
}

void add_linear_layer(GRBModel &model, const StateVars &input,
                      const StateVars &output, const std::string &prefix) {
    LinearAux aux;
    for (int i = 0; i < NXOR; ++i)
        for (int bit = 0; bit < BITS; ++bit) {
            aux.generated[i][bit] = model.addVar(
                0, 1, 0, GRB_BINARY, prefix + "_x" + std::to_string(ROWS + i) +
                    "_b" + std::to_string(bit));
        }

    for (int i = 0; i < NXOR; ++i)
        for (int bit = 0; bit < BITS; ++bit) {
            const int rotated_bit = (bit - m[i] + BITS) % BITS;
            const GRBVar left = word_bit(input, aux, n1[i], bit);
            const GRBVar right = word_bit(input, aux, n2[i], rotated_bit);
            add_xor(model, left, right, aux.generated[i][bit],
                    prefix + "_xor" + std::to_string(i) + "_b" +
                        std::to_string(bit));
        }

    for (int row = 0; row < ROWS; ++row)
        for (int bit = 0; bit < BITS; ++bit) {
            const int source_bit = (bit - r[row] + BITS) % BITS;
            model.addConstr(output[row][bit] ==
                                aux.generated[8 + row][source_bit],
                            prefix + "_rho_r" + std::to_string(row) + "_b" +
                                std::to_string(bit));
        }
}

GRBLinExpr add_chi_layer(GRBModel &model, const StateVars &input,
                   const StateVars &output, const std::vector<ChiRow> &rows,
                   const std::string &prefix, GRBLinExpr &objective) {
    GRBLinExpr active_total = 0;
    const int estimated_cost[6] = {0, 2, 3, 4, 4, 4};
    for (int bit = 0; bit < CHI_ROWS; ++bit) {
        std::array<GRBVar, 6> onehot;
        for (int h = 0; h <= 5; ++h)
            onehot[h] = model.addVar(0, 1, 0, GRB_BINARY,
                                     prefix + "_weight" + std::to_string(bit) +
                                         "_h" + std::to_string(h));
        GRBLinExpr select_sum = 0;
        GRBLinExpr input_weight = 0;
        for (int h = 0; h <= 5; ++h) {
            select_sum += onehot[h];
            input_weight += h * onehot[h];
            objective += estimated_cost[h] * onehot[h];
        }
        active_total += 1 - onehot[0];
        model.addConstr(select_sum == 1,
                        prefix + "_weight_onehot_" + std::to_string(bit));
        GRBLinExpr actual_weight = 0;
        for (int row = 0; row < ROWS; ++row) actual_weight += input[row][bit];
        model.addConstr(actual_weight == input_weight,
                        prefix + "_weight_link_" + std::to_string(bit));

        for (const ChiRow &row : rows) {
            GRBLinExpr expression = row[10];
            for (int i = 0; i < CHI_BITS; ++i)
                expression += row[i] * input[i][bit];
            for (int i = 0; i < CHI_BITS; ++i)
                expression += row[CHI_BITS + i] * output[i][bit];
            model.addConstr(expression >= 0,
                            prefix + "_ineq" + std::to_string(bit) + "_" +
                                std::to_string(&row - rows.data()));
        }
    }
    return active_total;
}

GRBLinExpr add_input_weight_objective(GRBModel &model, const StateVars &input,
                                const std::string &prefix,
                                GRBLinExpr &objective) {
    GRBLinExpr active_total = 0;
    const int estimated_cost[6] = {0, 2, 3, 4, 4, 4};
    for (int bit = 0; bit < CHI_ROWS; ++bit) {
        std::array<GRBVar, 6> onehot;
        GRBLinExpr select_sum = 0;
        GRBLinExpr input_weight = 0;
        for (int h = 0; h <= 5; ++h) {
            onehot[h] = model.addVar(0, 1, 0, GRB_BINARY,
                                     prefix + "_weight" + std::to_string(bit) +
                                         "_h" + std::to_string(h));
            select_sum += onehot[h];
            input_weight += h * onehot[h];
            objective += estimated_cost[h] * onehot[h];
        }
        active_total += 1 - onehot[0];
        model.addConstr(select_sum == 1,
                        prefix + "_weight_onehot_" + std::to_string(bit));
        GRBLinExpr actual_weight = 0;
        for (int row = 0; row < ROWS; ++row) actual_weight += input[row][bit];
        model.addConstr(actual_weight == input_weight,
                        prefix + "_weight_link_" + std::to_string(bit));
    }
    return active_total;
}

uint32_t chi5(uint32_t x) {
    uint32_t result = 0;
    for (int i = 0; i < 5; ++i) {
        const uint32_t xi = (x >> i) & 1U;
        const uint32_t x1 = (x >> ((i + 1) % 5)) & 1U;
        const uint32_t x2 = (x >> ((i + 2) % 5)) & 1U;
        result |= (xi ^ ((1U ^ x1) & x2)) << i;
    }
    return result;
}

int ddt_count(int alpha, int beta) {
    int count = 0;
    for (int x = 0; x < 32; ++x)
        if ((chi5(static_cast<uint32_t>(x)) ^
             chi5(static_cast<uint32_t>(x ^ alpha))) ==
            static_cast<uint32_t>(beta))
            ++count;
    return count;
}

int exact_weight(int alpha, int beta) {
    const int count = ddt_count(alpha, beta);
    if (count == 0) return -1;
    if (count == 32) return 0;
    if (count == 8) return 2;
    if (count == 4) return 3;
    if (count == 2) return 4;
    throw std::runtime_error("χ5 DDT 出现未预期的非零计数");
}

int exact_input_weight(int alpha) {
    if (alpha == 0) return 0;
    int first_count = 0;
    for (int beta = 0; beta < 32; ++beta) {
        const int count = ddt_count(alpha, beta);
        if (count != 0) {
            if (first_count == 0) first_count = count;
            else if (first_count != count)
                throw std::runtime_error("χ5 固定输入差分的非零 DDT 计数不一致");
        }
    }
    if (first_count == 0) throw std::runtime_error("χ5 非零输入差分没有合法输出");
    if (first_count == 8) return 2;
    if (first_count == 4) return 3;
    if (first_count == 2) return 4;
    throw std::runtime_error("χ5 DDT 输入行出现未预期的计数");
}

StateValue read_state(const StateVars &state) {
    StateValue result{};
    for (int row = 0; row < ROWS; ++row)
        for (int bit = 0; bit < BITS; ++bit)
            if (state[row][bit].get(GRB_DoubleAttr_X) > 0.5)
                result[row] |= uint64_t(1) << bit;
    return result;
}

void print_state(const char *label, const StateValue &state) {
    std::printf("%s", label);
    for (int row = 0; row < ROWS; ++row)
        std::printf(" %016llx", static_cast<unsigned long long>(state[row]));
    std::printf("\n");
}

void print_chi_report(const char *label, const StateValue &input,
                      const StateValue &output) {
    int estimated = 0;
    int exact = 0;
    int active = 0;
    std::printf("%s\n", label);
    for (int bit = 0; bit < BITS; ++bit) {
        int alpha = 0;
        int beta = 0;
        for (int row = 0; row < ROWS; ++row) {
            alpha |= static_cast<int>((input[row] >> bit) & 1U) << row;
            beta |= static_cast<int>((output[row] >> bit) & 1U) << row;
        }
        const int h = __builtin_popcount(static_cast<unsigned>(alpha));
        const int estimate = h == 0 ? 0 : (h == 1 ? 2 : (h == 2 ? 3 : 4));
        const int weight = exact_weight(alpha, beta);
        if (h != 0) ++active;
        estimated += estimate;
        if (weight >= 0) exact += weight;
        else
            throw std::runtime_error("MILP 解中出现不属于精确 χ5 支集的局部转移");
        if (h != 0)
            std::printf("  z=%2d alpha=%02x beta=%02x h=%d estimate=%d exact=%d\n",
                        bit, alpha, beta, h, estimate, weight);
    }
    std::printf("  活跃 χ5 数=%d 估计重量=%d 精确重量=%d\n",
                active, estimated, exact);
}

void print_input_report(const char *label, const StateValue &input) {
    int estimated = 0;
    int exact = 0;
    int active = 0;
    std::printf("%s\n", label);
    for (int bit = 0; bit < BITS; ++bit) {
        int alpha = 0;
        for (int row = 0; row < ROWS; ++row)
            alpha |= static_cast<int>((input[row] >> bit) & 1U) << row;
        const int h = __builtin_popcount(static_cast<unsigned>(alpha));
        const int estimate = h == 0 ? 0 : (h == 1 ? 2 : (h == 2 ? 3 : 4));
        const int weight = exact_input_weight(alpha);
        if (h != 0) ++active;
        estimated += estimate;
        exact += weight;
        if (h != 0)
            std::printf("  z=%2d alpha=%02x h=%d estimate=%d exact=%d\n",
                        bit, alpha, h, estimate, weight);
    }
    std::printf("  活跃 χ5 数=%d 估计重量=%d 精确重量=%d\n",
                active, estimated, exact);
}

void print_status(GRBModel &model) {
    std::printf("状态码=%d 求解时间=%.3f 秒 解数=%d\n",
                model.get(GRB_IntAttr_Status), model.get(GRB_DoubleAttr_Runtime),
                model.get(GRB_IntAttr_SolCount));
}

} // namespace

int main(int argc, char **argv) {
    std::signal(SIGINT, handle_sigint);
    try {
        const Options options = parse_options(argc, argv);
        const std::vector<ChiRow> chi_rows = load_chi_rows(options.model_path);
        std::printf("χ5 模型: %s，行数=%zu\n", options.model_path.c_str(),
                    chi_rows.size());

        GRBEnv env;
        GRBModel model(env);
        ActiveModelGuard guard(model);
        model.set(GRB_StringAttr_ModelName, "Hsilu_two_chi_three_linear_trail");

        std::array<StateVars, 6> states;
        for (int layer = 0; layer < 6; ++layer)
            states[layer] = add_state(model, "s" + std::to_string(layer));

        add_linear_layer(model, states[0], states[1], "L0");
        GRBLinExpr objective = 0;
        const GRBLinExpr active0 =
            add_chi_layer(model, states[1], states[2], chi_rows, "chi0", objective);
        add_linear_layer(model, states[2], states[3], "L1");
        const GRBLinExpr active1 =
            add_chi_layer(model, states[3], states[4], chi_rows, "chi1", objective);
        add_linear_layer(model, states[4], states[5], "L2");
        // The third chi layer is intentionally not modeled.  Its input is
        // state5, and its estimated/exact row weight is still counted.
        const GRBLinExpr active2 =
            add_input_weight_objective(model, states[5], "chi2_input", objective);
        // Reuse h=0 indicators: a box is active exactly when onehot[0]=0.
        model.addConstr(active0 == 6, "chi0_active_count");
        model.addConstr(active1 == 6, "chi1_active_count");
        model.addConstr(active2 == 27, "chi2_active_count");
        std::printf("三层活跃 χ5 数固定为 6, 6, 27；各层重量由目标函数优化。\n");

        // Exclude the trivial all-zero differential path.
        GRBLinExpr nonzero_input = 0;
        for (int row = 0; row < ROWS; ++row)
            for (int bit = 0; bit < BITS; ++bit) nonzero_input += states[0][row][bit];
        model.addConstr(nonzero_input >= 1, "nonzero_input_difference");
        model.setObjective(objective, GRB_MINIMIZE);

        model.update();
        model.write(options.prefix + ".lp");
        model.optimize();
        print_status(model);
        if (model.get(GRB_IntAttr_SolCount) == 0) return 0;

        model.write(options.prefix + ".sol");
        const int status = model.get(GRB_IntAttr_Status);
        std::printf("当前可行估计目标值=%.10g", model.get(GRB_DoubleAttr_ObjVal));
        if (status == GRB_OPTIMAL) std::printf("（已证明该 MILP 目标最优）");
        std::printf(" 最优界=%.10g MIPGap=%.6g\n",
                    model.get(GRB_DoubleAttr_ObjBound),
                    model.get(GRB_DoubleAttr_MIPGap));

        std::array<StateValue, 6> values;
        for (int layer = 0; layer < 6; ++layer) {
            values[layer] = read_state(states[layer]);
            print_state(("state" + std::to_string(layer) + ":").c_str(),
                        values[layer]);
        }
        print_chi_report("chi0 局部转移:", values[1], values[2]);
        print_chi_report("chi1 局部转移:", values[3], values[4]);
        print_input_report("chi2 输入（L2 输出；未建模第三个 χ5）:", values[5]);

        int total_estimated = 0;
        int total_exact = 0;
        for (int bit = 0; bit < BITS; ++bit) {
            int a0 = 0, a1 = 0, a2 = 0;
            for (int row = 0; row < ROWS; ++row) {
                a0 |= static_cast<int>((values[1][row] >> bit) & 1U) << row;
                a1 |= static_cast<int>((values[3][row] >> bit) & 1U) << row;
                a2 |= static_cast<int>((values[5][row] >> bit) & 1U) << row;
            }
            const int h0 = __builtin_popcount(static_cast<unsigned>(a0));
            const int h1 = __builtin_popcount(static_cast<unsigned>(a1));
            const int h2 = __builtin_popcount(static_cast<unsigned>(a2));
            total_estimated += (h0 == 0 ? 0 : (h0 == 1 ? 2 : (h0 == 2 ? 3 : 4)));
            total_estimated += (h1 == 0 ? 0 : (h1 == 1 ? 2 : (h1 == 2 ? 3 : 4)));
            total_estimated += (h2 == 0 ? 0 : (h2 == 1 ? 2 : (h2 == 2 ? 3 : 4)));
            total_exact += exact_input_weight(a0) + exact_input_weight(a1) +
                           exact_input_weight(a2);
        }
        std::printf("三层 χ 输入的估计总差分重量=%d，按完整 χ5 DDT 重算的精确总差分重量=%d\n",
                    total_estimated, total_exact);

        std::printf("模型文件: %s.lp; 解文件: %s.sol\n",
                    options.prefix.c_str(), options.prefix.c_str());
    } catch (const GRBException &e) {
        std::fprintf(stderr, "Gurobi 错误 %d: %s\n", e.getErrorCode(),
                     e.getMessage().c_str());
        return 1;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "错误: %s\n", e.what());
        return 1;
    }
    return 0;
}

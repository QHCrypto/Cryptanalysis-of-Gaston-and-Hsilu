// Gaston, Algorithm 1 and Table 6 of ePrint 2023/799.
// Bit j is the coefficient of 2^j. rol(v, s) moves bit j to j+s mod 64.
#include <gurobi_c++.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
constexpr int R = 5, N = 64, B = R * N;
constexpr std::array<int, R> EAST{0, 60, 22, 27, 4};
constexpr std::array<int, R> WEST{0, 56, 31, 46, 43};
constexpr std::array<int, R> TWIN{25, 32, 52, 60, 63};
using State = std::array<uint64_t, R>;
using WordVars = std::array<GRBVar, N>;
using StateVars = std::array<WordVars, R>;
using Matrix = std::array<std::bitset<B>, B>;

GRBModel* active_model = nullptr;
void handle_sigint(int) {
    if (active_model != nullptr) active_model->terminate();
}
class ActiveModelGuard {
public:
    explicit ActiveModelGuard(GRBModel& model) { active_model = &model; }
    ~ActiveModelGuard() { active_model = nullptr; }
    ActiveModelGuard(const ActiveModelGuard&) = delete;
    ActiveModelGuard& operator=(const ActiveModelGuard&) = delete;
};

int mod64(int x) { return (x % N + N) % N; }
uint64_t rol(uint64_t x, int s) {
    const int n = mod64(s);
    return n == 0 ? x : (x << n) | (x >> (N - n));
}
// Equation (2), Definition 5: E is evaluated on A=rho_east(input), before
// the common ROL(E,23) addition. E=0 does not require P=Q=0.
uint64_t theta_effect(const State& input) {
    uint64_t p = 0, q = 0;
    for (int i = 0; i < R; ++i) {
        const uint64_t a = rol(input[i], EAST[i]);
        p ^= a;
        q ^= rol(a, TWIN[i]);
    }
    return p ^ rol(p, 1) ^ q ^ rol(q, 18);
}
State linear_layer(const State& x) {
    State a{}, y{};
    uint64_t p = 0, q = 0;
    for (int i = 0; i < R; ++i) {
        a[i] = rol(x[i], EAST[i]);
        p ^= a[i];
        q ^= rol(a[i], TWIN[i]);
    }
    const uint64_t effect = rol(p ^ rol(p, 1) ^ q ^ rol(q, 18), 23);
    for (int i = 0; i < R; ++i) y[i] = rol(a[i] ^ effect, WEST[i]);
    return y;
}
State transpose_layer(const State& beta) {
    State d{}, alpha{};
    uint64_t p = 0;
    for (int i = 0; i < R; ++i) {
        d[i] = rol(beta[i], -WEST[i]);
        p ^= d[i];
    }
    for (int i = 0; i < R; ++i) {
        const uint64_t effect = rol(p, -23) ^ rol(p, -24)
            ^ rol(p, -23 - TWIN[i]) ^ rol(p, -41 - TWIN[i]);
        alpha[i] = rol(d[i] ^ effect, -EAST[i]);
    }
    return alpha;
}
uint64_t active_columns(const State& s) {
    uint64_t result = 0;
    for (uint64_t row : s) result |= row;
    return result;
}
int column_weight(const State& s) {
    return __builtin_popcountll(active_columns(s));
}
int bit_weight(const State& s) {
    int result = 0;
    for (uint64_t row : s) result += __builtin_popcountll(row);
    return result;
}
int bit(const State& s, int i, int j) { return int((s[i] >> j) & 1); }
std::bitset<B> flatten(const State& s) {
    std::bitset<B> result;
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < N; ++j) result[i * N + j] = bit(s, i, j);
    return result;
}
Matrix matrix_of(bool linear) {
    Matrix result{};
    for (int c = 0; c < B; ++c) {
        State x{};
        x[c / N] = uint64_t(1) << (c % N);
        const auto y = flatten(linear ? transpose_layer(x) : linear_layer(x));
        for (int r = 0; r < B; ++r) result[r][c] = y[r];
    }
    return result;
}
Matrix inverse(Matrix a) {
    Matrix inv{};
    for (int i = 0; i < B; ++i) inv[i][i] = 1;
    for (int c = 0; c < B; ++c) {
        int pivot = c;
        while (pivot < B && !a[pivot][c]) ++pivot;
        if (pivot == B) throw std::runtime_error("线性层矩阵不可逆");
        std::swap(a[c], a[pivot]);
        std::swap(inv[c], inv[pivot]);
        for (int r = 0; r < B; ++r) {
            if (r != c && a[r][c]) { a[r] ^= a[c]; inv[r] ^= inv[c]; }
        }
    }
    return inv;
}
State apply_matrix(const Matrix& a, const State& x) {
    const auto v = flatten(x);
    State result{};
    for (int r = 0; r < B; ++r)
        if ((a[r] & v).count() % 2) result[r / N] |= uint64_t(1) << (r % N);
    return result;
}

enum class EffectConstraint { Any, Zero, Nonzero };
const char* effect_constraint_name(EffectConstraint constraint) {
    switch (constraint) {
        case EffectConstraint::Zero: return "zero";
        case EffectConstraint::Nonzero: return "nonzero";
        default: return "any";
    }
}
bool accepts_effect(EffectConstraint constraint, uint64_t effect) {
    return constraint == EffectConstraint::Any
        || (constraint == EffectConstraint::Zero ? effect == 0 : effect != 0);
}
struct Options {
    bool linear = false, input_fixed = true, symmetry = true, export_only = false;
    bool effect_explicit = false;
    EffectConstraint effect_constraint = EffectConstraint::Any;
    int count = -1, threads = 0;
    std::string prefix;
};
int activity_objective(const Options& o, const State& input, const State& output) {
    return o.count < 0 ? column_weight(input) + column_weight(output)
                      : column_weight(o.input_fixed ? output : input);
}
const char* objective_side_name(const Options& o) {
    return o.count < 0 ? "input+output" : (o.input_fixed ? "output" : "input");
}
void usage(const char* program) {
    std::printf(
        "用法: %s [--mode difference|linear] [--input-active K | --output-active K] [选项]\n"
        "  省略两个数量参数：输入的 320 个比特之和 >= 1，最小化输入加输出活跃 S 盒列数\n"
        "  此时默认以输入第 0 列的五个比特之和 >= 1 减少循环移位对称性\n"
        "  --input-active K   恰好 K 个输入活跃 S 盒，最小化输出活跃数\n"
        "  --output-active K  恰好 K 个输出活跃 S 盒，最小化输入活跃数\n"
        "  --mode difference 差分/状态传播: output = lambda(input)，默认模式\n"
        "  --mode linear     线性掩码传播: input = lambda^T(output)\n"
        "  --theta-effect any|zero|nonzero  差分模式的 E 条件，默认 any\n"
        "                    zero: E=0（核内）；nonzero: E!=0（核外）\n"
        "                    E 按 A=rho_east(input) 计算；linear 不接受此参数\n"
        "  --prefix PATH     导出 PATH.lp；有解时导出 PATH.sol 和 PATH.json\n"
        "  --threads N       Gurobi 线程数；0 表示自动（默认）\n"
        "  --no-symmetry     不固定活跃列代表元\n"
        "  --export-only     只构建/导出 LP，不求解；须指定 --prefix\n"
        "  --help            显示帮助\n"
        "K 范围 0..64；活跃指同一列的 5 比特非零。没有求解时限；Ctrl+C 可保留当前解。\n",
        program);
}
int parse_int(const std::string& s, int upper) {
    std::size_t used = 0;
    int v;
    try { v = std::stoi(s, &used); }
    catch (const std::exception&) { throw std::invalid_argument("不是有效整数: " + s); }
    if (used != s.size() || v < 0 || v > upper)
        throw std::invalid_argument("整数超出允许范围: " + s);
    return v;
}
Options options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { usage(argv[0]); std::exit(0); }
        if (arg == "--no-symmetry") { o.symmetry = false; continue; }
        if (arg == "--export-only") { o.export_only = true; continue; }
        if (arg != "--mode" && arg != "--input-active" && arg != "--output-active"
            && arg != "--prefix" && arg != "--threads" && arg != "--theta-effect")
            throw std::invalid_argument("未知参数: " + arg);
        if (++i >= argc) throw std::invalid_argument(arg + " 缺少参数");
        const std::string value = argv[i];
        if (arg == "--mode") {
            if (value != "difference" && value != "linear")
                throw std::invalid_argument("--mode 必须是 difference 或 linear");
            o.linear = value == "linear";
        } else if (arg == "--input-active" || arg == "--output-active") {
            if (o.count != -1) throw std::invalid_argument("最多只能指定一个活跃数量约束");
            o.count = parse_int(value, N);
            o.input_fixed = arg == "--input-active";
        } else if (arg == "--threads") {
            o.threads = parse_int(value, 1024);
        } else if (arg == "--theta-effect") {
            if (o.effect_explicit) throw std::invalid_argument("--theta-effect 只能指定一次");
            o.effect_explicit = true;
            if (value == "any") o.effect_constraint = EffectConstraint::Any;
            else if (value == "zero") o.effect_constraint = EffectConstraint::Zero;
            else if (value == "nonzero") o.effect_constraint = EffectConstraint::Nonzero;
            else throw std::invalid_argument("--theta-effect 必须是 any、zero 或 nonzero");
        } else {
            if (value.empty()) throw std::invalid_argument("--prefix 不能为空");
            o.prefix = value;
        }
    }
    if (o.export_only && o.prefix.empty())
        throw std::invalid_argument("--export-only 需要 --prefix PATH");
    if (o.linear && o.effect_explicit)
        throw std::invalid_argument("--theta-effect 仅支持 difference 模式；linear 模式不能显式使用此参数（包括 any）");
    return o;
}
WordVars add_word(GRBModel& model, const std::string& prefix) {
    WordVars vars;
    for (int j = 0; j < N; ++j)
        vars[j] = model.addVar(0, 1, 0, GRB_BINARY, prefix + "_" + std::to_string(j));
    return vars;
}
StateVars add_state(GRBModel& model, const std::string& prefix) {
    StateVars vars;
    for (int i = 0; i < R; ++i) vars[i] = add_word(model, prefix + "_" + std::to_string(i));
    return vars;
}
// Jeroslow facets for even parity of (inputs, output). For each odd subset S:
// sum(v in S) - sum(v outside S) <= |S|-1. No integer slack variables.
// Direct hulls are used only for arity 2 and 3 (4 and 8 inequalities).
void add_parity_hull(GRBModel& model, const std::vector<GRBVar>& inputs,
                     const GRBVar& output, const std::string& name) {
    if (inputs.empty()) { model.addConstr(output == 0, name + "_zero"); return; }
    if (inputs.size() == 1) { model.addConstr(output == inputs[0], name + "_copy"); return; }
    if (inputs.size() > 3) throw std::logic_error("直接异或凸包最多支持 3 个输入");
    auto vars = inputs;
    vars.push_back(output);
    for (unsigned mask = 1; mask < (1u << vars.size()); ++mask) {
        const int cardinality = __builtin_popcount(mask);
        if (cardinality % 2 == 0) continue;
        GRBLinExpr lhs = 0;
        for (std::size_t j = 0; j < vars.size(); ++j)
            lhs += ((mask >> j) & 1u) ? vars[j] : -vars[j];
        model.addConstr(lhs <= cardinality - 1,
                        name + "_xor" + std::to_string(inputs.size())
                        + "_facet_" + std::to_string(mask));
    }
}
void add_parity(GRBModel& model, const std::vector<GRBVar>& inputs,
                const GRBVar& output, const std::string& name) {
    if (inputs.size() <= 3) { add_parity_hull(model, inputs, output, name); return; }
    // Four inputs use xor3+xor2: 12 constraints and one auxiliary variable.
    // Balanced ternary reduction: xor3 gates and, for even arity, one final
    // xor2. Do not split xor3 into two xor2 gates with an extra variable.
    auto current = inputs;
    int gate = 0;
    while (current.size() > 3) {
        std::vector<GRBVar> next;
        std::size_t j = 0;
        for (; j + 2 < current.size(); j += 3) {
            const std::string node = name + "_tree_" + std::to_string(gate++);
            const auto aux = model.addVar(0, 1, 0, GRB_BINARY, node);
            add_parity_hull(model, {current[j], current[j + 1], current[j + 2]}, aux, node);
            next.push_back(aux);
        }
        for (; j < current.size(); ++j) next.push_back(current[j]);
        current = std::move(next);
    }
    add_parity_hull(model, current, output, name + "_root");
}

// CSE across inverse-matrix rows. In each pass, select frequent, disjoint
// pairs and replace every occurrence by ONE shared xor2 binary variable.
// Counts are recomputed after substitution; ties use wire IDs for determinism.
// xor2 here is the shared-subexpression exception, not a per-row xor2 chain.
void extract_shared_pairs(GRBModel& model, std::vector<GRBVar>& wires,
                          std::vector<std::vector<int>>& expressions) {
    while (true) {
        std::unordered_map<uint64_t, int> counts;
        for (const auto& expr : expressions)
            for (std::size_t i = 0; i < expr.size(); ++i)
                for (std::size_t j = i + 1; j < expr.size(); ++j)
                    ++counts[(uint64_t(expr[i]) << 32) | uint64_t(expr[j])];
        std::vector<std::pair<uint64_t, int>> candidates;
        for (const auto& pair : counts)
            if (pair.second >= 2) candidates.push_back(pair);
        if (candidates.empty()) break;
        std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
            return a.second != b.second ? a.second > b.second : a.first < b.first;
        });
        std::vector<int> partner(wires.size(), -1), replacement(wires.size(), -1);
        for (const auto& pair : candidates) {
            const int a = int(pair.first >> 32), b = int(pair.first & 0xffffffffu);
            if (partner[a] != -1 || partner[b] != -1) continue;
            const int id = int(wires.size());
            const std::string name = "inverse_cse_" + std::to_string(id);
            const auto aux = model.addVar(0, 1, 0, GRB_BINARY, name);
            add_parity_hull(model, {wires[a], wires[b]}, aux, name);
            wires.push_back(aux);
            partner[a] = b;
            partner[b] = a;
            replacement[a] = replacement[b] = id;
        }
        for (auto& expr : expressions) {
            std::vector<int> next;
            for (int a : expr) {
                const int b = partner[a];
                if (b >= 0 && std::binary_search(expr.begin(), expr.end(), b)) {
                    if (a < b) next.push_back(replacement[a]);
                } else {
                    next.push_back(a);
                }
            }
            std::sort(next.begin(), next.end());
            expr = std::move(next);
        }
    }
}
void add_difference(GRBModel& model, const StateVars& x, const StateVars& y,
                    EffectConstraint constraint = EffectConstraint::Any) {
    const auto p = add_word(model, "p");
    const auto q = add_word(model, "q");
    const auto e = add_word(model, "effect");
    for (int j = 0; j < N; ++j) {
        std::vector<GRBVar> pbits, qbits;
        for (int i = 0; i < R; ++i) {
            pbits.push_back(x[i][mod64(j - EAST[i])]);
            qbits.push_back(x[i][mod64(j - EAST[i] - TWIN[i])]);
        }
        add_parity(model, pbits, p[j], "p_eq_" + std::to_string(j));
        add_parity(model, qbits, q[j], "q_eq_" + std::to_string(j));
        add_parity(model, {p[j], p[mod64(j - 1)], q[j], q[mod64(j - 18)]},
                   e[j], "effect_eq_" + std::to_string(j));
        for (int i = 0; i < R; ++i)
            add_parity(model, {x[i][mod64(j - WEST[i] - EAST[i])],
                              e[mod64(j - WEST[i] - 23)]},
                       y[i][j], "y_eq_" + std::to_string(i) + "_" + std::to_string(j));
    }
    if (constraint == EffectConstraint::Zero) {
        for (int j = 0; j < N; ++j)
            model.addConstr(e[j] == 0, "theta_effect_zero_" + std::to_string(j));
    } else if (constraint == EffectConstraint::Nonzero) {
        GRBLinExpr sum = 0;
        for (const auto& v : e) sum += v;
        model.addConstr(sum >= 1, "theta_effect_nonzero");
    }
}
void add_linear(GRBModel& model, const StateVars& alpha, const StateVars& beta) {
    const auto p = add_word(model, "p");
    const auto f = add_word(model, "f");
    const auto g = add_word(model, "g");
    for (int j = 0; j < N; ++j) {
        std::vector<GRBVar> pbits;
        for (int i = 0; i < R; ++i) pbits.push_back(beta[i][mod64(j + WEST[i])]);
        add_parity(model, pbits, p[j], "p_eq_" + std::to_string(j));
        add_parity(model, {p[mod64(j + 23)], p[mod64(j + 24)]}, f[j],
                   "f_eq_" + std::to_string(j));
        add_parity(model, {p[mod64(j + 23)], p[mod64(j + 41)]}, g[j],
                   "g_eq_" + std::to_string(j));
        for (int i = 0; i < R; ++i)
            add_parity(model, {beta[i][mod64(j + EAST[i] + WEST[i])],
                              f[mod64(j + EAST[i])], g[mod64(j + EAST[i] + TWIN[i])]},
                       alpha[i][j], "alpha_eq_" + std::to_string(i) + "_" + std::to_string(j));
    }
}
// When fixing the opposite end, expose its bits directly via the exact inverse
// over GF(2). This avoids asking the MILP presolver to perform binary elimination.
void add_inverse(GRBModel& model, const StateVars& fixed, const StateVars& other,
                 bool linear) {
    const auto inv = inverse(matrix_of(linear));
    std::vector<GRBVar> wires;
    for (const auto& row : fixed)
        for (const auto& v : row) wires.push_back(v);
    std::vector<std::vector<int>> expressions(B);
    for (int r = 0; r < B; ++r)
        for (int c = 0; c < B; ++c)
            if (inv[r][c]) expressions[r].push_back(c);
    extract_shared_pairs(model, wires, expressions);
    for (int r = 0; r < B; ++r) {
        std::vector<GRBVar> inputs;
        for (int id : expressions[r]) inputs.push_back(wires[id]);
        add_parity(model, inputs, other[r / N][r % N], "inverse_eq_" + std::to_string(r));
    }
}
WordVars add_activity(GRBModel& model, const StateVars& state, const std::string& name) {
    const auto active = add_word(model, name);
    for (int j = 0; j < N; ++j) {
        GRBLinExpr sum = 0;
        for (int i = 0; i < R; ++i) {
            model.addConstr(state[i][j] <= active[j]);
            sum += state[i][j];
        }
        model.addConstr(active[j] <= sum);
    }
    return active;
}
// Feasible MIP start only: never fixes a state or adds an objective bound.
std::optional<std::pair<State, State>> initial_pair(const Options& o) {
    const Matrix forward = matrix_of(o.linear);
    // Without a fixed count, sample one active column on either side and
    // compare total activity. These samples never constrain the search space.
    const int sample_count = o.count < 0 ? 1 : o.count;
    std::mt19937_64 rng(20261014);
    std::optional<std::pair<State, State>> best;
    int best_weight = 2 * N + 1;
    for (int side = 0; side < (o.count < 0 ? 2 : 1); ++side) {
        const bool sample_input = o.count < 0 ? side == 0 : o.input_fixed;
        const bool use_direct = sample_input != o.linear;
        const Matrix map = use_direct ? forward : inverse(forward);
        // For one active column, include all 31 nonzero column values.
        const int trials = sample_count == 1 ? 31 : 128;
        for (int trial = 0; trial < trials; ++trial) {
            State fixed{};
            for (int j = 0; j < sample_count; ++j) {
                const unsigned value = sample_count == 1 ? unsigned(trial + 1)
                                                   : unsigned(rng() % 31 + 1);
                for (int i = 0; i < R; ++i)
                    if ((value >> i) & 1) fixed[i] |= uint64_t(1) << j;
            }
            const State other = apply_matrix(map, fixed);
            auto pair = sample_input ? std::make_pair(fixed, other) : std::make_pair(other, fixed);
            if (o.count < 0 && o.symmetry) {
                const uint64_t columns = active_columns(pair.first);
                if (columns == 0) continue;
                const int shift = -__builtin_ctzll(columns);
                for (auto& row : pair.first) row = rol(row, shift);
                for (auto& row : pair.second) row = rol(row, shift);
            }
            if (!o.linear && !accepts_effect(o.effect_constraint,
                                            theta_effect(pair.first))) continue;
            const int weight = activity_objective(o, pair.first, pair.second);
            if (weight < best_weight) { best = pair; best_weight = weight; }
        }
    }
    return best;
}
void start_state(StateVars& vars, WordVars& activity, const State& value) {
    const uint64_t columns = active_columns(value);
    for (int j = 0; j < N; ++j) {
        activity[j].set(GRB_DoubleAttr_Start, int((columns >> j) & 1));
        for (int i = 0; i < R; ++i) vars[i][j].set(GRB_DoubleAttr_Start, bit(value, i, j));
    }
}
State read_state(const StateVars& vars) {
    State s{};
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < N; ++j)
            if (vars[i][j].get(GRB_DoubleAttr_X) > 0.5) s[i] |= uint64_t(1) << j;
    return s;
}
void verify_activity(const State& s, const WordVars& vars) {
    const uint64_t columns = active_columns(s);
    for (int j = 0; j < N; ++j)
        if (int((columns >> j) & 1) != int(vars[j].get(GRB_DoubleAttr_X) > 0.5))
            throw std::runtime_error("解的活跃标记与实际比特不一致");
}
void print_state(const char* label, const State& s) {
    std::printf("%s: active=%d; bits=%d\n  columns:", label, column_weight(s), bit_weight(s));
    const uint64_t columns = active_columns(s);
    for (int j = 0; j < N; ++j) if ((columns >> j) & 1) std::printf(" %d", j);
    if (!columns) std::printf(" (none)");
    std::printf("\n");
    for (int i = 0; i < R; ++i) {
        std::printf("  row %d = 0x%016llx  bits[0..63]: ", i, (unsigned long long)s[i]);
        for (int j = 0; j < N; ++j) std::printf("%c", bit(s, i, j) ? '1' : '.');
        std::printf("\n");
    }
}
void write_json(const Options& o, int status, double objective, double bound,
                double gap, double runtime, const State& input, const State& output) {
    std::ofstream file(o.prefix + ".json");
    if (!file) throw std::runtime_error("无法写入 JSON: " + o.prefix + ".json");
    file.precision(17);
    file << "{\n  \"mode\": \"" << (o.linear ? "linear" : "difference") << '"';
    if (o.count < 0) {
        file << ",\n  \"fixed_side\": null,\n  \"fixed_count\": null"
             << ",\n  \"input_bit_weight_min\": 1";
    } else {
        file << ",\n  \"fixed_side\": \"" << (o.input_fixed ? "input" : "output")
             << "\",\n  \"fixed_count\": " << o.count;
    }
    file << ",\n  \"objective_side\": \"" << objective_side_name(o) << '"'
         << ",\n  \"total_active\": " << column_weight(input) + column_weight(output)
         << ",\n  \"status\": " << status
         << ",\n  \"proven_optimal\": " << (status == GRB_OPTIMAL ? "true" : "false")
         << ",\n  \"objective\": " << objective << ",\n  \"bound\": ";
    if (std::isfinite(bound) && std::abs(bound) < GRB_INFINITY) file << bound;
    else file << "null";
    file << ",\n  \"gap\": ";
    if (std::isfinite(gap) && gap < GRB_INFINITY) file << gap;
    else file << "null";
    file << ",\n  \"runtime_seconds\": " << runtime
         << ",\n  \"verified\": true,\n  \"bit_index\": \"LSB=0; rol moves j to j+s mod 64\"";
    if (!o.linear) {
        const uint64_t effect = theta_effect(input);
        char hex[19];
        std::snprintf(hex, sizeof(hex), "0x%016llx", (unsigned long long)effect);
        file << ",\n  \"theta_effect_constraint\": \"" << effect_constraint_name(o.effect_constraint) << '"'
             << ",\n  \"theta_effect_hex\": \"" << hex << '"'
             << ",\n  \"theta_effect_weight\": " << __builtin_popcountll(effect)
             << ",\n  \"in_theta_kernel\": " << (effect == 0 ? "true" : "false");
    }
    for (int side = 0; side < 2; ++side) {
        const State& s = side == 0 ? input : output;
        const char* name = side == 0 ? "input" : "output";
        file << ",\n  \"" << name << "_active\": " << column_weight(s)
             << ",\n  \"" << name << "_bit_weight\": " << bit_weight(s)
             << ",\n  \"" << name << "_columns\": [";
        bool first = true;
        for (int j = 0; j < N; ++j) if ((active_columns(s) >> j) & 1) {
            if (!first) file << ", ";
            file << j;
            first = false;
        }
        file << "],\n  \"" << name << "_words\": [";
        for (int i = 0; i < R; ++i) {
            char hex[19];
            std::snprintf(hex, sizeof(hex), "0x%016llx", (unsigned long long)s[i]);
            if (i) file << ", ";
            file << '"' << hex << '"';
        }
        file << ']';
    }
    file << "\n}\n";
    file.close();
    if (!file) throw std::runtime_error("JSON 写入失败");
}
}  // namespace

int main(int argc, char** argv) {
    std::signal(SIGINT, handle_sigint);
    try {
        const Options o = options(argc, argv);
        GRBEnv env;
        env.set(GRB_IntParam_OutputFlag, 1);
        env.set(GRB_IntParam_LogToConsole, 1);
        GRBModel model(env);
        ActiveModelGuard guard(model);
        model.set(GRB_StringAttr_ModelName, o.linear ? "Gaston_linear_mask" : "Gaston_difference");
        model.set(GRB_IntParam_Threads, o.threads);
        model.set(GRB_DoubleParam_MIPGap, 0);
        model.set(GRB_DoubleParam_MIPGapAbs, 0);
        auto input = add_state(model, "input");
        auto output = add_state(model, "output");
        // E is explicit in the sparse forward circuit. Use that circuit for
        // either fixed side when restricting E; both express the same Y=L(X).
        if (o.count < 0) {
            if (o.linear) add_linear(model, input, output);
            else add_difference(model, input, output, o.effect_constraint);
        } else if (!o.linear && o.effect_constraint != EffectConstraint::Any)
            add_difference(model, input, output, o.effect_constraint);
        else if (o.input_fixed == o.linear)
            add_inverse(model, o.input_fixed ? input : output,
                        o.input_fixed ? output : input, o.linear);
        else if (o.linear) add_linear(model, input, output);
        else add_difference(model, input, output);
        auto ain = add_activity(model, input, "active_input");
        auto aout = add_activity(model, output, "active_output");
        GRBLinExpr nin = 0, nout = 0;
        for (int j = 0; j < N; ++j) { nin += ain[j]; nout += aout[j]; }
        if (!o.linear && o.effect_constraint == EffectConstraint::Nonzero) {
            // E!=0 implies X!=0 and, because lambda is invertible, Y!=0.
            // These valid redundant constraints make K=0 immediately infeasible.
            model.addConstr(nin >= 1, "theta_nonzero_requires_input");
            model.addConstr(nout >= 1, "theta_nonzero_requires_output");
        }
        if (o.count < 0) {
            GRBLinExpr input_bits = 0;
            for (const auto& row : input)
                for (const auto& v : row) input_bits += v;
            model.addConstr(input_bits >= 1, "input_nonzero");
            // Both lambda and lambda^T are invertible, so a nonzero input
            // also has a nonzero output. This is a valid redundant bound.
            model.addConstr(nout >= 1, "nonzero_input_requires_output");
        } else {
            model.addConstr((o.input_fixed ? nin : nout) == o.count, "fixed_count");
        }
        if (o.count < 0) model.setObjective(nin + nout, GRB_MINIMIZE);
        else model.setObjective(o.input_fixed ? nout : nin, GRB_MINIMIZE);
        // All operations commute with simultaneous cyclic rotation of all rows.
        // Every nonzero feasible pair therefore has a representative with column 0 active.
        if (o.symmetry) {
            if (o.count < 0) {
                GRBLinExpr first_column = 0;
                for (int i = 0; i < R; ++i) first_column += input[i][0];
                model.addConstr(first_column >= 1, "rotation_representative");
            } else if (o.count > 0) {
                model.addConstr((o.input_fixed ? ain[0] : aout[0]) == 1, "rotation_representative");
            }
        }
        auto& fixed = o.input_fixed ? input : output;
        for (auto& row : fixed)
            for (auto& v : row) v.set(GRB_IntAttr_BranchPriority, 10);
        if (!o.prefix.empty()) model.write(o.prefix + ".lp");
        if (o.export_only) return 0;
        const auto seed = initial_pair(o);
        if (seed) {
            if ((o.linear ? transpose_layer(seed->second) != seed->first
                          : linear_layer(seed->first) != seed->second)
                || (o.count < 0 ? bit_weight(seed->first) < 1
                               : column_weight(o.input_fixed ? seed->first : seed->second) != o.count)
                || (o.count < 0 && o.symmetry && !(active_columns(seed->first) & 1))
                || (!o.linear && !accepts_effect(o.effect_constraint, theta_effect(seed->first))))
                throw std::runtime_error("初始可行解构造失败");
            start_state(input, ain, seed->first);
            start_state(output, aout, seed->second);
        }
        if (o.count < 0) {
            std::printf("Mode=%s; input bit weight >= 1; minimize input + output active columns; "
                        "column index=LSB 0..63\n", o.linear ? "linear" : "difference");
            if (o.symmetry) std::printf("Symmetry: sum(input[row][0], row=0..4) >= 1\n");
        } else {
            std::printf("Mode=%s; fixed %s=%d; minimize %s; column index=LSB 0..63\n",
                        o.linear ? "linear" : "difference", o.input_fixed ? "input" : "output",
                        o.count, o.input_fixed ? "output" : "input");
        }
        if (!o.linear)
            std::printf("Theta-effect constraint=%s; E evaluated at A=rho_east(input)\n",
                        effect_constraint_name(o.effect_constraint));
        model.optimize();
        if (model.get(GRB_IntAttr_SolCount) == 0) return 0;

        const State x = read_state(input), y = read_state(output);
        if (o.linear ? transpose_layer(y) != x : linear_layer(x) != y)
            throw std::runtime_error("解未通过独立字级传播复核");
        if (!o.linear && !accepts_effect(o.effect_constraint, theta_effect(x)))
            throw std::runtime_error("解未通过 theta-effect 条件复核");
        verify_activity(x, ain);
        verify_activity(y, aout);
        if (o.count < 0 ? bit_weight(x) < 1
                       : column_weight(o.input_fixed ? x : y) != o.count)
            throw std::runtime_error("解未满足活跃数量或输入非零条件");
        if (o.count < 0 && o.symmetry && !(active_columns(x) & 1))
            throw std::runtime_error("解的输入第 0 列必须非零");
        const int status = model.get(GRB_IntAttr_Status);
        const double objective = model.get(GRB_DoubleAttr_ObjVal);
        if (std::abs(objective - activity_objective(o, x, y)) > 1e-5)
            throw std::runtime_error("解的目标值与实际活跃数量不一致");
        const double bound = model.get(GRB_DoubleAttr_ObjBound);
        const double gap = model.get(GRB_DoubleAttr_MIPGap);
        const double runtime = model.get(GRB_DoubleAttr_Runtime);
        std::printf("\nStatus: %d (%s)\nObjective: %.10g\nBound: %.10g\nMIPGap: %.10g\nRuntime: %.3f s\n",
                    status, status == GRB_OPTIMAL ? "已证明最优" : "当前可行解，尚未证明最优",
                    objective, bound, gap, runtime);
        print_state(o.linear ? "Input mask alpha" : "Input difference X", x);
        print_state(o.linear ? "Output mask beta" : "Output difference Y", y);
        std::printf("Active columns: input=%d; output=%d; total=%d\n",
                    column_weight(x), column_weight(y), column_weight(x) + column_weight(y));
        if (!o.linear) {
            const uint64_t effect = theta_effect(x);
            std::printf("Theta effect E=0x%016llx; weight=%d; in kernel=%s\n",
                        (unsigned long long)effect, __builtin_popcountll(effect), effect == 0 ? "yes" : "no");
        }
        std::printf("Verification: PASS (word-level propagation, activity counts, objective%s)\n",
                    o.linear ? "" : ", theta-effect condition");
        if (!o.prefix.empty()) {
            model.write(o.prefix + ".sol");
            write_json(o, status, objective, bound, gap, runtime, x, y);
            std::printf("Saved: %s.lp, %s.sol, %s.json\n", o.prefix.c_str(), o.prefix.c_str(), o.prefix.c_str());
        }
    } catch (const GRBException& e) {
        std::fprintf(stderr, "Gurobi 错误 %d: %s\n", e.getErrorCode(), e.getMessage().c_str());
        return 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "错误: %s\n", e.what());
        return 1;
    }
    return 0;
}

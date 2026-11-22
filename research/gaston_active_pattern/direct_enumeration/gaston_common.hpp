#pragma once
// Exact enumeration of Gaston active S-box columns. No solver dependency.
// Bit j has integer weight 2^j; ROL moves j to (j+s) mod 64.
#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr int ROWS = 5, COLS = 64, BITS = 320;
constexpr std::array<int, ROWS> EAST{0, 60, 22, 27, 4};
constexpr std::array<int, ROWS> WEST{0, 56, 31, 46, 43};
constexpr std::array<int, ROWS> TWIN{25, 32, 52, 60, 63};
using State = std::array<uint64_t, ROWS>;
using Matrix = std::array<std::bitset<BITS>, BITS>;
using Table = std::array<std::array<State, 32>, COLS>;
using Clock = std::chrono::steady_clock;
volatile std::sig_atomic_t stop_requested = 0;
void handle_sigint(int) { stop_requested = 1; }
double seconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}
uint64_t rol(uint64_t x, int s) {
    s = (s % COLS + COLS) % COLS;
    return s == 0 ? x : (x << s) | (x >> (COLS - s));
}
State xor_states(const State& x, const State& y) {
    return {x[0] ^ y[0], x[1] ^ y[1], x[2] ^ y[2], x[3] ^ y[3], x[4] ^ y[4]};
}
uint64_t active_columns(const State& x) { return x[0] | x[1] | x[2] | x[3] | x[4]; }
int weight(const State& x) { return __builtin_popcountll(active_columns(x)); }
int bit_weight(const State& x) {
    int result = 0;
    for (auto word : x) result += __builtin_popcountll(word);
    return result;
}
State linear_layer(const State& x) {
    State a{}, y{};
    uint64_t p = 0, q = 0;
    for (int i = 0; i < ROWS; ++i) {
        a[i] = rol(x[i], EAST[i]);
        p ^= a[i];
        q ^= rol(a[i], TWIN[i]);
    }
    const uint64_t effect = rol(p ^ rol(p, 1) ^ q ^ rol(q, 18), 23);
    for (int i = 0; i < ROWS; ++i) y[i] = rol(a[i] ^ effect, WEST[i]);
    return y;
}
State transpose_layer(const State& beta) {
    State d{}, alpha{};
    uint64_t p = 0;
    for (int i = 0; i < ROWS; ++i) {
        d[i] = rol(beta[i], -WEST[i]);
        p ^= d[i];
    }
    for (int i = 0; i < ROWS; ++i)
        alpha[i] = rol(d[i] ^ rol(p, -23) ^ rol(p, -24)
                      ^ rol(p, -23 - TWIN[i]) ^ rol(p, -41 - TWIN[i]), -EAST[i]);
    return alpha;
}
Matrix forward_matrix() {
    Matrix result{};
    for (int c = 0; c < BITS; ++c) {
        State unit{};
        unit[c / COLS] = uint64_t(1) << (c % COLS);
        const State image = linear_layer(unit);
        for (int r = 0; r < BITS; ++r)
            result[r][c] = (image[r / COLS] >> (r % COLS)) & 1;
    }
    return result;
}
Matrix transpose(const Matrix& matrix) {
    Matrix result{};
    for (int r = 0; r < BITS; ++r)
        for (int c = 0; c < BITS; ++c) result[r][c] = matrix[c][r];
    return result;
}
Matrix inverse(Matrix matrix) {
    Matrix result{};
    for (int i = 0; i < BITS; ++i) result[i][i] = 1;
    for (int c = 0; c < BITS; ++c) {
        int pivot = c;
        while (pivot < BITS && !matrix[pivot][c]) ++pivot;
        if (pivot == BITS) throw std::runtime_error("线性层矩阵不可逆");
        std::swap(matrix[c], matrix[pivot]);
        std::swap(result[c], result[pivot]);
        for (int r = 0; r < BITS; ++r) if (r != c && matrix[r][c]) {
            matrix[r] ^= matrix[c];
            result[r] ^= result[c];
        }
    }
    return result;
}

struct Options {
    bool linear = false, input_fixed = false, symmetry = true;
    int count = -1;
    std::string prefix;
};
const char* mode_name(const Options& o) { return o.linear ? "linear" : "difference"; }
const char* fixed_name(const Options& o) { return o.input_fixed ? "input" : "output"; }
const char* mapping_name(const Options& o) {
    if (o.linear) return o.input_fixed ? "lambda^(-T)" : "lambda^T";
    return o.input_fixed ? "lambda" : "lambda^(-1)";
}
void usage(const char* program) {
    std::printf(
        "用法: %s [--mode difference|linear] (--input-active K | --output-active K) [选项]\n"
        "  K 指恰好 K 个非零的 5 比特列，目标是最小化另一端的非零列数。\n"
        "  --mode difference  差分 Y=lambda(X)，默认模式\n"
        "  --mode linear      掩码 alpha=lambda^T(beta)\n"
        "  --no-symmetry      枚举全部位置组合；默认固定端第 0 列活跃\n"
        "  --prefix PATH      保存 PATH.json；自动创建父目录\n"
        "  --help             显示帮助\n"
        "K 范围 0..64。适合小 K；候选数随 K 迅速增长。Ctrl+C 可报告当前最优候选。\n",
        program);
}
[[maybe_unused]] Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { usage(argv[0]); std::exit(0); }
        if (arg == "--no-symmetry") { o.symmetry = false; continue; }
        if (arg != "--mode" && arg != "--input-active" && arg != "--output-active" && arg != "--prefix")
            throw std::invalid_argument("未知参数: " + arg);
        if (++i == argc) throw std::invalid_argument(arg + " 缺少参数");
        const std::string value = argv[i];
        if (arg == "--mode") {
            if (value != "difference" && value != "linear")
                throw std::invalid_argument("--mode 必须是 difference 或 linear");
            o.linear = value == "linear";
        } else if (arg == "--prefix") {
            if (value.empty()) throw std::invalid_argument("--prefix 不能为空");
            o.prefix = value;
        } else {
            if (o.count >= 0) throw std::invalid_argument("只能指定一个活跃数量约束");
            std::size_t used = 0;
            const int count = std::stoi(value, &used);
            if (used != value.size() || count < 0 || count > COLS)
                throw std::invalid_argument("活跃列数量必须在 0..64 内");
            o.count = count;
            o.input_fixed = arg == "--input-active";
        }
    }
    if (o.count < 0) throw std::invalid_argument("需要 --input-active K 或 --output-active K");
    return o;
}
long double choose(int n, int k) {
    k = std::min(k, n - k);
    long double result = 1;
    for (int i = 1; i <= k; ++i) result = result * (n - k + i) / i;
    return result;
}
long double planned_candidates(const Options& o) {
    if (o.count == 0) return 1;
    long double result = o.symmetry ? choose(COLS - 1, o.count - 1) : choose(COLS, o.count);
    for (int i = 0; i < o.count; ++i) result *= 31;
    return result;
}
Table make_table(const Options& o) {
    Matrix matrix = forward_matrix();
    if (o.linear) matrix = transpose(matrix);
    if (o.input_fixed == o.linear) matrix = inverse(matrix);
    std::array<State, BITS> basis{};
    for (int c = 0; c < BITS; ++c) {
        for (int r = 0; r < BITS; ++r)
            if (matrix[r][c]) basis[c][r / COLS] |= uint64_t(1) << (r % COLS);
        State unit{};
        unit[c / COLS] = uint64_t(1) << (c % COLS);
        const State x = o.input_fixed ? unit : basis[c];
        const State y = o.input_fixed ? basis[c] : unit;
        if (o.linear ? transpose_layer(y) != x : linear_layer(x) != y)
            throw std::runtime_error("矩阵未通过全部 320 个基向量的字级复核");
    }
    Table table{};
    for (int j = 0; j < COLS; ++j)
        for (unsigned value = 1; value < 32; ++value) {
            const unsigned row = unsigned(__builtin_ctz(value));
            table[j][value] = xor_states(table[j][value & (value - 1)], basis[row * COLS + j]);
        }
    return table;
}

struct Result {
    uint64_t visited = 0;
    int minimum = COLS + 1;
    State input{}, output{};
    bool complete = false;
};
void print_state(const char* name, const State& state) {
    std::printf("%s: active_columns=%d, active_bits=%d\n  columns (index:value):",
                name, weight(state), bit_weight(state));
    for (int j = 0; j < COLS; ++j) {
        unsigned value = 0;
        for (int i = 0; i < ROWS; ++i) value |= unsigned((state[i] >> j) & 1) << i;
        if (value) std::printf(" %d:0x%02x", j, value);
    }
    std::printf("\n");
    for (int i = 0; i < ROWS; ++i)
        std::printf("  row %d = 0x%016llx\n", i, static_cast<unsigned long long>(state[i]));
}
struct ExecutionInfo {
    std::string backend, device;
    int threads;
    unsigned position_sets, completed_sets, cpu_sets, gpu_sets;
    uint64_t planned;
};
[[maybe_unused]] void save_result(const Options& o, const Result& r, double setup, double search,
                 const ExecutionInfo* parallel = nullptr) {
    if (o.prefix.empty()) return;
    const std::filesystem::path prefix(o.prefix);
    if (!prefix.parent_path().empty()) std::filesystem::create_directories(prefix.parent_path());
    std::ofstream file(o.prefix + ".json");
    if (!file) throw std::runtime_error("无法写入结果 JSON");
    file.precision(18);
    file << "{\n  \"method\": \"exact_enumeration\",\n  \"mode\": \"" << mode_name(o)
         << "\",\n  \"mapping\": \"" << mapping_name(o)
         << "\",\n  \"fixed_metric\": \"active_sbox_columns\",\n  \"fixed_side\": \"" << fixed_name(o)
         << "\",\n  \"fixed_count\": " << o.count
         << ",\n  \"symmetry\": " << (o.symmetry ? "true" : "false")
         << ",\n  \"candidate_count_estimate\": "
         << (parallel ? static_cast<long double>(parallel->planned) : planned_candidates(o))
         << ",\n  \"visited\": " << r.visited
         << ",\n  \"enumeration_complete\": " << (r.complete ? "true" : "false")
         << ",\n  \"proven_optimal\": " << (r.complete ? "true" : "false")
         << ",\n  \"objective\": " << r.minimum
         << ",\n  \"setup_seconds\": " << setup
         << ",\n  \"search_seconds\": " << search
         << ",\n  \"runtime_seconds\": " << setup + search
         << ",\n  \"verified\": true,\n  \"bit_index\": \"LSB=0; column value=sum(row_bit<<row)\"";
    for (int side = 0; side < 2; ++side) {
        const State& state = side == 0 ? r.input : r.output;
        const char* name = side == 0 ? "input" : "output";
        file << ",\n  \"" << name << "_active\": " << weight(state)
             << ",\n  \"" << name << "_bit_weight\": " << bit_weight(state)
             << ",\n  \"" << name << "_words\": [";
        for (int i = 0; i < ROWS; ++i) {
            char hex[19];
            std::snprintf(hex, sizeof(hex), "0x%016llx", static_cast<unsigned long long>(state[i]));
            if (i) file << ", ";
            file << '"' << hex << '"';
        }
        file << "],\n  \"" << name << "_columns\": [";
        bool first = true;
        for (int j = 0; j < COLS; ++j) if ((active_columns(state) >> j) & 1) {
            if (!first) file << ", ";
            file << j;
            first = false;
        }
        file << ']';
    }
    if (parallel) {
        file << ",\n  \"backend\": \"" << parallel->backend << '"'
             << ",\n  \"device\": \"" << parallel->device << '"'
             << ",\n  \"threads\": " << parallel->threads
             << ",\n  \"symmetry_method\": \"" << (o.symmetry ? "position_rotation_orbits" : "none") << '"'
             << ",\n  \"position_sets\": " << parallel->position_sets
             << ",\n  \"completed_position_sets\": " << parallel->completed_sets
             << ",\n  \"cpu_position_sets\": " << parallel->cpu_sets
             << ",\n  \"gpu_position_sets\": " << parallel->gpu_sets;
    }
    file << "\n}\n";
    file.close();
    if (!file) throw std::runtime_error("JSON 写入失败");
}
}  // namespace

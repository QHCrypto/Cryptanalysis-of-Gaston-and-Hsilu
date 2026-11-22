#include "gaston_common.hpp"

namespace {
class Enumerator {
public:
    Enumerator(const Options& options, const Table& images)
        : o(options), table(images), positions(o.count), values(o.count) {}
    Result run() {
        if (o.symmetry && o.count > 0) { positions[0] = 0; choose_positions(1, 1); }
        else choose_positions(0, 0);
        result.complete = !stop_requested;
        return result;
    }
private:
    const Options& o;
    const Table& table;
    std::vector<int> positions, values;
    Result result;
    void choose_positions(int depth, int next) {
        if (stop_requested) return;
        if (depth == o.count) { enumerate_values(0, State{}); return; }
        for (int j = next; j <= COLS - (o.count - depth) && !stop_requested; ++j) {
            positions[depth] = j;
            choose_positions(depth + 1, j + 1);
        }
    }
    void enumerate_values(int depth, const State& partial) {
        if (stop_requested) return;
        if (depth < o.count) {
            const auto& images = table[positions[depth]];
            for (int value = 1; value < 32 && !stop_requested; ++value) {
                values[depth] = value;
                enumerate_values(depth + 1, xor_states(partial, images[value]));
            }
            return;
        }
        if (result.visited == std::numeric_limits<uint64_t>::max())
            throw std::overflow_error("枚举计数超过 uint64_t 范围");
        ++result.visited;
        const int objective = weight(partial);
        if (objective >= result.minimum) return;
        result.minimum = objective;
        State fixed{};
        for (int d = 0; d < o.count; ++d)
            for (int row = 0; row < ROWS; ++row)
                if ((values[d] >> row) & 1) fixed[row] |= uint64_t(1) << positions[d];
        result.input = o.input_fixed ? fixed : partial;
        result.output = o.input_fixed ? partial : fixed;
        std::printf("Improved: minimum=%d, visited=%llu\n", objective,
                    static_cast<unsigned long long>(result.visited));
        std::fflush(stdout);
    }
};

}  // namespace

int main(int argc, char** argv) {
    std::signal(SIGINT, handle_sigint);
    try {
        const Options o = parse_options(argc, argv);
        const auto begin = Clock::now();
        std::printf("Mode=%s; fixed %s active columns=%d; map=%s\n",
                    mode_name(o), fixed_name(o), o.count, mapping_name(o));
        std::printf("Candidates planned: %.12Lg; symmetry=%s; all 31 nonzero values per active column\n",
                    planned_candidates(o), o.symmetry ? "column 0 active" : "off");
        std::fflush(stdout);
        const Table table = make_table(o);
        const auto ready = Clock::now();
        const Result result = Enumerator(o, table).run();
        const auto end = Clock::now();
        if (result.minimum > COLS) {
            std::printf("Interrupted before the first candidate; no solution report.\n");
            return 130;
        }
        if ((o.linear ? transpose_layer(result.output) != result.input
                      : linear_layer(result.input) != result.output)
            || weight(o.input_fixed ? result.input : result.output) != o.count
            || weight(o.input_fixed ? result.output : result.input) != result.minimum)
            throw std::runtime_error("结果未通过字级传播及活跃数量复核");
        const double setup = seconds(begin, ready), search = seconds(ready, end);
        std::printf("\nStatus: %s\nMinimum %s active columns: %d\nVisited: %llu\n"
                    "Setup: %.6f s\nEnumeration: %.6f s\nTotal: %.6f s\n",
                    result.complete ? "已完成枚举，已证明最优" : "已中断，仅为当前最优候选",
                    o.input_fixed ? "output" : "input", result.minimum,
                    static_cast<unsigned long long>(result.visited), setup, search, setup + search);
        print_state(o.linear ? "Input mask alpha" : "Input difference X", result.input);
        print_state(o.linear ? "Output mask beta" : "Output difference Y", result.output);
        std::printf("Verification: PASS (320 basis vectors and final word-level propagation)\n");
        save_result(o, result, setup, search);
        if (!o.prefix.empty()) std::printf("Saved: %s.json\n", o.prefix.c_str());
        return result.complete ? 0 : 130;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "错误: %s\n", error.what());
        return 1;
    }
}

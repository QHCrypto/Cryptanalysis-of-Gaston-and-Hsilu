#include "experiment.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <tuple>

namespace fs = std::filesystem;
using namespace gdl;
using Options = std::map<std::string, std::string>;

static std::string get(const Options &o, const std::string &key, const std::string &fallback) {
    auto it = o.find(key);
    return it == o.end() ? fallback : it->second;
}
static uint64_t integer(const std::string &s) {
    if (s.empty() || s[0] == '-') throw std::runtime_error("expected nonnegative integer: " + s);
    size_t end = 0;
    uint64_t v = std::stoull(s, &end, 0);
    if (end != s.size()) throw std::runtime_error("invalid integer: " + s);
    return v;
}
static int bounded(const std::string &s, int lo, int hi) {
    const uint64_t v = integer(s);
    if (v < uint64_t(lo) || v > uint64_t(hi)) throw std::runtime_error("integer outside allowed range: " + s);
    return int(v);
}
static std::vector<std::string> split(const std::string &s) {
    std::vector<std::string> out;
    std::stringstream in(s);
    std::string part;
    while (std::getline(in, part, ',')) out.push_back(part);
    return out;
}
static std::vector<int> int_list(const std::string &s, int lo, int hi) {
    std::set<int> values;
    for (const auto &x : split(s)) values.insert(bounded(x, lo, hi));
    if (values.empty()) throw std::runtime_error("empty integer list");
    return {values.begin(), values.end()};
}
static Options parse(int argc, char **argv, const std::set<std::string> &allowed) {
    Options o;
    for (int i = 2; i < argc; i += 2) {
        std::string key = argv[i];
        if (!allowed.count(key)) throw std::runtime_error("unknown option: " + key);
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + key);
        if (!o.emplace(key, argv[i + 1]).second) throw std::runtime_error("duplicate option: " + key);
    }
    return o;
}
static std::string json(const std::string &s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c == '\t') out += "\\t";
        else out += c;
    }
    return out + '"';
}
static std::ofstream output_file(const fs::path &p) {
    std::ofstream f(p);
    if (!f) throw std::runtime_error("cannot write " + p.string());
    f.exceptions(std::ios::badbit | std::ios::failbit);
    f << std::setprecision(std::numeric_limits<Real>::max_digits10);
    return f;
}
static void prepare(const fs::path &dir, std::initializer_list<std::string> files) {
    for (const auto &f : files)
        if (fs::exists(dir / f)) throw std::runtime_error("refusing to overwrite " + (dir / f).string() + "; choose a new --output directory");
    fs::create_directories(dir);
}

struct Candidate {
    int rounds = 0, rank = 0, dc = 0, dp = 0, mc = 0, mp = 0;
    Real rho = 0;
    bool underflow = false;
};
static bool better(const Candidate &a, const Candidate &b) {
    if (std::abs(a.rho) != std::abs(b.rho)) return std::abs(a.rho) > std::abs(b.rho);
    return std::tie(a.dc, a.dp, a.mc, a.mp) < std::tie(b.dc, b.dp, b.mc, b.mp);
}
static void log_value(std::ostream &f, Real rho) {
    if (rho != 0) f << std::log2(std::abs(rho)); // blank means no finite logarithm
}
static void candidate_header(std::ostream &f) {
    f << "rounds,rank,input_col,input_pattern,output_col,output_pattern,model_rho,model_log2_abs,status,input_words,output_words\n";
}
static void candidate_row(std::ostream &f, const Candidate &c) {
    f << c.rounds << ',' << c.rank << ',' << c.dc << ',' << c.dp << ',' << c.mc << ',' << c.mp << ',' << c.rho << ',';
    log_value(f, c.rho);
    f << ',' << (c.underflow ? "underflow_detected" : c.rho == 0 ? "numerical_zero" : "finite")
      << ',' << hex_state(embed(c.dc, c.dp)) << ',' << hex_state(embed(c.mc, c.mp)) << '\n';
}

static void search(const Options &o) {
    const std::string scope = get(o, "--scope", "bit");
    if (scope != "bit" && scope != "box") throw std::runtime_error("--scope must be bit or box");
    const auto rounds = int_list(get(o, "--rounds", "2,3"), 1, 12);
    const int top = bounded(get(o, "--top", "16"), 1, 10000);
    const std::string symmetry_opt = get(o, "--symmetry", "on");
    if (symmetry_opt != "on" && symmetry_opt != "off") throw std::runtime_error("--symmetry must be on or off");
    const bool symmetry = symmetry_opt == "on";
    const fs::path dir = get(o, "--output", "results/search");
    prepare(dir, {"scores.csv", "top.csv", "search.json"});
    auto scores = output_file(dir / "scores.csv");
    candidate_header(scores);
    const auto began = std::chrono::steady_clock::now();
    Estimator estimator;
    std::vector<int> patterns;
    for (int a = 1; a < 32; ++a) if (scope == "box" || __builtin_popcount(unsigned(a)) == 1) patterns.push_back(a);
    std::map<int, std::vector<Candidate>> best;
    std::map<int, uint64_t> nonzero;
    uint64_t evaluated = 0;
    bool underflow = false;
    for (int dc = 0; dc < (symmetry ? 1 : 64); ++dc) {
        for (int a : patterns) {
            const auto e = estimator.estimate(dc, a, rounds.back());
            underflow |= e.underflow;
            for (int r : rounds) {
                for (int mc = 0; mc < 64; ++mc)
                    for (int b : patterns) {
                        Candidate c{r, 0, dc, a, mc, b, e.rounds[r][mc][b], e.underflow};
                        if (!std::isfinite(c.rho) || std::abs(c.rho) > 1 + 1e-15L)
                            throw std::runtime_error("invalid correlation produced by estimator");
                        candidate_row(scores, c);
                        nonzero[r] += c.rho != 0;
                        best[r].push_back(c);
                        ++evaluated;
                    }
                std::sort(best[r].begin(), best[r].end(), better);
                if (best[r].size() > size_t(top)) best[r].resize(top);
            }
            std::cerr << "searched input (column=" << dc << ", pattern=" << a << ")\n";
        }
    }
    auto top_file = output_file(dir / "top.csv");
    candidate_header(top_file);
    for (auto &[r, list] : best) {
        if (list.front().rho == 0) {
            std::cout << "R=" << r << " no nonzero candidate in this model; no endpoint selected\n";
            continue;
        }
        int rank = 0;
        for (auto &c : list) { if (c.rho == 0) break; c.rank = ++rank; candidate_row(top_file, c); }
        const auto &c = list.front();
        std::cout << "R=" << r << " delta=(" << c.dc << ',' << c.dp << ") mask=(" << c.mc << ',' << c.mp
                  << ") rho=" << std::setprecision(18) << c.rho << " log2_abs=";
        log_value(std::cout, c.rho);
        std::cout << '\n';
    }
    auto meta = output_file(dir / "search.json");
    const uint64_t endpoints = 64 * patterns.size();
    meta << "{\n  \"model\": \"HNW-zero-value-mask-rd1-rl0\",\n"
         << "  \"scope\": " << json(scope) << ",\n  \"rounds\": [";
    for (size_t i = 0; i < rounds.size(); ++i) { if (i) meta << ','; meta << rounds[i]; }
    meta << "],\n  \"omit_first_linear\": true,\n  \"input_boundary\": \"input of first chi\",\n"
         << "  \"output_boundary\": \"output of last chi\",\n  \"rd\": 1,\n  \"rl\": 0,\n"
         << "  \"constants_in_model\": \"ignored by zero-value-mask assumption\",\n"
         << "  \"rotation_reduction\": " << (symmetry ? "true" : "false") << ",\n"
         << "  \"orbit_size\": " << (symmetry ? 64 : 1) << ",\n"
         << "  \"endpoint_pairs_per_round\": " << endpoints * endpoints << ",\n"
         << "  \"evaluated_representatives_per_round\": " << evaluated / rounds.size() << ",\n"
         << "  \"top_per_round\": " << top << ",\n  \"numeric_type\": \"long double\",\n"
         << "  \"mantissa_bits\": " << std::numeric_limits<Real>::digits << ",\n"
         << "  \"underflow_detected\": " << (underflow ? "true" : "false") << ",\n"
         << "  \"elapsed_seconds\": " << std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count() << ",\n"
         << "  \"method_reference\": \"https://eprint.iacr.org/2026/358\",\n"
         << "  \"reference_code_commit\": \"783e8edc560b848d0f2ba63f11838aecd23da40f\",\n"
         << "  \"round_summaries\": [";
    for (size_t i = 0; i < rounds.size(); ++i) {
        if (i) meta << ',';
        const int r = rounds[i];
        meta << "{\"rounds\":" << r << ",\"nonzero_representatives\":" << nonzero[r]
             << ",\"status\":" << json(nonzero[r] ? "candidate_found" : "no_nonzero_model_candidate") << '}';
    }
    meta << "]\n}\n";
    if (underflow) std::cerr << "WARNING: numerical underflow occurred; zero estimates are not certified zeros.\n";
}

static std::vector<Candidate> read_candidates(const fs::path &p, int max_rank) {
    std::ifstream f(p);
    if (!f) throw std::runtime_error("cannot read candidates: " + p.string());
    std::string line;
    std::getline(f, line);
    const auto header = split(line);
    const std::vector<std::string> expected = {"rounds", "rank", "input_col", "input_pattern", "output_col", "output_pattern", "model_rho", "model_log2_abs", "status", "input_words", "output_words"};
    if (header != expected) throw std::runtime_error("candidate CSV header mismatch");
    std::vector<Candidate> out;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        const auto a = split(line);
        if (a.size() != expected.size()) throw std::runtime_error("malformed candidate CSV row");
        Candidate c;
        c.rounds = bounded(a[0], 1, 12); c.rank = bounded(a[1], 0, 10000000);
        if (c.rank == 0 || c.rank > max_rank) continue;
        c.dc = bounded(a[2], 0, 63); c.dp = bounded(a[3], 1, 31);
        c.mc = bounded(a[4], 0, 63); c.mp = bounded(a[5], 1, 31);
        c.rho = std::stold(a[6]);
        if (!std::isfinite(c.rho) || std::abs(c.rho) > 1) throw std::runtime_error("invalid candidate correlation");
        c.underflow = a[8] == "underflow_detected";
        out.push_back(c);
    }
    if (out.empty()) throw std::runtime_error("no ranked candidates selected");
    return out;
}

static void verify(const Options &o) {
    const auto candidates_path = get(o, "--candidates", "results/search/top.csv");
    const auto candidates = read_candidates(candidates_path, bounded(get(o, "--top-per-round", "1"), 1, 10000));
    const auto rotations = int_list(get(o, "--rotations", "0"), 0, 63);
    const int threads = bounded(get(o, "--threads", "4"), 1, 256);
    const uint64_t samples = integer(get(o, "--samples", "4194304"));
    if (!samples || samples > (uint64_t(1) << 60)) throw std::runtime_error("invalid sample count");
    const uint64_t seed = integer(get(o, "--seed", "20261015"));
    const std::string constants_opt = get(o, "--constants", "reference");
    if (constants_opt != "reference" && constants_opt != "zero") throw std::runtime_error("--constants must be reference or zero");
    const bool constants = constants_opt == "reference";
    const std::string start_opt = get(o, "--start-round", "last");
    for (const auto &c : candidates) {
        int start = start_opt == "last" ? 12 - c.rounds : bounded(start_opt, 0, 11);
        if (start + c.rounds > 12) throw std::runtime_error("round interval exceeds the reference's 12 rounds");
    }
    const fs::path dir = get(o, "--output", "results/verification");
    prepare(dir, {"experiments.csv", "verification.json"});
    auto f = output_file(dir / "experiments.csv");
    f << "rounds,rank,rotation,input_col,input_pattern,output_col,output_pattern,model_rho,model_log2_abs,samples,n0,n1,experimental_rho,experimental_log2_abs,ci95_low,ci95_high,significant,model_in_ci95,start_round,constants,seed,seconds,input_words,output_words\n";
    const auto began = std::chrono::steady_clock::now();
    for (const auto &c : candidates) {
        for (int rotation : rotations) {
            const int dc = (c.dc + rotation) % 64, mc = (c.mc + rotation) % 64;
            const int start = start_opt == "last" ? 12 - c.rounds : int(integer(start_opt));
            const uint64_t run_seed = splitmix_final(seed ^ (uint64_t(c.rounds) << 40) ^ (uint64_t(dc) << 32)
                                                   ^ (uint64_t(c.dp) << 24) ^ (uint64_t(mc) << 8) ^ uint64_t(c.mp));
            const auto t0 = std::chrono::steady_clock::now();
            std::cerr << "sampling R=" << c.rounds << " rank=" << c.rank << " rotation=" << rotation << " N=" << samples << '\n';
            const auto m = measure(dc, c.dp, mc, c.mp, c.rounds, start, constants, samples, run_seed, threads);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            f << c.rounds << ',' << c.rank << ',' << rotation << ',' << dc << ',' << c.dp << ',' << mc << ',' << c.mp << ',' << c.rho << ',';
            log_value(f, c.rho);
            f << ',' << samples << ',' << m.n0 << ',' << m.n1 << ',' << m.correlation << ',';
            log_value(f, m.correlation);
            const bool significant = m.low > 0 || m.high < 0;
            const bool in_ci = c.rho >= m.low && c.rho <= m.high;
            f << ',' << m.low << ',' << m.high << ',' << significant << ',' << in_ci << ',' << start
              << ',' << constants_opt << ',' << run_seed << ',' << seconds << ',' << hex_state(embed(dc, c.dp)) << ',' << hex_state(embed(mc, c.mp)) << '\n';
            f.flush();
            std::cout << "R=" << c.rounds << " rank=" << c.rank << " rotation=" << rotation
                      << " rho=" << std::setprecision(12) << m.correlation << " CI95=[" << m.low << ',' << m.high
                      << "] model_in_ci=" << in_ci << '\n';
        }
    }
    auto meta = output_file(dir / "verification.json");
    meta << "{\n  \"candidates\": " << json(fs::absolute(candidates_path).string())
         << ",\n  \"samples_per_candidate\": " << samples << ",\n  \"base_seed\": " << json(std::to_string(seed))
         << ",\n  \"threads\": " << threads << ",\n  \"start_round\": " << json(start_opt)
         << ",\n  \"constants\": " << json(constants_opt) << ",\n  \"rotations\": [";
    for (size_t i = 0; i < rotations.size(); ++i) { if (i) meta << ','; meta << rotations[i]; }
    meta << "],\n  \"generator\": \"counter-indexed SplitMix64; 5 outputs per state\",\n"
         << "  \"interval\": \"two-sided 95% Wilson, rho=2p-1; pointwise, not simultaneous\",\n"
         << "  \"first_chi_constant\": \"absorbed into the uniform input\",\n"
         << "  \"elapsed_seconds\": " << std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count() << "\n}\n";
}

int main(int argc, char **argv) {
    try {
        if (argc < 2 || std::string(argv[1]) == "--help") {
            std::cout << "Gaston differential-linear endpoint search (first linear layer omitted)\n"
                         "search --scope bit|box --rounds 2,3 --top 16 --symmetry on|off --output DIR\n"
                         "verify --candidates top.csv --top-per-round 1 --rotations 0 --samples 4194304\n"
                         "       --seed 20261015 --threads 4 --start-round last|NUMBER\n"
                         "       --constants reference|zero --output DIR\n"
                         "Patterns: bit i of the 5-bit integer selects lane i; column 0 is the uint64 LSB.\n"
                         "Search is model-based. Verification uses the real fixed-constant reduced permutation.\n";
            return 0;
        }
        const std::string command = argv[1];
        if (command == "search") search(parse(argc, argv, {"--scope", "--rounds", "--top", "--symmetry", "--output"}));
        else if (command == "verify") verify(parse(argc, argv, {"--candidates", "--top-per-round", "--rotations", "--samples", "--seed", "--threads", "--start-round", "--constants", "--output"}));
        else throw std::runtime_error("unknown command: " + command);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}

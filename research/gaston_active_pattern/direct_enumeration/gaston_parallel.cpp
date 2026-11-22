// Exact CPU/GPU enumeration, with one representative per position rotation orbit.
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <omp.h>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <numeric>
#include <thread>
#include "gaston_common.hpp"

namespace {
using Positions = std::array<int, 5>;
constexpr uint64_t NONE = UINT64_MAX;
struct ParallelOptions {
    Options common;
    std::string backend = "cpu";
    int threads = omp_get_max_threads();
    unsigned batch = 64;
    unsigned limit = 0; // Explicit benchmark subset; never proves optimality.
};
ParallelOptions parse_parallel(int argc, char** argv) {
    ParallelOptions p;
    std::vector<char*> serial_args{argv[0]};
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            usage(argv[0]);
            std::printf("并行版支持 K=0..5，默认按位置集合的循环移位轨道去重。\n"
                        "  --backend cpu|opencl|hybrid  CPU / GPU / 协同搜索，默认 cpu\n"
                        "  --threads N                 CPU 工作线程数，默认可用逻辑核数\n"
                        "  --batch-positions N         每次 GPU 批处理的位置数，默认 64\n"
                        "  --limit-positions N         仅测试前 N 个位置集合，不证明全局最优\n");
            std::exit(0);
        }
        if (a == "--backend" || a == "--threads" || a == "--batch-positions" || a == "--limit-positions") {
            if (++i == argc) throw std::invalid_argument(a + " 缺少参数");
            std::string v = argv[i];
            if (a == "--backend") {
                if (v != "cpu" && v != "opencl" && v != "hybrid")
                    throw std::invalid_argument("backend 必须是 cpu/opencl/hybrid");
                p.backend = v;
            } else {
                std::size_t used = 0;
                int n = std::stoi(v, &used);
                if (used != v.size() || n < 1) throw std::invalid_argument(a + " 必须是正整数");
                if (a == "--threads") p.threads = n;
                if (a == "--batch-positions") p.batch = n;
                if (a == "--limit-positions") p.limit = n;
            }
        } else serial_args.push_back(argv[i]);
    }
    p.common = parse_options(int(serial_args.size()), serial_args.data());
    if (p.common.count > 5) throw std::invalid_argument("并行版支持 K=0..5；更大 K 可使用原串行程序");
    if (p.batch > 65536 || p.threads > 1024) throw std::invalid_argument("批量或线程数过大");
    return p;
}
uint64_t power31(int k) {
    uint64_t n = 1;
    while (k--) n *= 31;
    return n;
}
uint64_t position_mask(const Positions& p, int k) {
    uint64_t m = 0;
    for (int i = 0; i < k; ++i) m |= uint64_t(1) << p[i];
    return m;
}
std::vector<Positions> make_positions(const Options& o) {
    std::vector<Positions> positions;
    Positions p{};
    std::function<void(int,int)> visit = [&](int depth, int next) {
        if (depth == o.count) {
            if (o.symmetry && o.count > 0) {
                const uint64_t m = position_mask(p, o.count);
                // Any minimum nonempty rotation has bit 0 set. Only rotations
                // that move an active position to 0 need to be compared.
                for (int i = 1; i < o.count; ++i)
                    if (rol(m, -p[i]) < m) return;
            }
            positions.push_back(p);
            return;
        }
        for (int j = next; j <= COLS - (o.count - depth); ++j) {
            p[depth] = j;
            visit(depth + 1, j + 1);
        }
    };
    if (o.symmetry && o.count > 0) { p[0] = 0; visit(1, 1); }
    else visit(0, 0);
    // Independent Burnside count includes periodic position sets (K=2,4).
    uint64_t expected = 0;
    if (o.symmetry) {
        for (int shift = 0; shift < 64; ++shift) {
            int cycles = std::gcd(64, shift), length = 64 / cycles;
            if (o.count % length == 0)
                expected += uint64_t(choose(cycles, o.count / length) + 0.5L);
        }
        expected /= 64;
    } else expected = uint64_t(choose(64, o.count) + 0.5L);
    if (positions.size() != expected) throw std::runtime_error("位置轨道计数未通过 Burnside 复核");
    return positions;
}
void check_rotation(const Table& t) {
    for (int col = 0; col < 64; ++col)
        for (int v = 1; v < 32; ++v)
            for (int row = 0; row < 5; ++row)
                if (t[col][v][row] != rol(t[0][v][row], col))
                    throw std::runtime_error("映射未通过循环移位可交换性复核");
}

// The two last values form a 961-entry table in the worker's private L1 cache.
// No atomics, signals, allocation or recursion in the candidate loop.
uint64_t cpu_position(const Table& table, const Positions& p, int k) {
    if (k == 0) return 0;
    uint64_t best = NONE;
    if (k == 1) {
        for (unsigned v = 1; v < 32; ++v)
            best = std::min(best, (uint64_t(weight(table[p[0]][v])) << 32) | (v - 1));
        return best;
    }
    alignas(64) uint64_t tail[5][961];
    for (unsigned b = 0; b < 31; ++b)
        for (unsigned a = 0; a < 31; ++a)
            for (int r = 0; r < 5; ++r)
                tail[r][31*b+a] = table[p[k-2]][a+1][r] ^ table[p[k-1]][b+1][r];
    const unsigned prefixes = unsigned(power31(k - 2));
    int minimum = 65;
    unsigned best_code = 0;
    for (unsigned prefix = 0; prefix < prefixes; ++prefix) {
        State s{};
        unsigned digits = prefix;
        for (int i = 0; i < k - 2; ++i) {
            s = xor_states(s, table[p[i]][digits % 31 + 1]);
            digits /= 31;
        }
        for (unsigned j = 0; j < 961; ++j) {
            const uint64_t active = (s[0] ^ tail[0][j]) | (s[1] ^ tail[1][j])
                                  | (s[2] ^ tail[2][j]) | (s[3] ^ tail[3][j])
                                  | (s[4] ^ tail[4][j]);
            const int w = __builtin_popcountll(active);
            const unsigned code = prefix + prefixes * j;
            if (w < minimum || (w == minimum && code < best_code)) {
                minimum = w;
                best_code = code;
            }
        }
    }
    return (uint64_t(minimum) << 32) | best_code;
}

const char* kernel_source = R"CLC(
// K, PREFIXES, GROUPS and WG are build-time constants.
__kernel void enumerate(__global const ulong* table,
                        __global const int* positions,
                        const uint start,
                        __global ulong* result) {
    const uint lane = get_local_id(0);
    const uint group = get_group_id(0);
    const uint pos = start + group / GROUPS;
    const uint part = group % GROUPS;
    __local ulong images[5 * 32 * 5];
    __local ulong reduction[WG];
    for (uint i = lane; i < K * 32 * 5; i += WG) {
        uint d = i / (32 * 5);
        uint offset = i % (32 * 5);
        images[i] = table[positions[pos*5+d] * 32 * 5 + offset];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    ulong best = (ulong)(-1);
    for (uint prefix = part*WG+lane; prefix < PREFIXES; prefix += GROUPS*WG) {
        ulong s0=0, s1=0, s2=0, s3=0, s4=0;
        uint digits = prefix;
        #pragma unroll
        for (uint d=0; d<K-2; ++d) {
            uint i = (d*32 + digits%31+1)*5;
            digits /= 31;
            s0 ^= images[i]; s1 ^= images[i+1]; s2 ^= images[i+2];
            s3 ^= images[i+3]; s4 ^= images[i+4];
        }
        for (uint b=0; b<31; ++b) {
            uint ib = ((K-1)*32+b+1)*5;
            ulong t0=s0^images[ib], t1=s1^images[ib+1], t2=s2^images[ib+2];
            ulong t3=s3^images[ib+3], t4=s4^images[ib+4];
            for (uint a=0; a<31; ++a) {
                uint ia = ((K-2)*32+a+1)*5;
                ulong active=(t0^images[ia]) | (t1^images[ia+1]) | (t2^images[ia+2])
                            | (t3^images[ia+3]) | (t4^images[ia+4]);
                uint w = popcount((uint)active) + popcount((uint)(active>>32));
                uint code = prefix + PREFIXES*(31*b+a);
                best = min(best, ((ulong)w<<32) | code);
            }
        }
    }
    reduction[lane]=best;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (uint stride=WG/2; stride; stride/=2) {
        if (lane<stride) reduction[lane]=min(reduction[lane],reduction[lane+stride]);
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (lane==0) result[group]=reduction[0];
}
)CLC";
void cl_check(cl_int code, const char* operation) {
    if (code != CL_SUCCESS) throw std::runtime_error(std::string(operation) + ": OpenCL error " + std::to_string(code));
}
class GPU {
    cl_context context = nullptr;
    cl_command_queue queue = nullptr;
    cl_program program = nullptr;
    cl_kernel kernel = nullptr;
    cl_mem tables = nullptr, positions = nullptr, results = nullptr;
    unsigned groups = 8, wg = 128;
public:
    std::string name;
    GPU(const Table& table, const std::vector<Positions>& p, int k, unsigned batch) {
        cl_uint n = 0;
        cl_check(clGetPlatformIDs(0, nullptr, &n), "clGetPlatformIDs");
        std::vector<cl_platform_id> platforms(n);
        cl_check(clGetPlatformIDs(n, platforms.data(), nullptr), "clGetPlatformIDs");
        cl_device_id device = nullptr;
        for (auto platform : platforms) {
            cl_int code = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
            if (code == CL_DEVICE_NOT_FOUND) continue;
            cl_check(code, "clGetDeviceIDs");
            break;
        }
        if (!device) throw std::runtime_error("没有可用的 OpenCL GPU");
        char buffer[1024]{};
        cl_check(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(buffer), buffer, nullptr), "clGetDeviceInfo");
        name = buffer;
        cl_int code = 0;
        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &code);
        cl_check(code, "clCreateContext");
        queue = clCreateCommandQueue(context, device, 0, &code);
        cl_check(code, "clCreateCommandQueue");
        program = clCreateProgramWithSource(context, 1, &kernel_source, nullptr, &code);
        cl_check(code, "clCreateProgramWithSource");
        // A single prefix for K=2 would waste most lanes; small K uses CPU.
        std::string flags = "-cl-std=CL1.2 -D K=" + std::to_string(k)
            + " -D PREFIXES=" + std::to_string(power31(k - 2))
            + " -D GROUPS=" + std::to_string(groups) + " -D WG=" + std::to_string(wg);
        code = clBuildProgram(program, 1, &device, flags.c_str(), nullptr, nullptr);
        if (code != CL_SUCCESS) {
            std::size_t length = 0;
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &length);
            std::vector<char> log(length + 1);
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, length, log.data(), nullptr);
            throw std::runtime_error("OpenCL build failed: " + std::string(log.data()));
        }
        kernel = clCreateKernel(program, "enumerate", &code);
        cl_check(code, "clCreateKernel");
        tables = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                sizeof(table), const_cast<Table*>(&table), &code);
        cl_check(code, "clCreateBuffer(table)");
        positions = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   sizeof(Positions)*p.size(), const_cast<Positions*>(p.data()), &code);
        cl_check(code, "clCreateBuffer(positions)");
        results = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ulong)*batch*groups, nullptr, &code);
        cl_check(code, "clCreateBuffer(results)");
        cl_check(clSetKernelArg(kernel, 0, sizeof(tables), &tables), "clSetKernelArg");
        cl_check(clSetKernelArg(kernel, 1, sizeof(positions), &positions), "clSetKernelArg");
        cl_check(clSetKernelArg(kernel, 3, sizeof(results), &results), "clSetKernelArg");
    }
    ~GPU() {
        if (results) clReleaseMemObject(results);
        if (positions) clReleaseMemObject(positions);
        if (tables) clReleaseMemObject(tables);
        if (kernel) clReleaseKernel(kernel);
        if (program) clReleaseProgram(program);
        if (queue) clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }
    std::vector<uint64_t> run(unsigned start, unsigned count) {
        cl_check(clSetKernelArg(kernel, 2, sizeof(start), &start), "clSetKernelArg(start)");
        const std::size_t local = wg, global = std::size_t(count)*groups*wg;
        cl_check(clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, &local, 0, nullptr, nullptr),
                 "clEnqueueNDRangeKernel");
        std::vector<uint64_t> raw(count*groups), out(count, NONE);
        cl_check(clEnqueueReadBuffer(queue, results, CL_TRUE, 0, raw.size()*sizeof(uint64_t),
                                    raw.data(), 0, nullptr, nullptr), "clEnqueueReadBuffer");
        for (unsigned i=0; i<count; ++i)
            for (unsigned g=0; g<groups; ++g) out[i]=std::min(out[i], raw[i*groups+g]);
        return out;
    }
};

State fixed_state(const Positions& p, int k, unsigned code) {
    State s{};
    for (int d=0; d<k; ++d) {
        unsigned v=code%31+1;
        code/=31;
        for (int r=0; r<5; ++r) if ((v>>r)&1) s[r] |= uint64_t(1)<<p[d];
    }
    return s;
}
State mapped_state(const Table& table, const Positions& p, int k, unsigned code) {
    State s{};
    for (int d=0; d<k; ++d) {
        s=xor_states(s,table[p[d]][code%31+1]);
        code/=31;
    }
    return s;
}
} // namespace

int main(int argc, char** argv) {
    std::signal(SIGINT, handle_sigint);
    try {
        ParallelOptions p = parse_parallel(argc, argv);
        const Options& o = p.common;
        const auto begin = Clock::now();
        const Table table = make_table(o);
        check_rotation(table);
        const auto positions = make_positions(o);
        const uint64_t per_position = power31(o.count);
        const unsigned total = unsigned(positions.size());
        const unsigned target = p.limit ? std::min(total, p.limit) : total;
        const uint64_t planned = uint64_t(total)*per_position;
        std::unique_ptr<GPU> gpu;
        if (p.backend != "cpu" && o.count >= 2)
            gpu = std::make_unique<GPU>(table, positions, o.count, p.batch);
        const bool cpu_enabled = p.backend != "opencl" || !gpu;
        std::printf("Mode=%s; fixed %s=%d; map=%s; backend=%s; CPU threads=%d\n",
                    mode_name(o), fixed_name(o), o.count, mapping_name(o), p.backend.c_str(),
                    cpu_enabled ? p.threads : 0);
        if (gpu) std::printf("GPU: %s\n", gpu->name.c_str());
        std::printf("Position sets=%u; target=%u; values per set=%llu; candidates=%llu; rotation orbits=%s\n",
                    total, target, (unsigned long long)per_position, (unsigned long long)planned,
                    o.symmetry ? "on (Burnside checked)" : "off");
        std::fflush(stdout);
        const auto ready = Clock::now();
        auto last_report = ready;
        std::vector<uint64_t> minima(total, NONE);
        std::atomic<unsigned> next{0};
        std::atomic<bool> failed{false};
        std::mutex mutex;
        std::exception_ptr error;
        unsigned done = 0, cpu_done = 0, gpu_done = 0;
        Result result;
        auto record = [&](unsigned index, uint64_t packed, bool on_gpu) {
            if ((packed>>32)>64 || uint32_t(packed)>=per_position)
                throw std::runtime_error("计算后端返回非法候选");
            minima[index] = packed; // Unique owner for every position set.
            std::lock_guard<std::mutex> lock(mutex);
            ++done;
            if (on_gpu) ++gpu_done; else ++cpu_done;
            result.visited += per_position;
            const int w = int(packed>>32);
            if (w < result.minimum) {
                State fixed = fixed_state(positions[index], o.count, uint32_t(packed));
                State image = mapped_state(table, positions[index], o.count, uint32_t(packed));
                if (weight(image) != w) throw std::runtime_error("后端最小值与 CPU 字级复核不一致");
                result.minimum = w;
                result.input = o.input_fixed ? fixed : image;
                result.output = o.input_fixed ? image : fixed;
                std::printf("Improved: minimum=%d; completed=%u/%u; elapsed=%.3f s\n",
                            w, done, total, seconds(ready, Clock::now()));
                std::fflush(stdout);
            }
            const auto now = Clock::now();
            if (seconds(last_report, now) >= 10) {
                const double elapsed = seconds(ready, now);
                std::printf("Progress: %.3f%%; sets=%u/%u (CPU=%u GPU=%u); minimum=%d; "
                            "elapsed=%.1f s; ETA=%.1f s\n", 100.0*done/total, done, total,
                            cpu_done, gpu_done, result.minimum, elapsed, elapsed*(total-done)/done);
                std::fflush(stdout);
                last_report = now;
            }
        };
        auto fail = [&] {
            std::lock_guard<std::mutex> lock(mutex);
            if (!error) error=std::current_exception();
            failed.store(true);
        };
        auto gpu_work = [&] {
            try {
                while (!stop_requested && !failed.load()) {
                    unsigned start = next.fetch_add(p.batch);
                    if (start >= target) break;
                    unsigned count = std::min(p.batch, target-start);
                    auto batch = gpu->run(start, count);
                    for (unsigned i=0; i<count; ++i) record(start+i, batch[i], true);
                }
            } catch (...) { fail(); }
        };
        std::thread device_worker;
        if (gpu) device_worker = std::thread(gpu_work);
        if (cpu_enabled) {
            #pragma omp parallel num_threads(p.threads)
            {
                try {
                    while (!stop_requested && !failed.load()) {
                        unsigned i=next.fetch_add(1);
                        if (i>=target) break;
                        record(i, cpu_position(table, positions[i], o.count), false);
                    }
                } catch (...) { fail(); }
            }
        }
        if (device_worker.joinable()) device_worker.join();
        if (error) std::rethrow_exception(error);
        const auto end = Clock::now();
        result.complete = done == total;
        if (result.minimum > 64) { std::puts("No completed position set."); return 130; }
        if ((o.linear ? transpose_layer(result.output) != result.input : linear_layer(result.input) != result.output)
            || weight(o.input_fixed ? result.input : result.output) != o.count
            || weight(o.input_fixed ? result.output : result.input) != result.minimum)
            throw std::runtime_error("最终状态未通过传播关系复核");
        std::printf("\nStatus: %s\nMinimum %s active columns: %d\nVisited: %llu\n"
                    "Setup: %.6f s\nEnumeration: %.6f s\nCPU sets: %u; GPU sets: %u\n",
                    result.complete ? "已完成枚举，已证明最优" : "部分搜索，未证明最优",
                    o.input_fixed ? "output" : "input", result.minimum,
                    (unsigned long long)result.visited, seconds(begin,ready), seconds(ready,end), cpu_done,gpu_done);
        print_state(o.linear ? "Input mask alpha" : "Input difference X", result.input);
        print_state(o.linear ? "Output mask beta" : "Output difference Y", result.output);
        ExecutionInfo info{p.backend, gpu ? gpu->name : "", cpu_enabled ? p.threads : 0,
                           total, done, cpu_done, gpu_done, planned};
        save_result(o, result, seconds(begin,ready), seconds(ready,end), &info);
        if (!o.prefix.empty()) {
            std::ofstream file(o.prefix+".positions.tsv");
            if (!file) throw std::runtime_error("无法写入逐位置枚举记录");
            file << "position_index\tposition_mask_hex\tminimum\tbase31_value_code\n";
            for (unsigned i=0; i<total; ++i) if (minima[i]!=NONE) {
                file << i << '\t' << std::hex << position_mask(positions[i],o.count) << std::dec
                     << '\t' << (minima[i]>>32) << '\t' << uint32_t(minima[i]) << '\n';
            }
            file.close();
            if (!file) throw std::runtime_error("逐位置记录写入失败");
            std::printf("Saved: %s.json, %s.positions.tsv\n", o.prefix.c_str(),o.prefix.c_str());
        }
        return result.complete ? 0 : 130;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"错误: %s\n",e.what());
        return 1;
    }
}

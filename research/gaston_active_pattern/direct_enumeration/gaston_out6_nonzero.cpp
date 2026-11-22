// Exact bounded enumeration for six active output columns and E != 0.
// A state with fewer than B active input columns has a zero pair among
// any B disjoint input-column pairs. Each selected pair gives an invertible
// 10x10 system determining the last two output-column values uniquely.
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <memory>
#include <numeric>
#include "gaston_common.hpp"

namespace {
using Positions = std::array<int, 6>;
using Pair = std::array<unsigned, 2>;
constexpr uint64_t NONE = UINT64_MAX;
constexpr unsigned PREFIXES = 31 * 31 * 31 * 31;
constexpr uint64_t PER_POSITION = uint64_t(PREFIXES) * 31 * 31;
struct Config {
    std::string backend = "opencl", prefix;
    unsigned cutoff = 26, start = 0, limit = 0, batch = 64;
};
Config parse(int argc, char** argv) {
    Config c;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help") {
            std::puts("差分输出恰好 6 个活跃列，要求 E!=0，最小化输入活跃列数。\n"
                      "  --backend opencl|cpu|brute  默认 opencl；brute 用于小样本独立对照\n"
                      "  --cutoff N                  搜索输入活跃列数 < N，默认 26，范围 1..26\n"
                      "  --start-position N          起始位置代表编号，默认 0\n"
                      "  --limit-positions N         仅搜索 N 个位置代表，默认全部\n"
                      "  --batch-positions N         GPU 批量位置数，默认 64\n"
                      "  --prefix PATH               保存 JSON 和逐位置证明记录\n"
                      "遍历全部位置且找到小于 cutoff 的解，才报告全局最优。");
            std::exit(0);
        }
        if (++i == argc) throw std::invalid_argument(a + " 缺少参数");
        const std::string v = argv[i];
        if (a == "--backend") {
            if (v != "opencl" && v != "cpu" && v != "brute")
                throw std::invalid_argument("backend 必须是 opencl/cpu/brute");
            c.backend = v;
        } else if (a == "--prefix") c.prefix = v;
        else {
            if (a != "--cutoff" && a != "--start-position" && a != "--limit-positions"
                && a != "--batch-positions") throw std::invalid_argument("未知参数: " + a);
            std::size_t used = 0;
            const auto n = std::stoull(v, &used);
            if (v.empty() || v[0] == '-' || used != v.size() || n > UINT32_MAX)
                throw std::invalid_argument("非法整数: " + v);
            if (a == "--cutoff") c.cutoff = unsigned(n);
            if (a == "--start-position") c.start = unsigned(n);
            if (a == "--limit-positions") c.limit = unsigned(n);
            if (a == "--batch-positions") c.batch = unsigned(n);
        }
    }
    if (c.cutoff < 1 || c.cutoff > 26 || c.batch < 1 || c.batch > 4096)
        throw std::invalid_argument("cutoff 范围 1..26，batch 范围 1..4096");
    return c;
}
uint64_t mask_of(const Positions& p) {
    uint64_t m = 0;
    for (int j : p) m |= uint64_t(1) << j;
    return m;
}
std::vector<Positions> position_orbits() {
    std::vector<Positions> out;
    for (int a = 1; a < 60; ++a)
    for (int b = a+1; b < 61; ++b)
    for (int c = b+1; c < 62; ++c)
    for (int d = c+1; d < 63; ++d)
    for (int e = d+1; e < 64; ++e) {
        Positions p{0,a,b,c,d,e};
        uint64_t m = mask_of(p);
        bool keep = true;
        for (int j : p) if (rol(m, -j) < m) { keep = false; break; }
        if (keep) out.push_back(p);
    }
    uint64_t expected = 0;
    for (int s=0; s<64; ++s) {
        const int g = std::gcd(64,s), len = 64/g;
        if (6 % len == 0) expected += uint64_t(choose(g,6/len)+0.5L);
    }
    if (out.size() != expected/64) throw std::runtime_error("位置轨道计数错误");
    return out;
}
unsigned column(const State& s, unsigned j) {
    unsigned v = 0;
    for (int r=0; r<5; ++r) v |= unsigned((s[r] >> j) & 1) << r;
    return v;
}
uint64_t effect(const State& x) {
    uint64_t p=0,q=0;
    for (int r=0;r<5;++r) { const auto a=rol(x[r],EAST[r]);p^=a;q^=rol(a,TWIN[r]); }
    return p ^ rol(p,1) ^ q ^ rol(q,18);
}
struct PairTables {
    unsigned count;
    std::vector<Pair> pairs;
    std::vector<uint16_t> solve;
    // A linear 10-bit map is two 5-bit lookup tables, XORed together.
    PairTables(const Table& table, unsigned n): count(n), pairs(64*n), solve(64*n*64) {
        for (unsigned gap=1; gap<64; ++gap) {
            uint64_t used=0;
            unsigned found=0;
            for (unsigned a=0; a<64 && found<n; ++a) {
                if ((used>>a)&1) continue;
                for (unsigned b=a+1; b<64; ++b) {
                    if ((used>>b)&1) continue;
                    std::array<unsigned,10> rows{};
                    for (unsigned j=0; j<10; ++j) {
                        const State& image=table[j<5 ? 0 : gap][1u<<(j%5)];
                        unsigned rhs=column(image,a)|(column(image,b)<<5);
                        for (unsigned r=0; r<10; ++r) rows[r]|=((rhs>>r)&1)<<j;
                    }
                    for (unsigned r=0; r<10; ++r) rows[r]|=1u<<(10+r);
                    bool invertible=true;
                    for (unsigned j=0; j<10; ++j) {
                        unsigned pivot=j;
                        while (pivot<10 && !((rows[pivot]>>j)&1)) ++pivot;
                        if (pivot==10) { invertible=false;break; }
                        std::swap(rows[j],rows[pivot]);
                        for (unsigned r=0;r<10;++r) if (r!=j && ((rows[r]>>j)&1)) rows[r]^=rows[j];
                    }
                    if (!invertible) continue;
                    const unsigned index=gap*n+found;
                    pairs[index]={a,b};
                    for (unsigned part=0;part<64;++part) {
                        const unsigned rhs=part<32 ? part : (part-32)<<5;
                        unsigned values=0;
                        for (unsigned r=0;r<10;++r)
                            values|=unsigned(__builtin_parity((rows[r]>>10)&rhs))<<r;
                        solve[index*64+part]=uint16_t(values);
                    }
                    for (unsigned rhs=0;rhs<1024;++rhs) {
                        const unsigned values=solve[index*64+(rhs&31)]^solve[index*64+32+(rhs>>5)];
                        const auto image=xor_states(table[0][values&31],table[gap][values>>5]);
                        if ((column(image,a)|(column(image,b)<<5))!=rhs)
                            throw std::runtime_error("10x10 反解表未通过全部 1024 个值的复核");
                    }
                    used|=(uint64_t(1)<<a)|(uint64_t(1)<<b);
                    ++found;
                    break;
                }
            }
            if (found!=n || __builtin_popcountll(used)!=int(2*n))
                throw std::runtime_error("未找到足够的互不重叠可逆输入列对");
        }
    }
};
struct Answer { uint64_t best=NONE, tested=0; };
State fixed_state(const Positions& p, unsigned code) {
    State y{};
    for (int d=0;d<6;++d) {
        const unsigned v=code%31+1;code/=31;
        for (int r=0;r<5;++r) if ((v>>r)&1) y[r]|=uint64_t(1)<<p[d];
    }
    return y;
}
State mapped_state(const Table& t,const Positions& p,unsigned code) {
    State x{};
    for (int d=0;d<6;++d) { x=xor_states(x,t[p[d]][code%31+1]);code/=31; }
    return x;
}
Answer cpu_position(const Table& t,const PairTables& pt,const Positions& p,unsigned bound,bool brute) {
    Answer out;
    const unsigned origin=p[4],gap=p[5]-p[4];
    for (unsigned prefix=0;prefix<PREFIXES;++prefix) {
        State s{};
        uint64_t y0=0;
        unsigned digits=prefix;
        for (unsigned d=0;d<4;++d) {
            const unsigned v=digits%31+1,j=(p[d]+64-origin)%64;digits/=31;
            s=xor_states(s,t[j][v]);
            if (v&1) y0|=uint64_t(1)<<j;
        }
        auto consider=[&](unsigned a,unsigned b) {
            ++out.tested;
            const auto x=xor_states(xor_states(s,t[0][a]),t[gap][b]);
            const unsigned w=weight(x);
            const unsigned code=prefix+PREFIXES*((b-1)*31+a-1);
            const uint64_t packed=(uint64_t(w)<<32)|code;
            const uint64_t output0=y0|(a&1)|(uint64_t(b&1)<<gap);
            if (w<bound && packed<out.best && x[0]!=output0) out.best=packed;
        };
        if (brute) {
            for (unsigned b=1;b<32;++b) for (unsigned a=1;a<32;++a) consider(a,b);
        } else {
            for (unsigned q=0;q<bound;++q) {
                const auto pair=pt.pairs[gap*pt.count+q];
                const unsigned rhs=column(s,pair[0])|(column(s,pair[1])<<5);
                const unsigned index=(gap*pt.count+q)*64;
                const unsigned v=pt.solve[index+(rhs&31)]^pt.solve[index+32+(rhs>>5)];
                if ((v&31) && (v>>5)) consider(v&31,v>>5);
            }
        }
    }
    return out;
}
const char* source=R"CLC(
uint col(ulong s0,ulong s1,ulong s2,ulong s3,ulong s4,uint j) {
    return ((s0>>j)&1) | (((s1>>j)&1)<<1) | (((s2>>j)&1)<<2)
        | (((s3>>j)&1)<<3) | (((s4>>j)&1)<<4);
}
__kernel void search(__global const ulong* table,__global const int* positions,
                     __global const uint2* pairs,__global const ushort* solve,
                     uint start,uint bound,__global ulong2* results) {
    uint lane=get_local_id(0),group=get_group_id(0),pos=start+group/GROUPS,part=group%GROUPS;
    uint origin=positions[pos*6+4],gap=positions[pos*6+5]-origin;
    __local ulong images[6*32*5];
    __local ushort lut[CUTOFF*64];
    __local uint2 zero_pairs[CUTOFF];
    __local ulong minima[WG],counts[WG];
    for (uint i=lane;i<6*32*5;i+=WG) {
        uint d=i/(32*5),off=i%(32*5);
        uint j=(positions[pos*6+d]+64-origin)%64;
        images[i]=table[j*32*5+off];
    }
    for (uint i=lane;i<CUTOFF*64;i+=WG) lut[i]=solve[gap*CUTOFF*64+i];
    if (lane<CUTOFF) zero_pairs[lane]=pairs[gap*CUTOFF+lane];
    barrier(CLK_LOCAL_MEM_FENCE);
    ulong best=(ulong)(-1),tested=0;
    for (uint prefix=part*WG+lane;prefix<PREFIXES;prefix+=GROUPS*WG) {
        ulong s0=0,s1=0,s2=0,s3=0,s4=0,y0=0;
        uint digits=prefix;
        #pragma unroll
        for (uint d=0;d<4;++d) {
            uint v=digits%31+1,i=(d*32+v)*5;digits/=31;
            s0^=images[i];s1^=images[i+1];s2^=images[i+2];s3^=images[i+3];s4^=images[i+4];
            y0|=((ulong)(v&1))<<((positions[pos*6+d]+64-origin)%64);
        }
        for (uint q=0;q<bound;++q) {
            uint2 pair=zero_pairs[q];
            uint rhs=col(s0,s1,s2,s3,s4,pair.x)|(col(s0,s1,s2,s3,s4,pair.y)<<5);
            uint v=lut[q*64+(rhs&31)]^lut[q*64+32+(rhs>>5)],a=v&31,b=v>>5;
            if (!a || !b) continue;
            ++tested;
            uint ia=(4*32+a)*5,ib=(5*32+b)*5;
            ulong x0=s0^images[ia]^images[ib];
            ulong active=x0|(s1^images[ia+1]^images[ib+1])|(s2^images[ia+2]^images[ib+2])
                |(s3^images[ia+3]^images[ib+3])|(s4^images[ia+4]^images[ib+4]);
            uint w=popcount((uint)active)+popcount((uint)(active>>32));
            uint code=prefix+PREFIXES*((b-1)*31+a-1);
            ulong packed=((ulong)w<<32)|code;
            // EAST[0]=WEST[0]=0 => X[0] XOR Y[0]=ROL(E,23).
            if (w<bound && packed<best && x0!=(y0|(ulong)(a&1)|((ulong)(b&1)<<gap))) best=packed;
        }
    }
    minima[lane]=best;counts[lane]=tested;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (uint stride=WG/2;stride;stride/=2) {
        if (lane<stride) {minima[lane]=min(minima[lane],minima[lane+stride]);counts[lane]+=counts[lane+stride];}
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (lane==0) results[group]=(ulong2)(minima[0],counts[0]);
}
)CLC";
void check(cl_int e,const char* op) {
    if (e!=CL_SUCCESS) throw std::runtime_error(std::string(op)+": OpenCL "+std::to_string(e));
}
class GPU {
    cl_context context=nullptr;
    cl_command_queue queue=nullptr;
    cl_program program=nullptr;
    cl_kernel kernel=nullptr;
    std::array<cl_mem,5> buffers{};
    static constexpr unsigned groups=8,wg=128;
public:
    std::string name;
    GPU(const Table& t,const PairTables& pt,const std::vector<Positions>& p,unsigned batch) {
        cl_uint n=0;check(clGetPlatformIDs(0,nullptr,&n),"platform count");
        std::vector<cl_platform_id> platforms(n);check(clGetPlatformIDs(n,platforms.data(),nullptr),"platforms");
        cl_device_id device=nullptr;
        for (auto platform:platforms) {
            auto e=clGetDeviceIDs(platform,CL_DEVICE_TYPE_GPU,1,&device,nullptr);
            if (e==CL_DEVICE_NOT_FOUND) continue;
            check(e,"device");break;
        }
        if (!device) throw std::runtime_error("没有 OpenCL GPU");
        char text[1024]{};check(clGetDeviceInfo(device,CL_DEVICE_NAME,sizeof(text),text,nullptr),"device name");name=text;
        cl_int e=0;
        context=clCreateContext(nullptr,1,&device,nullptr,nullptr,&e);check(e,"context");
        queue=clCreateCommandQueue(context,device,0,&e);check(e,"queue");
        program=clCreateProgramWithSource(context,1,&source,nullptr,&e);check(e,"program");
        const auto flags="-cl-std=CL1.2 -D CUTOFF="+std::to_string(pt.count)+" -D PREFIXES="
            +std::to_string(PREFIXES)+" -D GROUPS="+std::to_string(groups)+" -D WG="+std::to_string(wg);
        e=clBuildProgram(program,1,&device,flags.c_str(),nullptr,nullptr);
        if (e!=CL_SUCCESS) {
            std::size_t len=0;clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,0,nullptr,&len);
            std::vector<char> log(len+1);clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,len,log.data(),nullptr);
            throw std::runtime_error(log.data());
        }
        kernel=clCreateKernel(program,"search",&e);check(e,"kernel");
        const std::array<std::pair<const void*,std::size_t>,4> inputs{{
            {&t,sizeof(t)},{p.data(),p.size()*sizeof(Positions)},
            {pt.pairs.data(),pt.pairs.size()*sizeof(Pair)},{pt.solve.data(),pt.solve.size()*sizeof(uint16_t)}}};
        for (unsigned i=0;i<4;++i) {
            buffers[i]=clCreateBuffer(context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,inputs[i].second,
                                      const_cast<void*>(inputs[i].first),&e);check(e,"input buffer");
            check(clSetKernelArg(kernel,i,sizeof(cl_mem),&buffers[i]),"input argument");
        }
        buffers[4]=clCreateBuffer(context,CL_MEM_WRITE_ONLY,batch*groups*2*sizeof(uint64_t),nullptr,&e);check(e,"output buffer");
        check(clSetKernelArg(kernel,6,sizeof(cl_mem),&buffers[4]),"output argument");
    }
    ~GPU() {
        for (auto b:buffers) if (b) clReleaseMemObject(b);
        if (kernel) clReleaseKernel(kernel);
        if (program) clReleaseProgram(program);
        if (queue) clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }
    std::vector<Answer> run(unsigned start,unsigned count,unsigned bound) {
        check(clSetKernelArg(kernel,4,sizeof(start),&start),"start");
        check(clSetKernelArg(kernel,5,sizeof(bound),&bound),"bound");
        const std::size_t global=std::size_t(count)*groups*wg,local=wg;
        check(clEnqueueNDRangeKernel(queue,kernel,1,nullptr,&global,&local,0,nullptr,nullptr),"launch");
        std::vector<uint64_t> raw(count*groups*2);
        check(clEnqueueReadBuffer(queue,buffers[4],CL_TRUE,0,raw.size()*sizeof(uint64_t),raw.data(),0,nullptr,nullptr),"read");
        std::vector<Answer> out(count);
        for (unsigned i=0;i<count;++i) for (unsigned g=0;g<groups;++g) {
            const unsigned index=(i*groups+g)*2;
            out[i].best=std::min(out[i].best,raw[index]);out[i].tested+=raw[index+1];
        }
        return out;
    }
};
void write_state(std::ostream& f,const char* name,const State& s) {
    f<<",\n  \""<<name<<"_active\": "<<weight(s)<<",\n  \""<<name<<"_words\": [";
    for (int r=0;r<5;++r) {
        char hex[19];std::snprintf(hex,sizeof(hex),"0x%016llx",(unsigned long long)s[r]);
        if (r) f<<", ";
        f<<'"'<<hex<<'"';
    }
    f<<']';
}
}

int main(int argc,char** argv) {
    std::signal(SIGINT,handle_sigint);
    try {
        const auto c=parse(argc,argv);
        const auto begin=Clock::now();
        Options o;o.count=6;
        const auto table=make_table(o);
        // Confirm the normalization used in the GPU on all column values.
        for (int j=0;j<64;++j) for (int v=0;v<32;++v) for (int r=0;r<5;++r)
            if (table[j][v][r]!=rol(table[0][v][r],j)) throw std::runtime_error("移位复核失败");
        const PairTables pt(table,c.cutoff);
        const auto positions=position_orbits();
        const unsigned total=unsigned(positions.size());
        if (c.start>=total) throw std::invalid_argument("起始位置编号超出范围");
        const unsigned end=c.limit ? c.start+std::min(c.limit,total-c.start) : total;
        std::unique_ptr<GPU> gpu;
        if (c.backend=="opencl") gpu=std::make_unique<GPU>(table,pt,positions,c.batch);
        std::ofstream records;
        if (!c.prefix.empty()) {
            const auto parent=std::filesystem::path(c.prefix).parent_path();
            if (!parent.empty()) std::filesystem::create_directories(parent);
            records.open(c.prefix+".positions.tsv");
            if (!records) throw std::runtime_error("无法写入位置记录");
            records<<"position_index\tposition_mask_hex\tsearch_below\tminimum_below_cutoff\tbase31_value_code\ttested\n";
        }
        std::printf("Output active=6; theta-effect=nonzero; cutoff=%u; backend=%s; device=%s\n"
                    "Position orbits=%u; range=[%u,%u); full output assignments=%llu\n"
                    "All disjoint column pairs and 1024-entry inverse tables verified.\n",
                    c.cutoff,c.backend.c_str(),gpu?gpu->name.c_str():"CPU",total,c.start,end,
                    (unsigned long long)(uint64_t(total)*PER_POSITION));
        std::fflush(stdout);
        const auto ready=Clock::now();auto last=ready;
        uint64_t best=NONE,tested=0;
        unsigned best_index=0,done=0,bound=c.cutoff;
        for (unsigned start=c.start;start<end && !stop_requested;) {
            const unsigned n=std::min(gpu?c.batch:1u,end-start),search_bound=bound;
            std::vector<Answer> answers;
            if (gpu) answers=gpu->run(start,n,search_bound);
            else answers.push_back(cpu_position(table,pt,positions[start],search_bound,c.backend=="brute"));
            for (unsigned j=0;j<n;++j) {
                const auto& a=answers[j];const unsigned index=start+j;
                ++done;tested+=a.tested;
                if (a.best!=NONE) {
                    const unsigned code=uint32_t(a.best),w=unsigned(a.best>>32);
                    if (code>=PER_POSITION || w>=search_bound) throw std::runtime_error("后端候选非法");
                    const auto x=mapped_state(table,positions[index],code),y=fixed_state(positions[index],code);
                    if (weight(x)!=int(w) || linear_layer(x)!=y || effect(x)==0 || x[0]==y[0])
                        throw std::runtime_error("后端候选未通过独立传播、重量及 E 复核");
                    if (a.best<best) {
                        best=a.best;best_index=index;bound=std::min(bound,w);
                        std::printf("Improved: input=%u; position=%u; done=%u/%u; elapsed=%.3f s\n",
                                    w,index,done,total,seconds(ready,Clock::now()));
                        std::fflush(stdout);
                    }
                }
                if (records.is_open()) {
                    records<<index<<'\t'<<std::hex<<mask_of(positions[index])<<std::dec<<'\t'<<search_bound<<'\t';
                    if (a.best==NONE) records<<"null\tnull";
                    else records<<(a.best>>32)<<'\t'<<uint32_t(a.best);
                    records<<'\t'<<a.tested<<'\n';
                }
            }
            start+=n;
            const auto now=Clock::now();
            if (seconds(last,now)>=10) {
                std::printf("Progress %.3f%%; completed=%u/%u; best=%d; elapsed=%.1f s; ETA=%.1f s\n",
                            100.0*done/(end-c.start),done,end-c.start,best==NONE?-1:int(best>>32),
                            seconds(ready,now),seconds(ready,now)*(end-c.start-done)/done);
                std::fflush(stdout);if (records.is_open()) records.flush();last=now;
            }
        }
        const auto finish=Clock::now();
        const bool complete=c.start==0 && done==total,proven=complete && best!=NONE;
        if (records.is_open()) {records.close();if (!records) throw std::runtime_error("位置记录写入失败");}
        State x{},y{};
        if (best!=NONE) {
            x=mapped_state(table,positions[best_index],uint32_t(best));
            y=fixed_state(positions[best_index],uint32_t(best));
            std::printf("\nMinimum input active columns: %llu; proven_optimal=%s\n",
                        (unsigned long long)(best>>32),proven?"true":"false");
            print_state("Input X",x);print_state("Output Y",y);
            std::printf("Theta effect E=0x%016llx\n",(unsigned long long)effect(x));
        } else std::printf("\nNo state below cutoff %u in the completed position sets.\n",c.cutoff);
        std::printf("Completed positions=%u/%u; tested=%llu; covered output assignments=%llu\n"
                    "Setup=%.6f s; search=%.6f s\n",done,total,(unsigned long long)tested,
                    (unsigned long long)(uint64_t(done)*PER_POSITION),seconds(begin,ready),seconds(ready,finish));
        if (!c.prefix.empty()) {
            std::ofstream f(c.prefix+".json");if (!f) throw std::runtime_error("无法写入 JSON");
            f.precision(18);
            f<<"{\n  \"method\": \"exact_enumeration_with_zero_pair_pruning\",\n  \"mode\": \"difference\","
              "\n  \"fixed_side\": \"output\",\n  \"fixed_count\": 6,\n  \"theta_effect_constraint\": \"nonzero\","
              "\n  \"cutoff\": "<<c.cutoff<<",\n  \"position_sets\": "<<total
             <<",\n  \"start_position\": "<<c.start<<",\n  \"completed_position_sets\": "<<done
             <<",\n  \"tested_candidates_with_duplicates\": "<<tested
             <<",\n  \"covered_output_states\": "<<uint64_t(done)*PER_POSITION
             <<",\n  \"enumeration_complete\": "<<(complete?"true":"false")
             <<",\n  \"proven_optimal\": "<<(proven?"true":"false")
             <<",\n  \"objective\": ";
            if (best==NONE) f<<"null";else f<<(best>>32);
            f<<",\n  \"backend\": \""<<c.backend<<"\",\n  \"device\": \""<<(gpu?gpu->name:"CPU")<<"\""
             <<",\n  \"setup_seconds\": "<<seconds(begin,ready)<<",\n  \"search_seconds\": "<<seconds(ready,finish);
            if (best!=NONE) {
                write_state(f,"input",x);write_state(f,"output",y);
                char hex[19];std::snprintf(hex,sizeof(hex),"0x%016llx",(unsigned long long)effect(x));
                f<<",\n  \"theta_effect_hex\": \""<<hex<<"\",\n  \"verified\": true";
            }
            f<<"\n}\n";f.close();if (!f) throw std::runtime_error("JSON 写入失败");
        }
        return proven?0:(complete?2:130);
    } catch (const std::exception& e) {std::fprintf(stderr,"错误: %s\n",e.what());return 1;}
}

// Gaston-320, three chi layers, exact linear layers, Hsilu Table 14 lossy chi.
// Free endpoints: a0 -> chi -> b0 -> L -> a1 -> chi -> b1 -> L -> a2 -> chi -> b2.
// The initial invertible L is eliminated. Optional explicit kernel/count subproblems.
#include <gurobi_c++.h>
#include <array>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "src/activity_profile_data.hpp"
#include "src/kernel_type_data.hpp"

namespace fs = std::filesystem;
using Word = std::array<GRBVar, 64>;
using State = std::array<Word, 5>;
using Value = std::array<uint64_t, 5>;
constexpr std::array<int, 5> EAST{0, 60, 22, 27, 4};
constexpr std::array<int, 5> WEST{0, 56, 31, 46, 43};
constexpr std::array<int, 5> TWIN{25, 32, 52, 60, 63};
constexpr std::array<int, 6> COST{0, 2, 3, 3, 4, 4};
GRBModel* active_model = nullptr;
void interrupt(int) { if (active_model) active_model->terminate(); }
struct ModelGuard {
    explicit ModelGuard(GRBModel& m) { active_model = &m; }
    ~ModelGuard() { active_model = nullptr; }
};
int mod64(int z) { return (z % 64 + 64) % 64; }
uint64_t rol(uint64_t x, int n) { return (x << n) | (x >> ((64 - n) & 63)); }
std::string bitname(const std::string& s, int z) { return s + "_" + std::to_string(z); }
Word word(GRBModel& m, const std::string& s) {
    Word v;
    for (int z = 0; z < 64; ++z) v[z] = m.addVar(0, 1, 0, GRB_BINARY, bitname(s, z));
    return v;
}
State state(GRBModel& m, const std::string& s) {
    State v;
    for (int i = 0; i < 5; ++i) v[i] = word(m, s + "_" + std::to_string(i));
    return v;
}
struct Gate { std::string name; std::vector<GRBVar> in; GRBVar out; };
std::vector<Gate> gates;
struct KernelCheck { std::string name; std::vector<GRBVar> bits; };
std::vector<KernelCheck> kernel_checks;
void parity(GRBModel& m, const std::vector<GRBVar>& in, GRBVar out,
            const std::string& name) {
    if (in.size() != 2 && in.size() != 3) throw std::logic_error("XOR arity");
    std::vector<GRBVar> v = in;
    v.push_back(out);
    for (unsigned mask = 1; mask < (1U << v.size()); ++mask) {
        const int k = __builtin_popcount(mask);
        if (!(k & 1)) continue;
        GRBLinExpr e = 0;
        for (unsigned j = 0; j < v.size(); ++j) e += ((mask >> j) & 1) ? v[j] : -v[j];
        m.addConstr(e <= k - 1, name + "_f" + std::to_string(mask));
    }
    gates.push_back({name, in, out});
}
Word linear(GRBModel& m, const State& x, const State& y, const std::string& s) {
    State a;
    for (int i = 0; i < 5; ++i)
        for (int z = 0; z < 64; ++z) a[i][z] = x[i][mod64(z - EAST[i])];
    auto p3 = word(m, s + "_p3"), p = word(m, s + "_p");
    auto q3 = word(m, s + "_q3"), q = word(m, s + "_q");
    auto e3 = word(m, s + "_e3"), e = word(m, s + "_e");
    for (int z = 0; z < 64; ++z) {
        parity(m, {a[0][z], a[1][z], a[2][z]}, p3[z], bitname(s + "_P3", z));
        parity(m, {p3[z], a[3][z], a[4][z]}, p[z], bitname(s + "_P", z));
        parity(m, {a[0][mod64(z-TWIN[0])], a[1][mod64(z-TWIN[1])],
                   a[2][mod64(z-TWIN[2])]}, q3[z], bitname(s + "_Q3", z));
        parity(m, {q3[z], a[3][mod64(z-TWIN[3])], a[4][mod64(z-TWIN[4])]},
               q[z], bitname(s + "_Q", z));
    }
    for (int z = 0; z < 64; ++z) {
        parity(m, {p[z], p[mod64(z-1)], q[z]}, e3[z], bitname(s + "_E3", z));
        parity(m, {e3[z], q[mod64(z-18)]}, e[z], bitname(s + "_E", z));
    }
    for (int i = 0; i < 5; ++i)
        for (int z = 0; z < 64; ++z)
            parity(m, {a[i][z], e[mod64(z-23)]}, y[i][mod64(z+WEST[i])],
                   bitname(s + "_Y" + std::to_string(i), z));
    return e;
}
void kernel_linear(GRBModel& m, const State& x, State& y, const std::string& s) {
    State a;
    for(int i=0;i<5;++i) for(int z=0;z<64;++z) {
        a[i][z]=x[i][mod64(z-EAST[i])];
        y[i][z]=x[i][mod64(z-EAST[i]-WEST[i])];
    }
    auto p3=word(m,s+"_p3"), p=word(m,s+"_p");
    auto q3=word(m,s+"_q3"), q=word(m,s+"_q");
    for(int z=0;z<64;++z) {
        parity(m,{a[0][z],a[1][z],a[2][z]},p3[z],bitname(s+"_P3",z));
        parity(m,{p3[z],a[3][z],a[4][z]},p[z],bitname(s+"_P",z));
        parity(m,{a[0][mod64(z-TWIN[0])],a[1][mod64(z-TWIN[1])],a[2][mod64(z-TWIN[2])]},q3[z],bitname(s+"_Q3",z));
        parity(m,{q3[z],a[3][mod64(z-TWIN[3])],a[4][mod64(z-TWIN[4])]},q[z],bitname(s+"_Q",z));
    }
    for(int z=0;z<64;++z) {
        std::vector<GRBVar> v{p[z],p[mod64(z-1)],q[z],q[mod64(z-18)]};
        auto name=bitname(s+"_kernel",z);
        for(unsigned mask=1;mask<16;++mask) if(__builtin_popcount(mask)&1) {
            GRBLinExpr e=0;
            for(int j=0;j<4;++j) e+=((mask>>j)&1)?v[j]:-v[j];
            m.addConstr(e<=__builtin_popcount(mask)-1,name+"_f"+std::to_string(mask));
        }
        kernel_checks.push_back({name,v});
    }
}
struct Chi { Word active; Word extra; };
Chi lossy(GRBModel& m, const State& x, const State& y, int r, GRBLinExpr& obj, bool or_hull) {
    const auto s = "chi" + std::to_string(r);
    Chi c;
    c.active = word(m, s + "_active");
    for (int z = 0; z < 64; ++z) {
        c.extra[z] = m.addVar(0, 2, 0, GRB_INTEGER, bitname(s + "_extra", z));
        GRBLinExpr h = 0, k = 0;
        for (int i = 0; i < 5; ++i) { h += x[i][z]; k += y[i][z]; }
        const auto p = bitname(s, z);
        // Four lossy support inequalities: nonzero input iff nonzero output.
        m.addConstr(h >= c.active[z], p + "_input_min");
        m.addConstr(h <= 5*c.active[z], p + "_input_max");
        m.addConstr(k >= c.active[z], p + "_output_min");
        m.addConstr(k <= 5*c.active[z], p + "_output_max");
        // Two cost inequalities. Under minimization, extra=ceil((h-active)/2).
        // This integer is a cost variable, never an XOR parity slack.
        m.addConstr(h <= c.active[z] + 2*c.extra[z], p + "_cost_min");
        m.addConstr(c.extra[z] <= 2*c.active[z], p + "_cost_max");
        if (or_hull) for (int i=0; i<5; ++i) {
            // Same binary lossy relation; stronger fractional activity relaxation.
            m.addConstr(x[i][z] <= c.active[z], p + "_input_bit" + std::to_string(i));
            m.addConstr(y[i][z] <= c.active[z], p + "_output_bit" + std::to_string(i));
        }
        obj += 2*c.active[z] + c.extra[z];
    }
    return c;
}
void activity_cuts(GRBModel& m, const std::array<Chi,3>& chis,
                   const std::array<Word,2>& effects, bool kernel_only=false) {
    std::array<Word,3> select;
    std::array<GRBLinExpr,3> count, forward, reverse;
    for (int r=0; r<3; ++r) {
        select[r] = word(m,"count"+std::to_string(r));
        GRBLinExpr one=0, weighted=0;
        count[r]=0; forward[r]=0; reverse[r]=0;
        for (int z=0; z<64; ++z) {
            // Selector index z represents exactly z+1 active columns.
            one += select[r][z];
            weighted += (z+1)*select[r][z];
            count[r] += chis[r].active[z];
            forward[r] += FORWARD_ACTIVITY_LOWER[z+1]*select[r][z];
            reverse[r] += REVERSE_ACTIVITY_LOWER[z+1]*select[r][z];
        }
        m.addConstr(one==1,"profile_onehot_"+std::to_string(r));
        m.addConstr(count[r]==weighted,"profile_count_"+std::to_string(r));
    }
    GRBLinExpr total=count[0]+count[1]+count[2];
    m.addConstr(total>=18,"profile_three_round_total");
    GRBLinExpr middle_total=0;
    for(int k=1;k<=64;++k) middle_total+=MIDDLE_TOTAL_ACTIVITY_LOWER[k]*select[1][k-1];
    m.addConstr(total>=middle_total,"profile_middle_total");
    for (int r=0; r<2; ++r) {
        const auto p="profile_L"+std::to_string(r);
        if(kernel_only) {
            m.addConstr(count[r+1]>=forward[r],p+"_forward");
            m.addConstr(count[r]>=reverse[r+1],p+"_reverse");
            m.addConstr(count[r]+count[r+1]>=12,p+"_branch");
            m.addConstr(count[r+1]<=5*count[r],p+"_kernel_forward");
            m.addConstr(count[r]<=5*count[r+1],p+"_kernel_reverse");
            continue;
        }
        GRBVar n=m.addVar(0,1,0,GRB_BINARY,"L"+std::to_string(r)+"_nonzero");
        GRBLinExpr effect_weight=0;
        for (int z=0; z<64; ++z) {
            effect_weight += effects[r][z];
            m.addConstr(effects[r][z]<=n,p+"_effect_bit"+std::to_string(z));
        }
        m.addConstr(n<=effect_weight,p+"_effect_or");
        // E has even parity, so a nonzero effect contains at least two bits.
        m.addConstr(effect_weight>=2*n,p+"_effect_even_min");
        m.addConstr(count[r+1]>=forward[r],p+"_forward");
        m.addConstr(count[r]>=reverse[r+1],p+"_reverse");
        m.addConstr(count[r]+count[r+1]>=12+4*n,p+"_branch");
        // Output=6 -> input>=13 only outside the kernel. No unconditional 13 cut.
        m.addConstr(count[r]>=reverse[r+1]+7*(select[r+1][5]+n-1),p+"_outside_output6");
        // In the kernel, L only moves existing bits within rows.
        m.addConstr(count[r+1]<=5*count[r]+59*n,p+"_kernel_forward");
        m.addConstr(count[r]<=5*count[r+1]+59*n,p+"_kernel_reverse");
        // Exhaustive 1..64 count-triple combination of the imported profiles.
        m.addConstr(total>=18+4*n,p+"_three_round_total");
    }
}
void kernel_type_cuts(GRBModel& m, const std::array<State,6>& states,
                      const std::array<Chi,3>& chis, const GRBLinExpr& obj, int middle) {
    std::array<GRBLinExpr,3> counts;
    std::array<GRBLinExpr,2> weights,extras;
    for(int r=0;r<3;++r) {
        counts[r]=0;
        for(int z=0;z<64;++z) counts[r]+=chis[r].active[z];
    }
    for(int r=0;r<2;++r) {
        weights[r]=0;extras[r]=0;
        for(int z=0;z<64;++z) {
            extras[r]+=chis[r+1].extra[z];
            for(int i=0;i<5;++i) weights[r]+=states[2*(r+1)][i][z];
        }
    }
    GRBLinExpr one=0,floor=0;
    std::array<GRBLinExpr,3> tc{0,0,0};
    std::array<GRBLinExpr,2> tw{0,0},te{0,0};
    for(std::size_t j=0;j<KERNEL_TYPES.size();++j) {
        const auto& t=KERNEL_TYPES[j];
        if(t.c1!=middle) continue;
        auto q=m.addVar(0,1,0,GRB_BINARY,"kernel_type_"+std::to_string(j));
        one+=q; floor+=t.weight_floor*q;
        tc[0]+=t.c0*q;tc[1]+=t.c1*q;tc[2]+=t.c2*q;
        tw[0]+=t.h1*q;tw[1]+=t.h2*q;
        te[0]+=((t.h1-t.c1+1)/2)*q;te[1]+=((t.h2-t.c2+1)/2)*q;
    }
    m.addConstr(one==1,"kernel_type_onehot");
    for(int r=0;r<3;++r) m.addConstr(counts[r]==tc[r],"kernel_type_count"+std::to_string(r));
    for(int r=0;r<2;++r) {
        m.addConstr(weights[r]==tw[r],"kernel_type_bits"+std::to_string(r));
        m.addConstr(extras[r]>=te[r],"kernel_type_extra"+std::to_string(r));
    }
    m.addConstr(obj>=floor,"kernel_type_weight_floor");
}
int chi5(int x) {
    int y = 0;
    for (int i=0; i<5; ++i)
        y |= (((x>>i)&1) ^ ((1^((x>>((i+1)%5))&1)) & ((x>>((i+2)%5))&1))) << i;
    return y;
}
int ddt(int a, int b) {
    int n=0;
    for (int x=0; x<32; ++x) n += ((chi5(x)^chi5(x^a))==b);
    return n;
}
int col(const Value& x, int z) {
    int a=0;
    for (int i=0; i<5; ++i) a |= ((x[i]>>z)&1) << i;
    return a;
}
Value linear_value(const Value& x) {
    Value a, y;
    uint64_t p=0, q=0;
    for (int i=0; i<5; ++i) { a[i]=rol(x[i],EAST[i]); p^=a[i]; q^=rol(a[i],TWIN[i]); }
    uint64_t e=rol(p^rol(p,1)^q^rol(q,18),23);
    for (int i=0; i<5; ++i) y[i]=rol(a[i]^e,WEST[i]);
    return y;
}
void atomic_write(const fs::path& path, const std::string& text) {
    auto tmp=path; tmp += ".tmp";
    { std::ofstream f(tmp); if (!f) throw std::runtime_error("Cannot write " + tmp.string());
      f << text; if (!f) throw std::runtime_error("Write failed " + tmp.string()); }
    fs::rename(tmp,path);
}
std::string number(double x) {
    if (!std::isfinite(x) || std::abs(x)>=GRB_INFINITY/2) return "null";
    std::ostringstream s; s << std::setprecision(16) << x; return s.str();
}
using Values = std::array<Value,6>;
std::string solution_json(const Values& v, double obj) {
    if (linear_value(v[1])!=v[2] || linear_value(v[3])!=v[4])
        throw std::runtime_error("Incumbent violates exact linear propagation");
    std::ostringstream s;
    s << "{\n  \"scope\": \"lossy MILP incumbent, not automatically a valid DDT trail\",\n"
      << "  \"solver_objective\": " << number(obj) << ",\n  \"states\": {\n";
    const std::array<std::string,6> labels{"a0","b0","a1","b1","a2","b2"};
    for (int t=0; t<6; ++t) {
        s << "    \"" << labels[t] << "\": [";
        for (int i=0; i<5; ++i) {
            char h[32]; std::snprintf(h,sizeof(h),"\"%016llx\"",static_cast<unsigned long long>(v[t][i]));
            s << (i ? ", " : "") << h;
        }
        s << "]" << (t==5 ? "\n" : ",\n");
    }
    s << "  },\n  \"canonical_lossy_round_weights\": [";
    int canonical=0, exact=0, invalid=0;
    std::array<int,3> act{};
    for (int r=0; r<3; ++r) {
        int w=0;
        for (int z=0; z<64; ++z) {
            int a=col(v[2*r],z), b=col(v[2*r+1],z);
            if ((a==0)!=(b==0)) throw std::runtime_error("Incumbent violates lossy activity");
            w += COST[__builtin_popcount(static_cast<unsigned>(a))];
            act[r] += a!=0;
            int n=ddt(a,b);
            if (!n) ++invalid; else exact += 5-__builtin_ctz(static_cast<unsigned>(n));
        }
        canonical += w;
        s << (r ? ", " : "") << w;
    }
    if (!act[0] || canonical > obj+1e-5) throw std::runtime_error("Incumbent objective/nonzero mismatch");
    s << "],\n  \"canonical_lossy_weight\": " << canonical
      << ",\n  \"active_columns\": [" << act[0] << ", " << act[1] << ", " << act[2]
      << "],\n  \"invalid_ddt_columns\": " << invalid
      << ",\n  \"exact_ddt_weight\": " << (invalid ? "null" : std::to_string(exact)) << "\n}\n";
    return s.str();
}
class Progress : public GRBCallback {
    fs::path output;
    std::array<State,6>& states;
    double last=-1;
public:
    std::string failure;
    Progress(fs::path p, std::array<State,6>& v) : output(std::move(p)), states(v) {}
protected:
    void callback() override {
        try {
            if (where==GRB_CB_MIP) {
                double now=getDoubleInfo(GRB_CB_RUNTIME);
                if (last>=0 && now-last<10) return;
                last=now;
                std::ostringstream s;
                s << "{\"status\":\"running\",\"runtime_seconds\":" << number(now)
                  << ",\"lossy_incumbent\":" << number(getDoubleInfo(GRB_CB_MIP_OBJBST))
                  << ",\"lossy_best_bound\":" << number(getDoubleInfo(GRB_CB_MIP_OBJBND))
                  << ",\"nodes\":" << number(getDoubleInfo(GRB_CB_MIP_NODCNT)) << "}\n";
                atomic_write(output/"progress.json",s.str());
            } else if (where==GRB_CB_MIPSOL) {
                Values v{};
                for (int t=0;t<6;++t) for(int i=0;i<5;++i) for(int z=0;z<64;++z)
                    if (getSolution(states[t][i][z])>0.5) v[t][i] |= uint64_t(1)<<z;
                atomic_write(output/"incumbent.json",solution_json(v,getDoubleInfo(GRB_CB_MIPSOL_OBJ)));
            }
        } catch (const GRBException& e) {
            failure="Callback Gurobi error: "+e.getMessage(); abort();
        } catch (const std::exception& e) { failure=e.what(); abort(); }
    }
};
void encoding_json(const fs::path& output, int vars, int constraints, const std::array<State,6>& states) {
    std::ostringstream s;
    s << "{\"variable_count\":" << vars << ",\"constraint_count\":" << constraints << ",\"gates\":[\n";
    for (std::size_t k=0;k<gates.size();++k) {
        const auto& g=gates[k];
        s << "{\"name\":\"" << g.name << "\",\"inputs\":[";
        for(std::size_t j=0;j<g.in.size();++j)
            s << (j ? "," : "") << "\"" << g.in[j].get(GRB_StringAttr_VarName) << "\"";
        s << "],\"output\":\"" << g.out.get(GRB_StringAttr_VarName) << "\"}" << (k+1==gates.size()?"\n":",\n");
    }
    s << "],\"kernel_checks\":[";
    for(std::size_t k=0;k<kernel_checks.size();++k) {
        const auto& check=kernel_checks[k];
        s << (k?",":"") << "{\"name\":\"" << check.name << "\",\"bits\":[";
        for(std::size_t j=0;j<check.bits.size();++j)
            s << (j?",":"") << "\"" << check.bits[j].get(GRB_StringAttr_VarName) << "\"";
        s << "]}";
    }
    s << "],\"state_variables\":[";
    for(int t=0;t<6;++t) {
        s << (t?",":"") << "[";
        for(int i=0;i<5;++i) {
            s << (i?",":"") << "[";
            for(int z=0;z<64;++z) s << (z?",":"") << "\"" << states[t][i][z].get(GRB_StringAttr_VarName) << "\"";
            s << "]";
        }
        s << "]";
    }
    s << "]}\n"; atomic_write(output/"encoding.json",s.str());
}
int main(int argc,char** argv) {
    std::signal(SIGINT,interrupt);
    std::signal(SIGTERM,interrupt);
    try {
        fs::path output;
        int threads=16, middle_active=0, weight_cap=-1;
        bool symmetry=true, export_only=false, use_activity=true, or_hull=true, kernel_both=false, kernel_reduced=false;
        for (int i=1;i<argc;++i) {
            std::string a=argv[i];
            if (a=="--output" && i+1<argc) output=argv[++i];
            else if (a=="--threads" && i+1<argc) threads=std::stoi(argv[++i]);
            else if (a=="--no-symmetry") symmetry=false;
            else if (a=="--export-only") export_only=true;
            else if (a=="--no-activity-cuts") use_activity=false;
            else if (a=="--no-or-hull") or_hull=false;
            else if (a=="--kernel-both") kernel_both=true;
            else if (a=="--kernel-reduced") kernel_reduced=true;
            else if (a=="--middle-active" && i+1<argc) middle_active=std::stoi(argv[++i]);
            else if (a=="--weight-cap" && i+1<argc) weight_cap=std::stoi(argv[++i]);
            else if (a=="--help") {
                std::printf("%s --output NEW_DIR [--threads 1..16] [--no-symmetry] [--export-only] [--no-activity-cuts] [--no-or-hull] [--kernel-both] [--kernel-reduced] [--middle-active 1..64] [--weight-cap NONNEGATIVE]\n",argv[0]); return 0;
            } else throw std::invalid_argument("Unknown or incomplete argument: "+a);
        }
        if (output.empty() || threads<1 || threads>16) throw std::invalid_argument("Require --output and 1<=threads<=16");
        if (middle_active<0 || middle_active>64 || weight_cap < -1)
            throw std::invalid_argument("Invalid middle active count or weight cap");
        if(kernel_reduced && (!kernel_both || middle_active<6 || middle_active>8 || weight_cap!=43 || !use_activity))
            throw std::invalid_argument("Reduced kernel mode requires --kernel-both --middle-active 6..8 --weight-cap 43 and activity cuts");
        if (fs::exists(output)) throw std::invalid_argument("Output directory already exists: "+output.string());
        fs::create_directories(output);
        GRBEnv env;
        env.set(GRB_IntParam_OutputFlag,1);
        env.set(GRB_IntParam_LogToConsole,1);
        GRBModel model(env);
        ModelGuard guard(model);
        model.set(GRB_StringAttr_ModelName,"Gaston_3round_Hsilu_lossy");
        model.set(GRB_IntParam_Threads,threads);
        model.set(GRB_IntParam_MIPFocus,3);
        model.set(GRB_DoubleParam_MIPGap,0);
        model.set(GRB_DoubleParam_MIPGapAbs,0);
        std::array<State,6> states;
        const std::array<std::string,6> labels{"a0","b0","a1","b1","a2","b2"};
        for(int t=0;t<6;++t) if(!kernel_reduced || (t!=2 && t!=4)) states[t]=state(model,labels[t]);
        std::array<Word,2> effects;
        if(kernel_reduced) {
            kernel_linear(model,states[1],states[2],"L0");
            kernel_linear(model,states[3],states[4],"L1");
        } else {
            effects[0]=linear(model,states[1],states[2],"L0");
            effects[1]=linear(model,states[3],states[4],"L1");
        }
        GRBLinExpr obj=0;
        std::array<Chi,3> chis;
        for(int r=0;r<3;++r) chis[r]=lossy(model,states[2*r],states[2*r+1],r,obj,or_hull);
        GRBLinExpr active=0;
        for(int z=0;z<64;++z) active+=chis[0].active[z];
        model.addConstr(active>=1,"nonzero_trail");
        if(symmetry) model.addConstr(chis[0].active[0]==1,"column_shift_representative");
        if(use_activity) activity_cuts(model,chis,effects,kernel_reduced);
        if(kernel_reduced) kernel_type_cuts(model,states,chis,obj,middle_active);
        if(kernel_both && !kernel_reduced) for(int r=0;r<2;++r) for(int z=0;z<64;++z)
            model.addConstr(effects[r][z]==0,"subproblem_kernel_L"+std::to_string(r)+"_"+std::to_string(z));
        if(middle_active) {
            GRBLinExpr middle=0;
            for(int z=0;z<64;++z) middle+=chis[1].active[z];
            model.addConstr(middle==middle_active,"subproblem_middle_active");
        }
        if(weight_cap>=0) model.addConstr(obj<=weight_cap,"subproblem_weight_cap");
        model.setObjective(obj,GRB_MINIMIZE);
        // write() updates pending changes before exporting the auditable model.
        model.write((output/"model.lp").string());
        model.write((output/"model.mps").string());
        model.write((output/"parameters.prm").string());
        encoding_json(output,model.get(GRB_IntAttr_NumVars),model.get(GRB_IntAttr_NumConstrs),states);
        std::ostringstream settings;
        settings << "{\"algorithm\":\"Gaston-320\",\"rounds\":3,\"threads\":" << threads
                 << ",\"time_limit\":null,\"MIPFocus\":3,\"MIPGap\":0,\"MIPGapAbs\":0,\"column_symmetry\":"
                 << (symmetry?"true":"false")
                 << ",\"endpoints\":\"free nonzero XOR trail\",\"sbox_model\":\"Hsilu section 8.1.1 Table 14\","
                 << "\"cost_by_input_hw\":[0,2,3,3,4,4],\"linear_model\":\"binary parity-polytope xor2/xor3\","
                 << "\"extra_branch_cuts\":" << (use_activity?"true":"false")
                 << ",\"or_hull\":" << (or_hull?"true":"false")
                 << ",\"kernel_both\":" << (kernel_both?"true":"false")
                 << ",\"kernel_reduced\":" << (kernel_reduced?"true":"false")
                 << ",\"kernel_audit_sha256\":\"" << KERNEL_AUDIT_SHA256 << "\""
                 << ",\"middle_active\":" << (middle_active?std::to_string(middle_active):"null")
                 << ",\"weight_cap\":" << (weight_cap>=0?std::to_string(weight_cap):"null")
                 << ",\"activity_profile_sha256\":\"" << ACTIVITY_PROFILE_SHA256
                 << "\",\"independent_unsat_certificate\":false}\n";
        atomic_write(output/"settings.json",settings.str());
        std::printf("Model exported: %d variables, %d constraints. Threads=%d; no solve limit.\n",
                    model.get(GRB_IntAttr_NumVars),model.get(GRB_IntAttr_NumConstrs),threads);
        if(export_only) return 0;
        Progress progress(output,states);
        model.setCallback(&progress);
        model.optimize();
        if(!progress.failure.empty()) throw std::runtime_error(progress.failure);
        const int status=model.get(GRB_IntAttr_Status), count=model.get(GRB_IntAttr_SolCount);
        const double runtime=model.get(GRB_DoubleAttr_Runtime);
        std::ostringstream result;
        result << "{\"status_code\":" << status << ",\"solution_count\":" << count
               << ",\"runtime_seconds\":" << number(runtime) << ",\"optimal\":" << (status==GRB_OPTIMAL?"true":"false");
        result << ",\"subproblem_infeasible\":" << (status==GRB_INFEASIBLE?"true":"false")
               << ",\"kernel_both\":" << (kernel_both?"true":"false")
               << ",\"kernel_reduced\":" << (kernel_reduced?"true":"false")
               << ",\"middle_active\":" << (middle_active?std::to_string(middle_active):"null")
               << ",\"weight_cap\":" << (weight_cap>=0?std::to_string(weight_cap):"null");
        if(count>0) {
            const double val=model.get(GRB_DoubleAttr_ObjVal), bound=model.get(GRB_DoubleAttr_ObjBound);
            Values v{};
            for(int t=0;t<6;++t) for(int i=0;i<5;++i) for(int z=0;z<64;++z)
                if(states[t][i][z].get(GRB_DoubleAttr_X)>0.5) v[t][i] |= uint64_t(1)<<z;
            atomic_write(output/"incumbent.json",solution_json(v,val));
            model.write((output/"incumbent.sol").string());
            std::printf("Status=%d (%s), lossy objective=%.10g, best bound=%.10g, gap=%.8g, runtime=%.3f s\n",
                        status,status==GRB_OPTIMAL?"optimal lossy model":"incumbent only",val,bound,model.get(GRB_DoubleAttr_MIPGap),runtime);
            for(int t=0;t<6;++t) {
                std::printf("%s:",labels[t].c_str());
                for(auto w:v[t]) std::printf(" %016llx",static_cast<unsigned long long>(w));
                std::printf("\n");
            }
            result << ",\"lossy_objective\":" << number(val) << ",\"lossy_best_bound\":" << number(bound)
                   << ",\"gap\":" << number(model.get(GRB_DoubleAttr_MIPGap));
        } else std::printf("Status=%d, no incumbent, runtime=%.3f s\n",status,runtime);
        result << "}\n";
        atomic_write(output/"result.json",result.str());
        atomic_write(output/"progress.json",result.str());
        return 0;
    } catch(const GRBException& e) {
        std::fprintf(stderr,"Gurobi error %d: %s\n",e.getErrorCode(),e.getMessage().c_str()); return 1;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Error: %s\n",e.what()); return 1;
    }
}

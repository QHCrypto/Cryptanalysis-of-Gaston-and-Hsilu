// Reconstructed, value-sensitive structural MILP for the two-round Hsilu path.
// C/L/N tags are conservative degree bounds, not polynomial coefficients.
// The model does NOT certify condition rank, output rank, or sponge reachability.
#include <gurobi_c++.h>
#include <algorithm>
#include <array>
#include <csignal>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int W = 64, R = 5;
constexpr std::array<int,13> N1{0,3,1,0,5,5,7,11,12,12,13,14,15};
constexpr std::array<int,13> N2{4,4,2,1,2,3,6,8,9,6,10,11,9};
constexpr std::array<int,13> ROT{47,51,54,30,4,55,0,14,59,4,10,45,53};
constexpr std::array<int,5> SHIFT{0,57,13,31,20};
constexpr std::array<unsigned,2> RC{0x5a,0x4b};
using Grid = std::array<std::array<int,W>,R>;
GRBModel* active_model = nullptr;
void interrupt(int) { if (active_model) active_model->terminate(); }
struct ActiveGuard {
    explicit ActiveGuard(GRBModel& m) { active_model = &m; }
    ~ActiveGuard() { active_model = nullptr; }
};

struct Path {
    std::array<int,W> free{}, selected{};
    std::array<std::array<Grid,2>,2> condition{}; // [round][value][row][bit]
};
struct Options {
    int digest = 128, threads = 0;
    std::string start = "hsilu_2round_best_attack.txt";
    std::string output = "hsilu_2round_rebuilt";
    std::string lp;
    bool use_start = true, fixed = false;
};

std::string field(const std::string& text, const std::string& key) {
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines,line))
        if (line.rfind(key+": ",0)==0) return line.substr(key.size()+2);
    throw std::runtime_error("Missing field: "+key);
}
std::vector<int> numbers(const std::string& s) {
    std::vector<int> result;
    const std::regex number("[0-9]+");
    for (std::sregex_iterator it(s.begin(),s.end(),number), end; it!=end; ++it)
        result.push_back(std::stoi(it->str()));
    return result;
}
void read_rows(const std::string& text, Grid& out) {
    const std::regex row("([0-9]+)\\s*:\\s*\\[([^\\]]*)\\]");
    for (std::sregex_iterator it(text.begin(),text.end(),row), end; it!=end; ++it) {
        const int i = std::stoi((*it)[1]);
        if (i<0 || i>=R) throw std::runtime_error("Row index outside 0..4");
        for (int j:numbers((*it)[2])) {
            if (j<0 || j>=W || out[i][j]) throw std::runtime_error("Invalid/duplicate bit index");
            out[i][j]=1;
        }
    }
}
Path load_path(const std::string& filename) {
    std::ifstream f(filename);
    if (!f) throw std::runtime_error("Cannot read start path: "+filename);
    std::string text((std::istreambuf_iterator<char>(f)),{});
    if (std::stoi(field(text,"rounds"))!=2 || field(text,"algorithm").find("Hsilu")!=0)
        throw std::runtime_error("Only two-round Hsilu paths are supported");
    if (numbers(field(text,"rate_lanes"))!=std::vector<int>{0})
        throw std::runtime_error("This model requires rate_lanes=[0]");
    Path p;
    Grid free{}, selected{};
    read_rows(field(text,"free_variables_by_lane"),free);
    read_rows(field(text,"selected_outputs_by_lane"),selected);
    for (int i=1;i<R;++i) for(int j=0;j<W;++j)
        if (free[i][j] || selected[i][j]) throw std::runtime_error("Only row 0 can be active/selected");
    p.free=free[0]; p.selected=selected[0];
    const std::regex stage("'X([12])'\\s*:\\s*\\{([^}]*)\\}");
    for (int value=0;value<2;++value) {
        const std::string s=field(text,"conditions_eq_"+std::to_string(value));
        for (std::sregex_iterator it(s.begin(),s.end(),stage), end;it!=end;++it)
            read_rows((*it)[2],p.condition[std::stoi((*it)[1])-1][value]);
    }
    return p;
}
int count(const std::array<int,W>& a) {
    int result=0; for(int v:a) result+=v; return result;
}
int conditions(const Path& p,int round=-1) {
    int result=0;
    for(int r=0;r<2;++r) if(round<0 || round==r)
        for(int b=0;b<2;++b) for(const auto& row:p.condition[r][b]) result+=count(row);
    return result;
}
bool rc_bit(int round,int row,int j) {
    return row==0 && j>=56 && ((RC[round]>>(63-j))&1U);
}
// Path numbering is MSB=0. A left rotation reads input column j+s.
Grid mix(const Grid& input) {
    std::vector<std::array<int,W>> x(input.begin(),input.end());
    for(int k=0;k<13;++k) {
        std::array<int,W> row{};
        for(int j=0;j<W;++j) row[j]=std::max(x[N1[k]][j],x[N2[k]][(j+ROT[k])%W]);
        x.push_back(row);
    }
    Grid out{};
    for(int i=0;i<R;++i) for(int j=0;j<W;++j) out[i][j]=x[13+i][(j+SHIFT[i])%W];
    return out;
}
struct Replay { std::array<Grid,3> x{}; std::array<Grid,2> y{}; };
Replay replay(const Path& p) {
    Replay result;
    result.x[0][0]=p.free;
    if(count(p.free)!=count(p.selected) || count(p.free)<1 || count(p.free)+conditions(p)>64)
        throw std::runtime_error("Path violates selected=free>=1 or free+conditions<=64");
    for(int r=0;r<2;++r) {
        for(int i=0;i<R;++i) for(int j=0;j<W;++j) {
            const int q=p.condition[r][0][i][j]+p.condition[r][1][i][j];
            if(q>1 || (q && result.x[r][i][j]!=0))
                throw std::runtime_error("A value condition is allowed only on a constant tag");
        }
        for(int i=0;i<R;++i) for(int j=0;j<W;++j) {
            const int b=(i+1)%R,c=(i+2)%R;
            const bool b_one=p.condition[r][rc_bit(r,b,j)?0:1][b][j];
            const bool c_zero=p.condition[r][rc_bit(r,c,j)?1:0][c][j];
            const int product=(b_one||c_zero)?0:std::min(2,result.x[r][b][j]+result.x[r][c][j]);
            result.y[r][i][j]=std::max(result.x[r][i][j],product);
        }
        result.x[r+1]=mix(result.y[r]);
    }
    for(int j=0;j<W;++j) if(p.selected[j] && result.x[2][0][j]!=1)
        throw std::runtime_error("Selected output is not tagged linear");
    return result;
}

struct Atom { GRBVar var; int seed=-1; int fixed=-1; };
struct Literal { Atom a; bool neg=false; };
struct Node { Atom active, nonlinear; };
using ModelRow=std::array<Node,W>;
using ModelGrid=std::array<ModelRow,R>;
using ConditionGrid=std::array<std::array<Atom,W>,R>;
class Builder {
public:
    GRBModel& model;
    bool seed_enabled;
    Atom zero,one;
    explicit Builder(GRBModel& m,bool use_seed):model(m),seed_enabled(use_seed) {
        zero=atom("zero",0,0,0); one=atom("one",1,1,1);
    }
    Atom atom(const std::string& name,int seed=-1,int lb=0,int ub=1) {
        Atom a{model.addVar(lb,ub,0,GRB_BINARY,name),seed,lb==ub?lb:-1};
        if(seed_enabled && seed>=0) a.var.set(GRB_DoubleAttr_Start,seed);
        return a;
    }
    Atom gate(const std::string& name,std::vector<Literal> inputs,bool conjunction) {
        std::vector<Literal> args;
        for(const auto& l:inputs) {
            if(l.a.fixed>=0) {
                const int v=l.a.fixed ^ int(l.neg);
                if(v==int(!conjunction)) return conjunction?zero:one;
            } else args.push_back(l);
        }
        if(args.empty()) return conjunction?one:zero;
        if(args.size()==1 && !args[0].neg) return args[0].a;
        int seed=conjunction?1:0;
        for(const auto& l:args) {
            if(l.a.seed<0) {seed=-1;break;}
            const int v=l.a.seed ^ int(l.neg);
            seed=conjunction?(seed&v):(seed|v);
        }
        Atom out=atom(name,seed);
        GRBLinExpr sum=0;
        for(size_t k=0;k<args.size();++k) {
            GRBLinExpr e=args[k].neg?1-GRBLinExpr(args[k].a.var):GRBLinExpr(args[k].a.var);
            sum+=e;
            if(conjunction) model.addConstr(out.var<=e,name+"_le_"+std::to_string(k));
            else model.addConstr(out.var>=e,name+"_ge_"+std::to_string(k));
        }
        if(conjunction) model.addConstr(out.var>=sum-int(args.size())+1,name+"_lower");
        else model.addConstr(out.var<=sum,name+"_upper");
        return out;
    }
    Node xor_nodes(const Node& a,const Node& b,const std::string& n) {
        return {gate(n+"_A",{{a.active},{b.active}},false),
                gate(n+"_N",{{a.nonlinear},{b.nonlinear}},false)};
    }
    ModelGrid linear(const ModelGrid& input,int round) {
        std::vector<ModelRow> x(input.begin(),input.end());
        for(int k=0;k<13;++k) {
            ModelRow row;
            for(int j=0;j<W;++j)
                row[j]=xor_nodes(x[N1[k]][j],x[N2[k]][(j+ROT[k])%W],
                    "M"+std::to_string(round)+"_"+std::to_string(k+5)+"_"+std::to_string(j));
            x.push_back(row);
        }
        ModelGrid out;
        for(int i=0;i<R;++i) for(int j=0;j<W;++j) out[i][j]=x[13+i][(j+SHIFT[i])%W];
        return out;
    }
};
int value(const Atom& a) {return a.var.get(GRB_DoubleAttr_X)>0.5;}
int degree(const Node& n) {return value(n.active)+value(n.nonlinear);}

std::string bit_list(const std::array<int,W>& a) {
    std::ostringstream s; s<<"["; bool first=true;
    for(int j=0;j<W;++j) if(a[j]) {if(!first)s<<", ";s<<j;first=false;}
    s<<"]"; return s.str();
}
std::string condition_map(const Path& p,int value) {
    std::ostringstream s;s<<"{";bool first_stage=true;
    for(int r=0;r<2;++r) {
        int total=0;for(const auto& row:p.condition[r][value]) total+=count(row);
        if(!total)continue;
        if(!first_stage)s<<", ";
        first_stage=false;
        s<<"'X"<<r+1<<"': {";bool first_row=true;
        for(int i=0;i<R;++i) if(count(p.condition[r][value][i])) {
            if(!first_row)s<<", ";
            first_row=false;
            s<<i<<": "<<bit_list(p.condition[r][value][i]);
        }
        s<<"}";
    }
    s<<"}";return s.str();
}
std::string indexed(const std::array<int,13>& a) {
    std::ostringstream s;s<<"{";
    for(int k=0;k<13;++k) {if(k)s<<", ";s<<k+5<<": "<<a[k];}
    s<<"}";return s.str();
}
void export_text(const Options& opt,const Path& p,const Replay& r,GRBModel& m) {
    std::ofstream f(opt.output+".txt");
    if(!f)throw std::runtime_error("Cannot write result file");
    const int F=count(p.free),L=count(p.selected),C=conditions(p),status=m.get(GRB_IntAttr_Status);
    std::array<int,W> guessed{},linear{};
    for(int j=0;j<W;++j) {guessed[j]=1-p.free[j];linear[j]=r.x[2][0][j]==1;}
    f<<"algorithm: Hsilu\nrounds: 2\nmode: "<<(opt.fixed?"check-reference":"best")
     <<"\nstatus: "<<status<<"\nsolver_optimal: "<<(status==GRB_OPTIMAL?"True":"False")
     <<"\noptimality_scope: "<<(opt.fixed?"fixed supplied decisions only":"reconstructed structural model only")
     <<"\nruntime_seconds: "<<m.get(GRB_DoubleAttr_Runtime)
     <<"\nsolution_count: "<<m.get(GRB_IntAttr_SolCount)
     <<"\ncomplexity_exponent: "<<opt.digest-L
     <<"\nobjective: "<<m.get(GRB_DoubleAttr_ObjVal)
     <<"\nobjective_bound: "<<m.get(GRB_DoubleAttr_ObjBound)
     <<"\ndigest_bits: "<<opt.digest<<"\nrate_lanes: [0]\noutput_bits: 64"
     <<"\nguessed_count: "<<W-F<<"\nfree_variable_count: "<<F
     <<"\ncondition_count: "<<C<<"\ncondition_counts: {'X1': "<<conditions(p,0)<<", 'X2': "<<conditions(p,1)<<"}"
     <<"\ncondition_weights: {'X1': 1.0, 'X2': 1.0}"
     <<"\nweighted_condition_count: "<<C<<"\nresource_usage: "<<F+C<<"\nresource_bound: 64"
     <<"\nselected_count: "<<L<<"\nguessed_by_lane: {0: "<<bit_list(guessed)<<"}"
     <<"\nfree_variables_by_lane: {0: "<<bit_list(p.free)<<"}"
     <<"\nselected_outputs_by_lane: {0: "<<bit_list(p.selected)<<"}"
     <<"\nfinal_linearizable_x3_count: "<<count(linear)
     <<"\nfinal_linearizable_x3_by_lane: {0: "<<bit_list(linear)<<"}"
     <<"\nconditions_eq_0: "<<condition_map(p,0)<<"\nconditions_eq_1: "<<condition_map(p,1)
     <<"\nmrm_parameters: {'n1': "<<indexed(N1)<<", 'n2': "<<indexed(N2)<<", 'm': "<<indexed(ROT)<<", 's': [0, 57, 13, 31, 20]}"
     <<"\nbit_numbering: MSB=0; path bit j is native integer bit 63-j"
     <<"\nround_constants: [90, 75]"
     <<"\ncondition_eligibility: constant tags only; actual consistency and independence need algebraic validation"
     <<"\nvalidation: {'structural_replay': 'passed', 'equation_rank': 'not_checked', 'concrete_condition_consistency': 'not_checked', 'real_preimage': 'not_computed'}\n";
    if(!f)throw std::runtime_error("Result write failed");
}
void export_svg(const Options& opt,const Path& p,const Replay& r) {
    std::ofstream f(opt.output+".svg");
    if(!f)throw std::runtime_error("Cannot write SVG");
    f<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"980\" height=\"590\" viewBox=\"0 0 980 590\">\n"
     <<"<rect width=\"100%\" height=\"100%\" fill=\"white\"/><style>text{font-family:monospace;font-size:10px;fill:#222}</style>\n"
     <<"<text x=\"20\" y=\"24\" style=\"font-size:17px\">Hsilu 2-round structural path; "<<opt.digest<<"-bit target</text>\n"
     <<"<text x=\"20\" y=\"43\">MSB=0 | gray: constant | yellow: affine | white: nonlinear | red/green: condition 0/1 | blue: selected</text>\n";
    for(int panel=0;panel<5;++panel) {
        const bool is_x=panel%2==0;const int rindex=panel/2;
        const Grid& grid=is_x?r.x[rindex]:r.y[rindex];
        const std::string label=std::string(is_x?"X":"Y")+std::to_string(rindex+1);
        const int top=78+panel*98;
        f<<"<text x=\"20\" y=\""<<top<<"\" style=\"font-size:15px\">"<<label<<"</text>\n";
        for(int j=0;j<W;++j)if(j%4==0)f<<"<text x=\""<<70+14*j<<"\" y=\""<<top-5<<"\">"<<j<<"</text>\n";
        for(int i=0;i<R;++i) {
            f<<"<text x=\"48\" y=\""<<top+13*i+11<<"\">"<<i<<"</text>\n";
            for(int j=0;j<W;++j) {
                std::string color=grid[i][j]==0?"#b8b8b8":grid[i][j]==1?"#ffd84d":"#ffffff";
                if(is_x && rindex<2) {
                    if(p.condition[rindex][0][i][j])color="#d73a49";
                    if(p.condition[rindex][1][i][j])color="#2da44e";
                }
                if(panel==4 && i==0 && p.selected[j])color="#8ec5ff";
                f<<"<rect x=\""<<70+14*j<<"\" y=\""<<top+13*i<<"\" width=\"13\" height=\"12\" fill=\""<<color<<"\" stroke=\"#555\" stroke-width=\"0.4\"><title>"<<label<<"["<<i<<"]["<<j<<"]</title></rect>\n";
            }
        }
    }
    f<<"<text x=\"20\" y=\"582\">free="<<count(p.free)<<", conditions="<<conditions(p)<<", selected="<<count(p.selected)<<"; rank and full XOF feasibility not certified</text></svg>\n";
}
void help() {
    std::printf("Usage: hsilu_2round_milp [options]\n"
        "  --digest-bits N         target length, default 128 (N >= 64)\n"
        "  --start FILE            optional path used only as a MIP start\n"
        "  --no-start              search without the supplied incumbent\n"
        "  --check-reference FILE  fix only path decisions and verify the model\n"
        "  --output PREFIX         write PREFIX.txt and PREFIX.svg\n"
        "  --threads N             Gurobi thread count, otherwise solver default\n"
        "  --write-lp FILE         optionally export the model\n"
        "No default time limit. Ctrl+C requests termination and exports an incumbent.\n");
}
Options options(int argc,char** argv) {
    Options o;
    for(int i=1;i<argc;++i) {
        const std::string a=argv[i];
        if(a=="--no-start") {o.use_start=false;continue;}
        if(i+1>=argc)throw std::runtime_error("Missing option value: "+a);
        const std::string b=argv[++i];
        if(a=="--digest-bits")o.digest=std::stoi(b);
        else if(a=="--start") {o.start=b;o.use_start=true;}
        else if(a=="--check-reference") {o.start=b;o.use_start=true;o.fixed=true;}
        else if(a=="--output")o.output=b;
        else if(a=="--threads")o.threads=std::stoi(b);
        else if(a=="--write-lp")o.lp=b;
        else throw std::runtime_error("Unknown option: "+a);
    }
    if(o.digest<64 || o.threads<0)throw std::runtime_error("Invalid digest length or thread count");
    if(o.fixed && !o.use_start)throw std::runtime_error("--check-reference needs its input file");
    if(std::filesystem::path(o.output).filename()=="hsilu_2round_best_attack")
        throw std::runtime_error("Refusing to overwrite the supplied original path");
    if(o.use_start && std::filesystem::weakly_canonical(o.output+".txt")==std::filesystem::weakly_canonical(o.start))
        throw std::runtime_error("Output must differ from the input path");
    return o;
}
} // namespace

int main(int argc,char** argv) {
    std::signal(SIGINT,interrupt);
    try {
        if(argc==2 && std::string(argv[1])=="--help") {help();return 0;}
        const Options opt=options(argc,argv);
        Path seed;
        if(opt.use_start) {seed=load_path(opt.start);replay(seed);}
        GRBEnv env;
        env.set(GRB_IntParam_OutputFlag,1);
        env.set(GRB_IntParam_LogToConsole,1);
        GRBModel model(env);
        if(opt.threads)model.set(GRB_IntParam_Threads,opt.threads);
        Builder b(model,opt.use_start);
        std::array<ModelGrid,3> x;
        std::array<ModelGrid,2> y;
        std::array<std::array<ConditionGrid,2>,2> q;
        std::array<Atom,W> chosen;
        GRBLinExpr free_sum=0,condition_sum=0,selected_sum=0;
        GRBLinExpr condition_sum_x1=0,condition_sum_x2=0;
        for(int i=0;i<R;++i)for(int j=0;j<W;++j) {
            const Atom a=i==0?b.atom("free_"+std::to_string(j),opt.use_start?seed.free[j]:-1):b.zero;
            x[0][i][j]={a,b.zero};
            if(i==0) {
                free_sum+=a.var;
                if(opt.fixed)model.addConstr(a.var==seed.free[j]);
            }
        }
        for(int round=0;round<2;++round) {
            for(int i=0;i<R;++i)for(int j=0;j<W;++j) {
                for(int v=0;v<2;++v) {
                    const std::string name="q"+std::to_string(round+1)+"_"+std::to_string(i)+"_"+std::to_string(j)+"_"+std::to_string(v);
                    q[round][v][i][j]=b.atom(name,opt.use_start?seed.condition[round][v][i][j]:-1);
                    condition_sum+=q[round][v][i][j].var;
                    if(round==0)condition_sum_x1+=q[round][v][i][j].var;
                    else condition_sum_x2+=q[round][v][i][j].var;
                    if(opt.fixed)model.addConstr(q[round][v][i][j].var==seed.condition[round][v][i][j]);
                }
                model.addConstr(q[round][0][i][j].var+q[round][1][i][j].var+x[round][i][j].active.var<=1,
                    "condition_constant_"+std::to_string(round)+"_"+std::to_string(i)+"_"+std::to_string(j));
            }
            for(int i=0;i<R;++i)for(int j=0;j<W;++j) {
                const int ib=(i+1)%R,ic=(i+2)%R;
                const Node& a=x[round][i][j];const Node& nb=x[round][ib][j];const Node& nc=x[round][ic][j];
                const Atom cb_one=q[round][rc_bit(round,ib,j)?0:1][ib][j];
                const Atom cc_zero=q[round][rc_bit(round,ic,j)?1:0][ic][j];
                const std::string n="Y"+std::to_string(round+1)+"_"+std::to_string(i)+"_"+std::to_string(j);
                const Atom kill=b.gate(n+"_kill",{{cb_one},{cc_zero}},false);
                const Atom any=b.gate(n+"_any",{{nb.active},{nc.active}},false);
                const Atom both=b.gate(n+"_both",{{nb.active},{nc.active}},true);
                const Atom raw_bad=b.gate(n+"_rawbad",{{nb.nonlinear},{nc.nonlinear},{both}},false);
                const Atom pa=b.gate(n+"_productA",{{any},{kill,true}},true);
                const Atom pn=b.gate(n+"_productN",{{raw_bad},{kill,true}},true);
                y[round][i][j]={b.gate(n+"_A",{{a.active},{pa}},false),b.gate(n+"_N",{{a.nonlinear},{pn}},false)};
            }
            x[round+1]=b.linear(y[round],round+1);
        }
        for(int j=0;j<W;++j) {
            chosen[j]=b.atom("selected_"+std::to_string(j),opt.use_start?seed.selected[j]:-1);
            selected_sum+=chosen[j].var;
            model.addConstr(chosen[j].var<=x[2][0][j].active.var);
            model.addConstr(chosen[j].var+x[2][0][j].nonlinear.var<=1);
            if(opt.fixed)model.addConstr(chosen[j].var==seed.selected[j]);
        }
        model.addConstr(selected_sum==free_sum,"selected_equals_free");
        model.addConstr(free_sum>=1,"nonempty_free_space");
        model.addConstr(free_sum+condition_sum_x2<=64,"condition_resource_budget2");
        model.addConstr(free_sum+condition_sum_x1<=64,"condition_resource_budget1");
        model.setObjective(opt.digest-selected_sum,GRB_MINIMIZE);
        if(!opt.lp.empty())model.write(opt.lp);
        std::printf("Model: two rounds, MSB=0, constant-only value conditions, selected=free, free+conditions<=64.\n");
        std::printf("Target=%d bits; %s. No rank/full-XOF complexity certificate is implied.\n",opt.digest,opt.fixed?"fixed-path feasibility check":"unrestricted structural search");
        {
            ActiveGuard guard(model);
            model.optimize();
        }
        if(model.get(GRB_IntAttr_SolCount)==0)return 2;
        Path found;
        for(int j=0;j<W;++j) {found.free[j]=value(x[0][0][j].active);found.selected[j]=value(chosen[j]);}
        for(int r=0;r<2;++r)for(int v=0;v<2;++v)for(int i=0;i<R;++i)for(int j=0;j<W;++j)
            found.condition[r][v][i][j]=value(q[r][v][i][j]);
        const Replay checked=replay(found);
        for(int r=0;r<3;++r)for(int i=0;i<R;++i)for(int j=0;j<W;++j)
            if(degree(x[r][i][j])!=checked.x[r][i][j])throw std::runtime_error("X state replay mismatch");
        for(int r=0;r<2;++r)for(int i=0;i<R;++i)for(int j=0;j<W;++j)
            if(degree(y[r][i][j])!=checked.y[r][i][j])throw std::runtime_error("Y state replay mismatch");
        export_text(opt,found,checked,model);export_svg(opt,found,checked);
        const int status=model.get(GRB_IntAttr_Status);
        std::printf("\nStatus: %d; optimal: %s; scope: %s\n",status,status==GRB_OPTIMAL?"yes":"no",opt.fixed?"fixed input decisions":"this structural model");
        std::printf("Objective exponent: %.10g; bound: %.10g; gap: %.6g; runtime: %.3f s\n",model.get(GRB_DoubleAttr_ObjVal),model.get(GRB_DoubleAttr_ObjBound),model.get(GRB_DoubleAttr_MIPGap),model.get(GRB_DoubleAttr_Runtime));
        std::printf("Free=%d; guessed=%d; selected=%d; conditions=%d (X1=%d, X2=%d); resource=%d/64\n",count(found.free),64-count(found.free),count(found.selected),conditions(found),conditions(found,0),conditions(found,1),count(found.free)+conditions(found));
        std::printf("Free positions: %s\nSelected outputs: %s\n",bit_list(found.free).c_str(),bit_list(found.selected).c_str());
        std::printf("Conditions = 0: %s\nConditions = 1: %s\n",condition_map(found,0).c_str(),condition_map(found,1).c_str());
        std::printf("Structural replay passed. Results: %s.txt and %s.svg\n",opt.output.c_str(),opt.output.c_str());
        return 0;
    } catch(const GRBException& e) {
        std::fprintf(stderr,"Gurobi error %d: %s\n",e.getErrorCode(),e.getMessage().c_str());return 1;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Error: %s\n",e.what());return 1;
    }
}

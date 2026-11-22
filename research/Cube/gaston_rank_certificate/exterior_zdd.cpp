// Exact highest-degree propagation in the square-zero polynomial algebra.
// A ZDD stores the support of a homogeneous GF(2) polynomial. Addition is
// symmetric difference; products containing repeated variables vanish.
// Resource-limit exits produce no coefficient or rank claim.
#include <algorithm>
#include <array>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using ID=uint32_t;
using Word=uint64_t;
struct Node { ID lo,hi; unsigned var,degree; };
struct Triple {
    ID lo,hi;unsigned var;
    bool operator==(const Triple& b)const{return lo==b.lo&&hi==b.hi&&var==b.var;}
};
struct Hash {
    size_t operator()(const Triple& k)const {
        Word x=(Word(k.lo)<<32)|k.hi;
        x^=Word(k.var)*0x9e3779b97f4a7c15ULL;
        x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;
        return x^(x>>31);
    }
};
class ZDD {
public:
    unsigned variables;
    size_t node_limit,cache_limit;
    double seconds_limit;
    std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    std::vector<Node> nodes{{0,0,100,0},{0,0,100,0}};
    std::unordered_map<Triple,ID,Hash> unique;
    std::unordered_map<Word,ID> sums,products;
    Word calls=0;
    Word collections=0;
    size_t peak_nodes=2;
    ZDD(unsigned n,size_t limit,double sec):variables(n),node_limit(limit),
       cache_limit(limit*2),seconds_limit(sec) {}
    double elapsed()const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    void tick() {
        if((++calls&65535)==0&&elapsed()>seconds_limit)throw std::runtime_error("time_limit");
    }
    ID make(unsigned v,ID lo,ID hi) {
        if(!hi)return lo;
        Triple key{lo,hi,v};auto it=unique.find(key);
        if(it!=unique.end())return it->second;
        if(nodes.size()>=node_limit)throw std::runtime_error("node_limit");
        unsigned d=nodes[hi].degree+1;
        if(lo&&nodes[lo].degree!=d)throw std::runtime_error("inhomogeneous_node");
        ID id=nodes.size();nodes.push_back({lo,hi,v,d});unique.emplace(key,id);
        peak_nodes=std::max(peak_nodes,nodes.size());
        return id;
    }
    void collect(const std::vector<ID*>& roots) {
        std::vector<ID> mapping(nodes.size(),UINT32_MAX);mapping[0]=0;mapping[1]=1;
        std::vector<Node> compact{nodes[0],nodes[1]};
        std::function<ID(ID)> copy=[&](ID a)->ID {
            if(mapping[a]!=UINT32_MAX)return mapping[a];
            Node x=nodes[a];x.lo=copy(x.lo);x.hi=copy(x.hi);
            ID b=compact.size();compact.push_back(x);mapping[a]=b;return b;
        };
        for(auto p:roots)*p=copy(*p);
        nodes.swap(compact);unique.clear();sums.clear();products.clear();
        for(ID a=2;a<nodes.size();++a)unique.emplace(Triple{nodes[a].lo,nodes[a].hi,nodes[a].var},a);
        ++collections;
    }
    ID add(ID a,ID b) {
        tick();if(a==b)return 0;if(!a)return b;if(!b)return a;
        if(nodes[a].degree!=nodes[b].degree)throw std::runtime_error("inhomogeneous_add");
        if(a>b)std::swap(a,b);
        Word key=(Word(a)<<32)|b;auto it=sums.find(key);if(it!=sums.end())return it->second;
        Node x=nodes[a],y=nodes[b];unsigned v=std::min(x.var,y.var);
        ID a0=x.var==v?x.lo:a,a1=x.var==v?x.hi:0;
        ID b0=y.var==v?y.lo:b,b1=y.var==v?y.hi:0;
        ID lo=add(a0,b0),hi=add(a1,b1),out=make(v,lo,hi);
        if(sums.size()>=cache_limit)sums.clear();
        sums.emplace(key,out);return out;
    }
    ID multiply(ID a,ID b) {
        tick();if(!a||!b)return 0;if(a==1)return b;if(b==1)return a;
        if(a==b)return 0; // Every positive-degree element squares to zero.
        if(a>b)std::swap(a,b);
        Node x=nodes[a],y=nodes[b];unsigned v=std::min(x.var,y.var);
        if(x.degree+y.degree>variables-v)return 0;
        Word key=(Word(a)<<32)|b;
        auto it=products.find(key);if(it!=products.end())return it->second;
        ID a0=x.var==v?x.lo:a,a1=x.var==v?x.hi:0;
        ID b0=y.var==v?y.lo:b,b1=y.var==v?y.hi:0;
        ID lo=multiply(a0,b0),u=multiply(a0,b1),w=multiply(a1,b0);
        ID out=make(v,lo,add(u,w));
        if(products.size()>=cache_limit)products.clear();
        products.emplace(key,out);return out;
    }
    unsigned parity_with_key_variable(ID root) {
        std::unordered_map<ID,unsigned> memo;
        std::function<unsigned(ID)> f=[&](ID a)->unsigned {
            if(a<2)return a;
            auto it=memo.find(a);if(it!=memo.end())return it->second;
            Node x=nodes[a];
            unsigned bit=f(x.lo)^f(x.hi);
            memo.emplace(a,bit);return bit;
        };
        return root>=2&&nodes[root].var==0?f(nodes[root].hi):0;
    }
};
constexpr unsigned E[5]={0,60,22,27,4},W[5]={0,56,31,46,43},T[5]={25,32,52,60,63};
using State=std::array<std::vector<ID>,5>;
State linear(ZDD& z,const State& s,unsigned n) {
    State a,b;for(auto& x:a)x.resize(n);for(auto& x:b)x.resize(n);
    std::vector<ID> p(n),q(n),effect(n);
    auto back=[&](unsigned j,unsigned shift){return (j+n-shift%n)%n;};
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<n;++j)a[i][j]=s[i][back(j,E[i])];
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<n;++j) {
        p[j]=z.add(p[j],a[i][j]);q[j]=z.add(q[j],a[i][back(j,T[i])]);
    }
    for(unsigned j=0;j<n;++j)
        effect[j]=z.add(z.add(p[back(j,23)],p[back(j,24)]),
                        z.add(q[back(j,23)],q[back(j,41)]));
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<n;++j)
        b[i][j]=z.add(a[i][back(j,W[i])],effect[back(j,W[i])]);
    return b;
}
int main(int argc,char** argv) {
    if(argc!=7) {
        std::fprintf(stderr,"usage: %s rounds key_word node_limit seconds low_terms output.json\n",argv[0]);
        return 2;
    }
    unsigned rounds=std::stoul(argv[1]),key_word=std::stoul(argv[2]),m=std::stoul(argv[5]);
    if(rounds<1||rounds>6||key_word>1||m<1||m>(1u<<rounds))return 2;
    unsigned n=1u<<rounds,completed=0;std::string stage="initial",failure;
    ZDD z(n+1,std::stoull(argv[3]),std::stod(argv[4]));
    std::vector<std::array<unsigned,2>> bits;
    try {
        State s;for(auto& x:s)x.resize(n);
        // Put the key first so the common key-free polynomial is shared.
        for(unsigned j=0;j<n;++j)s[3][j]=z.make(j+1,0,1);
        s[1+key_word][0]=z.make(0,0,1);
        for(unsigned r=0;r<rounds;++r) {
            stage="linear_"+std::to_string(r+1);State a=linear(z,s,n);
            std::vector<ID*> roots;for(auto& word:a)for(auto& v:word)roots.push_back(&v);
            z.collect(roots);
            std::printf("%s nodes=%zu seconds=%.3f\n",stage.c_str(),z.nodes.size(),z.elapsed());std::fflush(stdout);
            stage="quadratic_"+std::to_string(r+1);
            for(auto& word:s)std::fill(word.begin(),word.end(),0);
            for(auto& word:s)for(auto& v:word)roots.push_back(&v);
            for(unsigned i=0;i<(r+1==rounds?2u:5u);++i)for(unsigned j=0;j<n;++j) {
                s[i][j]=z.multiply(a[(i+1)%5][j],a[(i+2)%5][j]);
                if(j%8==7)z.collect(roots);
            }
            completed=r+1;
            std::printf("%s nodes=%zu seconds=%.3f\n",stage.c_str(),z.nodes.size(),z.elapsed());std::fflush(stdout);
        }
        stage="projection";
        for(unsigned i=0;i<m;++i) {
            std::array<unsigned,2> row{};
            for(unsigned w=0;w<2;++w) {
                ID sum=0;for(unsigned j=0;j<n;++j)if((j&i)==i)sum=z.add(sum,s[w][j]);
                // Degree n in n+1 variables: monomials containing the key
                // have exactly n-1 nonce variables. Summing their coefficients
                // is the restriction to the even-parity nonce hyperplane.
                row[w]=z.parity_with_key_variable(sum);
            }
            bits.push_back(row);
        }
    }catch(const std::exception& e){failure=e.what();}
    FILE* out=std::fopen(argv[6],"w");if(!out)return 3;
    std::fprintf(out,"{\n  \"rounds\": %u,\n  \"word_bits\": %u,\n  \"key_word\": %u,\n"
                 "  \"key_variable_order\": \"first\",\n"
                 "  \"completed_rounds\": %u,\n  \"stage\": \"%s\",\n  \"nodes\": %zu,\n"
                 "  \"peak_nodes\": %zu,\n  \"collections\": %" PRIu64 ",\n"
                 "  \"recursive_calls\": %" PRIu64 ",\n  \"seconds\": %.9f,\n"
                 "  \"coefficient_certificate\": %s,\n  \"failure\": \"%s\",\n  \"low_coefficients\": ",
                 rounds,n,key_word,completed,stage.c_str(),z.nodes.size(),z.peak_nodes,z.collections,z.calls,z.elapsed(),
                 failure.empty()?"true":"false",failure.c_str());
    if(!failure.empty())std::fprintf(out,"null");
    else {
        std::fprintf(out,"[");
        for(unsigned i=0;i<bits.size();++i)std::fprintf(out,"%s[%u,%u]",i?",":"",bits[i][0],bits[i][1]);
        std::fprintf(out,"]");
    }
    std::fprintf(out,"\n}\n");std::fclose(out);
    std::printf("certificate=%s reason=%s\n",failure.empty()?"yes":"no",failure.c_str());
}

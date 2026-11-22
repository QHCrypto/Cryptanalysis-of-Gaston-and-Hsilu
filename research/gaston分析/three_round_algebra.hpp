#pragma once
#include <array>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gaston3 {
using Word = std::uint64_t;
using State = std::array<Word, 5>;
constexpr std::array<int, 5> east = {0,60,22,27,4};
constexpr std::array<int, 5> west = {0,56,31,46,43};
constexpr std::array<int, 5> twin = {25,32,52,60,63};
constexpr std::array<Word, 3> rc = {0x69,0x5a,0x4b};
inline Word rot(Word v, int n) { n &= 63; return n ? (v<<n)|(v>>(64-n)) : v; }
inline bool parity(Word x) { return __builtin_parityll(x); }
inline int index(int j) { return (-j)&63; }
inline bool bit(Word x,int j) { return (x>>index(j))&1; }
inline Word urow(int j) {
    Word r=0; for (int d : {23,24,48,2}) r ^= Word(1)<<index(j+d); return r;
}
inline State linear(State x) {
    State d{}; Word p=0,q=0;
    for(int i=0;i<5;++i) { d[i]=rot(x[i],east[i]); p^=d[i]; q^=rot(d[i],twin[i]); }
    Word shared=rot(p^rot(p,1)^q^rot(q,18),23);
    for(int i=0;i<5;++i) x[i]=rot(d[i]^shared,west[i]);
    return x;
}
inline State chi(const State& x) {
    State z{}; for(int i=0;i<5;++i) z[i]=x[i]^((~x[(i+1)%5])&x[(i+2)%5]); return z;
}
inline State before(State x,int r) { x=linear(x); x[0]^=rc.at(r); return x; }
inline State perm(State x) { for(int r=0;r<3;++r) x=chi(before(x,r)); return x; }
inline State input(Word m,const State& s) { State x=s; x[0]^=m; return x; }
struct Trace { State y1,y2,y3,out; };
inline Trace trace(Word m,const State& s) {
    Trace t; t.y1=before(input(m,s),0); t.y2=before(chi(t.y1),1);
    t.y3=before(chi(t.y2),2); t.out=chi(t.y3); return t;
}
struct Equation { Word row=0; bool rhs=false; };
struct Space {
    bool consistent=true; int rank=0; Word base=0;
    std::vector<Word> basis;
};
inline Space eliminate(std::vector<Equation> rows,int width=64) {
    Space s; std::vector<int> pivots;
    for(int col=0;col<width;++col) {
        int k=s.rank;
        while(k<int(rows.size()) && !(rows[k].row>>col&1)) ++k;
        if(k==int(rows.size())) continue;
        std::swap(rows[k],rows[s.rank]);
        for(int i=0;i<int(rows.size());++i) if(i!=s.rank && (rows[i].row>>col&1)) {
            rows[i].row^=rows[s.rank].row; rows[i].rhs^=rows[s.rank].rhs;
        }
        pivots.push_back(col); ++s.rank;
    }
    for(const auto& e:rows) if(!e.row && e.rhs) { s.consistent=false; return s; }
    Word pivot_mask=0;
    for(int k=0;k<s.rank;++k) {
        pivot_mask|=Word(1)<<pivots[k]; if(rows[k].rhs) s.base|=Word(1)<<pivots[k];
    }
    for(int col=0;col<width;++col) if(!(pivot_mask>>col&1)) {
        Word v=Word(1)<<col;
        for(int k=0;k<s.rank;++k) if(rows[k].row>>col&1) v|=Word(1)<<pivots[k];
        s.basis.push_back(v);
    }
    return s;
}
inline Word point(const Space& s,Word z) {
    Word m=s.base;
    while(z) { int i=__builtin_ctzll(z); m^=s.basis.at(i); z&=z-1; }
    return m;
}
inline Space restrict_space(const Space& outer,const Space& inner) {
    Space out; out.consistent=outer.consistent&&inner.consistent;
    if(!out.consistent) return out;
    out.base=point(outer,inner.base); out.rank=outer.rank+inner.rank;
    for(Word v:inner.basis) out.basis.push_back(point(outer,v)^outer.base);
    return out;
}
inline std::vector<State> full_anf(std::vector<State> values,int dimension) {
    if(values.size()!=(std::size_t(1)<<dimension)) throw std::runtime_error("ANF table size mismatch");
    for(int k=0;k<dimension;++k) for(std::size_t start=0;start<values.size();start+=(std::size_t(2)<<k))
        for(std::size_t j=0;j<(std::size_t(1)<<k);++j) for(int row=0;row<5;++row)
            values[start+(std::size_t(1)<<k)+j][row]^=values[start+j][row];
    return values;
}
inline Equation affine_equation(const std::vector<State>& anf,int dimension,int row,int j,bool value) {
    Equation e; e.rhs=value^bit(anf[0][row],j);
    for(std::size_t k=1;k<anf.size();++k) if(bit(anf[k][row],j)) {
        if(__builtin_popcountll(k)>1) throw std::runtime_error("Selected second-round condition is nonlinear");
        e.row|=Word(1)<<__builtin_ctzll(k);
    }
    if(dimension<64 && (e.row>>dimension)) throw std::runtime_error("Bad affine coefficient");
    return e;
}
inline bool is_affine(const std::vector<State>& anf,int row,int j) {
    for(std::size_t k=1;k<anf.size();++k) if(__builtin_popcountll(k)>1 && bit(anf[k][row],j)) return false;
    return true;
}
} // namespace gaston3

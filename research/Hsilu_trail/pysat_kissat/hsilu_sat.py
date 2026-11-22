"""PySAT CNF construction + external Kissat for the supplied MRM trail interfaces.

Reference support data derives from MRM (Apache-2.0): original copyright
2022 Rusydi H. Makarim and Raghvendra Rohit; MRM modifications 2024.
New backend and API implementation, 2026-10-11. See LICENSE and documentation.
All words use row 0 first and bit 0 as the least significant bit.
"""
from __future__ import annotations
import argparse
from contextlib import ExitStack
from functools import lru_cache
import json
import math
from pathlib import Path
import subprocess
import tempfile
import time
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE/'reference_data.json').read_text())
DEFAULT_KISSAT = HERE.parent/'tools/sat/bin/kissat'
MASK = (1 << 64) - 1
CARD_ENCODINGS = {name: getattr(EncType, name) for name in
                  ('seqcounter', 'totalizer', 'mtotalizer', 'kmtotalizer',
                   'sortnetwrk', 'cardnetwrk')}


# BEGIN MANUAL SBOX CNF
# Copy this entire block into hsilu_sat.py between the same markers.
# Literal 1..5: input row 0..4; 6..10: output row 0..4.
# Positive = variable, negative = NOT variable; no trailing DIMACS 0.
SBOX_TABLE = (0, 9, 18, 11, 5, 12, 22, 15, 10, 3, 24, 1, 13, 4, 30, 7, 20, 21, 6, 23, 17, 16, 2, 19, 26, 27, 8, 25, 29, 28, 14, 31)
# 切换时只启用一组差分子句；无需设置额外的排除输入配置。
SBOX_CNF = {
    'diff': (
        # ===== 完整差分关系：46 条，全部注释（317 个合法差分对） =====
        (1, 2, 3, -6),
        (-1, 2, 3, 6),
        (-2, 3, 4, 7),
        (-3, 4, 5, 8),
        (1, -4, 5, 9),
        (1, 2, -5, 10),
        (3, 4, 5, -8),
        (2, 5, 6, -10),
        (2, 4, -7, 8),
        (1, 4, -9, 10),
        (-1, 3, 6, 7),
        (-1, -2, 3, -4, -6, -7),
        (1, -4, 9, 10),
        (1, -2, -4, -5, -9, -10),
        (-2, 4, 7, 8),
        (-2, -3, 4, -5, -7, -8),
        (-3, 5, 8, 9),
        (-1, -3, -4, 5, -8, -9),
        (2, -5, 6, 10),
        (-1, 2, -3, -5, -6, -10),
        (1, 3, -6, 7),
        (1, -2, 3, -4, 6, -7),
        (2, -3, 4, -5, 7, -8),
        (1, 4, 5, -9),
        (3, 5, -8, 9),
        (-1, 3, -4, 5, 8, -9),
        (-1, 2, -3, 5, -6, 10),
        (1, -2, 4, -5, 9, -10),
        (2, 3, 4, -7),
        (1, 2, 5, -10),
        (-1, 6, 7, 8, 9, 10),
        (-1, -3, -6, 7, 8, 9, -10),
        (-1, -4, 6, 7, -8, -9, 10),
        (-2, -4, -6, -7, 8, 9, 10),
        (-1, -2, -4, -6, -7, -8, -9, 10),
        (-1, -3, -4, -6, 7, -8, -9, -10),
        (-1, -2, -3, -4, 6, -7, -8, -9, -10),
        (-1, -2, -3, -4, 6, -7, 8, 9, -10),
        (-2, -5, 6, 7, 8, -9, -10),
        (-3, -5, 6, -7, -8, 9, 10),
        (-1, -3, -5, -6, -7, -8, 9, -10),
        (-1, -2, -3, -5, -6, 7, 8, -9, 10),
        (-2, -4, -5, -6, -7, 8, -9, -10),
        (-1, -2, -4, -5, 6, 7, -8, 9, -10),
        (-1, -3, -4, -5, 6, -7, 8, -9, 10),
        (-2, -3, -4, -5, -6, 7, -8, 9, 10),


        # (-1, 6, 7, 8, 9, 10),
        # (-2, 6, 7, 8, 9, 10),
        # (-3, 6, 7, 8, 9, 10),
        # (-4, 6, 7, 8, 9, 10),
        # (-5, 6, 7, 8, 9, 10),



        # ===== 排除输入 11111：31 条，当前启用（301 个合法差分对） =====
        (-1, -2, -3, -4, -5),  # 禁止输入 11111
        # (-1, 2, 3, 6),
        # (1, 2, 3, -6),
        # (-2, 3, 4, 7),
        # (-3, 4, 5, 8),
        # (1, -4, 5, 9),
        # (1, 2, -5, 10),
        # (3, 4, 5, -8),
        # (2, 5, 6, -10),
        # (2, 4, -7, 8),
        # (1, 4, -9, 10),
        # (-1, 3, 6, 7),
        # (1, -4, 9, 10),
        # (-1, -2, 3, -4, -6, -7),
        # (1, -2, -4, -5, -9, -10),
        # (-2, 4, 7, 8),
        # (-3, 5, 8, 9),
        # (2, -5, 6, 10),
        # (-2, -3, 4, -5, -7, -8),
        # (-1, -3, -4, 5, -8, -9),
        # (-1, 2, -3, -5, -6, -10),
        # (1, 3, -6, 7),
        # (1, 4, 5, -9),
        # (3, 5, -8, 9),
        # (2, -3, 4, -5, 7, -8),
        # (-1, 3, -4, 5, 8, -9),
        # (2, 3, 4, -7),
        # (1, 2, 5, -10),
        # (1, -2, 3, -4, 6, -7),
        # (-1, 2, -3, 5, -6, 10),
        # (1, -2, 4, -5, 9, -10),
    ),
    'linear': (
        (-3, 6, 7, 8),
        (3, 6, 7, -8),
        (2, 6, -7, 10),
        (1, -6, 9, 10),
        (5, 8, 9, -10),
        (4, 7, 8, -9),
        (-1, 6, 9, 10),
        (-4, 7, 8, 9),
        (-2, 6, 7, 10),
        (-5, 8, 9, 10),
        (1, 3, 5, 9, -10),
        (2, 5, 6, -7, 8),
        (1, 3, 4, 7, -9),
        (1, 2, 4, -7, 10),
        (2, -3, -5, 6, -7, -8),
        (-1, 3, -4, 7, -8, -9),
        (1, -2, -4, -6, -7, 10),
        (-1, -3, 5, -6, 9, -10),
        (-2, 4, -5, 8, -9, -10),
        (2, 3, 5, 6, -8),
        (1, -3, -5, -7, -8, 9, -10),
        (1, 3, 7, -8, 9),
        (1, 4, -6, 7, 10),
        (2, 4, 5, 8, -10),
        (2, 4, 8, -9, 10),
        (-2, -4, 5, -6, -7, 8, -10),
        (-1, -3, 4, -6, 7, -9, -10),
        (-2, 3, -5, 6, -8, -9, -10),
        (-1, 2, -4, -7, -8, -9, 10),
        (-1, -3, -6, 7, -8, 9, -10),
        (-1, -4, -6, 7, -8, -9, 10),
        (-2, -4, -6, -7, 8, -9, 10),
        (-2, -5, 6, -7, 8, -9, -10),
        (-1, -2, 4, -6, -7, -8, -9, 10),
        (1, -3, -4, -6, 7, -8, -9, -10),
        (-1, 3, -5, -6, -7, -8, 9, -10),
        (-2, -3, 5, 6, -7, -8, -9, -10),
        (2, -4, -5, -6, -7, 8, -9, -10),
        (-2, 5, 6, 7, 8),
        (1, -3, 7, 8, 9),
        (-1, 4, 6, 7, 10),
        (3, -5, 6, 9, 10),
        (3, 5, 6, 9, -10),
        (2, -4, 8, 9, 10),
        (-1, 3, -6, 7, 8, 9, -10),
        (1, -4, 6, 7, -8, -9, 10),
        (-2, 4, -6, -7, 8, 9, 10),
        (2, -5, 6, 7, 8, -9, -10),
        (-3, 5, 6, -7, -8, 9, 10),
        (-3, -5, 6, -7, -8, 9, -10),
        (1, 2, 3, 4, 5, -10),
        (1, -2, 3, -4, 5, -6, -7),
        (1, 2, -3, 4, -5, -7, -8),
        (-1, 2, 3, -4, 5, -8, -9),
        (-1, 2, -3, 4, 5, -6, -10),
        (1, -2, 3, 4, -5, -9, -10),
        (-1, 2, -3, -4, -5, -7, -8, -9),
        (-1, -2, -3, 4, -5, -6, -9, -10),
        (-1, -2, 3, -4, -5, -8, -9, -10),
        (-1, -2, 3, 4, 5, -6, -7, -8, -9),
        (-1, 2, 3, 4, -5, -6, -7, -8, -10),
        (1, 2, 3, -4, -5, -6, -7, -9, -10),
        (1, 2, -3, -4, 5, -6, -8, -9, -10),
        (1, -2, -3, -4, -5, -6, -8, -9, -10),
        (1, -2, -3, 4, 5, -7, -8, -9, -10),
        (-1, -2, -3, -4, 5, -7, -8, -9, -10),
    ),
}
# END MANUAL SBOX CNF


def rotate(w, c, direction):
    c %= 64
    if direction == 'left':
        return ((w << c) | (w >> ((64-c) % 64))) & MASK
    return ((w >> c) | (w << ((64-c) % 64))) & MASK


def chi5(x):
    return sum((((x >> i) & 1) ^ ((1 ^ ((x >> ((i+1)%5)) & 1)) &
                ((x >> ((i+2)%5)) & 1))) << i for i in range(5))


@lru_cache(None)
def exact_table(mode):
    """Table for weight evaluation only; never synthesizes or checks the pasted CNF."""
    if mode == 'diff':
        table = [[0]*32 for _ in range(32)]
        for a in range(32):
            for x in range(32):
                table[a][SBOX_TABLE[x] ^ SBOX_TABLE[x ^ a]] += 1
    else:
        table = [[sum(1 if ((a & x).bit_count() + (b & SBOX_TABLE[x]).bit_count()) % 2 == 0
                      else -1 for x in range(32)) for b in range(32)] for a in range(32)]
    return table


def entry_weight(count, factor):
    log_count=count.bit_length()-1 if count & (count-1)==0 else math.log2(count)
    return factor*(5-log_count)


def diff_xor_plan(params, max_inputs=3):
    """Inline internal words with exactly one consumer, keeping shared words.

    Terms are (word index, rotation amount), in params['direction'].
    Final words remain explicit to preserve the existing layer interface.
    """
    n=len(params['n1'])
    def parity(terms):
        odd=set()
        for word,shift in terms:
            term=(word,shift%64)
            if term in odd: odd.remove(term)
            else: odd.add(term)
        return tuple(sorted(odd))
    plan={j+5:parity(((a,0),(b,m)))
          for j,(a,b,m) in enumerate(zip(params['n1'],params['n2'],params['m']))}
    removed=[]
    for word in range(5,n):  # n..n+4 are layer outputs, never inline them.
        uses=[(parent,shift) for parent,terms in plan.items()
              for child,shift in terms if child==word]
        if len(uses)!=1: continue
        parent,shift=uses[0]
        expanded=parity([t for t in plan[parent] if t!=(word,shift)]+
                        [(child,rot+shift) for child,rot in plan[word]])
        if len(expanded)>max_inputs: continue
        plan[parent]=expanded
        del plan[word]
        removed.append(word)
    return plan,removed


def linear_value(words, params, mode='diff', relation='exact'):
    n = len(params['n1'])
    if mode == 'diff':
        z = list(words)
        for a,b,m in zip(params['n1'], params['n2'], params['m']):
            z.append(z[a] ^ rotate(z[b], m, params['direction']))
        return [rotate(z[n+i], params['rot'][i], params['direction']) for i in range(5)]
    # Reference linear_smt explicitly uses right rotations in the backward
    # accumulator. Exact mode uses the adjoint of the selected forward layer.
    direction = 'right' if relation == 'reference' else ('left' if params['direction']=='right' else 'right')
    z = [0]*(n+5)
    for i in range(5):
        z[n+i] = rotate(words[i], params['rot'][i], direction)
    for j in range(n-1,-1,-1):
        v = z[j+5]
        z[params['n1'][j]] ^= v
        z[params['n2'][j]] ^= rotate(v, params['m'][j], direction)
    return z[:5]


class Trail:
    def __init__(self, states, labels, edges, rounds, mode):
        self._states, self.labels, self.edges = states, labels, edges
        self._rounds, self.mode = rounds, mode
    def trail(self): return self._states
    def nrounds(self): return self._rounds
    def input_values(self): return self._states[0]
    def output_values(self): return self._states[-1]
    def __repr__(self): return f'{self._rounds}-round {self.mode} trail core'
    def __str__(self):
        blocks = []
        for i,(label,state) in enumerate(zip(self.labels,self._states)):
            blocks.append(label + '\n' + '\n'.join(f'{w:064b}'.replace('0','.') for w in state))
            if i < len(self.edges): blocks.append('-'*64+' '+self.edges[i])
        return '\n'.join(blocks)


class HsiluBase:
    """Default: paper Hsilu, exact support, (chi,L) repeated nrounds-1 times.

    S-box activity is counted at y0, x1, ..., x[nrounds-1].
    Bit counts retain their input-state convention (x0, or y0 in reference).
    The final x state is counted as the input to the omitted final chi.
    remove_outer_layer=False selects full rounds. For the original MRM topology,
    explicitly use topology='reference'; it also omits the initial chi.
    """
    word_size, nwords, size = 64, 5, 320
    mode = 'diff'
    def __init__(self, nrounds, remove_outer_layer=True, *, profile='hsilu',
                 relation='exact', topology=None, kissat=None, output_prefix=None,
                 xor_encoding='binary', xor_max_inputs=3,
                 card_encoding='totalizer', shared_activity=True,
                 merge_equalities=True, preprocess=False, column_symmetry=False,
                 first_nonzero=True):
        if type(nrounds) is not int or nrounds < 1: raise ValueError('nrounds must be positive')
        topology = topology or ('core' if remove_outer_layer else 'full')
        if topology not in ('core','reference','full'): raise ValueError('unknown topology')
        if topology != 'full' and nrounds < 2: raise ValueError('truncated core requires at least 2 rounds')
        if profile not in DATA['profiles']: raise ValueError('unknown profile')
        if relation not in ('reference','exact'): raise ValueError('unknown relation')
        if xor_encoding not in ('binary','inline'): raise ValueError('unknown XOR encoding')
        if card_encoding not in CARD_ENCODINGS: raise ValueError('unknown cardinality encoding')
        self.card_encoding = card_encoding
        self.shared_activity = shared_activity
        self.merge_equalities = merge_equalities
        self.preprocess = preprocess
        self.column_symmetry = column_symmetry
        self.first_nonzero = first_nonzero
        self.solve_stats = {}
        if type(xor_max_inputs) is not int or not 2<=xor_max_inputs<=8:
            raise ValueError('xor_max_inputs must be in 2..8')
        if self.mode!='diff' and xor_encoding!='binary':
            raise ValueError('inline XOR encoding currently supports differential propagation only')
        self.xor_encoding,self.xor_max_inputs=xor_encoding,xor_max_inputs
        self.linear_encoding_stats=[]
        self.nrounds, self.profile, self.relation = nrounds, profile, relation
        self.topology = topology
        self._removed = topology != 'full'
        self.params = DATA['profiles'][profile]
        if relation=='exact':
            self.sbox_clauses = SBOX_CNF[self.mode]
        else:
            # Legacy snapshot only; no CNF synthesis during model construction.
            reference_rows=json.loads((HERE/'reference_cnf.json').read_text())
            key=profile+':'+self.mode
            if key not in reference_rows:
                raise ValueError(f'No supplied reference implementation for {key}; use relation="exact"')
            self.sbox_clauses = reference_rows[key]
        if any(type(v) is not int or not 1<=abs(v)<=10 for row in self.sbox_clauses for v in row):
            raise ValueError('S-box CNF literals must be in -10..10, excluding 0')
        self.kissat = str(Path(kissat or DEFAULT_KISSAT).resolve())
        self.output_prefix = output_prefix
        self.pool, self.cnf = IDPool(), CNF()
        self._states, self._activity = {}, {}
        self._model, self._solving_time, self.status = None, None, None
        self.additional_constraints = []  # extra CNF clauses, not STP expressions
        self._checks = []
        self._false = self.pool.id('constant_false')
        self.cnf.append([-self._false])
        self._linear_rounds = list(range(nrounds-1 if self._removed else nrounds))
        self._chi_rounds = list(range(1,nrounds-1) if topology=='reference' else range(nrounds-1 if self._removed else nrounds))
        for r in self._linear_rounds: self._add_linear(r)
        for r in self._chi_rounds: self._add_chi(r)
        self._model = None

    @property
    def solver(self): return self.cnf
    def is_outer_layer_removed(self): return self._removed
    @staticmethod
    def rotate_right(w,c): return rotate(w,c,'right')
    @staticmethod
    def rotate_left(w,c): return rotate(w,c,'left')
    def _round(self,r):
        if type(r) is not int or not 0 <= r < self.nrounds: raise ValueError('round index out of range')
    def vars(self,name,r):
        if name not in ('x','y') or type(r) is not int or not 0 <= r <= self.nrounds:
            raise ValueError('state must be x/y and round in 0..nrounds')
        key = (name,r)
        if key not in self._states:
            self._states[key] = [[self.pool.id((name,r,i,z)) for z in range(64)] for i in range(5)]
        return self._states[key]
    def input_nonlinear_vars(self,r): return self.vars('x',r)
    def input_linear_vars(self,r): return self.vars('y',r)
    def output_nonlinear_vars(self,r): return self.vars('y',r)
    def output_linear_vars(self,r): return self.vars('x',r+1)
    def input_vars(self): return self.vars('y' if self.topology=='reference' else 'x',0)
    def output_vars(self): return self.vars('x',self.nrounds-1 if self._removed else self.nrounds)
    def _add(self,clause):
        self.cnf.append(clause)
        self._model = None
    def _xor(self,a,b):
        if a == b: return self._false
        if a == self._false: return b
        if b == self._false: return a
        c = self.pool.id()
        for clause in ([a,b,-c],[a,-b,c],[-a,b,c],[-a,-b,-c]): self._add(clause)
        return c
    def _word_xor(self,a,b): return [self._xor(x,y) for x,y in zip(a,b)]
    def _xor_many(self,inputs):
        odd=set()
        for v in inputs:
            if v==self._false: continue
            if v in odd: odd.remove(v)
            else: odd.add(v)
        terms=sorted(odd)
        if not terms: return self._false
        if len(terms)==1: return terms[0]
        if len(terms)==2: return self._xor(*terms)
        c=self.pool.id(); variables=terms+[c]
        # XOR(inputs,c)=0: forbid every odd-parity assignment.
        for mask in range(1<<len(variables)):
            if mask.bit_count()%2:
                self._add([-v if mask>>i&1 else v for i,v in enumerate(variables)])
        return c
    def _rot(self,a,c,direction):
        return [a[(i-c if direction=='left' else i+c)%64] for i in range(64)]
    def _equal_words(self,a,b):
        for x,y in zip(a,b):
            if x != y:
                self._add([-x,y]); self._add([x,-y])
    def _add_linear(self,r):
        # Allocate only independent endpoint variables when using aliases.
        x = None if self.merge_equalities and self.mode!='diff' else self.input_linear_vars(r)
        y = None if self.merge_equalities and self.mode=='diff' else self.output_linear_vars(r)
        p=self.params; n=len(p['n1'])
        if self.mode == 'diff':
            if self.xor_encoding=='binary':
                z=list(x)
                for a,b,m in zip(p['n1'],p['n2'],p['m']):
                    z.append(self._word_xor(z[a],self._rot(z[b],m,p['direction'])))
            else:
                plan,removed=diff_xor_plan(p,self.xor_max_inputs)
                z=dict(enumerate(x))
                for word,terms in sorted(plan.items()):
                    operands=[self._rot(z[child],shift,p['direction']) for child,shift in terms]
                    z[word]=[self._xor_many([operand[bit] for operand in operands]) for bit in range(64)]
                histogram={str(k):sum(len(terms)==k for terms in plan.values())
                           for k in sorted({len(terms) for terms in plan.values()})}
                self.linear_encoding_stats.append(dict(round=r,removed_words=removed,
                                                       word_xors_by_input_count=histogram))
            output=[self._rot(z[n+i],p['rot'][i],p['direction']) for i in range(5)]
            if self.merge_equalities: self._states[('x',r+1)]=output
            else:
                for i in range(5): self._equal_words(y[i],output[i])
        else:
            direction='right' if self.relation=='reference' else ('left' if p['direction']=='right' else 'right')
            z=[[self._false]*64 for _ in range(n+5)]
            for i in range(5): z[n+i]=self._rot(y[i],p['rot'][i],direction)
            for j in range(n-1,-1,-1):
                v=z[j+5]; a,b=p['n1'][j],p['n2'][j]
                z[a]=self._word_xor(z[a],v)
                z[b]=self._word_xor(z[b],self._rot(v,p['m'][j],direction))
            if self.merge_equalities: self._states[('y',r)]=z[:5]
            else:
                for i in range(5): self._equal_words(x[i],z[i])
    def _add_chi(self,r):
        x,y=self.input_nonlinear_vars(r),self.output_nonlinear_vars(r)
        rows=self.sbox_clauses  # manually pasted local CNF, instantiated at each column
        for z in range(64):
            vs=[w[z] for w in x+y]
            for row in rows: self._add([vs[abs(v)-1] if v>0 else -vs[abs(v)-1] for v in row])
    def _count_key(self,r,bits=False):
        self._round(r)
        return ('y',0) if r==0 and (not bits or self.topology=='reference') else ('x',r)
    def _count_state(self,r,bits=False):
        return self.vars(*self._count_key(r,bits))
    def nactive_sboxes_var(self,r):
        """Return 64 Boolean activity literals (not a bit-vector integer)."""
        if r not in self._activity:
            state=self._count_state(r); active=[]
            sides=[state]
            if self.shared_activity and r in self._chi_rounds:
                sides=[self.vars('x',r),self.vars('y',r)]
            for z in range(64):
                a=self.pool.id(('active',r,z))
                for side in sides:
                    bits=[w[z] for w in side]
                    self._add([-a]+bits)
                    for b in bits: self._add([-b,a])
                active.append(a)
            self._activity[r]=active
        return self._activity[r]
    def nactive_sboxes_vars(self): return [self.nactive_sboxes_var(r) for r in range(self.nrounds)]
    def nactive_bits_var(self,r): return [v for w in self._count_state(r,bits=True) for v in w]
    def _card(self,lits,k,kind):
        if type(k) is not int or not 0 <= k <= len(lits): raise ValueError('invalid count')
        enc=getattr(CardEnc,kind)(lits=lits,bound=k,vpool=self.pool,encoding=CARD_ENCODINGS[self.card_encoding])
        self.cnf.extend(enc.clauses); self._model=None
    def _limit(self,r,k,kind,bits=False):
        self._round(r)
        lits=self.nactive_bits_var(r) if bits else self.nactive_sboxes_var(r)
        self._card(lits,k,kind); self._checks.append(('bits' if bits else 'active',r,kind,k))
    def set_nactive_sboxes_at_round(self,r,w): self._limit(r,w,'equals')
    def set_min_nactive_sboxes_at_round(self,r,w): self._limit(r,w,'atleast')
    def set_max_nactive_sboxes_at_round(self,r,w): self._limit(r,w,'atmost')
    def set_nactive_bits_at_round(self,r,w): self._limit(r,w,'equals',True)
    def set_min_nactive_bits_at_round(self,r,w): self._limit(r,w,'atleast',True)
    def set_max_nactive_bits_at_round(self,r,w): self._limit(r,w,'atmost',True)
    def set_active_sboxes_at_round(self,r,active_sboxes):
        indices=set(active_sboxes)
        if any(type(i) is not int or not 0 <= i < 64 for i in indices): raise ValueError('invalid active position')
        for i,a in enumerate(self.nactive_sboxes_var(r)): self._add([a if i in indices else -a])
        self._checks.append(('positions',r,sorted(indices)))
    def _fix(self,name,r,words):
        if len(words)!=5 or any(type(w) is not int or not 0<=w<=MASK for w in words): raise ValueError('expected five uint64 words')
        for vs,value in zip(self.vars(name,r),words):
            for z,v in enumerate(vs): self._add([v if value>>z&1 else -v])
        self._checks.append(('fixed',name,r,list(words)))
    def set_input_nonlinear_state_at_round(self,r,*words): self._fix('x',r,words)
    def set_output_nonlinear_state_at_round(self,r,*words): self._fix('y',r,words)
    def set_input_linear_state_at_round(self,r,*words): self._fix('y',r,words)
    def set_output_linear_state_at_round(self,r,*words): self._fix('x',r+1,words)
    def set_input_round_state(self,r,*words): self._fix('x',r,words)
    def set_output_round_state(self,r,*words): self._fix('x',r+1,words)
    def require_nonzero(self): self._add([v for w in self._count_state(0,bits=True) for v in w])
    def _formula(self,total=None):
        symmetry=[]
        positive_count=any(c[0]=='active' and c[1]==0 and c[2] in ('equals','atleast') and c[3]>0
                           for c in self._checks)
        if self.column_symmetry:
            # Absolute positions/custom clauses may destroy cyclic symmetry.
            if self.additional_constraints or any(c[0] in ('fixed','positions') for c in self._checks):
                raise ValueError('column symmetry requires count-only constraints (no fixed states, positions or custom clauses)')
            active=self.nactive_sboxes_var(0)
            nonzero=self.first_nonzero or positive_count
            # Keep the zero orbit when nonzero activity is not required.
            symmetry=[[active[0]]] if nonzero else [[-a,active[0]] for a in active[1:]]
        if total is not None:
            lits=[v for row in self.nactive_sboxes_vars() for v in row]
        cnf=CNF(from_clauses=self.cnf.clauses)
        cnf.extend(self.additional_constraints)
        cnf.extend(symmetry)
        if self.first_nonzero and not positive_count and not self.column_symmetry:
            cnf.append([v for w in self._count_state(0) for v in w])
        if total is not None:
            if type(total) is not int or not 0 <= total <= 64*self.nrounds: raise ValueError('invalid total count')
            # temporary pool: repeated queries must not permanently constrain total
            enc=CardEnc.equals(lits=lits,bound=total,top_id=max(self.pool.top,cnf.nv),encoding=CARD_ENCODINGS[self.card_encoding])
            cnf.extend(enc.clauses)
        return cnf
    def export_cnf(self,path,nactive_sboxes=None):
        cnf=self._formula(nactive_sboxes); cnf.to_file(str(path)); return cnf.nv,len(cnf.clauses)
    def has_trail(self,nactive_sboxes=None,timeout=None):
        if timeout is not None and timeout<=0: raise ValueError('timeout must be positive')
        self._model=None; self.status=None
        cnf=self._formula(nactive_sboxes)
        original_cnf=cnf
        self.solve_stats=dict(shared_activity=self.shared_activity,
            merge_equalities=self.merge_equalities,preprocess=self.preprocess,
            column_symmetry=self.column_symmetry,
            first_nonzero=self.first_nonzero,
            original_variables=cnf.nv,original_clauses=len(cnf.clauses),
            preprocessing_seconds=0.0,restore_seconds=0.0)
        with ExitStack() as stack, tempfile.TemporaryDirectory(prefix='hsilu-kissat-') as tmp:
            prefix=Path(self.output_prefix) if self.output_prefix else Path(tmp)/'trail'
            prefix.parent.mkdir(parents=True,exist_ok=True)
            cnfpath=Path(str(prefix)+'.cnf'); logpath=Path(str(prefix)+'.log')
            resultpath=Path(str(prefix)+'.json')
            if self.output_prefix:
                resultpath.write_text(json.dumps({'status':'RUNNING'})+'\n')
            processor=None
            if self.preprocess:
                from pysat.process import Processor
                if self.output_prefix: original_cnf.to_file(str(prefix)+'.original.cnf')
                before=time.monotonic()
                processor=stack.enter_context(Processor(bootstrap_with=original_cnf))
                cnf=processor.process(rounds=1)
                self.solve_stats['preprocessing_seconds']=time.monotonic()-before
            self.solve_stats.update(solver_variables=cnf.nv,
                solver_used_variables=len({abs(v) for c in cnf.clauses for v in c}),
                solver_clauses=len(cnf.clauses))
            cnf.to_file(str(cnfpath)); start=time.monotonic()
            with logpath.open('w') as log:
                try:
                    proc=subprocess.Popen([self.kissat,str(cnfpath)],stdout=log,stderr=subprocess.STDOUT)
                except OSError as exc:
                    self.status='error'
                    if self.output_prefix:
                        resultpath.write_text(json.dumps({'status':'error','message':str(exc)})+'\n')
                    raise
                try: code=proc.wait(timeout=timeout)
                except (subprocess.TimeoutExpired,KeyboardInterrupt) as exc:
                    proc.terminate()
                    try: proc.wait(timeout=5)
                    except subprocess.TimeoutExpired: proc.kill(); proc.wait()
                    self._solving_time=time.monotonic()-start
                    self.status='timeout' if isinstance(exc,subprocess.TimeoutExpired) else 'interrupted'
                    if self.output_prefix:
                        resultpath.write_text(json.dumps({'status':self.status,'seconds':self._solving_time,'solve_stats':self.solve_stats})+'\n')
                    return self.status
            self._solving_time=time.monotonic()-start
            text=logpath.read_text()
            if code==20 and 's UNSATISFIABLE' in text:
                self.status='UNSAT'
                if self.output_prefix:
                    resultpath.write_text(json.dumps({'status':'UNSAT','seconds':self._solving_time,'solve_stats':self.solve_stats})+'\n')
                return False
            if code!=10 or 's SATISFIABLE' not in text:
                self.status='error'
                if self.output_prefix:
                    resultpath.write_text(json.dumps({'status':'error','exit_code':code})+'\n')
                raise RuntimeError(f'Kissat exit {code}: {text[-2000:]}')
            values={}
            for line in text.splitlines():
                if line.startswith('v '):
                    for token in line.split()[1:]:
                        lit=int(token)
                        if lit: values[abs(lit)]=lit>0
            if not all(any(values.get(abs(v))==(v>0) for v in clause) for clause in cnf.clauses):
                raise RuntimeError('Kissat assignment does not satisfy exported CNF')
            if processor is not None:
                before=time.monotonic()
                restored=processor.restore([v if truth else -v for v,truth in values.items()])
                values={abs(v):v>0 for v in restored}
                self.solve_stats['restore_seconds']=time.monotonic()-before
                if not all(any(values.get(abs(v))==(v>0) for v in c) for c in original_cnf.clauses):
                    raise RuntimeError('Restored assignment does not satisfy original CNF')
            self._model=values; self.status='SAT'
            self.validate(nactive_sboxes)
            if self.output_prefix:
                Path(str(prefix)+'.json').write_text(json.dumps(self.result(),indent=2)+'\n')
            return True
    def solving_time(self):
        if self._solving_time is None: raise RuntimeError('run has_trail first')
        return self._solving_time
    def _values(self,name,r):
        if self._model is None: raise RuntimeError('no current SAT model')
        # Omitted outer states have no propagation semantics, even if constrained.
        vs=self._states.get((name,r))
        if vs is None: raise ValueError(f'{name}{r} is omitted/unallocated; use trail() for connected core')
        return [sum((1<<z) for z,v in enumerate(w) if self._model.get(v,False)) for w in vs]
    def input_nonlinear_values(self,r): return self._values('x',r)
    def input_linear_values(self,r): return self._values('y',r)
    def output_nonlinear_values(self,r): return self._values('y',r)
    def output_linear_values(self,r): return self._values('x',r+1)
    def output_values(self): return self._values('x',self.nrounds-1 if self._removed else self.nrounds)
    def input_values(self): return self._values('y' if self.topology=='reference' else 'x',0)
    def _count_values(self,r,bits=False):
        return self._values(*self._count_key(r,bits))
    def active_sboxes_at_round(self,r):
        words=self._count_values(r); union=0
        for w in words: union |= w
        return tuple(z for z in range(64) if union>>z&1)
    def nactive_sboxes_at_round(self,r): return len(self.active_sboxes_at_round(r))
    def nactive_bits_at_round(self,r): return sum(w.bit_count() for w in self._count_values(r,bits=True))
    def nactive_sboxes(self): return sum(self.nactive_sboxes_at_round(r) for r in range(self.nrounds))
    def trail_at_round(self,r): return [self._values('x',r),self._values('y',r)]
    def trail(self):
        if self.topology=='reference':
            keys=[('y',0),('x',1)]; edges=['pL']
            for r in range(1,self.nrounds-1): keys.extend([('y',r),('x',r+1)]); edges.extend(['pS','pL'])
        else:
            keys=[('x',0)]; edges=[]
            for r in self._chi_rounds: keys.extend([('y',r),('x',r+1)]); edges.extend(['pS','pL'])
        return Trail([self._values(*k) for k in keys],[f'{k[0]}{k[1]}' for k in keys],edges,self.nrounds,self.mode)
    def validate(self,total=None):
        if self.first_nonzero and not self.active_sboxes_at_round(0):
            raise RuntimeError('first round must be nonzero')
        if self.column_symmetry:
            positions=self.active_sboxes_at_round(0)
            if positions and 0 not in positions: raise RuntimeError('column symmetry anchor mismatch')
        for r in self._linear_rounds:
            x,y=self._values('y',r),self._values('x',r+1)
            actual=linear_value(x if self.mode=='diff' else y,self.params,self.mode,self.relation)
            if actual != (y if self.mode=='diff' else x): raise RuntimeError('linear layer mismatch')
        for check in self._checks:
            if check[0] in ('active','bits'):
                kind,r,op,k=check; val=self.nactive_sboxes_at_round(r) if kind=='active' else self.nactive_bits_at_round(r)
                ok={'equals':val==k,'atleast':val>=k,'atmost':val<=k}[op]
            elif check[0]=='positions': ok=list(self.active_sboxes_at_round(check[1]))==check[2]
            else: ok=self._values(check[1],check[2])==check[3]
            if not ok: raise RuntimeError(f'constraint mismatch: {check}')
        if total is not None and self.nactive_sboxes()!=total: raise RuntimeError('total activity mismatch')
        return True
    def result(self):
        trail=self.trail(); table=exact_table(self.mode); local=[]
        factor=1 if self.mode=='diff' else 2
        for r in self._chi_rounds:
            x,y=self._values('x',r),self._values('y',r); invalid=[]; weights=[]
            for z in range(64):
                a=sum(((w>>z)&1)<<i for i,w in enumerate(x)); b=sum(((w>>z)&1)<<i for i,w in enumerate(y))
                c=abs(table[a][b])
                if not c: invalid.append({'z':z,'input':a,'output':b})
                else: weights.append(entry_weight(c,factor))
            local.append({'round':r,'standard_chi_invalid':invalid,'weight':sum(weights) if not invalid else None})
        terminal=None
        if self._removed:
            words=self._values('x',self.nrounds-1); output=[0]*5; weight=0; uniform=True
            for z in range(64):
                a=sum(((w>>z)&1)<<i for i,w in enumerate(words))
                b=max(range(32),key=lambda b:abs(table[a][b]))
                c=abs(table[a][b])
                weight+=entry_weight(c,factor)
                uniform &= len({abs(v) for v in table[a] if v})==1
                for i in range(5): output[i] |= ((b>>i)&1)<<z
            terminal=dict(round=self.nrounds-1,weight=weight,
                          selection='maximum absolute DDT/LAT entry for each fixed input',
                          all_feasible_outputs_have_same_weight=uniform,
                          output=[f'{w:016x}' for w in output])
        modeled_total=sum(v['weight'] for v in local) if all(v['weight'] is not None for v in local) else None
        completed_total=None
        if modeled_total is not None and self.topology!='reference':
            completed_total=modeled_total+(terminal['weight'] if terminal else 0)
        return dict(status=self.status,mode=self.mode,profile=self.profile,relation=self.relation,
                    card_encoding=self.card_encoding,
                    solve_stats=self.solve_stats,
                    xor_encoding=self.xor_encoding,xor_max_inputs=self.xor_max_inputs,
                    linear_encoding_stats=self.linear_encoding_stats,
                    nrounds=self.nrounds,topology=self.topology,remove_outer_layer=self._removed,seconds=self.solving_time(),
                    active=[self.nactive_sboxes_at_round(r) for r in range(self.nrounds)],
                    active_count_states=[''.join(map(str,self._count_key(r))) for r in range(self.nrounds)],
                    active_bits=[self.nactive_bits_at_round(r) for r in range(self.nrounds)],
                    bit_count_states=[''.join(map(str,self._count_key(r,bits=True))) for r in range(self.nrounds)],
                    states={label:[f'{w:016x}' for w in s] for label,s in zip(trail.labels,trail.trail())},
                    modeled_sbox_weights=local,
                    modeled_sbox_total_weight=modeled_total,
                    terminal_sbox_completion=terminal,
                    completed_trail_weight=completed_total,
                    weight_convention='diff: -log2(DDT/32); linear: -log2((Walsh/32)^2). Terminal completion is postprocessing, not a SAT constraint or an optimality claim. Round-constant signs are not reported.')


class HsiluDifferential(HsiluBase): mode='diff'
class HsiluLinear(HsiluBase): mode='linear'
AsconDifferential=HsiluDifferential
AsconLinear=HsiluLinear


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['diff','linear'],default='diff')
    p.add_argument('--profile',choices=sorted(DATA['profiles']),default='hsilu')
    p.add_argument('--relation',choices=['reference','exact'],default='exact')
    p.add_argument('--rounds',type=int,default=3)
    p.add_argument('--full',action='store_true',help='alias for --topology full')
    p.add_argument('--topology',choices=['core','reference','full'],default='core')
    p.add_argument('--active',type=int,nargs='+',help='exact active columns at y0, x1, x2, ...')
    p.add_argument('--min-active',type=int,nargs='+')
    p.add_argument('--max-active',type=int,nargs='+')
    p.add_argument('--bits',type=int,nargs='+')
    p.add_argument('--min-bits',type=int,nargs='+')
    p.add_argument('--max-bits',type=int,nargs='+')
    p.add_argument('--total-active',type=int)
    p.add_argument('--active-positions',action='append',default=[],metavar='ROUND:Z,Z,...')
    p.add_argument('--fix',action='append',default=[],metavar='x0=WORD,WORD,WORD,WORD,WORD')
    p.add_argument('--nonzero',action='store_true')
    p.add_argument('--allow-zero',action='store_true',help='disable default nonzero first-round activity')
    p.add_argument('--xor-encoding',choices=['binary','inline'],default='binary')
    p.add_argument('--xor-max-inputs',type=int,default=3,help='max inputs per XOR after inlining (2..8)')
    p.add_argument('--card-encoding',choices=list(CARD_ENCODINGS),default='totalizer',
                   help='encoding for all activity/bit/total cardinality constraints')
    p.add_argument('--shared-activity',action=argparse.BooleanOptionalAction,default=True,
                   help='connect activity to both sides of each modeled chi')
    p.add_argument('--merge-equalities',action=argparse.BooleanOptionalAction,default=True,
                   help='reuse linear endpoint variables instead of equality clauses')
    p.add_argument('--preprocess',action='store_true',help='one round of PySAT Processor, restore original SAT assignment')
    p.add_argument('--column-symmetry',action='store_true',
                   help='anchor nonzero round-0 activity at column 0; count-only queries')
    p.add_argument('--timeout',type=float)
    p.add_argument('--kissat',type=Path,default=DEFAULT_KISSAT)
    p.add_argument('--prefix',default='hsilu_sat_run')
    p.add_argument('--export-only',action='store_true')
    o=p.parse_args()
    try:
        A=(HsiluDifferential if o.mode=='diff' else HsiluLinear)(o.rounds,not o.full,
            profile=o.profile,relation=o.relation,topology='full' if o.full else o.topology,kissat=o.kissat,output_prefix=o.prefix,
            xor_encoding=o.xor_encoding,xor_max_inputs=o.xor_max_inputs,
            card_encoding=o.card_encoding,shared_activity=o.shared_activity,
            merge_equalities=o.merge_equalities,preprocess=o.preprocess,
            column_symmetry=o.column_symmetry,first_nonzero=not o.allow_zero)
        for vals,fn in [(o.active,A.set_nactive_sboxes_at_round),(o.min_active,A.set_min_nactive_sboxes_at_round),
                        (o.max_active,A.set_max_nactive_sboxes_at_round),(o.bits,A.set_nactive_bits_at_round),
                        (o.min_bits,A.set_min_nactive_bits_at_round),(o.max_bits,A.set_max_nactive_bits_at_round)]:
            if vals is not None:
                if len(vals)!=o.rounds: raise ValueError('counts must contain one value per round')
                for r,v in enumerate(vals): fn(r,v)
        for spec in o.active_positions:
            r,positions=spec.split(':'); A.set_active_sboxes_at_round(int(r),[int(i) for i in positions.split(',') if i])
        for spec in o.fix:
            name,words=spec.split('='); A._fix(name[0],int(name[1:]),[int(w,0) for w in words.split(',')])
        if o.nonzero: A.require_nonzero()
        print(f'{o.mode} profile={o.profile} relation={o.relation}; S layers={len(A._chi_rounds)}, L layers={len(A._linear_rounds)}',flush=True)
        print('counted states:', [''.join(map(str,A._count_key(r))) for r in range(o.rounds)],flush=True)
        print('XOR encoding:',A.xor_encoding,flush=True)
        print('cardinality encoding:',A.card_encoding,flush=True)
        if A.linear_encoding_stats: print('linear encoding:',A.linear_encoding_stats,flush=True)
        if o.export_only:
            Path(o.prefix).parent.mkdir(parents=True,exist_ok=True)
            print('CNF variables, clauses:',A.export_cnf(o.prefix+'.cnf',o.total_active)); return
        print('Kissat log:',str(Path(o.prefix+'.log').resolve()),flush=True)
        answer=A.has_trail(o.total_active,o.timeout)
        print(answer, 'seconds=',A.solving_time())
        if answer is True:
            print(A.trail()); print(json.dumps(A.result(),indent=2))
        elif answer is not False:
            raise SystemExit(2)
    except (ValueError,RuntimeError,OSError) as exc:
        p.exit(2,f'Error: {exc}\n')

if __name__=='__main__': main()

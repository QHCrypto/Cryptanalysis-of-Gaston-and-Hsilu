from sage.all import *
from sage.crypto.sboxes import SBox
from copy import deepcopy
from sage.crypto.boolean_function import BooleanFunction


# The linear layer of Ascon
def LT(state_lt):
    def ROR(v, offset):
        l = list(v)
        return vector(l[-offset:]+l[:-offset])
    state_lt[0] += ROR(state_lt[0], 19) + ROR(state_lt[0], 28)
    state_lt[1] += ROR(state_lt[1], 61) + ROR(state_lt[1], 39)
    state_lt[2] += ROR(state_lt[2], 1) + ROR(state_lt[2], 6)
    state_lt[3] += ROR(state_lt[3], 10) + ROR(state_lt[3], 17)
    state_lt[4] += ROR(state_lt[4], 7) + ROR(state_lt[4], 41)
    return state_lt


# Round constants addition, may omitted.
def PC(state_pc, Round):
    rc = RoundConstants[Round]
    rc = bin(rc)[2:].zfill(64)
    rc = list(map(lambda x: int(x), rc))
    rc = vector(GF(2), rc)
    state_pc[2] += rc
    return state_pc


def PA(state_pa, r):
    state_pa = PC(state_pa, r)
    state_pa[0] += state_pa[4]
    state_pa[4] += state_pa[3]
    state_pa[2] += state_pa[1]
    return state_pa


def Psn(state_pb):
    state = deepcopy(state_pb)

    def comprod(v1, v2):
        v = [v1[i] * v2[i] for i in range(len(v1))]
        return vector(v)
    state_pb[0] = state[0] + comprod(state[2], state[1]+vector(GF(2), [1]*64))
    state_pb[1] = state[1] + comprod(state[3], state[2]+vector(GF(2), [1]*64))
    state_pb[2] = state[2] + comprod(state[4], state[3]+vector(GF(2), [1]*64))
    state_pb[3] = state[3] + comprod(state[0], state[4]+vector(GF(2), [1]*64))
    state_pb[4] = state[4] + comprod(state[1], state[0]+vector(GF(2), [1]*64))
    state_pb[1] += state_pb[0]
    state_pb[3] += state_pb[2]
    state_pb[0] += state_pb[4]
    state_pb[2] += vector(GF(2), [1]*64)
    return state_pb


def PB(state_pb):
    state_pb = Psn(state_pb)
    state_pb = LT(state_pb)
    return state_pb


RoundConstants = [0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5]
S = [4, 11, 31, 20, 26, 21, 9, 2, 27, 5, 8, 18, 29, 3, 6, 28,
     30, 19, 7, 14, 0, 13, 17, 24, 16, 12, 1, 25, 22, 10, 15, 23]
anf = [0] * 5
for j in range(5):
    anf[j] = SBox(S).component_function(1 << (j)).algebraic_normal_form()
sbox_anf = deepcopy(anf)


def IsolatedVariables(f):
    return list(filter(lambda x: f/x == 1, f.variables()))


# Calculate the bias of a boolean polynomial by calculate its hamming weight
def fDis(f):
    if f == 1:
        return -1/2
    if f == 0:
        return 1/2
    if len(IsolatedVariables(f)) > 0:
        return 0
    varf = f.variables()
    B = BooleanPolynomialRing(len(varf), varf)
    h = B(f)
    h = BooleanFunction(h)
    onesf = sum(h.truth_table())
    return 1/2 - onesf/2**h.nvariables()


# Separate a boolean polynomial as a sum of polynomials whose variable sets are mutually
def Separate(f):
    if f.nvariables() <= 1:
        return [f]
    varf = f.variables()
    fp = (f/varf[0]) * varf[0]
    fq = f.subs({varf[0]: 0})
    varintersec = set(fp.variables()).intersection(fq.variables())
    while len(varintersec) > 0:
        for var in varintersec:
            fp += (fq/var) * var
            fq = fq.subs({var: 0})
        varintersec = set(fq.variables()).intersection(fp.variables())
    if fq.is_constant():
        return [fp+fq]
    else:
        return [fp] + Separate(fq)


# Calculate the bias of a boolean polynomial using Shannon decomposition and piling-up lemma
def Bias(f, depth=0):
    assert (depth < 20)
    if f == 1:
        return -1/2
    if f == 0:
        return 1/2
    if len(IsolatedVariables(f)) > 0:
        return 0
    # When the number of variables is small, calculate the bias by hamming weight directly
    if f.nvariables() <= 27:
        return fDis(f)
    fl = Separate(f)
    epsilon = 1/2
    for fi in fl:
        if fi.nvariables() <= 27:
            epsiloni = fDis(fi)
        else:
            # For complex polynomials, select a variable v that minimizing
            # the maximum caridinality of the variable sets of the polynomials
            # in Separate(f_{v=0}) and Separate(f_{v=1})
            maxnvar = fi.nvariables()
            maxvar = 0
            for var in fi.variables():
                f0nvar = max(map(lambda x: x.nvariables(),
                             Separate(fi.subs({var: 0}))))
                f1nvar = max(map(lambda x: x.nvariables(),
                             Separate(fi.subs({var: 1}))))
                if max(f0nvar, f1nvar) < maxnvar:
                    maxnvar = max(f0nvar, f1nvar)
                    maxvar = var
            # print(maxvar, maxnvar)

        # piling-up lemma
            epsiloni = Bias(fi.subs({maxvar: 0}), depth+1) + \
                Bias(fi.subs({maxvar: 1}), depth+1)
            epsiloni *= 1/2
        epsilon *= 2*epsiloni
        if epsilon == 0:
            return 0
    return epsilon


fullround = 5
diff_in = 63 - 63
BPR = BooleanPolynomialRing(128, "k")
BPR = BooleanPolynomialRing(
    BPR.ngens() + 128, BPR.gens()+BooleanPolynomialRing(128, "n").gens())
for r in range(1, fullround):
    BPR = BooleanPolynomialRing(
        BPR.ngens() + 320, BPR.gens()+BooleanPolynomialRing(320, "r"+str(r)+"a").gens())
    BPR = BooleanPolynomialRing(
        BPR.ngens() + 320, BPR.gens()+BooleanPolynomialRing(320, "r"+str(r)+"b").gens())
BPR = BooleanPolynomialRing(
    BPR.ngens() + 1, BPR.gens()+BooleanPolynomialRing(1, "diff").gens())
x = BPR.gens()[-1]
IV = vector(GF(2), [bin(0x80400c0600000000)[2:][i] for i in range(64)])
state = Matrix([[BPR(0)]*64] * 5)
state[0] = IV
for i in range(1, 5):
    state[i] = vector(BPR.gens()[64*(i-1):64*i][::-1])

# k63=0
state[1] = vector([state[1][i] if i != diff_in else BPR(0) for i in range(64)])
# k127=0
state[2] = vector([state[2][i] if i != diff_in else BPR(0) for i in range(64)])
state[3] += vector([x if i == diff_in else BPR(0) for i in range(64)])
# n63=n127
state[4] += vector([state[4][i]+state[3][i] if i ==
                   diff_in else BPR(0) for i in range(64)])

state = PA(state, 0)
state = PB(state)
state = PA(state, 1)
Conditions = set()

# When calculating the bias that all the 9 polynomials equals 0, comment the following conditions
Conditions.add(BPR.gens()[256+127])
Conditions.add(BPR.gens()[256+300] + 1)
Conditions.add(BPR.gens()[256+108])
Conditions.add(BPR.gens()[256+99])

I = BPR.ideal(list(Conditions))
alphaQ = {}
Q = {}
for r in range(1, fullround):
    for i in range(5):
        tmps = [0] * 64
        for j in range(64):
            alpha = state[i][j].subs({x: 0})
            alpha = alpha.reduce(I)
            beta = state[i][j]/x
            beta = beta.reduce(I)
            varalpha = BPR.gens()[r*640+(i * 64 + 63-j) - 384]
            varbeta = BPR.gens()[r*640+(i * 64 + 63-j) - 64]
            tmps[j] = BPR(0)
            if alpha.nvariables() > 1:
                tmps[j] += varalpha
                alphaQ[varalpha] = alpha
            else:
                tmps[j] += alpha
            if beta in Q.values():
                beta = list(filter(lambda x: Q[x] == beta, Q.keys()))[0]
            elif beta + 1 in Q.values():
                beta = list(filter(lambda x: Q[x] == beta+1, Q.keys()))[0] + 1

            # When calculating the bias that all the 9 polynomials equals 0, uncomment the following 4 lines
            # if r <= 2 and not beta.is_constant():
            #     Conditions.add(beta)
            #     I = BPR.ideal(list(Conditions))
            #     continue
            if beta.nvariables() > 1:
                tmps[j] += varbeta * x
                Q[varbeta] = beta
            else:
                tmps[j] += beta * x
        state[i] = vector(tmps)
    if r == fullround-1:
        state = Psn(state)
        break
    state = PB(state)
    state = PA(state, r+1)
print("Conditions:", Conditions)
print("")
for condition in Conditions:
    print(condition, "=", condition.subs(alphaQ).subs(alphaQ))
print("")
f = state[0][63-16]
f = f/x
f = f.subs(Q).reduce(I).subs(Q).reduce(I).subs(Q).reduce(I)
biasf = Bias(f)
print("Bias:", biasf, log(abs(biasf), 2).n())

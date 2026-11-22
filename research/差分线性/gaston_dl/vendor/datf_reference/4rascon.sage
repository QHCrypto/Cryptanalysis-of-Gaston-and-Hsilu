from sage.all import *
from sage.crypto.sboxes import SBox
from copy import deepcopy
from sage.crypto.boolean_function import BooleanFunction


RoundConstants = [0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5]
S = [4, 11, 31, 20, 26, 21, 9, 2, 27, 5, 8, 18, 29, 3, 6, 28,
     30, 19, 7, 14, 0, 13, 17, 24, 16, 12, 1, 25, 22, 10, 15, 23]
anf = [0] * 5
for j in range(5):
    anf[j] = SBox(S).component_function(1 << (j)).algebraic_normal_form()
# Calculate the Algebraic Normal Form (ANF) of the S-box of Ascon
sbox_anf = deepcopy(anf)


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


# The non-linear layer of Ascon
def NLT(state_nlt):
    state_nlt = state_nlt.transpose()
    x = sbox_anf[0].parent().gens()
    for i in range(64):
        input_anf = {}
        for j in range(5):
            input_anf[x[4-j]] = state_nlt[i][j]
        output_anf = [0] * 5
        for j in range(5):
            output_anf[j] = sbox_anf[4-j].subs(input_anf)
        state_nlt[i] = vector(output_anf)
    return state_nlt.transpose()


# Round constants addition, may omitted.
def PC(state_pc, Round):
    rc = RoundConstants[Round]
    rc = bin(rc)[2:].zfill(64)
    rc = list(map(lambda x: int(x), rc))
    rc = vector(GF(2), rc)
    state_pc[2] += rc
    return state_pc


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


startround = 0
fullround = 4

# The array use little endian format, opposite to big endian in paper.
diff_in = 63 - 63
mask_out = 63 - 9

# Construct a boolean polynomial ring
BPR = BooleanPolynomialRing(128, "k")
BPR = BooleanPolynomialRing(
    BPR.ngens() + 128, BPR.gens()+BooleanPolynomialRing(128, "n").gens())
for r in range(1, fullround):
    # For each round, the alpha variables are denoted as r*a*
    BPR = BooleanPolynomialRing(
        BPR.ngens() + 320, BPR.gens()+BooleanPolynomialRing(320, "r"+str(r)+"a").gens())
    # And the beta variables are denoted as r*b*
    BPR = BooleanPolynomialRing(
        BPR.ngens() + 320, BPR.gens()+BooleanPolynomialRing(320, "r"+str(r)+"b").gens())
BPR = BooleanPolynomialRing(
    BPR.ngens() + 1, BPR.gens()+BooleanPolynomialRing(1, "diff").gens())
x = BPR.gens()[-1]
IV = vector(GF(2), [bin(0x80400c0600000000)[2:][i] for i in range(64)])
state = Matrix([[BPR(0)]*64] * 5)
# Initialize states with IVs, keys, and nonces
state[0] = IV
for i in range(1, 5):
    state[i] = vector(BPR.gens()[64*(i-1):64*i][::-1])
# state[1] = vector([state[1][i] if i != diff_in else BPR(0) for i in range(64)])
# state[2] = vector([state[2][i] if i != diff_in else BPR(1) for i in range(64)])

# Set difference on x_3 and x_4
state[3] += vector([x if i == diff_in else BPR(0) for i in range(64)])
state[4] += vector([x if i == diff_in else BPR(0) for i in range(64)])
state = PC(state, 0)
state = NLT(state)
state = LT(state)

# betaQ to store DATF of beta variables
betaQ = {}
# alphaQ to store DATF of alpha variables
alphaQ = {}
for r in range(startround + 1, startround + fullround):
    for i in range(5):
        tmps = [0] * 64
        for j in range(64):
            alpha = state[i][j].subs({x: 0})
            beta = state[i][j]/x
            varalpha = BPR.gens()[r*640+(i * 64 + 63-j) - 384]
            varbeta = BPR.gens()[r*640+(i * 64 + 63-j) - 64]
            tmps[j] = BPR(0)
            # Transitional Rule
            if alpha.nvariables() > 1:
                tmps[j] += varalpha
                alphaQ[varalpha] = alpha
            else:
                tmps[j] += alpha
            # If beta or beta + 1 is already in the value of betaQ, then use the same beta variable
            if beta in betaQ.values():
                beta = list(
                    filter(lambda x: betaQ[x] == beta, betaQ.keys()))[0]
            elif beta + 1 in betaQ.values():
                beta = list(
                    filter(lambda x: betaQ[x] == beta + 1, betaQ.keys()))[0] + 1
            if beta.nvariables() > 1:
                tmps[j] += varbeta * x
                betaQ[varbeta] = beta
            else:
                tmps[j] += beta * x
        state[i] = vector(tmps)
    state = PC(state, r)
    state = NLT(state)
    if r == startround + fullround - 1:
        break
    state = LT(state)

f = state[0][mask_out]
f = f/x
print("Original DATF:", f)
g = f.subs(betaQ).subs(betaQ).subs(betaQ)
print("Expanded DATF:", g)
print("Bias: ", Bias(g), log(Bias(g), 2).n())

# For verify the independence, output the ANF of variables in g
#  and separate the isolated variables

# for varg in g.variables()[1:]:
#     varganf = varg.subs(alphaQ).subs(alphaQ).subs(alphaQ)
#     lvarg = set(IsolatedVariables(BPR(varganf)))
#     print(varg,str(varg)[1],int(str(varg)[3:])//64,int(str(varg)[3:])%64)
#     print(lvarg,set(varganf.variables())-lvarg)

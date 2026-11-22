from sage.all import *
from sage.crypto.sboxes import SBox
from time import perf_counter
from copy import deepcopy
from sage.crypto.boolean_function import BooleanFunction

S = [[3, 8, 15, 1, 10, 6, 5, 11, 14, 13, 4, 2, 7, 0, 9, 12],  # S0
     [15, 12, 2, 7, 9, 0, 5, 10, 1, 11, 14, 8, 6, 13, 3, 4],  # S1
     [8, 6, 7, 9, 3, 12, 10, 15, 13, 1, 14, 4, 0, 11, 5, 2],  # S2
     [0, 15, 11, 8, 12, 9, 6, 3, 13, 1, 2, 4, 10, 7, 5, 14],  # S3
     [1, 15, 8, 3, 12, 0, 11, 6, 2, 5, 4, 10, 9, 14, 7, 13],  # S4
     [15, 5, 2, 11, 4, 10, 9, 12, 0, 3, 14, 8, 13, 6, 7, 1],  # S5
     [7, 2, 12, 5, 8, 4, 6, 11, 14, 9, 1, 15, 13, 3, 10, 0],  # S6
     [1, 13, 15, 0, 14, 8, 2, 11, 7, 4, 12, 10, 9, 3, 5, 6],  # S7
     ]

# Calculate the Algebraic Normal Form (ANF) of the S-box of each round
sbox_anf = [0] * 8
for i in range(8):
    anf = [0] * 4
    for j in range(4):
        anf[j] = SBox(S[i]).component_function(1 << j).algebraic_normal_form()
    sbox_anf[i] = deepcopy(anf)


# The linear layer of Serpent
def LT(state_lt):
    def lll(v, offset):
        l = list(v)
        return vector(l[offset:]+l[:offset])
    state_lt[0] = lll(state_lt[0], 13)
    state_lt[2] = lll(state_lt[2], 3)
    state_lt[1] += state_lt[0] + state_lt[2]
    state_lt[3] += state_lt[2] + vector(list(state_lt[0][3:]) + [0] * 3)
    state_lt[1] = lll(state_lt[1], 1)
    state_lt[3] = lll(state_lt[3], 7)
    state_lt[0] += state_lt[1] + state_lt[3]
    state_lt[2] += state_lt[3] + vector(list(state_lt[1][7:]) + [0] * 7)
    state_lt[0] = lll(state_lt[0], 5)
    state_lt[2] = lll(state_lt[2], 22)
    return state_lt


# The non-linear layer of Serpent
def NLT(state_nlt, roundi):
    state_nlt = state_nlt.transpose()
    x = sbox_anf[roundi][0].parent().gens()
    for i in range(32):
        input_anf = {}
        for j in range(4):
            input_anf[x[j]] = state_nlt[i][j]
        output_anf = [0] * 4
        for j in range(4):
            output_anf[j] = sbox_anf[roundi][j].subs(input_anf)
        state_nlt[i] = vector(output_anf)
    return state_nlt.transpose()


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


startround = 3
fullround = 3

# the array use little endian format, opposite to big endian in paper.
# For convenience, X_i corresponds to [32*i : 32*i + 32]
diff_in = 66
mask_out = 120
# Set the output mask to 124 will get a better distinguisher (The running time will increase to about 1 minutes)

# There are no beta variables of startround
BPR = BooleanPolynomialRing(128, 'r'+str(startround)+'a')
for r in range(startround+1, startround+fullround+1):
    # For each round, the alpha variables are denoted as r*a*
    BPR = BooleanPolynomialRing(
        BPR.ngens() + 128, BPR.gens() + BooleanPolynomialRing(128, 'r'+str(r)+'a').gens())
    # And the beta variables are denoted as r*b*
    BPR = BooleanPolynomialRing(
        BPR.ngens() + 128, BPR.gens() + BooleanPolynomialRing(128, 'r'+str(r)+'b').gens())
BPR = BooleanPolynomialRing(
    BPR.ngens() + 1, BPR.gens() + BooleanPolynomialRing(1, 'diff').gens())
x = BPR.gens()[-1]

# betaQ to store DATF of beta variables
betaQ = {}
# alphaQ to store DATF of alpha variables
alphaQ = {}


# Initialize states and set difference
state = Matrix([BPR.gens()[i*32:i*32+32][::-1] for i in range(4)])
state[diff_in//32] += vector([x if i == diff_in %
                             32 else BPR(0) for i in range(32)])
state = NLT(state, startround)
state = LT(state)
starttime = perf_counter()
for r in range(startround + 1, startround + fullround):
    for i in range(4):
        tmps = [0] * 32
        for j in range(32):
            alpha = state[i][j].subs({x: 0})
            beta = state[i][j]/x
            tmps[j] = BPR(0)
            # Transitional Rule
            if alpha.nvariables() > 1:
                varalpha = BPR.gens()[(r-startround)*256-128 + (i * 32 + 31-j)]
                tmps[j] += varalpha
                alphaQ[varalpha] = alpha
            # If beta or beta + 1 is already in the value of betaQ, then use the same beta variable
            if beta in betaQ.values():
                beta = list(
                    filter(lambda x: betaQ[x] == beta, betaQ.keys()))[0]
            elif beta + 1 in betaQ.values():
                beta = list(
                    filter(lambda x: betaQ[x] == beta+1, betaQ.keys()))[0] + 1
            if beta.nvariables() > 1:
                varbeta = BPR.gens()[(r-startround)*256 + (i * 32 + j)]
                tmps[j] += varbeta * x
                betaQ[varbeta] = beta
            else:
                tmps[j] += beta * x
        state[i] = vector(tmps)
    state = NLT(state, r)
    state = LT(state)

f = state[mask_out//32][mask_out % 32]
f = f/x
print("Original DATF:", f)
g = f.subs(betaQ).subs(betaQ)
print("Variables in the expanded DATF:")
print(g.nvariables(), g.variables())

# The following codes are used to simplify the DATF,
# the main goal is to reduce the number of variables
r5beta = BPR.gens()[(fullround-1)*256:(fullround-1)*256+128]
fr5beta = set(f.variables()).intersection(r5beta)
r4beta = BPR.gens()[(fullround-2)*256:(fullround-2)*256+128]
r4beta = set(betaQ.keys()).intersection(r4beta)

# Count the number of occurrences of r4 alpha variables
ReferenceCount = {}
for var in fr5beta:
    g = betaQ[var]
    for r4betavar in r4beta:
        for r4avar in (g/r4betavar).variables():
            if r4avar in ReferenceCount.keys():
                ReferenceCount[r4avar] += 1
            else:
                ReferenceCount[r4avar] = 1
lessvar = True
while lessvar:
    lessvar = False
    for var in fr5beta:
        g = betaQ[var]
        newQ = BPR(0)
        for r4betavar in r4beta:
            if (g/r4betavar) == 0:
                continue
            if (g/r4betavar).nvariables() == 1:
                newQ += g/r4betavar*r4betavar
                continue
            # If a alpha variable is isolated variable, and it only occurs in this expression, then this expression is balanced
            lv = list(
                filter(lambda x: ReferenceCount[x] == 1, IsolatedVariables(g/r4betavar)))
            if len(lv) > 0:
                newQ += r4betavar * lv[0]
                lessvar = True
            # If some alpha variables only occurs in this expression and the occured terms are balaced, then the whole expression is balanced
            elif sum([ReferenceCount[x] for x in (g/r4betavar).variables()]) == (g/r4betavar).nvariables() and Bias(g/r4betavar) == 0:
                newQ += r4betavar * (g/r4betavar).variable(0)
                lessvar = True
            else:
                newQ += g/r4betavar*r4betavar
        betaQ[var] = newQ
    # Update the counters for further simplify
    ReferenceCount = {}
    for var in fr5beta:
        g = betaQ[var]
        for r4betavar in r4beta:
            for r4avar in (g/r4betavar).variables():
                if r4avar in ReferenceCount.keys():
                    ReferenceCount[r4avar] += 1
                else:
                    ReferenceCount[r4avar] = 1

g = f.subs(betaQ).subs(betaQ)
print("Variables in the equivalent DATF:")
print(g.nvariables(), g.variables())
biasg = Bias(g)
print("Bias", biasg, log(biasg, 2).n())
print(perf_counter()-starttime, "s")

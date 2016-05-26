from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public
from sympy.combinatorics.free_group import FreeGroupElm


@public
def fp_group(fr_grp, relators=[]):
    _fp_group = FpGroup(fr_grp, relators)
    return (_fp_group,) + tuple(_fp_group._generators)

@public
def xfp_group(fr_grp, relators=[]):
    _fp_group = FpGroup(fr_grp, relators)
    return (_fp_group, _fp_group._generators)

@public
def vfp_group(fr_grpm, relators):
    _fp_group = FpGroup(symbols, relators)
    pollute([sym.name for sym in _fp_group.symbols], _fp_group.generators)
    return _fp_group


def _parse_relators(rels):
    """Parse the passed relators."""
    return rels

_fp_group_cache = {}

class FpGroup(Basic):
    """
    The FpGroup would take a FreeGroup and a list/tuple of relators, the
    relators would be specified in such a way that each of them be equal to the
    identity of the provided free group.
    """
    is_FpGroup = True

    def __new__(cls, fr_grp, relators):
        relators = _parse_relators(relators)
        # return the corresponding FreeGroup if no relators
        # is specified
        if not relators:
            return fr_grp
        obj = Basic.__new__(cls, fr_grp, relators)
        obj._free_group = fr_grp
        obj._relators = relators
        obj.generators = obj._generators()
        obj.dytype = type("FpGroupElm", (FpGroupElm,), {"group": obj})
        return obj

    @property
    def free_group(self):
        return self._free_group

    def relators(self):
        return tuple(self._relators)

    def _generators(self):
        """Returns the generators of the associated free group."""
        return self.free_group.generators

    def subgroup(self, gens):
        # gens is list of generators of self.
        pass

    def coset_enumeration_c(G, Y):
        Y_copy = []
        for relator in Y:
            Y_copy.append(relator.identity_cyclic_reduction())
        Y = Y_copy

    def coset_table(self, H):
        pass


# sets the upper limit on the number of cosets generated during
# Coset Enumeration. "M" from Derek Holt's. It is supposed to be
# user definable.
DefaultMaxLimit = 500

# fp_grp: Finitely Presented Group with < X|R > as presentation.
# H: subgroup of fp_grp.
# NOTE: We start with H as being only a list of words in generators
#       of "fp_grp". Since `.subgroup` method has not been implemented.

def coset_table(fp_grp, H):
    # list of generators and their inverses, should be made a tuple.
    A = []
    for gen in fp_grp:
        A.extend([gen, gen.inverse()])
    r = len(A)
    # start with table of 1 row (can change) and 2*(no. of gens) coloumns (fixed)
    table = [r*[-1]]

# coset_table: Mathematically a coset table
#               represented using a list of lists
# alpha: Mathematically a coset (precisely, a live coset)
#       represented by an integer between i with 1 <= i <= n
#       α ∈ c
# x: Mathematically an element of "A" (set of generators and
#   their inverses), represented using "FpGroupElm"

def define(C, alpha, x):
    # now C.n is now obtainable
#    if C.n == DefaultMaxLimit:
#        # abort the further generation of cosets
#        return
#    C.n += 1
    # beta is the new coset generated
    beta = n
    C.p[beta] = beta
    C[alpha][C.A.index(x)] = beta
    C[alpha][C.A.index(x.inverse())] = alpha


def scan(C, alpha, word, A):
    """
    >>> F, x, y = free_group("x, y")
    >>> A = [x, x**-1, y, y**-1]

    # Example 5.1, Pg 150
    >>> c = [[0, 0, 1, 2], [-1, -1, -1, 0], [-1, -1, 0, -1]]
    >>> scan(c, 1, y**3, A)
    # this lead to deductions 2^y = 3, 3^(y^-1) = 2 (1, 2, 3 represent cosets in this comment)
    >>> c
    [[0, 0, 1, 2], [-1, -1, 2, 0], [-1, -1, 0, 1]]
    >>> scan(c, 0, x**-1*y**-1*x*y, A)
    >>> c
    [[0, 0, 1, 2], [-1, -1, 2, 0], [2, 2, 0, 1]]
    >>> scan(c, 1, x**-1*y**-1*x*y, A)
    >>> c
    [[0, 0, 1, 2], [1, 1, 2, 0], [2, 2, 0, 1]]

    # Example 5.2, Pg 154
    >>> c1 = [[1, 1, -1, -1], [0, 0, -1, -1]]
    >>> scan(c1, 0, x*y, A)
    >>> c1
    [[1, 1, -1, 1], [0, 0, 0, -1]]

    >>> c2 = [[1, 1, 2, 1], [0, 0, 0, -1], [-1, -1, -1, 0]]
    >>> scan(c2, 1, y**3, A)
    >>> c2
    [[1, 1, 2, 1], [0, 0, 0, 2], [-1, -1, 1, 0]]

    >>> c3 = [[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0], [2, 2, -1, -1]]
    >>> scan(c4, 2, (x*y)**3, A)
    >>> c3
    [[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0], [2, 2, 3, 3]]

    """
    # alpha is an integer representing a "coset"
    # since scanning can be in two cases
    # 1. for alpha=0 and w in Y (i.e generating set of H)
    # 2. alpha in sigma (set of live cosets), w in R (relators)
    f = alpha
    i = 0
    r = len(word)
    # list of union of generators and their inverses
    # A = coset_table.A
    while i < r and C[f][A.index(word.subword(i, i+1))] != -1:
        f = C[f][A.index(word.subword(i, i+1))]
        i += 1
    if i >= r:
        if f != alpha:
            # implement the "coincidence" routine on Pg 158 of Handbook.
            coincidence(C, f, alpha)
            return
    b = alpha
    j = r
    while j >= i and C[b][A.index(word.subword(j-1, j).inverse())] != -1:
        b = C[b][A.index(word.subword(j-1, j).inverse())]
        j -= 1
    if j < i:
        # run the "coincidence" routine
        coincidence(C, f, b)
    elif j == i + 1:
        # deduction process
        C[f][A.index(word.subword(i, i+1))] = b
        C[b][A.index(word.subword(i, i+1).inverse())] = f

# here alpha, beta represent the pair of cosets where
# coicidence occurs

def merge(k, lamda, p, q, l):
    phi = rep(k, p)
    chi = rep(lamda, p)
    if phi != chi:
        mu = min(phi, chi)
        v = max(phi, chi)
        p[v] = mu
        l += 1
        q[l] = v

def rep(k, p):
    lamda = k
    pho = p[lamda]
    while p != lamda:
        lamda = pho
        pho = p[lamda]
    mu = k
    pho = p[mu]
    while pho != lamda:
        p[mu] = lamda
        mu = pho
        pho = p[mu]
    return lamda

def coincidence(C, alpha, beta):
    l = 0
    q = []
    merge(alpha, beta, p, q, l)
    i = 1
    while i <= l:
        gamma = q[i]
        i += 1
        del C[gamma]
        # commenting this out
        # sigma = sigma - set(gamma)
        for ind in range(len(A)):
            if C[gamma][ind] != -1:
                delta = C[gamma][ind]
                # index of x**-1 in A is index(x)+1
                C[delta][ind + 1] = -1
                mu = rep(gamma, p)
                v = rep(delta, p)
                if C[mu][ind] == -1:
                    merge(v, C[mu][ind], p, q, l)
                elif C[v][ind + 1] != -1:
                    merge(mu, C[v][ind + 1], p, q, l)
                else:
                    C[mu][ind] = v
                    C[v][ind + 1] = mu


# method used in the HLT strategy
def scan_and_fill(C, alpha, word):
    f = alpha
    i = 0
    b = alpha
    j = len(word)
    while True:
        # do the forward scanning
        while i < r and C[f][A.index(word.subword(i, i+1))] != -1:
            f = C[f][A.index(word.subword(i, i+1))]
            i += 1
        if i >= r:
            if f != alpha:
                coincidence(C, f, alpha)
                return
        # forward scan was incomplete, scan backwards
        while j >= i and C[b][A.index(word.subword(i, i+1))] != -1:
            b = C[b][A.index(word.subword(i, i+1))]
            j -= 1
        if j < i:
            coincidence(C, f, b)
        elif j == i + 1:
            C[f][A.index(word.subword(i, i+1))] = b
            C[b][A.index(word.subword(i, i+1).inverse())] = f
        else:
            define(C, f, word.subword(i, i+1))


# Pg. 166
def process_deductions(C, R):
    while deduction_stack is not empty:
        if deduction_stack is full:
            lookahead(C, R)
            empty_deduction_stack
        else:
            deduction_stack.pop((alpha, x))
            if p[alpha] == alpha:
                for w in R_c:
                    scan(C, alpha, w)
                    if p[alpha] < alpha:
                        break
        beta = C[alpha][x]
        if p[beta] == beta:
            for w in R_c_x_inv:
                scan(C, beta, w)
                if p[beta] < beta:
                    break

def standardize(C):
    gamma = 2
    n = C.n
    for a in range(n):
        x in A
        beta = C[alpha][A.index(x)]*gamma
        if beta >= gamma:
            if beta > gamma:
                switch(C, gamma, beta)*gamma
                gamma += 1
                if gamma == n:
                    return


# Pg. 166
def coset_enumeration_c(fp_grp, Y):
    X = fp_grp.generators
    R = fp_grp.relators
    for i in range(len(R)):
        R[i] = R[i].identity_cyclic_reduction()
    R_union_R_inverse = R + [rel.inverse() for rel in R]
    p = []
    p[0] = 1
    n = 0
    for w in Y:
        scan_and_fill(C, 0, w)
    process_deductions(C, R)
    for alpha in range(len(C)):
        for x in range(len(A)):
            if C[alpha][x] == -1:
                define(C, alpha, x)
                process_deductions(C, R)
    return C

# Pg. 167 5.2.3
def compress(C):
    gamma = 0
    for alpha in range(C.n):
        if p[alpha] == alpha:
            gamma += 1
            if gamma != alpha:
                for x in A:
                    beta = C[alpha][A.index(x)]
                    if beta == alpha:
                        beta = gamma
                    C[gamma][A.index(x)] = beta
                    C[beta][A.index(x**-1)] == gamma
    C.n = gamma
    for alpha in range(C.n):
        p[alpha] = alpha

FpGroupElm = FreeGroupElm

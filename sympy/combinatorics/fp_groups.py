# -*- coding: utf-8 -*-
from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public
from sympy.utilities.iterables import flatten
from sympy.combinatorics.free_group import FreeGroupElm
from itertools import chain


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

    def coset_table(self, Y):
        R = self.relators


DefaultMaxLimit = 2000

class CosetTable(Basic):
    # coset_table: Mathematically a coset table
    #               represented using a list of lists
    # alpha: Mathematically a coset (precisely, a live coset)
    #       represented by an integer between i with 1 <= i <= n
    #       α ∈ c
    # x: Mathematically an element of "A" (set of generators and
    #   their inverses), represented using "FpGroupElm"
    # fp_grp: Finitely Presented Group with < X|R > as presentation.
    # H: subgroup of fp_grp.
    # NOTE: We start with H as being only a list of words in generators
    #       of "fp_grp". Since `.subgroup` method has not been implemented.

    # sets the upper limit on the number of cosets generated during
    # Coset Enumeration. "M" from Derek Holt's. It is supposed to be
    # user definable.

    """

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of Computational Group Theory"

    [2] John J. Cannon; Lucien A. Dimino; George Havas; Jane M. Watson
    Mathematics of Computation, Vol. 27, No. 123. (Jul., 1973), pp. 463-490.
    "Implementation and Analysis of the Todd-Coxeter Algorithm"

    """

    def __init__(self, fp_grp, subgroup):
        self.fp_group = fp_grp
        self.subgroup = subgroup
        self.p = [0]
        self.A = list(chain.from_iterable((gen, gen**-1) \
                for gen in self.fp_group.generators))
        self.table = [[None]*len(self.A)]
        self.A_dict = {x: self.A.index(x) for x in self.A}

    @property
    def omega(self):
        return [coset for coset in self.p if self.p[coset] == coset]

    @property
    def n(self):
        """The number `n` represents the largest number that has been used so
        far for a live coset
        """
        return max(self.omega) + 1

    # checks whether the Coset Table is complete or not
    def is_complete(self):
        return not None in flatten(self.table)

    # Pg. 153
    def define(self, alpha, x):
        A = self.A
        if self.n == DefaultMaxLimit:
            # abort the further generation of cosets
            return
        self.table.append([None]*len(A))
        # beta is the new coset generated
        beta = len(self.table) - 1
        self.p.append(beta)
        self.table[alpha][A.index(x)] = beta
        self.table[beta][A.index(x.inverse())] = alpha

    def scan(self, alpha, word):
        """
        >>> from sympy.combinatorics.free_group import free_group
        >>> F, x, y = free_group("x, y")
        >>> f = free_group(F, [x**3, y**3, x**-1*y**-1*x*y])
        >>> c = CosetTable(f, [x])

        # Example 5.1, Pg 150
        >>> c = [[0, 0, 1, 2], [None, None, None, 0], [None, None, 0, None]]
        >>> scan(c, 0, y**3, A)
        # this lead to deductions 2^y = 3, 3^(y^-1) = 2 (1, 2, 3 represent cosets in this comment)
        >>> c
        [[0, 0, 1, 2], [None, None, 2, 0], [None, None, 0, 1]]
        >>> scan(c, 0, x**-1*y**-1*x*y, A)
        >>> c
        [[0, 0, 1, 2], [None, None, 2, 0], [2, 2, 0, 1]]
        >>> c.scan(1, x**-1*y**-1*x*y)
        >>> c
        [[0, 0, 1, 2], [1, 1, 2, 0], [2, 2, 0, 1]]

        # Example 5.2, Pg 154
        >>> c1 = [[1, 1, None, None], [0, 0, None, None]]
        >>> c1.scan(0, x*y)
        >>> c1
        [[1, 1, None, 1], [0, 0, 0, None]]

        >>> c2 = [[1, 1, 2, 1], [0, 0, 0, None], [None, None, None, 0]]
        >>> c2.scan(1, y**3)
        >>> c2
        [[1, 1, 2, 1], [0, 0, 0, 2], [None, None, 1, 0]]

        >>> c3 = [[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0], [2, 2, None, None]]
        >>> c3.scan(2, (x*y)**3)
        >>> c3
        [[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0], [2, 2, 3, 3]]

        # Example 5.3
        >>> c_1 = [[1, None, None, None], [2, 0, None, None],
                    [None, 1, 3, None], [None, None, None, 2]]
        >>> c_1.scan(0, x**2*y**2)
        >>> c_1
        [[1, None, None, 3],
         [2, 0, None, None],
         [None, 1, 3, None],
         [None, None, 0, 2]]
        >>> c_2 = [[1, None, None, 3], [2, 0, None, None],
                   [4, 1, 3, None], [None, None, 0, 2],
                   [None, 2, 5, None], [None, None, None, 4]]
        >>> c_2.scan(1, x**2*y**2)
        >>> c_2
        [[1, None, None, 3],
         [2, 0, None, 5],
         [4, 1, 3, None],
         [None, None, 0, 2],
         [None, 2, 5, None],
         [None, None, 1, 4]]
        >>> c_3 = c_2.copy()
        >>> c_3.scan(0, x**3*y**5)
        >>> c_3
        [[1, None, None, 3],
         [2, 0, 2, 5],
         [4, 1, 3, 1],
         [None, None, 0, 2],
         [None, 2, 5, None],
         [None, None, 1, 4]]
        >>> c_4 = c_3.copy()
        >>> c_4.scan(2, x**2*y**2)
        >>> c_4
        [[1, None, None, 3],
         [2, 0, 2, 5],
         [4, 1, 3, 1],
         [None, None, 0, 2],
         [5, 2, 5, None],
         [None, 4, 1, 4]]
        >>> c_4.scan(1, x**3*y**5)
        >>> c_4
        [[1, 3, 1, 3],
         [2, 0, 2, 0],
         [3, 1, 3, 1],
         [0, 2, 0, 2],
         [0, 2, None, None],
         [None, 4, 1, 4]]

        # Example 5.4
        >>> F, a, b, c, d, e = free_group("a, b, c, d, e")
        >>> A = [a, a**-1, b, b**-1, c, c**-1, d, d**-1, e, e**-1]
        >>> c1 = [[0, None, 1, None, None, None, None, None, None, None],
                   [None, None, None, 0, None, None, None, None, None, None],
                   [None, None, None, None, None, None, None, None, None, None]]
        >>> c1.scan(0, a*b*c**-1)
        >>> c1
        [[0, None, 1, None, 1, None, None, None, None, None],
         [None, None, None, 0, None, 0, None, None, None, None],
         [None, None, None, None, None, None, None, None, None, None]]
        >>> c1 = [[0, None, 1, None, 1, None, None, None, None, None],
                  [None, None, None, 0, 2, 0, None, None, None, None],
                  [None, None, None, None, None, 1, None, None, None, None]]
        >>> c1.scan(0, b*c*d**-1)
        >>> c1
        [[0, None, 1, None, 1, None, 2, None, None, None],
         [None, None, None, 0, 2, 0, None, None, None, None],
         [None, None, None, None, None, 1, None, 0, None, None]]
        >>> c1.scan(0, d*e*a**-1)
        >>> c1
        [[0, None, 1, None, 1, None, 2, None, None, 2],
         [None, None, None, 0, 2, 0, None, None, None, None],
         [None, None, None, None, None, 1, None, 0, 0, None]]
        >>> c1.scan(2, e*a*b**-1)
        >>> c1
        [[0, None, 1, 2, 1, None, 2, None, None, 2],
         [None, None, None, 0, 2, 0, None, None, None, None],
         [None, None, 0, None, None, 1, None, 0, 0, None]]
        >>> c1.scan(2, b*c*d**-1)
        >>> c1
        [[0, None, 1, 2, 1, None, 2, None, None, 2],
         [None, None, None, 0, 2, 0, None, 2, None, None],
         [None, None, 0, None, None, 1, 1, 0, 0, None]]
        >>> c1.scan(1, c*d*e**-1)
        >>> c1
        [[0, None, 1, 2, 1, None, 2, None, None, 2],
         [None, None, None, 0, 2, 0, None, 2, 1, 1],
         [None, None, 0, None, None, 1, 1, 0, 0, None]]
        >>> c1.scan(2, d*e*a**-1)
        >>> c1
        [[0, None, 1, 2, 1, None, 2, None, None, 2],
         [None, 2, None, 0, 2, 0, None, 2, 1, 1],
         [1, None, 0, None, None, 1, 1, 0, 0, None]]
        >>> c1.scan(0, e*a*b**-1)
        >>> c1
        [[0, None, 1, 2, 1, None, 2, None, 2, 2],
         [None, 2, None, 0, 2, 0, None, 2, 1, 1],
         [1, None, 0, None, None, 1, 1, 0, 0, 0]]
        >>> c1.scan(0, c*d*e**-1)
        # sudden-collapse
        # I don't know what this means!! (programmatically, mathematically clear)
        >>> c1
        [[0, None, 0, 0, 0, 0, 0, 0, 0, 0],
         [None, 2, None, 0, 2, 0, None, 2, 1, None],
         [0, None, 0, None, None, None, None, 0, 0, 0]]

        # Example from original Todd-Coxeter paper
        >>> F = FpGroup(F, [x**5, y**3, (x*y)**2])
        # NOTE: In comments we name cosets as 1, 2, 3, ... (starting with 1)
        # define H with Y = {x}, so H is subgroup generated by {x}
        # 1. Start with "definitions" 1^x = 1 (trivial), 1^y = 2, 2^y = 3.
        # here "^" operation represent the function χ (chi, from Derek Holt)
        >>> C = [[0, 0, 1, None], [None, None, 2, 0], [None, None, None, 1]]
        # scanning y**3 under 1
        >>> C.scan(0, y**3)
        # added 3^y = 1 (by deduction)
        >>> C
        [[0, 0, 1, 2], [None, None, 2, 0], [None, None, 0, 1]]
        # scanning (x*y)^2 under 1
        >>> C.scan(0, (x*y)**2)
        # added 2^x = 3 (by deduction)
        >>> C
        [[0, 0, 1, 2], [2, None, 2, 0], [None, 1, 0, 1]]
        # define 3^x = 4, 4^x = 5, 5^x = 6
        >>> C = [[0, 0, 1, 2], [2, None, 2, 0], [3, 1, 0, 1], [4, 2, None, None], [5, 3, None, None], [None, 4, None, None]]
        # scanning x**5 under 2
        >>> C.scan(1, x**5)
        # added 2^y = 3, 6^x = 2 (by deductions)
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, None, None],
         [5, 3, None, None],
         [1, 4, None, None]]

        # scanning (x*y)**3 under 3
        >>> C.scan(2, (x*y)**2)
        # added 4^y = 6 (by deductions)
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, None],
         [5, 3, None, None],
         [1, 4, None, 3]]

        # scanning y**5 under 4
        >>> C.scan(3, y**3)
        # added 7^y = 4 (by deductions)
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, None, None],
         [1, 4, 6, 3],
         [None, None, 3, 5]]

        >>> C = [[0, 0, 1, 2], [2, 5, 2, 0],
                 [3, 1, 0, 1], [4, 2, 5, 6],
                 [5, 3, 7, None], [1, 4, 6, 3],
                 [None, None, 3, 5], [None, None, None, 4]]
        >>> C.scan(3, (x*y)**2)
        # added 5^y = 8, 8^x = 7
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, None],
         [1, 4, 6, 3],
         [None, 7, 3, 5],
         [6, None, None, 4]]

        # define 8^y = 9
        >>> C = [[0, 0, 1, 2], [2, 5, 2, 0],
                 [3, 1, 0, 1], [4, 2, 5, 6],
                 [5, 3, 7, None], [1, 4, 6, 3],
                 [None, 7, 3, 5], [6, None, 8, 4],
                 [None, None, None, 7]]
        >>> C.scan(4, y**3)
        # added 9^y = 5
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, 8],
         [1, 4, 6, 3],
         [None, 7, 3, 5],
         [6, None, 8, 4],
         [None, None, 4, 7]]
        >>> C.scan(4, (x*y)**2)
        # added 7^x = 9
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, 8],
         [1, 4, 6, 3],
         [8, 7, 3, 5],
         [6, None, 8, 4],
         [None, 6, 4, 7]]

        >>> C = [[0, 0, 1, 2], [2, 5, 2, 0],
                 [3, 1, 0, 1], [4, 2, 5, 6],
                 [5, 3, 7, 8], [1, 4, 6, 3],
                 [8, 7, 3, 5], [6, None, 8, 4],
                 [9, 6, 4, 7], [10, 8, None, None],
                 [None, 9, None, None]]
        >>> C.scan(6, x**5)
        # added 11^x = 8
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, 8],
         [1, 4, 6, 3],
         [8, 7, 3, 5],
         [6, 10, 8, 4],
         [9, 6, 4, 7],
         [10, 8, None, None],
         [7, 9, None, None]]

        >>> C.scan(8, (x*y)**2)
        # added 10^y = 11
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, 8],
         [1, 4, 6, 3],
         [8, 7, 3, 5],
         [6, 10, 8, 4],
         [9, 6, 4, 7],
         [10, 8, 10, None],
         [7, 9, None, 9]]
        # define 11^y = 12
        >>> C = [[0, 0, 1, 2], [2, 5, 2, 0],
                 [2, 5, 2, 0], [4, 2, 5, 6],
                 [5, 3, 7, 8], [1, 4, 6, 3],
                 [8, 7, 3, 5], [6, 10, 8, 4],
                 [9, 6, 4, 7], [10, 8, 10, None],
                 [7, 9, 11, 9], [None, None, None, 10]]
        # added 12^y = 10
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, 8],
         [1, 4, 6, 3],
         [8, 7, 3, 5],
         [6, 10, 8, 4],
         [9, 6, 4, 7],
         [10, 8, 10, 11],
         [7, 9, 11, 9],
         [None, None, 9, 10]]

        >>> C.scan(9, (x*y)**2)
        # define 12^x = 12 (by deductions)
        >>> C
        [[0, 0, 1, 2],
         [2, 5, 2, 0],
         [3, 1, 0, 1],
         [4, 2, 5, 6],
         [5, 3, 7, 8],
         [1, 4, 6, 3],
         [8, 7, 3, 5],
         [6, 10, 8, 4],
         [9, 6, 4, 7],
         [10, 8, 10, 11],
         [7, 9, 11, 9],
         [11, 11, 9, 10]]

        """
        # alpha is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for alpha=0 and w in Y (i.e generating set of H)
        # 2. alpha in omega (set of live cosets), w in R (relators)
        f = alpha
        i = 0
        r = len(word)
        # list of union of generators and their inverses
        A_index = self.A_index
        while i < r and self.table[f][A_index[word.subword(i, i+1)]] is not None:
            f = self.table[f][A_index[word.subword(i, i+1)]]
            i += 1
        # can this be replaced with i == r ?
        if i >= r:
            if f != alpha:
                # implement the "coincidence" routine on Pg 158 of Handbook.
                self.coincidence(f, alpha)
                return
        b = alpha
        j = r - 1
        while j >= i and self.table[b][A_index[word.subword(j, j+1)**-1]] is not None:
            b = self.table[b][A_index[word.subword(j, j+1).inverse()]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # run the "coincidence" routine
            # line: We assume that before each call of COINCIDENCE, for α ∈ [1..n],
            # we have p[α]=α iff α ∈ Ω.
            self.p.extend(list(range(len(self))))
            self.coincidence(f, b)
        elif j == i:
            # deduction process
            self.table[f][A_index[word.subword(i, i+1)]] = b
            self.table[b][A.index[word.subword(i, i+1).inverse()]] = f
        # otherwise scan is incomplete and yields no information

    def merge(self, k, lamda, q):
        p = self.p
        phi = self.rep(k)
        chi = self.rep(lamda)
        if phi != chi:
            mu = min(phi, chi)
            v = max(phi, chi)
            p[v] = mu
            q.append(v)

    def rep(self, k):
        p = self.p
        lamda = k
        rho = p[lamda]
        while rho != lamda:
            lamda = rho
            rho = p[lamda]
        mu = k
        rho = p[mu]
        while rho != lamda:
            p[mu] = lamda
            mu = rho
            rho = p[mu]
        return lamda

    # alpha, beta coincide
    # here alpha, beta represent the pair of cosets where
    # coincidence occurs
    def coincidence(self, alpha, beta):
        """
        >>> F, x, y = free_group("x, y")
        # when coincidence occurs with 5 ~ 0, then found 4 ~ 3.
        >>> p = [0, 1, 2, 3, 4, 5]
        >>> C = [[1, None, None, 3], [2, 0, 2, 5],
                 [4, 1, 3, 1],[None, None, 0, 2],
                 [5, 2, 5, None], [None, 4, 1, 4]]
        >>> coincidence(C, 0, 5, [x, x**-1, y, y**-1], p)
        >>> C
        [[1, 3, 1, 3],
         [2, 0, 2, 0],
         [3, 1, 3, 1],
         [0, 2, 0, 2],
         [0, 2, None, None],
         [None, 4, 1, 4]]

        """
        A_dict = self.A_dict
        p = self.p
        l = 0
        # behaves as a queue
        q = []
        self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            # comment by Kalevi, this is already done by p[v] = mu
            # del C[gamma]
            # commenting this out
            # omega = omega - set(gamma)
            for x in A_dict:
                delta = self.table[gamma][A_dict[x]]
                if delta is not None:
                    self.table[delta][A_dict[x**-1]] = None
                    mu = self.rep(gamma)
                    nu = self.rep(delta)
                    if self.table[mu][A_dict[x]] is not None:
                        self.merge(nu, self.table[mu][A_dict[x]], q)
                    elif self.table[nu][A_dict[x**-1]] is not None:
                        self.merge(mu, self.table[nu][A_dict[x**-1]], q)
                    else:
                        self.table[mu][A_dict[x]] = nu
                        self.table[nu][A_dict[x**-1]] = mu

    # method used in the HLT strategy
    def scan_and_fill(self, alpha, word):
        f = alpha
        i = 0
        r = len(word)
        A_dict = self.A_dict
        l_A = len(A_dict)
        i_A = 0
        while i_A < l_A:
            # do the forward scanning
            while i < r and self.table[f][A_dict[word.subword(i, i+1)]] is not None:
                f = self.table[f][A_dict[word.subword(i, i+1)]]
                i += 1
            if i >= r:
                if f != alpha:
                    self.coincidence(f, alpha)
                    return
            # forward scan was incomplete, scan backwards
            b = alpha
            j = r - 1
            while j >= i and self.table[b][A_dict[word.subword(j, j+1)**-1]] is not None:
                b = self.table[b][A_dict[word.subword(j, j+1)**-1]]
                j -= 1
            if j < i:
                self.coincidence(f, b)
            elif j == i:
                self.table[f][A_dict[word.subword(i, i+1)]] = b
                self.table[b][A_dict[word.subword(i, i+1)**-1]] = f
            else:
                self.define(f, word.subword(i, i+1))
            # loop until it has filled the alpha row in the table.
            i_A += 1


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
                    C.scan(alpha, w)
                    if p[alpha] < alpha:
                        break
        beta = C.table[alpha][x]
        if p[beta] == beta:
            for w in R_c_x_inv:
                C.scan(beta, w)
                if p[beta] < beta:
                    break


def standardize(C):
    gamma = 2
    n = C.n
    for a in range(n):
        x in A
        beta = C.table[alpha][A_dict[x]]*gamma
        if beta >= gamma:
            if beta > gamma:
                C.switch(gamma, beta)*gamma
                gamma += 1
                if gamma == n:
                    return


# relator-based method
def coset_enumeration_r(fp_grp, Y):
    """
    # Example 5.1
    >>> from sympy.combinatorics.free_group import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup
    >>> F, x, y = free_group("x, y")

    # Example 5.1
    >>> f = FpGroup(F, [x**3, y**3, x**-1*y**-1*x*y])
    >>> C = coset_enumeration_r(f, [x])
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [0, 0, 1, 2]
    [1, 1, 2, 0]
    [2, 2, 0, 1]
    >>> C.p
    [0, 1, 2, 1, 2]

    # Example 5.2
    >>> f = FpGroup(F, [x**2, y**3, (x*y)**3])
    >>> Y = [x*y]
    >>> C = coset_enumeration_r(f, Y)
    >>> C.table
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [1, 1, 2, 1]
    [0, 0, 0, 2]
    [3, 3, 1, 0]
    [2, 2, 3, 3]

    # Example 5.3
    >>> f = FpGroup(F, [x**2*y**2, x**3*y**5])
    >>> Y = [x*y]
    >>> C = coset_enumeration_r(f, Y)
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [1, 3, 1, 3]
    [2, 0, 2, 0]
    [3, 1, 3, 1]
    [0, 2, 0, 2]

    # Example 5.4
    >>> F = free_group("a, b, c, d, e")
    >>> f = FpGroup(F, [a*b*c**-1, b*c*d**-1, c*d*e**-1, d*e*a**-1, e*a*b**-1])
    >>> Y = [a]
    >>> C = coset_enumeration_r(f, Y)
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    # Exercises Pg. 161, Q2.
    >>> f = FpGroup(F, [x**2*y**2, y**-1*x*y*x**-3])
    >>> Y = []
    >>> C = coset_enumeration_r(f, Y)
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [1, 5, 7, 3]
    [2, 0, 4, 8]
    [5, 1, 3, 7]
    [4, 8, 0, 2]
    [7, 3, 5, 1]
    [0, 2, 8, 4]
    [8, 4, 2, 0]
    [3, 7, 1, 5]

    # John J. Cannon; Lucien A. Dimino; George Havas; Jane M. Watson
    # Mathematics of Computation, Vol. 27, No. 123. (Jul., 1973), pp. 463-490
    # from 1973chwd.pdf
    # Table 1. Ex. 1
    >>> F, r, s, t = free_group("r, s, t")
    >>> E1 = FpGroup(F, [t**-1*r*t*r**-2, r**-1*s*r*s**-2, s**-1*t*s*t**-2])
    >>> C = coset_enumeration_r(E1, [x])
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [0, 0, 0, 0, 0, 0]

    Ex. 2
    >>> F, a, b = free_group("a, b")
    >>> Cox = FpGroup(F, [a**6, b**6, (a*b)**2, (a**2*b**2)**2, (a**3*b**3)**5])
    >>> C = coset_enumeration_r(Cox, [a])
    >>> index = 0
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    ...         index += 1
    >>> index
    500

    # Ex. 3
    >>> F, a, b = free_group("a, b")
    >>> B_2_4 = FpGroup(F, [a**4, b**4, (a*b)**4, (a**-1*b)**4, (a**2*b)**4, (a*b**2)**4, (a**2*b**2)**4, (a**-1*b*a*b)**4, (a*b**-1*a*b)**4])
    >>> C = coset_enumeration_r(B_2_4, [a])
    >>> index = 0
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    ...         index += 1
    >>> index
    1024

    """
    # 1. Initialize a coset table C for < X|R >
    C = CosetTable(fp_grp, Y)
    R = fp_grp.relators()
    A_dict = C.A_dict
    p = C.p
    for w in Y:
        C.scan_and_fill(0, w)
    alpha = 0
    while alpha < C.n:
        if p[alpha] == alpha:
            for w in R:
                C.scan_and_fill(alpha, w)
                if p[alpha] < alpha:
                    break
            if p[alpha] < alpha:
                continue
            for x in A_dict:
                if C.table[alpha][A_dict[x]] is None:
                    C.define(alpha, x)
        alpha += 1
    return C


# Pg. 166
# coset-table bases method
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
            if C[alpha][x] is None:
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
                    beta = C[alpha][A_dict[x]]
                    if beta == alpha:
                        beta = gamma
                    C[gamma][A.index(x)] = beta
                    C[beta][A.index(x**-1)] == gamma
    C.n = gamma
    for alpha in range(C.n):
        p[alpha] = alpha


FpGroupElm = FreeGroupElm

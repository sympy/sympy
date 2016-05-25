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

def define(coset_table, alpha, x):
    if coset_table.n == DefaultMaxLimit:
        # abort the further generation of cosets
        return
    coset_table.n += 1
    # beta is the new coset generated
    beta = n
    coset_table.p[beta] = beta
    coset_table[alpha][coset_table.A.index(x)] = beta
    coset_table[alpha][coset_table.A.index(x.inverse())] = alpha


def scan(coset_table, alpha, word, A):
    """
    >>> F, x, y = free_group("x, y")
    >>> A = [x, x**-1, y, y**-1]

    # Example 5.1, Pg 150
    >>> c = [[0, 0, 1, 2], [-1, -1, -1, 0], [-1, -1, 0, -1]]
    >>> scan(c, 1, y**3, A)
    # this lead to deductions 2^y = 3, 3^(y^-1) = 2 (1, 2, 3 represent cosets in this comment)
    >>> c
    [[0, 0, 1, 2], [-1, -1, 2, 0], [-1, -1, 0, 1]]

    # Example 5.2, Pg 154
    >>> c1 = [[1, 1, -1, -1], [0, 0, -1, -1]]
    >>> scan(c1, 0, x*y, A)
    >>> c1
    [[1, 1, -1, 1], [0, 0, 0, -1]]

    >>> c2 = [[1, 1, 2, 1], [0, 0, 0, -1]] + [[-1, -1, -1, 0]]
    >>> scan(c2, 1, y**3, A)
    >>> c2
    [[1, 1, 2, 1], [0, 0, 0, 2], [-1, -1, 1, 0]]

    >>> c3 = [[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0]] + [[2, 2, -1, -1]]
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
    while i < r and coset_table[f][A.index(word.subword(i, i+1))] != -1:
        f = coset_table[f][A.index(word.subword(i, i+1))]
        i += 1
    if i >= r:
        if f != alpha:
            # implement the "coincidence" routine on Pg 158 of Handbook.
            concidence(coset_table, f, alpha)
            return
    b = alpha
    j = r
    while j >= i and coset_table[b][A.index(word.subword(j-1, j).inverse())] != -1:
        b = coset_table[b][A.index(word.subword(j-1, j).inverse())]
        j -= 1
    if j < i:
        # run the "coincidence" routine
        concidence(C, f, b)
    elif j == i + 1:
        # deduction process
        coset_table[f][A.index(word.subword(i, i+1))] = b
        coset_table[b][A.index(word.subword(i, i+1).inverse())] = f


FpGroupElm = FreeGroupElm

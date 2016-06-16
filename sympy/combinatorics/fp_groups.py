# -*- coding: utf-8 -*-
from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public
from sympy.utilities.iterables import flatten
from sympy.combinatorics.free_group import FreeGroupElement
from itertools import chain
from bisect import bisect_left


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


###############################################################################
#                           FINITELY PRESENTED GROUPS                         #
###############################################################################


class FpGroup(DefaultPrinting):
    """
    The FpGroup would take a FreeGroup and a list/tuple of relators, the
    relators would be specified in such a way that each of them be equal to the
    identity of the provided free group.
    """
    is_group = True
    is_FpGroup = True
    is_PermutationGroup = False

    def __new__(cls, fr_grp, relators):
        relators = _parse_relators(relators)
        # return the corresponding FreeGroup if no relators are specified
        if not relators:
            return fr_grp
        obj = object.__new__(cls)
        obj._free_group = fr_grp
        obj._relators = relators
        obj.generators = obj._generators()
        obj.dytype = type("FpGroupElement", (FpGroupElement,), {"group": obj})
        return obj

    @property
    def free_group(self):
        return self._free_group

    def relators(self):
        return tuple(self._relators)

    def _generators(self):
        """Returns the generators of the associated free group."""
        return self.free_group.generators

    def __str__(self):
        if self.free_group.rank > 30:
            str_form = "<fp group with %s generators>" % self.free_group.rank
        else:
            str_form = "<fp group on the generators %s>" % str(self.generators)
        return str_form

    __repr__ = __str__


# sets the upper limit on the number of cosets generated during
# Coset Enumeration. "M" from Derek Holt's. It is supposed to be
# user definable.
CosetTableDefaultMaxLimit = 4096000
max_stack_size = 500

class CosetTable(DefaultPrinting):
    # coset_table: Mathematically a coset table
    #               represented using a list of lists
    # alpha: Mathematically a coset (precisely, a live coset)
    #       represented by an integer between i with 1 <= i <= n
    #       α ∈ c
    # x: Mathematically an element of "A" (set of generators and
    #   their inverses), represented using "FpGroupElement"
    # fp_grp: Finitely Presented Group with < X|R > as presentation.
    # H: subgroup of fp_grp.
    # NOTE: We start with H as being only a list of words in generators
    #       of "fp_grp". Since `.subgroup` method has not been implemented.

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
        # "p" is setup independent of Ω and n
        self.p = [0]
        self.A = list(chain.from_iterable((gen, gen**-1) \
                for gen in self.fp_group.generators))
        self.table = [[None]*len(self.A)]
        self.A_dict = {x: self.A.index(x) for x in self.A}
        self.A_dict_inv = {}
        for x, index in self.A_dict.items():
            if index % 2 == 0:
                self.A_dict_inv[x] = self.A_dict[x] + 1
            else:
                self.A_dict_inv[x] = self.A_dict[x] - 1
        self.deduction_stack = []

    @property
    def omega(self):
        """Set of live cosets
        """
        return [coset for coset in range(len(self.p)) if self.p[coset] == coset]

    def copy(self):
        self_copy = self.__class__(self.fp_group, self.subgroup)
        self_copy.table = [perm_rep.copy() for perm_rep in self.table]
        self_copy.p = self.p.copy()
        self_copy.deduction_stack = self.deduction_stack.copy()
        return self_copy

    def __str__(self):
        return "Coset Table on %s with %s as subgroup generators" \
                % (self.fp_group, self.subgroup)

    __repr__ = __str__

    @property
    def n(self):
        """The number 'n' represents the length of the sublist containing the
        live cosets.
        """
        return max(self.omega) + 1

    # Pg 152 [1]
    def is_complete(self):
        """
        The coset table is called complete if it has no undefined entries
        on the live cosets; that is, whether α^x is defined for all
        α ∈ Ω and x ∈ A.
        """
        return not any([None in self.table[coset] for coset in self.omega])

    # Pg. 153
    def define(self, alpha, x):
        A = self.A
        if len(self.table) == CosetTableDefaultMaxLimit:
            # abort the further generation of cosets
            return
        self.table.append([None]*len(A))
        # beta is the new coset generated
        beta = len(self.table) - 1
        self.p.append(beta)
        self.table[alpha][self.A_dict[x]] = beta
        self.table[beta][self.A_dict_inv[x]] = alpha

    def define_f(self, alpha, x):
        A = self.A
        if self.table == CosetTableDefaultMaxLimit:
            # abort the further generation of cosets
            return
        self.table.append([None]*len(A))
        # beta is the new coset generated
        beta = len(self.table) - 1
        self.p.append(beta)
        self.table[alpha][self.A_dict[x]] = beta
        self.table[beta][self.A_dict_inv[x]] = alpha
        # append to deduction stack
        self.deduction_stack.append((alpha, x))

    def scan_f(self, alpha, word):
        # alpha is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for alpha=0 and w in Y (i.e generating set of H)
        # 2. alpha in omega (set of live cosets), w in R (relators)
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        # list of union of generators and their inverses
        while i <= j and self.table[f][A_dict[word[i]]] is not None:
            f = self.table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            if f != b:
                self.coincidence_f(f, b)
            return
        while j >= i and self.table[b][A_dict_inv[word[j]]] is not None:
            b = self.table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # run the "coincidence" routine
            self.coincidence_f(f, b)
        elif j == i:
            # deduction process
            self.table[f][A_dict[word[i]]] = b
            self.table[b][A_dict_inv[word[i]]] = f
            self.deduction_stack.append((f, word[i]))
        # otherwise scan is incomplete and yields no information

    # α, β coincide, i.e. α, β represent the pair of cosets where
    # coincidence occurs
    def coincidence_f(self, alpha, beta):
        """
        """
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        p = self.p
        l = 0
        # behaves as a queue
        q = []
        self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            for x in A_dict:
                delta = self.table[gamma][A_dict[x]]
                if delta is not None:
                    self.table[delta][A_dict_inv[x]] = None
                    self.deduction_stack.append((delta, x**-1))
                    mu = self.rep(gamma)
                    nu = self.rep(delta)
                    if self.table[mu][A_dict[x]] is not None:
                        self.merge(nu, self.table[mu][A_dict[x]], q)
                    elif self.table[nu][A_dict_inv[x]] is not None:
                        self.merge(mu, self.table[nu][A_dict_inv[x]], q)
                    else:
                        self.table[mu][A_dict[x]] = nu
                        self.table[nu][A_dict_inv[x]] = mu

    def scan(self, alpha, word):
        # alpha is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for alpha=0 and w in Y (i.e generating set of H)
        # 2. alpha in omega (set of live cosets), w in R (relators)
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        while i <= j and self.table[f][A_dict[word[i]]] is not None:
            f = self.table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            if f != b:
                self.coincidence(f, b)
            return
        while j >= i and self.table[b][A_dict_inv[word[j]]] is not None:
            b = self.table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # run the "coincidence" routine
            self.coincidence(f, b)
        elif j == i:
            # deduction process
            self.table[f][A_dict[word[i]]] = b
            self.table[b][A_dict_inv[word[i]]] = f
        # otherwise scan is incomplete and yields no information

    def scan_check(self, alpha, word):
        """
        Another version of "scan" routine, it checks whether α scans correctly
        under w, it is a straightforward modification of "scan". "scan_check"
        return false (rather than calling "coincidence") if the scan completes
        incorrectly; otherwise it returns true.

        Examples
        ========

        >>> from sympy.combinatorics.fp_groups import FpGroup, scan_check
        >>> from sympy.combinatorics.free_group import free_group
        >>> F, x, y = free_group("x, y")
        >>> f = FpGroup(F, [x**2, y**3, (x*y)**4])
        >>> C = CosetTable(f, [x])

        # scans correctly for (0, x)
        >>> C.scan_check(0, x)
        True
        >>> C.scan(0, x)

        # scans correctly under (0, x**2)
        >>> C.scan(0, x); C.scan_check(0, x**2)
        True
        >>> C.scan(0, x**2); C.scan_check(0, y**3)
        True

        """
        # alpha is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for alpha=0 and w in Y (i.e generating set of H)
        # 2. alpha in omega (set of live cosets), w in R (relators)
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        while i <= j and self.table[f][A_dict[word[i]]] is not None:
            f = self.table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            return f == b
        while j >= i and self.table[b][A_dict_inv[word[j]]] is not None:
            b = self.table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # return False, instead of calling coincidence routine
            return False
        else:
            # return True otherwise
            return True

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

    # α, β coincide, i.e. α, β represent the pair of cosets where
    # coincidence occurs
    def coincidence(self, alpha, beta):
        """
        """
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        p = self.p
        l = 0
        # behaves as a queue
        q = []
        self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            for x in A_dict:
                delta = self.table[gamma][A_dict[x]]
                if delta is not None:
                    self.table[delta][A_dict_inv[x]] = None
                    mu = self.rep(gamma)
                    nu = self.rep(delta)
                    if self.table[mu][A_dict[x]] is not None:
                        self.merge(nu, self.table[mu][A_dict[x]], q)
                    elif self.table[nu][A_dict_inv[x]] is not None:
                        self.merge(mu, self.table[nu][A_dict_inv[x]], q)
                    else:
                        self.table[mu][A_dict[x]] = nu
                        self.table[nu][A_dict_inv[x]] = mu

    # method used in the HLT strategy
    def scan_and_fill(self, alpha, word):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        r = len(word)
        f = alpha
        i = 0
        b = alpha
        j = r - 1
        # loop until it has filled the α row in the table.
        while True:
            # do the forward scanning
            while i <= j and self.table[f][A_dict[word[i]]] is not None:
                f = self.table[f][A_dict[word[i]]]
                i += 1
            if i > j:
                if f != b:
                    self.coincidence(f, b)
                return
            # forward scan was incomplete, scan backwards
            while j >= i and self.table[b][A_dict_inv[word[j]]] is not None:
                b = self.table[b][A_dict_inv[word[j]]]
                j -= 1
            if j < i:
                self.coincidence(f, b)
            elif j == i:
                self.table[f][A_dict[word[i]]] = b
                self.table[b][A_dict_inv[word[i]]] = f
            else:
                self.define(f, word[i])

    def scan_and_fill_f(self, alpha, word):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        r = len(word)
        f = alpha
        i = 0
        b = alpha
        j = r - 1
        # loop until it has filled the α row in the table.
        while True:
            # do the forward scanning
            while i <= j and self.table[f][A_dict[word[i]]] is not None:
                f = self.table[f][A_dict[word[i]]]
                i += 1
            if i > j:
                if f != b:
                    self.coincidence_f(f, b)
                return
            # forward scan was incomplete, scan backwards
            while j >= i and self.table[b][A_dict_inv[word[j]]] is not None:
                b = self.table[b][A_dict_inv[word[j]]]
                j -= 1
            if j < i:
                self.coincidence_f(f, b)
            elif j == i:
                self.table[f][A_dict[word[i]]] = b
                self.table[b][A_dict_inv[word[i]]] = f
                self.deduction_stack.append((f, word[i]))
            else:
                self.define_f(f, word[i])

    # currently not used anywhere
    def look_ahead(self):
        R = self.fp_group.relators()
        p = self.p
        for beta in self.omega:
            # complete scan all relators under all cosets(obviously live)
            # without making new definitions
            for w in R:
                self.scan(beta, w)
                if p[beta] < beta:
                    break

    # Pg. 166
    def process_deductions(self, R_c_x, R_c_x_inv):
        p = self.p
        while len(self.deduction_stack) > 0:
            if len(self.deduction_stack) >= max_stack_size:
                self.lookahead()
                del self.deduction_stack[:]
            else:
                alpha, x = self.deduction_stack.pop()
                if p[alpha] == alpha:
                    for w in R_c_x:
                        self.scan_f(alpha, w)
                        if p[alpha] < alpha:
                            break
            beta = self.table[alpha][self.A_dict[x]]
            if beta is not None and p[beta] == beta:
                for w in R_c_x_inv:
                    self.scan_f(beta, w)
                    if p[beta] < beta:
                        break

    # Pg. 166
    def process_deductions_check(self, R_c_x, R_c_x_inv):
        """
        A variation of "process_deductions", this calls "scan_check" wherever
        "process_deductions" calls "scan".
        """
        p = self.p
        while len(self.deduction_stack) > 0:
            if len(self.deduction_stack) >= max_stack_size:
                self.lookahead()
                del self.deduction_stack[:]
            else:
                alpha, x = self.deduction_stack.pop()
                if p[alpha] == alpha:
                    for w in R_c_x:
                        self.scan_check(alpha, w)
                        if p[alpha] < alpha:
                            break
            beta = self.table[alpha][self.A_dict[x]]
            if beta is not None and p[beta] == beta:
                for w in R_c_x_inv:
                    self.scan_check(beta, w)
                    if p[beta] < beta:
                        break



    def switch(self, beta, gamma):
        """
        Switch the elements β, γ of Ω in C
        """
        A = self.A
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        X = self.fp_group.generators
        table = self.table
        for x in A:
            z = table[gamma][A_dict[x]]
            table[gamma][A_dict[x]] = table[beta][A_dict[x]]
            table[beta][A_dict[x]] = z
            for alpha in range(len(self.p)):
                if self.p[alpha] == alpha:
                    if table[alpha][A_dict[x]] == beta:
                        table[alpha][A_dict[x]] = gamma
                    elif table[alpha][A_dict[x]] == gamma:
                        table[alpha][A_dict[x]] = beta

    def standardize(self):
        """A coset table is standardized if when running through the cosets
        and within each coset through the generator images (ignoring generator
        inverses), the cosets appear in order of the integers 0, 1, 2, ... n.
        "Standardize" reorders the elements of \Omega such that, if we scan
        the coset table first by elements of \Omega and then by elements of A,
        then the cosets occur in ascending order. `standardize()` is used at
        the end of an enumeration to permute the cosets so that they occur in
        some sort of standard order.

        >>> from sympy.combinatorics.free_group import free_group
        >>> from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_r
        >>> F, x, y = free_group("x, y")

        # Example 5.3 from [1]
        >>> f = FpGroup(F, [x**2*y**2, x**3*y**5])
        >>> C = coset_enumeration_r(f, [])
        >>> C.compress()
        >>> C.table
        [[1, 3, 1, 3], [2, 0, 2, 0], [3, 1, 3, 1], [0, 2, 0, 2]]
        >>> C.standardize()
        >>> C.table
        [[1, 2, 1, 2], [3, 0, 3, 0], [0, 3, 0, 3], [2, 1, 2, 1]]

        """
        A = self.A
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        gamma = 1
        for alpha in range(self.n):
            for x in A:
                beta = self.table[alpha][A_dict[x]]
                if beta >= gamma:
                    if beta > gamma:
                        self.switch(gamma, beta)
                    gamma += 1
                    if gamma == self.n:
                        return

    # Compression of a Coset Table
    # Pg. 167 5.2.3
    def compress(self):
        """Removes the non-live cosets from the coset table
        """
        gamma = -1
        A = self.A
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        chi = tuple([i for i in range(len(self.p)) if self.p[i] != i])
        for alpha in self.omega:
            gamma += 1
            if gamma != alpha:
                # replace α by γ in coset table
                for x in A:
                    beta = self.table[alpha][A_dict[x]]
                    self.table[gamma][A_dict[x]] = beta
                    self.table[beta][A_dict_inv[x]] == gamma
        # all the cosets in the table are live cosets
        self.p = list(range(gamma + 1))
        # delete the useless coloumns
        del self.table[len(self.p):]
        # re-define values
        for row in self.table:
            for j in range(len(self.A)):
                row[j] -= bisect_left(chi, row[j])


# relator-based method
def coset_enumeration_r(fp_grp, Y, C):
    """
    >>> from sympy.combinatorics.free_group import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_r
    >>> F, x, y = free_group("x, y")

    # Example 5.1 from [1]
    >>> f = FpGroup(F, [x**3, y**3, x**-1*y**-1*x*y])
    >>> C = coset_enumeration_r(f, [x])
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [0, 0, 1, 2]
    [1, 1, 2, 0]
    [2, 2, 0, 1]
    >>> C.p
    [0, 1, 2, 1, 1]

    # Example from exercises Q2 [1]
    >>> f = FpGroup(F, [x**2*y**2, y**-1*x*y*x**-3])
    >>> C = coset_enumeration_r(f, [])
    >>> C.compress(); C.standardize()
    >>> C.table
    [[1, 2, 3, 4],
    [5, 0, 6, 7],
    [0, 5, 7, 6],
    [7, 6, 5, 0],
    [6, 7, 0, 5],
    [2, 1, 4, 3],
    [3, 4, 2, 1],
    [4, 3, 1, 2]]

    # Example 5.2
    >>> f = FpGroup(F, [x**2, y**3, (x*y)**3])
    >>> Y = [x*y]
    >>> C = coset_enumeration_r(f, Y)
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [1, 1, 2, 1]
    [0, 0, 0, 2]
    [3, 3, 1, 0]
    [2, 2, 3, 3]

    # Example 5.3
    >>> f = FpGroup(F, [x**2*y**2, x**3*y**5])
    >>> Y = []
    >>> C = coset_enumeration_r(f, Y)
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [1, 3, 1, 3]
    [2, 0, 2, 0]
    [3, 1, 3, 1]
    [0, 2, 0, 2]

    # Example 5.4
    >>> F, a, b, c, d, e = free_group("a, b, c, d, e")
    >>> f = FpGroup(F, [a*b*c**-1, b*c*d**-1, c*d*e**-1, d*e*a**-1, e*a*b**-1])
    >>> Y = [a]
    >>> C = coset_enumeration_r(f, Y)
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         print(C.table[i])
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    # example of "compress" method
    >>> C.compress()
    >>> C.table
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

    # Exercises Pg. 161, Q2.
    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**2*y**2, y**-1*x*y*x**-3])
    >>> Y = []
    >>> C = coset_enumeration_r(f, Y)
    >>> C.compress()
    >>> C.standardize()
    >>> C.table
    [[1, 2, 3, 4],
    [5, 0, 6, 7],
    [0, 5, 7, 6],
    [7, 6, 5, 0],
    [6, 7, 0, 5],
    [2, 1, 4, 3],
    [3, 4, 2, 1],
    [4, 3, 1, 2]]

    # John J. Cannon; Lucien A. Dimino; George Havas; Jane M. Watson
    # Mathematics of Computation, Vol. 27, No. 123. (Jul., 1973), pp. 463-490
    # from 1973chwd.pdf
    # Table 1. Ex. 1
    >>> F, r, s, t = free_group("r, s, t")
    >>> E1 = FpGroup(F, [t**-1*r*t*r**-2, r**-1*s*r*s**-2, s**-1*t*s*t**-2])
    >>> C = coset_enumeration_r(E1, [r])
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
    ...         index += 1
    >>> index
    500

    # Ex. 3
    >>> F, a, b = free_group("a, b")
    >>> B_2_4 = FpGroup(F, [a**4, b**4, (a*b)**4, (a**-1*b)**4, (a**2*b)**4, \
            (a*b**2)**4, (a**2*b**2)**4, (a**-1*b*a*b)**4, (a*b**-1*a*b)**4])
    >>> C = coset_enumeration_r(B_2_4, [a])
    >>> index = 0
    >>> for i in range(len(C.p)):
    ...     if C.p[i] == i:
    ...         index += 1
    >>> index
    1024

    """
    # 1. Initialize a coset table C for < X|R >
    if not C:
        C = CosetTable(fp_grp, Y)
    R = fp_grp.relators()
    A_dict = C.A_dict
    A_dict_inv = C.A_dict_inv
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
            if p[alpha] >= alpha:
                for x in A_dict:
                    if C.table[alpha][A_dict[x]] is None:
                        C.define(alpha, x)
        alpha += 1
    return C


# Pg. 166
# coset-table based method
def coset_enumeration_c(fp_grp, Y):
    """
    >>> from sympy.combinatorics.free_group import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_c
    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**3, y**3, x**-1*y**-1*x*y])
    >>> C = coset_enumeration_c(f, [x])
    >>> C.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [2, 2, 0, 1]]

    """
    # Initialize a coset table C for < X|R >
    C = CosetTable(fp_grp, Y)
    X = fp_grp.generators
    R = fp_grp.relators()
    A = C.A
    # replace all the elements by cyclic reductions
    R_cyc_red = [rel.identity_cyclic_reduction() for rel in R]
    R_c = list(chain.from_iterable((rel.cyclic_conjugates(), (rel**-1).cyclic_conjugates()) \
            for rel in R_cyc_red))
    R_set = set()
    for conjugate in R_c:
        R_set = R_set.union(conjugate)
    # a list of subsets of R_c whose words start with "x".
    R_c_list = []
    for x in C.A:
        r = set([word for word in R_set if word[0] == x])
        R_c_list.append(r)
        R_set.difference_update(r)
    for w in Y:
        C.scan_and_fill_f(0, w)
    for x in A:
        C.process_deductions(R_c_list[C.A_dict[x]], R_c_list[C.A_dict_inv[x]])
    i = 0
    while i < len(C.omega):
        alpha = C.omega[i]
        i += 1
        for x in C.A:
            if C.table[alpha][C.A_dict[x]] is None:
                C.define_f(alpha, x)
                C.process_deductions(R_c_list[C.A_dict[x]], R_c_list[C.A_dict_inv[x]])
    return C


def low_index_subgroups(G, N):
    """
    Implements the Low Index Subgroups algorithm, i.e find all subgroups of
    "G" upto a given index "N". This implements the method described in
    [Sim94]. This procedure involves a backtrack search over incomplete Coset
    Tables, rather than over forced coincidences.

    G: An FpGroup < X|R >
    N: positive integer, representing the maximun index value for subgroups

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of Computational Group Theory"
    Section 5.4

    [2] Marston Conder and Peter Dobcsanyi
    "Applications and Adaptions of the Low Index Subgroups Procedure"

    """
    C = CosetTable(G, [])
    R = G.relators()
    # length chosen for the length of the short relators
    len_short_rel = 5
    # elements of R2 only checked at the last step for complete
    # coset tables
    R2 = set([rel for rel in R if len(rel) > len_short_rel])
    # elements of R1 are used in inner parts of the process to prune
    # branches of the search tree,
    R1 = set([rel.identity_cyclic_reduction() for rel in set(R) - R2])
    R1_c = list(chain.from_iterable((rel.cyclic_conjugates(), \
            (rel**-1).cyclic_conjugates()) for rel in R1))
    R1_set = set()
    for conjugate in R1_c:
        R1_set = R1_set.union(conjugate)
    R1_c_list = []
    for x in C.A:
        r = set([word for word in R1_set if word[0] == x])
        R1_c_list.append(r)
        R1_set.difference_update(r)
    S = []
    for x in C.A:
        descendant_subgroups(S, C, R1_c_list[C.A_dict[x]], \
                R1_c_list[C.A_dict_inv[x]], R2, N)
        return C.table, C.deduction_stack
    return S


def descendant_subgroups(S, C, R1_x_c, R1_x_c_inv, R2, N):
    A_dict = C.A_dict
    A_dict_inv = C.A_dict_inv
    if C.is_complete():
        # check whether the relators in R2 are satisfied
        for w in R2:
            for alpha in C.omega:
                if not C.scan_check(alpha, w):
                    return
        # relators in R2 are satisfied, append the table to list
        S.append(C)
    else:
        j = 0
        for alpha in range(len(C.table)):
            for x in C.A:
                alpha_c_x = C.table[alpha][A_dict[x]]
                if alpha_c_x is None:
                    # this is "x" in pseudo-code (using "y" makes it clear)
                    y = x
                    j = 1
                    break
            if j == 1:
                break
        reach = C.omega + [C.n]
        for beta in reach:
            if beta == C.n:
                C.table.append([None]*len(C.A))
            if beta < N and C.table[beta][A_dict_inv[y]] is None:
                try_descendant(S, C, R1_x_c, R1_x_c_inv, R2, N, alpha, y, beta)


def try_descendant(S, C, R1_x_c, R1_x_c_inv, R2, N, alpha, x, beta):
    """
    It solves the problem of trying out each individual possibility for α^x.
    """
    D = C.copy()
    A_dict = D.A_dict
    if beta == D.n:
        D.table.append([None]*len(D.A))
        D.p.append(beta)
    D.table[alpha][D.A_dict[x]] = beta
    D.table[beta][D.A_dict_inv[x]] = alpha
    D.deduction_stack.append((alpha, x))
    D.process_deductions_check(R1_x_c, R1_x_c_inv)
    if D.n == 0:
        return
    if first_in_class(D):
        descendant_subgroups(S, D, R1_x_c, R1_x_c_inv, R2, N)


def first_in_class(C):
    """
    >>> from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_c, coset_enumeration_r, CosetTable, first_in_class
    >>> from sympy.combinatorics.free_group import FreeGroup, free_group
    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**2, y**3])
    >>> C = CosetTable(f, [])
    >>> C.table
    [[None, None, None, None]]

    # C1
    >>> C.table[0][0] = 0; C.table[0][1] = 0
    >>> C.table.append([None]*4)
    >>> C.table[0][2] = 1; C.table[1][3] = 0 # C12
    >>> C.table.append([None]*4)
    >>> C.table[0][3] = 2; C.table[2][2] = 0 #C121
    >>> C.scan(0, y**3) # makes deductions 2^y = 3; 3^(y^-1) = 2
    >>> D = C.copy() # equals C121
    >>> C.table[1][0] = 1; C.table[1][1] = 1 # C1211
    >>> C.table.append([None]*4)
    >>> C.table[2][0] = 3; C.table[3][1] = 2 # C12112
    >>> C.table[3][2] = 3; C.table[3][3] = 3
    >>> C.table[2][1] = C.table[2][0]
    >>> C.table[3][0] = C.table[3][1] # this definition, since n=4
    >>> C.p = [0, 1, 2, 3]
    >>> first_in_class(C)
    True

    >>> D.table[1][0] = 2; D.table[1][1] = 2 # C1212
    >>> D.table[2][0] = 1; D.table[2][1] = 1

    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**2, y**3])
    >>> C = CosetTable(f, [])
    >>> C1 = C.copy()
    >>> C1.table
    [[None, None, None, None]]
    >>> C1.table[0][0] = 0; C1.table[0][1] = 0
    >>> C1.table
    [[0, 0, None, None]]
    >>> C11 = C1.copy()
    >>> C11.table[0][2] = 0; C11.table[0][3] = 0
    >>> C11.table
    [[0, 0, 0, 0]]
    >>> first_in_class(C11)
    True
    >>> C12 = C1.copy()
    >>> C12.table.append([None]*4); C12.p = [0, 1]
    >>> C12.table
    [[0, 0, None, None], [None, None, None, None]]
    >>> C12.table[0][2] = 1; C12.table[1][3] = 0
    >>> C12.table
    [[0, 0, 1, None], [None, None, None, 0]]
    >>> C122 = C12.copy()
    >>> C122.table[0][3] = 1; C122.table[1][2] = 0
    >>> C122.table
    [[0, 0, 1, 1], [None, None, 0, 0]]
    >>> C122.scan_check(0, y**3)
    False
    >>> C121 = C12.copy()
    >>> C121.table.append([None]*4)
    >>> C121.p = [0, 1, 2]
    >>> C121.table[0][3] = 2; C121.table[2][2] = 0
    >>> C121.table
    [[0, 0, 1, 2], [None, None, None, 0], [None, None, 0, None]]
    >>> C121.scan(0, y**3)
    >>> C121.table
    [[0, 0, 1, 2], [None, None, 2, 0], [None, None, 0, 1]]
    >>> C1211 = C121.copy()
    >>> C1211.table[1][0] = 1; C1211.table[1][1] = 1
    >>> C1211.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [None, None, 0, 1]]
    >>> C12111 = C1211.copy()
    >>> C12111.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [None, None, 0, 1]]
    >>> C12111.table[2][0] = 2; C12111.table[2][1] = 2
    >>> C12111.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [2, 2, 0, 1]]
    >>> first_in_class(C12111)
    True
    >>> C12111.scan_check(0, (x*y)**4)
    False
    >>> C12112 = C1211.copy()
    >>> C12112.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [None, None, 0, 1]]
    >>> C12112.p
    [0, 1, 2]
    >>> C12112.table.append([None]*4)
    >>> C12112.p = [0, 1, 2, 3]
    >>> C12112.table[2][0] = 3; C12112.table[3][1] = 2
    >>> C12112.table[3][0] = 2; C12112.table[2][1] = 3
    >>> C12112.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [3, 3, 0, 1], [2, 2, None, None]]
    >>> C12112.table[3][2] = 3; C12112.table[3][3] = 3
    >>> C12112.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [3, 3, 0, 1], [2, 2, 3, 3]]
    >>> first_in_class(C12112) # giving the wrong-result
    False
    >>> C12112.table
    [[0, 0, 1, 2], [1, 1, 2, 0], [3, 3, 0, 1], [2, 2, 3, 3]]
    >>> C12112.scan_check(0, (x*y)**4)
    True
    >>> C1213 = C121.copy()
    >>> C1213.table
    [[0, 0, 1, 2], [None, None, 2, 0], [None, None, 0, 1]]
    >>> C1213.table.append([None]*4)
    >>> C1213.table
    [[0, 0, 1, 2],
     [None, None, 2, 0],
     [None, None, 0, 1],
     [None, None, None, None]]
    >>> C1213.table[1][0] = 3; C1213.table[3][1] = 1;
    >>> C1213.table[1][1] = 3; C1213.table[3][0] = 1;
    >>> C1213.table
    [[0, 0, 1, 2], [3, 3, 2, 0], [None, None, 0, 1], [1, 1, None, None]]
    >>> C1213.p = [0, 1, 2, 3]
    >>> C12131 = C1213.copy()
    >>> C12131.table
    [[0, 0, 1, 2], [3, 3, 2, 0], [None, None, 0, 1], [1, 1, None, None]]
    >>> C12131.scan_check(2, x**2)
    True
    >>> C12131.table[2][0] = 2; C12131.table[2][1] = 2;
    >>> C12131.table
    [[0, 0, 1, 2], [3, 3, 2, 0], [2, 2, 0, 1], [1, 1, None, None]]
    >>> first_in_class(C12131)
    False
    >>> C2 = C.copy()
    >>> C2.table
    [[None, None, None, None]]
    >>> C2.table.append([None]*4)
    >>> C2.p = [0, 1]
    >>> C2.table[0][0] = 1; C2.table[1][1] = 0
    >>> C2.table[0][1] = 1; C2.table[1][0] = 0
    >>> C2.table
    [[1, 1, None, None], [0, 0, None, None]]
    >>> C21 = C2.copy()
    >>> C21.table[0][2] = 0; C21.table[0][3] = 0
    >>> C21.table
    [[1, 1, 0, 0], [0, 0, None, None]]
    >>> C211 = C21.copy()
    >>> C211.table[1][2] = 1; C211.table[1][3] = 1
    >>> C211.table
    [[1, 1, 0, 0], [0, 0, 1, 1]]
    >>> first_in_class(C211)
    True
    >>> C211.scan_check(0, (x*y)**4)
    True
    >>> C212 = C21.copy()
    >>> C212.table.append([None]*4)
    >>> C212.p = [0, 1, 2]
    >>> C212.table
    [[1, 1, 0, 0], [0, 0, None, None], [None, None, None, None]]
    >>> C212.table[1][2] = 2; C212.table[2][3] = 1
    >>> C212.table
    [[1, 1, 0, 0], [0, 0, 2, None], [None, None, None, 1]]
    >>> new_C212 = C212.copy()
    >>> new_C212.table[1][3] = 2; new_C212.table[2][2] = 1
    >>> new_C212.scan_check(1, y**3)
    False
    >>> C2121 = C212.copy()
    >>> C2121.table.append([None]*4)
    >>> C2121.table
    >>> C2121.table
    [[1, 1, 0, 0], [0, 0, 2, 3], [None, None, None, 1], [None, None, 1, None]]
    >>> C2121.p = [0, 1, 2, 3]
    >>> C2121.scan(2, y**3)
    >>> C2121.table
    >>> C2121.table
    [[1, 1, 0, 0], [0, 0, 2, 3], [None, None, 3, 1], [None, None, 1, 2]]
    >>> new_C2121 = C2121.copy()
    >>> new_C2121.table[2][0] = 2; new_C2121.table[2][1] = 2
    >>> new_C2121.table
    >>> new_C2121.table
    [[1, 1, 0, 0], [0, 0, 2, 3], [2, 2, 3, 1], [None, None, 1, 2]]
    >>> first_in_class(new_C2121)
    False
    >>> C21211 = C2121.copy()
    >>> C21211.table[2][0] = 3; C21211.table[2][1] = 3;
    >>> C21211.table[3][0] = 2; C21211.table[3][1] = 2;
    >>> C21211.table
    [[1, 1, 0, 0], [0, 0, 2, 3], [3, 3, 3, 1], [2, 2, 1, 2]]
    >>> first_in_class(C21211)  # gives wrong-result
    False
    >>> C21211.scan_check(1, (x*y)**4)
    False
    >>> C22 = C2.copy()
    >>> C22.table[0][2] = 1; C22.table[1][3] = 0
    >>> C22.table
    [[1, 1, 1, None], [0, 0, None, 0]]
    >>> C221 = C22.copy()
    >>> C221.table.append([None]*4)
    >>> C221.p = [0, 1, 2]
    >>> C221.table[2][2] = 0; C221.table[0][3] = 2;
    >>> C221.table
    [[1, 1, 1, 2], [0, 0, None, 0], [None, None, 0, None]]
    >>> C221.scan(1, y**3)
    >>> C221.table
    >>> C221.scan(2, y**3)
    [[1, 1, 1, 2], [0, 0, 2, 0], [None, None, 0, 1]]
    >>> first_in_class(C221)
    True
    >>> C221.table
    [[1, 1, 1, 2], [0, 0, 2, 0], [None, None, 0, 1]]
    >>> new_C221 = C221.copy()
    >>> new_C221.table[2][0] = 2; new_C221.table[2][1] = 2
    >>> first_in_class(new_C221)
    False
    >>> C2211 = C221.copy()
    >>> C2211.table.append([None]*4)
    >>> C2211.p = [0, 1, 2, 3]
    >>> C2211.table[2][0] = 3; C2211.table[2][1] = 3
    >>> C2211.table
    [[1, 1, 1, 2], [0, 0, 2, 0], [3, 3, 0, 1], [None, None, None, None]]
    >>> C2211.table[3][0] = 2; C2211.table[3][1] = 2;
    >>> C2211.table
    [[1, 1, 1, 2], [0, 0, 2, 0], [3, 3, 0, 1], [2, 2, None, None]]
    >>> C2211.table[3][2] = 3; C2211.table[3][3] = 3
    >>> C2211.table
    [[1, 1, 1, 2], [0, 0, 2, 0], [3, 3, 0, 1], [2, 2, 3, 3]]
    >>> first_in_class(C2211)
    False
    >>> C23 = C2.copy()
    >>> C23.table
    [[1, 1, 2, None], [0, 0, None, None]]
    >>> C23.table.append([None]*4);
    >>> C23.p = [0, 1, 2]
    >>> C23.table[0][2] = 2; C23.table[2][3] = 0
    >>> C23.table
    [[1, 1, 2, None], [0, 0, None, None], [None, None, None, 0]]
    >>> C231 = C23.copy()
    >>> C231.table[0][3] = 1; C231.table[1][2] = 0
    >>> C231.table
    [[1, 1, 2, 1], [0, 0, 0, None], [None, None, None, 0]]
    >>> first_in_class(C231)
    False
    >>> C232 = C23.copy()
    >>> C232.table
    [[1, 1, 2, None], [0, 0, None, None], [None, None, None, 0]]
    >>> C232.p
    [0, 1, 2]
    >>> C232.table.append([None]*4)
    >>> C232.table
    [[1, 1, 2, None],
     [0, 0, None, None],
     [None, None, None, 0],
     [None, None, None, None]]
    >>> C232.table[0][3] = 3; C232.table[3][2] = 0
    >>> C232.table
    [[1, 1, 2, 3],
     [0, 0, None, None],
     [None, None, None, 0],
     [None, None, 0, None]]
    >>> C232.scan(2, y**3)
    >>> C232.table
    [[1, 1, 2, 3], [0, 0, None, None], [None, None, 3, 0], [None, None, 0, 2]]
    >>> C232.table[1][2] = 1; C232.table[1][3] = 1
    >>> C232.table
    [[1, 1, 2, 3], [0, 0, 1, 1], [None, None, 3, 0], [None, None, 0, 2]]
    >>> first_in_class(C232)
    False

    """
    n = C.n
    lamda = -1
    nu = [None]*n
    mu = []
    breaker = False
    for alpha in range(1, n):
        # reset ν to "None" after previous value of α
        for beta in range(lamda+1):
            nu[mu[beta]] = None
        # try α as the new point 0 in Ω_C_α
        try:
            mu[0] = alpha
        except IndexError:
            mu.append(alpha)
        nu[alpha] = 0
        # compare corresponding entries in C and C_α
        lamda = 0
        for beta in range(C.n):
            if breaker:
                break
            for x in C.A:
                if C.table[beta][C.A_dict[x]] is None or \
                        C.table[mu[beta]][C.A_dict[x]] is None:
                    # continue with α
                    breaker = True
                    break
                gamma = C.table[beta][C.A_dict[x]]
                delta = C.table[mu[beta]][C.A_dict[x]]
                if nu[delta] is None:
                    # delta becomes the next point in Ω_C_α
                    lamda += 1
                    nu[delta] = lamda
                    try:
                        mu[lamda] = delta
                    except:
                        mu.append(delta)
                if nu[delta] < gamma:
                    return False
                if nu[delta] > gamma:
                    # continue with α
                    breaker = True
                    break
    return True


FpGroupElement = FreeGroupElement

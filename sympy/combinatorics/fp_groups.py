# -*- coding: utf-8 -*-
from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public
from sympy.utilities.iterables import flatten
from sympy.combinatorics.free_group import FreeGroupElement
from itertools import chain, product
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
        self_copy.table = [list(perm_rep) for perm_rep in self.table]
        self_copy.p = list(self.p)
        self_copy.deduction_stack = list(self.deduction_stack)
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
        if not self.table:
            return 0
        return max(self.omega) + 1

    # Pg 152 [1]
    def is_complete(self):
        """
        The coset table is called complete if it has no undefined entries
        on the live cosets; that is, α^x is defined for all α ∈ Ω and x ∈ A.
        """
        return not any([None in self.table[coset] for coset in self.omega])

    # Pg. 153 [1]
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

    # used in the low-index subgroups algorithm
    def scan_check(self, alpha, word):
        """
        Another version of "scan" routine, it checks whether α scans correctly
        under w, it is a straightforward modification of "scan". "scan_check"
        return false (rather than calling "coincidence") if the scan completes
        incorrectly; otherwise it returns true.

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
        elif j == i:
            # deduction process
            self.table[f][A_dict[word[i]]] = b
            self.table[b][A_dict_inv[word[i]]] = f
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
        # complete scan all relators under all cosets(obviously live)
        # without making new definitions
        for beta in self.omega:
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

    def process_deductions_check(self, R_c_x, R_c_x_inv):
        """
        A variation of "process_deductions", this calls "scan_check" wherever
        "process_deductions" calls "scan".
        """
        p = self.p
        while len(self.deduction_stack) > 0:
            alpha, x = self.deduction_stack.pop()
            for w in R_c_x:
                if not self.scan_check(alpha, w):
                    return False
            beta = self.table[alpha][self.A_dict[x]]
            if beta is not None:
                for w in R_c_x_inv:
                    if not self.scan_check(beta, w):
                        return False
        return True

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
        for alpha, x in product(range(self.n), A):
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

    def conjugates(self, R):
        R_c = list(chain.from_iterable((rel.cyclic_conjugates(), \
                (rel**-1).cyclic_conjugates()) for rel in R))
        R_set = set()
        for conjugate in R_c:
            R_set = R_set.union(conjugate)
        R_c_list = []
        for x in self.A:
            r = set([word for word in R_set if word[0] == x])
            R_c_list.append(r)
            R_set.difference_update(r)
        return R_c_list

# relator-based method
def coset_enumeration_r(fp_grp, Y):
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


def low_index_subgroups(G, N, Y=[]):
    """
    Implements the Low Index Subgroups algorithm, i.e find all subgroups of
    "G" upto a given index "N". This implements the method described in
    [Sim94]. This procedure involves a backtrack search over incomplete Coset
    Tables, rather than over forced coincidences.

    G: An FpGroup < X|R >
    N: positive integer, representing the maximun index value for subgroups
    Y: (an optional argument) specifying a list of subgroup generators, such
    that each of the resulting subgroup contains the subgroup generated by Y.

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of Computational Group Theory"
    Section 5.4

    [2] Marston Conder and Peter Dobcsanyi
    "Applications and Adaptions of the Low Index Subgroups Procedure"

    Examples
    ========

    >>> from sympy.combinatorics.free_group import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup, low_index_subgroups
    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**2, y**3, (x*y)**4])
    >>> L = low_index_subgroups(f, 4)
    >>> for coset_table in L:
    ...     print(coset_table.table)
    [[0, 0, 0, 0]]
    [[0, 0, 1, 2], [1, 1, 2, 0], [3, 3, 0, 1], [2, 2, 3, 3]]
    [[0, 0, 1, 2], [2, 2, 2, 0], [1, 1, 0, 1]]
    [[1, 1, 0, 0], [0, 0, 1, 1]]

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
    R1_c_list = C.conjugates(R1)
    S = []
    descendant_subgroups(S, C, R1_c_list, C.A[0], R2, N, Y)
    return S


def descendant_subgroups(S, C, R1_c_list, x, R2, N, Y):
    A_dict = C.A_dict
    A_dict_inv = C.A_dict_inv
    if C.is_complete():
        # if C is complete then it only needs to test
        # whether the relators in R2 are satisfied
        for w, alpha in product(R2, C.omega):
            if not C.scan_check(alpha, w):
                return
        # relators in R2 are satisfied, append the table to list
        S.append(C)
    else:
        # find the first undefined entry in Coset Table
        for alpha, x in product(range(len(C.table)), C.A):
            if C.table[alpha][A_dict[x]] is None:
                # this is "x" in pseudo-code (using "y" makes it clear)
                undefined_coset, undefined_gen = alpha, x
                break
        # for filling up the undefine entry we try all possible values
        # of β ∈ Ω or β = n where β^(undefined_gen^-1) is undefined
        reach = C.omega + [C.n]
        for beta in reach:
            if beta < N:
                if beta == C.n or C.table[beta][A_dict_inv[undefined_gen]] is None:
                    try_descendant(S, C, R1_c_list, R2, N, undefined_coset, \
                            undefined_gen, beta, Y)


def try_descendant(S, C, R1_c_list, R2, N, alpha, x, beta, Y):
    """
    It solves the problem of trying out each individual possibility for α^x.
    """
    D = C.copy()
    A_dict = D.A_dict
    if beta == D.n and beta < N:
        D.table.append([None]*len(D.A))
        D.p.append(beta)
    D.table[alpha][D.A_dict[x]] = beta
    D.table[beta][D.A_dict_inv[x]] = alpha
    D.deduction_stack.append((alpha, x))
    if not D.process_deductions_check(R1_c_list[D.A_dict[x]], \
            R1_c_list[D.A_dict_inv[x]]):
        return
    for w in Y:
        if not D.scan_check(0, w):
            return
    if first_in_class(D, Y):
        descendant_subgroups(S, D, R1_c_list, x, R2, N, Y)


def first_in_class(C, Y=[]):
    """
    Checks whether the subgroup H=G1 corresponding to the Coset Table could
    possibly be the canonical representative of its conjugacy class.

    Parameters
    ==========

    C: CosetTable

    Returns
    =======

    bool: True/False

    If this returns False, then no descendant of C can have that property, and
    so we can abandon C. If it returns True, then we need to process further
    the node of the search tree corresponding to C, and so we call
    ``descendant_subgroups`` recursively on C.

    Examples
    ========

    >>> from sympy.combinatorics.free_group import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup, CosetTable, first_in_class
    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**2, y**3, (x*y)**4])
    >>> C = CosetTable(f, [])
    >>> C.table = [[0, 0, None, None]]
    >>> first_in_class(C)
    True
    >>> C.table = [[1, 1, 1, None], [0, 0, None, 1]]; C.p = [0, 1]
    >>> first_in_class(C)
    True
    >>> C.table = [[1, 1, 2, 1], [0, 0, 0, None], [None, None, None, 0]]
    >>> C.p = [0, 1, 2]
    >>> first_in_class(C)
    False
    >>> C.table = [[1, 1, 1, 2], [0, 0, 2, 0], [2, None, 0, 1]]
    >>> first_in_class(C)
    False

    # TODO:: Sims points out in [Sim94] that performance can be improved by
    # remembering some of the information computed by ``first_in_class``. If
    # the ``continue α`` statement is executed at line 14, then the same thing
    # will happen for that value of α in any descendant of the table C, and so
    # the values the values of α for which this occurs could profitably be
    # stored and passed through to the descendants of C. Of course this would
    # make the code more complicated.

    # The code below is taken directly from the function on page 208 of [Sim94]
    # ν[α]

    """
    n = C.n
    # lamda is the largest numbered point in Ω_c_α which is currently defined
    lamda = -1
    # for α ∈ Ω_c, ν[α] is the point in Ω_c_α corresponding to α
    nu = [None]*n
    # for α ∈ Ω_c_α, μ[α] is the point in Ω_c corresponding to α
    mu = [None]*n
    # mutually ν and μ are the mutually-inverse equivalence maps between
    # Ω_c_α and Ω_c
    next_alpha = False
    # For each 0≠α ∈ [0 .. nc-1], we start by constructing the equivalent
    # standardized coset table C_α corresponding to H_α
    for alpha in range(1, n):
        # reset ν to "None" after previous value of α
        for beta in range(lamda+1):
            nu[mu[beta]] = None
        # we only want to reject our current table in favour of a preceding
        # table in the ordering in which 1 is replaced by α, if the subgroup
        # G_α corresponding to this preceding table definitely contains the
        # given subgroup
        for w in Y:
            # TODO: this should support input of a list of general words
            # not just the words which are in "A" (i.e gen and gen^-1)
            if C.table[alpha][C.A_dict[w]] != alpha:
                # continue with α
                next_alpha = True
                break
        if next_alpha:
            next_alpha = False
            continue
        # try α as the new point 0 in Ω_C_α
        mu[0] = alpha
        nu[alpha] = 0
        # compare corresponding entries in C and C_α
        lamda = 0
        for beta in range(n):
            for x in C.A:
                gamma = C.table[beta][C.A_dict[x]]
                delta = C.table[mu[beta]][C.A_dict[x]]
                # if either of the entries is undefined,
                # we move with next α
                if gamma is None or delta is None:
                    # continue with α
                    next_alpha = True
                    break
                if nu[delta] is None:
                    # delta becomes the next point in Ω_C_α
                    lamda += 1
                    nu[delta] = lamda
                    mu[lamda] = delta
                if nu[delta] < gamma:
                    return False
                if nu[delta] > gamma:
                    # continue with α
                    next_alpha = True
                    break
            if next_alpha:
                next_alpha = False
                break
    return True


FpGroupElement = FreeGroupElement

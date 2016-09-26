# -*- coding: utf-8 -*-
"""Finitely Presented Groups and its algorithms. """

from __future__ import print_function, division
from sympy.core.basic import Basic
from sympy.core import Symbol, Mod
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public
from sympy.utilities.iterables import flatten
from sympy.combinatorics.free_groups import FreeGroupElement, free_group, zero_mul_simp

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
        obj.dtype = type("FpGroupElement", (FpGroupElement,), {"group": obj})

        # CosetTable instance on identity subgroup
        obj._coset_table = None
        # returns whether coset table on identity subgroup
        # has been standardized
        obj._is_standardized = False

        obj._order = None
        obj._center = None
        return obj

    @property
    def free_group(self):
        return self._free_group

    def coset_enumeration(self, H, strategy="relator_based"):
        """
        Return an instance of ``coset table``, when Todd-Coxeter algorithm is
        run over the ``self`` with ``H`` as subgroup, using ``strategy``
        argument as strategy. The returned coset table is compressed but not
        standardized.

        """
        if strategy == 'relator_based':
            C = coset_enumeration_r(self, H)
        else:
            C = coset_enumeration_c(self, H)
        C.compress()
        return C

    def standardize_coset_table(self):
        """
        Standardized the coset table ``self`` and makes the internal variable
        ``_is_standardized`` equal to ``True``.

        """
        self._coset_table.standardize()
        self._is_standardized = True

    def coset_table(self, H, strategy="relator_based"):
        """
        Return the mathematical coset table of ``self`` in ``H``.

        """
        if not H:
            if self._coset_table != None:
                if not self._is_standardized:
                    self.standardize_coset_table()
            else:
                C = self.coset_enumeration([], strategy)
                self._coset_table = C
                self.standardize_coset_table()
            return self._coset_table.table
        else:
            C = self.coset_enumeration(H, strategy)
            C.standardize()
            return C.table

    def order(self, strategy="relator_based"):
        """
        Returns the order of the finitely presented group ``self``. It uses
        the coset enumeration with identity group as subgroup, i.e ``H=[]``.

        Examples
        ========

        >>> from sympy.combinatorics.free_groups import free_group
        >>> from sympy.combinatorics.fp_groups import FpGroup
        >>> F, x, y = free_group("x, y")
        >>> f = FpGroup(F, [x, y**2])
        >>> f.order(strategy="coset_table_based")
        2

        """
        if self._order != None:
            return self._order
        if self._coset_table != None:
            self._order = len(self._coset_table.table)
        else:
            self._coset_table = self.coset_enumeration([], strategy)
            self._order = len(self._coset_table.table)
        return self._order

    def index(self, H, strategy="relator_based"):
        """
        Returns the index of subgroup ``H`` in group ``self``.

        Examples
        ========

        >>> from sympy.combinatorics.free_groups import free_group
        >>> from sympy.combinatorics.fp_groups import FpGroup
        >>> F, x, y = free_group("x, y")
        >>> f = FpGroup(F, [x**5, y**4, y*x*y**3*x**3])
        >>> f.index([x])
        4

        """
        # TODO: use |G:H| = |G|/|H| (currently H can't be made into a group)
        # when we know |G| and |H|

        if H == []:
            return self.order()
        else:
            C = self.coset_enumeration(H, strategy)
            return len(C.table)

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


###############################################################################
#                           COSET TABLE                                       #
###############################################################################

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

    r"""

    Properties
    ==========

    [1] `0 \in \Omega` and `\tau(1) = \epsilon`
    [2] `\alpha^x = \beta \Leftrightarrow \beta^{x^{-1}} = \alpha`
    [3] If `\alpha^x = \beta`, then `H \tau(\alpha)x = H \tau(\beta)`
    [4] `\forall \alpha \in \Omega, 1^{\tau(\alpha)} = \alpha`

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of Computational Group Theory"

    [2] John J. Cannon; Lucien A. Dimino; George Havas; Jane M. Watson
    Mathematics of Computation, Vol. 27, No. 123. (Jul., 1973), pp. 463-490.
    "Implementation and Analysis of the Todd-Coxeter Algorithm"

    """
    # default limit for the number of cosets allowed in a
    # coset enumeration.
    coset_table_max_limit = 4096000
    # maximum size of deduction stack above or equal to
    # which it is emptied
    max_stack_size = 500

    def __init__(self, fp_grp, subgroup):
        self.fp_group = fp_grp
        self.subgroup = subgroup
        # "p" is setup independent of Ω and n
        self.p = [0]
        # a list of the form `[gen_1, gen_1^{-1}, ... , gen_k, gen_k^{-1}]`
        self.A = list(chain.from_iterable((gen, gen**-1) \
                for gen in self.fp_group.generators))
        # the mathematical coset table which is a list of lists
        self.table = [[None]*len(self.A)]
        self.A_dict = {x: self.A.index(x) for x in self.A}
        self.A_dict_inv = {}
        for x, index in self.A_dict.items():
            if index % 2 == 0:
                self.A_dict_inv[x] = self.A_dict[x] + 1
            else:
                self.A_dict_inv[x] = self.A_dict[x] - 1
        # used in the coset-table based method of coset enumeration. Each of
        # the element is called a "deduction" which is the form (α, x) whenever
        # a value is assigned to α^x during a definition or "deduction process"
        self.deduction_stack = []

    @property
    def omega(self):
        """Set of live cosets. """
        return [coset for coset in range(len(self.p)) if self.p[coset] == coset]

    def copy(self):
        """
        Return a shallow copy of Coset Table instance ``self``.

        """
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
        """The number `n` represents the length of the sublist containing the
        live cosets.

        """
        if not self.table:
            return 0
        return max(self.omega) + 1

    # Pg. 152 [1]
    def is_complete(self):
        r"""
        The coset table is called complete if it has no undefined entries
        on the live cosets; that is, `\alpha^x` is defined for all
        `\alpha \in \Omega` and `x \in A`.

        """
        return not any([None in self.table[coset] for coset in self.omega])

    # Pg. 153 [1]
    def define(self, alpha, x):
        r"""
        This routine is used in the relator-based strategy of Todd-Coxeter
        algorithm if some `\alpha^x` is undefined. We check whether there is
        space available for defining a new coset. If there is enough space
        then we remedy this by adjoining a new coset `\beta` to `\Omega`
        (i.e to set of live cosets) and put that equal to `\alpha^x`, then
        make an assignment satisfying Property[1]. If there is not enough space
        then we halt the Coset Table creation. The maximum amount of space that
        can be used by Coset Table can be manipulated using the class variable
        ``CosetTable.coset_table_max_limit``.

        See Also
        ========
        define_c

        """
        A = self.A
        table = self.table
        len_table = len(table)
        if len_table == CosetTable.coset_table_max_limit:
            # abort the further generation of cosets
            raise ValueError("the coset enumeration has defined more than "
                    "%s cosets. Try with a greater value max number of cosets "
                    % CosetTable.coset_table_max_limit)
        table.append([None]*len(A))
        # beta is the new coset generated
        beta = len_table
        self.p.append(beta)
        table[alpha][self.A_dict[x]] = beta
        table[beta][self.A_dict_inv[x]] = alpha

    def define_c(self, alpha, x):
        r"""
        A variation of ``define`` routine, described on Pg. 165 [1], used in
        the coset table-based strategy of Todd-Coxeter algorithm. It differs
        from ``define`` routine in that for each definition it also adds the
        tuple `(\alpha, x)` to the deduction stack.

        See Also
        ========
        define

        """
        A = self.A
        table = self.table
        len_table = len(table)
        if len_table == CosetTable.coset_table_max_limit:
            # abort the further generation of cosets
            raise ValueError("the coset enumeration has defined more than "
                    "%s cosets. Try with a greater value max number of cosets "
                    % CosetTable.coset_table_max_limit)
        table.append([None]*len(A))
        # beta is the new coset generated
        beta = len_table
        self.p.append(beta)
        table[alpha][self.A_dict[x]] = beta
        table[beta][self.A_dict_inv[x]] = alpha
        # append to deduction stack
        self.deduction_stack.append((alpha, x))

    def scan_c(self, alpha, word):
        """
        A variation of ``scan`` routine, described on pg. 165 of [1], which
        puts at tuple, whenever a deduction occurs, to deduction stack.

        See Also
        ========
        scan, scan_check, scan_and_fill, scan_and_fill_c

        """
        # α is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for α=0 and w in Y (i.e generating set of H)
        # 2. α in Ω (set of live cosets), w in R (relators)
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        # list of union of generators and their inverses
        while i <= j and table[f][A_dict[word[i]]] is not None:
            f = table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            if f != b:
                self.coincidence_c(f, b)
            return
        while j >= i and table[b][A_dict_inv[word[j]]] is not None:
            b = table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # run the "coincidence" routine
            self.coincidence_c(f, b)
        elif j == i:
            # deduction process
            table[f][A_dict[word[i]]] = b
            table[b][A_dict_inv[word[i]]] = f
            self.deduction_stack.append((f, word[i]))
        # otherwise scan is incomplete and yields no information

    # α, β coincide, i.e. α, β represent the pair of cosets where
    # coincidence occurs
    def coincidence_c(self, alpha, beta):
        """
        A variation of ``coincidence`` routine used in the coset-table based
        method of coset enumeration. The only difference being on addition of
        a new coset in coset table(i.e new coset introduction), then it is
        appended to ``deduction_stack``.

        See Also
        ========
        coincidence

        """
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        p = self.p
        l = 0
        # behaves as a queue
        q = []
        self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            for x in A_dict:
                delta = table[gamma][A_dict[x]]
                if delta is not None:
                    table[delta][A_dict_inv[x]] = None
                    # only line of difference from ``coincidence`` routine
                    self.deduction_stack.append((delta, x**-1))
                    mu = self.rep(gamma)
                    nu = self.rep(delta)
                    if table[mu][A_dict[x]] is not None:
                        self.merge(nu, table[mu][A_dict[x]], q)
                    elif table[nu][A_dict_inv[x]] is not None:
                        self.merge(mu, table[nu][A_dict_inv[x]], q)
                    else:
                        table[mu][A_dict[x]] = nu
                        table[nu][A_dict_inv[x]] = mu

    def scan(self, alpha, word):
        r"""
        ``scan`` performs a scanning process on the input ``word``.
        It first locates the largest prefix ``s`` of ``word`` for which
        `\alpha^s` is defined (i.e is not ``None``), ``s`` may be empty. Let
        ``word=sv``, let ``t`` be the longest suffix of ``v`` for which
        `\alpha^{t^{-1}}` is defined, and let ``v=ut``. Then three
        possibilities are there:

        1. If ``t=v``, then we say that the scan completes, and if, in addition
        `\alpha^s = \alpha^{t^{-1}}`, then we say that the scan completes
        correctly.

        2. It can also happen that scan does not complete, but `|u|=1`; that
        is, the word ``u`` consists of a single generator `x \in A`. In that
        case, if `\alpha^s = \beta` and `\alpha^{t^{-1}} = \gamma`, then we can
        set `\beta^x = \gamma` and `\gamma^{x^{-1}} = \beta`. These assignments
        are known as deductions and enable the scan to complete correctly.

        3. See ``coicidence`` routine for explanation of third condition.

        Notes
        =====
        The code for the procedure of scanning `\alpha \in \Omega`
        under `w \in A*` is defined on pg. 155 [1]

        See Also
        ========
        scan_c, scan_check, scan_and_fill, scan_and_fill_c

        """
        # α is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for α=0 and w in Y (i.e generating set of H)
        # 2. α in Ω (set of live cosets), w in R (relators)
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        while i <= j and table[f][A_dict[word[i]]] is not None:
            f = table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            if f != b:
                self.coincidence(f, b)
            return
        while j >= i and table[b][A_dict_inv[word[j]]] is not None:
            b = table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # run the "coincidence" routine
            self.coincidence(f, b)
        elif j == i:
            # deduction process
            table[f][A_dict[word[i]]] = b
            table[b][A_dict_inv[word[i]]] = f
        # otherwise scan is incomplete and yields no information

    # used in the low-index subgroups algorithm
    def scan_check(self, alpha, word):
        r"""
        Another version of ``scan`` routine, described on, it checks whether
        `\alpha` scans correctly under `word`, it is a straightforward
        modification of ``scan``. ``scan_check`` returns ``False`` (rather than
        calling ``coincidence``) if the scan completes incorrectly; otherwise
        it returns ``True``.

        See Also
        ========
        scan, scan_c, scan_and_fill, scan_and_fill_c

        """
        # α is an integer representing a "coset"
        # since scanning can be in two cases
        # 1. for α=0 and w in Y (i.e generating set of H)
        # 2. α in Ω (set of live cosets), w in R (relators)
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        while i <= j and table[f][A_dict[word[i]]] is not None:
            f = table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            return f == b
        while j >= i and table[b][A_dict_inv[word[j]]] is not None:
            b = table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            # we have an incorrect completed scan with coincidence f ~ b
            # return False, instead of calling coincidence routine
            return False
        elif j == i:
            # deduction process
            table[f][A_dict[word[i]]] = b
            table[b][A_dict_inv[word[i]]] = f
        return True

    def merge(self, k, lamda, q):
        """
        Input: 'k', 'lamda' being the two class representatives to be merged.
        =====

        Merge two classes with representatives ``k`` and ``lamda``, described
        on Pg. 157 [1] (for pseudocode), start by putting ``p[k] = lamda``.
        It is more efficient to choose the new representative from the larger
        of the two classes being merged, i.e larger among ``k`` and ``lamda``.
        procedure ``merge`` performs the merging operation, adds the deleted
        class representative to the queue ``q``.

        Notes
        =====
        Pg. 86-87 [1] contains a description of this method.

        See Also
        ========
        coincidence, rep

        """
        p = self.p
        phi = self.rep(k)
        psi = self.rep(lamda)
        if phi != psi:
            mu = min(phi, psi)
            v = max(phi, psi)
            p[v] = mu
            q.append(v)

    def rep(self, k):
        r"""
        Input: `k \in [0 \ldots n-1]`, as for ``self`` only array ``p`` is used
        =====
        Output: Representative of the class containing ``k``.
        ======

        Returns the representative of `\sim` class containing ``k``, it also
        makes some modification to array ``p`` of ``self`` to ease further
        computations, described on Pg. 157 [1].

        The information on classes under `\sim` is stored in array `p` of
        ``self`` argument, which will always satisfy the property:

        `p[\alpha] \sim \alpha` and `p[\alpha]=\alpha \iff \alpha=rep(\alpha)`
        `\forall \in [0 \ldots n-1]`.

        So, for `\alpha \in [0 \ldots n-1]`, we find `rep(self, \alpha)` by
        continually replacing `\alpha` by `p[\alpha]` until it becomes
        constant (i.e satisfies `p[\alpha] = \alpha`):w

        To increase the efficiency of later ``rep`` calculations, whenever we
        find `rep(self, \alpha)=\beta`, we set
        `p[\gamma] = \beta \forall \gamma \in p-chain` from `\alpha` to `\beta`

        Notes
        =====
        ``rep`` routine is also described on Pg. 85-87 [1] in Atkinson's
        algorithm, this results from the fact that ``coincidence`` routine
        introduces functionality similar to that introduced by the
        ``minimal_block`` routine on Pg. 85-87 [1].

        See also
        ========
        coincidence, merge

        """
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

    # α, β coincide, i.e. α, β represent the pair of cosets
    # where coincidence occurs
    def coincidence(self, alpha, beta):
        r"""
        The third situation described in ``scan`` routine is handled by this
        routine, described on Pg. 156-161 [1].

        The unfortunate situation when the scan completes but not correctly,
        then ``coincidence`` routine is run. i.e when for some `i` with
        `1 \le i \le r+1`, we have `w=st` with `s=x_1*x_2 ... x_{i-1}`,
        `t=x_i*x_{i+1} ... x_r`, and `\beta = \alpha^s` and
        `\gamma = \alph^{t-1}` are defined but unequal. This means that
        `\beta` and `\gamma` represent the same coset of `H` in `G`. Described
        on Pg. 156 [1]. ``rep``

        See Also
        ========
        scan

        """
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        p = self.p
        l = 0
        # behaves as a queue
        q = []
        self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            for x in A_dict:
                delta = table[gamma][A_dict[x]]
                if delta is not None:
                    table[delta][A_dict_inv[x]] = None
                    mu = self.rep(gamma)
                    nu = self.rep(delta)
                    if table[mu][A_dict[x]] is not None:
                        self.merge(nu, table[mu][A_dict[x]], q)
                    elif table[nu][A_dict_inv[x]] is not None:
                        self.merge(mu, table[nu][A_dict_inv[x]], q)
                    else:
                        table[mu][A_dict[x]] = nu
                        table[nu][A_dict_inv[x]] = mu

    # method used in the HLT strategy
    def scan_and_fill(self, alpha, word):
        """
        A modified version of ``scan`` routine used in the relator-based
        method of coset enumeration, described on pg. 162-163 [1], which
        follows the idea that whenever the procedure is called and the scan
        is incomplete then it makes new definitions to enable the scan to
        complete; i.e it fills in the gaps in the scan of the relator or
        subgroup generator.

        """
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        r = len(word)
        f = alpha
        i = 0
        b = alpha
        j = r - 1
        # loop until it has filled the α row in the table.
        while True:
            # do the forward scanning
            while i <= j and table[f][A_dict[word[i]]] is not None:
                f = table[f][A_dict[word[i]]]
                i += 1
            if i > j:
                if f != b:
                    self.coincidence(f, b)
                return
            # forward scan was incomplete, scan backwards
            while j >= i and table[b][A_dict_inv[word[j]]] is not None:
                b = table[b][A_dict_inv[word[j]]]
                j -= 1
            if j < i:
                self.coincidence(f, b)
            elif j == i:
                table[f][A_dict[word[i]]] = b
                table[b][A_dict_inv[word[i]]] = f
            else:
                self.define(f, word[i])

    def scan_and_fill_c(self, alpha, word):
        """
        A modified version of ``scan`` routine, described on Pg. 165 second
        para. [1], with modification similar to that of ``scan_anf_fill`` the
        only difference being it calls the coincidence procedure used in the
        coset-table based method i.e. the routine ``coincidence_c`` is used.

        Also See
        ========
        scan, scan_and_fill

        """
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        r = len(word)
        f = alpha
        i = 0
        b = alpha
        j = r - 1
        # loop until it has filled the α row in the table.
        while True:
            # do the forward scanning
            while i <= j and table[f][A_dict[word[i]]] is not None:
                f = table[f][A_dict[word[i]]]
                i += 1
            if i > j:
                if f != b:
                    self.coincidence_c(f, b)
                return
            # forward scan was incomplete, scan backwards
            while j >= i and table[b][A_dict_inv[word[j]]] is not None:
                b = table[b][A_dict_inv[word[j]]]
                j -= 1
            if j < i:
                self.coincidence_c(f, b)
            elif j == i:
                table[f][A_dict[word[i]]] = b
                table[b][A_dict_inv[word[i]]] = f
                self.deduction_stack.append((f, word[i]))
            else:
                self.define_c(f, word[i])

    # method used in the HLT strategy
    def look_ahead(self):
        """
        When combined with the HLT method this is known as HLT+Lookahead
        method of coset enumeration, described on pg. 164 [1]. Whenever
        ``define`` aborts due to lack of space available this procedure is
        executed. This routine helps in recovering space resulting from
        "coincidence" of cosets.

        """
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
        """
        Processes the deductions that have been pushed onto ``deduction_stack``,
        described on Pg. 166 [1] and is used in coset-table based enumeration.

        See Also
        ========
        deduction_stack

        """
        p = self.p
        table = self.table
        while len(self.deduction_stack) > 0:
            if len(self.deduction_stack) >= CosetTable.max_stack_size:
                self.look_ahead()
                del self.deduction_stack[:]
            else:
                alpha, x = self.deduction_stack.pop()
                if p[alpha] == alpha:
                    for w in R_c_x:
                        self.scan_c(alpha, w)
                        if p[alpha] < alpha:
                            break
            beta = table[alpha][self.A_dict[x]]
            if beta is not None and p[beta] == beta:
                for w in R_c_x_inv:
                    self.scan_c(beta, w)
                    if p[beta] < beta:
                        break

    def process_deductions_check(self, R_c_x, R_c_x_inv):
        """
        A variation of ``process_deductions``, this calls ``scan_check``
        wherever ``process_deductions`` calls ``scan``, described on Pg. [1].

        See Also
        ========
        process_deductions

        """
        p = self.p
        table = self.table
        while len(self.deduction_stack) > 0:
            alpha, x = self.deduction_stack.pop()
            for w in R_c_x:
                if not self.scan_check(alpha, w):
                    return False
            beta = table[alpha][self.A_dict[x]]
            if beta is not None:
                for w in R_c_x_inv:
                    if not self.scan_check(beta, w):
                        return False
        return True

    def switch(self, beta, gamma):
        r"""Switch the elements `\beta, \gamma \in \Omega` of ``self``, used
        by the ``standardize`` procedure, described on Pg. 167 [1].

        See Also
        ========
        standardize

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
        """
        A coset table is standardized if when running through the cosets and
        within each coset through the generator images (ignoring generator
        inverses), the cosets appear in order of the integers
        `0, 1, , \ldots, n`. "Standardize" reorders the elements of `\Omega`
        such that, if we scan the coset table first by elements of `\Omega`
        and then by elements of A, then the cosets occur in ascending order.
        ``standardize()`` is used at the end of an enumeration to permute the
        cosets so that they occur in some sort of standard order.

        Notes
        =====
        procedure is described on pg. 167-168 [1], it also makes use of the
        ``switch`` routine to replace by smaller integer value.

        >>> from sympy.combinatorics.free_groups import free_group
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
    def compress(self):
        """Removes the non-live cosets from the coset table, described on
        pg. 167 [1].

        """
        gamma = -1
        A = self.A
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        chi = tuple([i for i in range(len(self.p)) if self.p[i] != i])
        for alpha in self.omega:
            gamma += 1
            if gamma != alpha:
                # replace α by γ in coset table
                for x in A:
                    beta = table[alpha][A_dict[x]]
                    table[gamma][A_dict[x]] = beta
                    table[beta][A_dict_inv[x]] == gamma
        # all the cosets in the table are live cosets
        self.p = list(range(gamma + 1))
        # delete the useless coloumns
        del table[len(self.p):]
        # re-define values
        for row in table:
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


###############################################################################
#                           COSET ENUMERATION                                 #
###############################################################################

# relator-based method
def coset_enumeration_r(fp_grp, Y):
    """
    This is easier of the two implemented methods of coset enumeration.
    and is often called the HLT method, after Hazelgrove, Leech, Trotter
    The idea is that we make use of ``scan_and_fill`` makes new definitions
    whenever the scan is incomplete to enable the scan to complete; this way
    we fill in the gaps in the scan of the relator or subgroup generator,
    that's why the name relator-based method.

    # TODO: complete the docstring

    See Also
    ========
    scan_and_fill,

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of computational group theory"

    >>> from sympy.combinatorics.free_groups import free_group
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
                # if α was eliminated during the scan then break
                if p[alpha] < alpha:
                    break
            if p[alpha] == alpha:
                for x in A_dict:
                    if C.table[alpha][A_dict[x]] is None:
                        C.define(alpha, x)
        alpha += 1
    return C


# Pg. 166
# coset-table based method
def coset_enumeration_c(fp_grp, Y):
    """
    >>> from sympy.combinatorics.free_groups import free_group
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
        C.scan_and_fill_c(0, w)
    for x in A:
        C.process_deductions(R_c_list[C.A_dict[x]], R_c_list[C.A_dict_inv[x]])
    alpha = 0
    while alpha < len(C.table):
        if C.p[alpha] == alpha:
            for x in C.A:
                if C.p[alpha] != alpha:
                    break
                if C.table[alpha][C.A_dict[x]] is None:
                    C.define_c(alpha, x)
                    C.process_deductions(R_c_list[C.A_dict[x]], R_c_list[C.A_dict_inv[x]])
        alpha += 1
    return C


###############################################################################
#                           LOW INDEX SUBGROUPS                               #
###############################################################################

def low_index_subgroups(G, N, Y=[]):
    """
    Implements the Low Index Subgroups algorithm, i.e find all subgroups of
    ``G`` upto a given index ``N``. This implements the method described in
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

    >>> from sympy.combinatorics.free_groups import free_group
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
    r"""
    Solves the problem of trying out each individual possibility
    for `\alpha^x.

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
    Checks whether the subgroup ``H=G1`` corresponding to the Coset Table
    could possibly be the canonical representative of its conjugacy class.

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

    >>> from sympy.combinatorics.free_groups import free_group
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


###############################################################################
#                           SUBGROUP PRESENTATIONS                            #
###############################################################################

# Pg 175 [1]
def define_schreier_generators(C):
    y = []
    gamma = 1
    f = C.fp_group
    X = f.generators
    C.P = [[None]*len(C.A) for i in range(C.n)]
    for alpha, x in product(C.omega, C.A):
        beta = C.table[alpha][C.A_dict[x]]
        if beta == gamma:
            C.P[alpha][C.A_dict[x]] = "<identity>"
            C.P[beta][C.A_dict_inv[x]] = "<identity>"
            gamma += 1
        elif x in X and C.P[alpha][C.A_dict[x]] is None:
            y_alpha_x = '%s_%s' % (x, alpha)
            y.append(y_alpha_x)
            C.P[alpha][C.A_dict[x]] = y_alpha_x
    grp_gens = list(free_group(', '.join(y)))
    C._schreier_free_group = grp_gens.pop(0)
    C._schreier_generators = grp_gens
    # replace all elements of P by, free group elements
    for i, j in product(range(len(C.P)), range(len(C.A))):
        # if equals "<identity>", replace by identity element
        if C.P[i][j] == "<identity>":
            C.P[i][j] = C._schreier_free_group.identity
        elif isinstance(C.P[i][j], str):
            r = C._schreier_generators[y.index(C.P[i][j])]
            C.P[i][j] = r
            beta = C.table[i][j]
            C.P[beta][j + 1] = r**-1


def reidemeister_relators(C):
    R = C.fp_group.relators()
    rels = [rewrite(C, coset, word) for word in R for coset in range(C.n)]
    identity = C._schreier_free_group.identity
    order_1_gens = set([i for i in rels if len(i) == 1])

    # remove all the order 1 generators from relators
    rels = list(filter(lambda rel: rel not in order_1_gens, rels))

    # replace order 1 generators by identity element in reidemeister relators
    for i in range(len(rels)):
        w = rels[i]
        for gen in order_1_gens:
            w = w.eliminate_word(gen, identity)
        rels[i] = w

    C._schreier_generators = [i for i in C._schreier_generators if i not in order_1_gens]

    # Tietze transformation 1 i.e TT_1
    # remove cyclic conjugate elements from relators
    i = 0
    while i < len(rels):
        w = rels[i]
        j = i + 1
        while j < len(rels):
            if w.is_cyclic_conjugate(rels[j]):
                del rels[j]
            else:
                j += 1
        i += 1

    C._reidemeister_relators = rels


def rewrite(C, alpha, w):
    """
    Parameters
    ----------

    C: CosetTable
    α: A live coset
    w: A word in `A*`

    Returns
    -------

    ρ(τ(α), w)

    Examples
    ========

    >>> from sympy.combinatorics.fp_groups import FpGroup, CosetTable, define_schreier_generators, rewrite
    >>> from sympy.combinatorics.free_groups import free_group
    >>> F, x, y = free_group("x ,y")
    >>> f = FpGroup(F, [x**2, y**3, (x*y)**6])
    >>> C = CosetTable(f, [])
    >>> C.table = [[1, 1, 2, 3], [0, 0, 4, 5], [4, 4, 3, 0], [5, 5, 0, 2], [2, 2, 5, 1], [3, 3, 1, 4]]
    >>> C.p = [0, 1, 2, 3, 4, 5]
    >>> define_schreier_generators(C)
    >>> rewrite(C, 0, (x*y)**6)
    x_4*y_2*x_3*x_1*x_2*y_4*x_5

    """
    v = C._schreier_free_group.identity
    for i in range(len(w)):
        x_i = w[i]
        v = v*C.P[alpha][C.A_dict[x_i]]
        alpha = C.table[alpha][C.A_dict[x_i]]
    return v


# Pg 350, section 2.5.1 from [2]
def elimination_technique_1(C):
    rels = C._reidemeister_relators
    # the shorter relators are examined first so that generators selected for
    # elimination will have shorter strings as equivalent
    rels.sort(reverse=True)
    gens = C._schreier_generators
    redundant_gens = {}
    contained_gens = []
    # examine each relator in relator list for any generator occuring exactly
    # once
    next_i = False
    for i in range(len(rels) -1, -1, -1):
        rel = rels[i]
        # don't look for a new generator occuring once in relator which
        # has already found to posses a
        for gen in redundant_gens:
            gen_sym = gen.array_form[0][0]
            if any([gen_sym == r[0] for r in rel.array_form]):
                next_i = True
                break
        if next_i:
            next_i = False
            continue
        for j in range(len(gens) - 1, -1, -1):
            gen = gens[j]
            if rel.generator_count(gen) == 1 and gen not in contained_gens:
                k = rel.exponent_sum(gen)
                gen_index = rel.index(gen**k)
                bk = rel.subword(gen_index + 1, len(rel))
                fw = rel.subword(0, gen_index)
                chi = (bk*fw).identity_cyclic_reduction()
                redundant_gens[gen] = chi**(-1*k)
                contained_gens.extend(chi.contains_generators())
                del rels[i]; del gens[j]
                break
    # eliminate the redundant generator from remaing relators
    for i, gen in product(range(len(rels)), redundant_gens):
        rels[i] = (rels[i].eliminate_word(gen, redundant_gens[gen])).identity_cyclic_reduction()
    rels.sort()
    try:
        rels.remove(C._schreier_free_group.identity)
    except ValueError:
        pass
    C._reidemeister_relators = rels
    C._schreier_generators = gens

# Pg 350, section 2.5.2 from [2]
def elimination_technique_2(C):
    """
    This technique eliminates one generator at a time. Heuristically this
    seems superior in that we may select for elimination the generator with
    shortest equivalent string at each stage.

    >>> from sympy.combinatorics.free_groups import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_r, \
            reidemeister_relators, define_schreier_generators, elimination_technique_2
    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**3, y**5, (x*y)**2]); H = [x*y, x**-1*y**-1*x*y*x]
    >>> C = coset_enumeration_r(f, H)
    >>> C.compress(); C.standardize()
    >>> define_schreier_generators(C)
    >>> reidemeister_relators(C)
    >>> elimination_technique_2(C)
    ([y_1, y_2], [y_2**-3, y_2*y_1*y_2*y_1*y_2*y_1, y_1**2])

    """
    rels = C._reidemeister_relators
    rels.sort(reverse=True)
    gens = C._schreier_generators
    for i in range(len(gens) - 1, -1, -1):
        rel = rels[i]
        for j in range(len(gens) - 1, -1, -1):
            gen = gens[j]
            if rel.generator_count(gen) == 1:
                k = rel.exponent_sum(gen)
                gen_index = rel.index(gen**k)
                bk = rel.subword(gen_index + 1, len(rel))
                fw = rel.subword(0, gen_index)
                rep_by = (bk*fw)**(-1*k)
                del rels[i]; del gens[j]
                for l in range(len(rels)):
                    rels[l] = rels[l].eliminate_word(gen, rep_by)
                break
    C._reidemeister_relators = rels
    C._schreier_generators = gens
    return C._schreier_generators, C._reidemeister_relators

def simplify_presentation(C):
    """Relies upon ``_simplification_technique_1`` for its functioning. """
    rels = C._reidemeister_relators
    rels_arr = _simplification_technique_1(rels)
    group = C._schreier_free_group

    # don't add "identity" element in relator list
    C._reidemeister_relators = [group.dtype(tuple(r)).identity_cyclic_reduction() for r in rels_arr if r]

def _simplification_technique_1(rels):
    """
    All relators are checked to see if they are of the form `gen^n`. If any
    such relators are found then all other relators are processed for strings
    in the `gen` known order.

    Examples
    ========

    >>> from sympy.combinatorics.free_groups import free_group
    >>> from sympy.combinatorics.fp_groups import _simplification_technique_1
    >>> F, x, y = free_group("x, y")
    >>> w1 = [x**2*y**4, x**3]
    >>> _simplification_technique_1(w1)
    [[(x, 3)], [(x, -1), (y, 4)]]

    >>> w2 = [x**2*y**-4*x**5, x**3, x**2*y**8, y**5]
    >>> _simplification_technique_1(w2)
    [[(x, 3)], [(y, 5)], [(x, -1), (y, -2)], [(x, -1), (y, 1), (x, -1)]]

    >>> w3 = [x**6*y**4, x**4]
    >>> _simplification_technique_1(w3)
    [[(x, 4)], [(x, 2), (y, 4)]]

    """
    rels = list(set(rels))
    rels.sort()
    l_rels = len(rels)

    # all syllables with single syllable
    one_syllable_rels = set()
    # since "nw" has a max size = l_rels, only identity element
    # removal can possibly happen
    nw = [None]*l_rels
    for i in range(l_rels):
        w = rels[i].identity_cyclic_reduction()
        if w.number_syllables() == 1:

            # replace one syllable relator with the corresponding inverse
            # element, for ex. x**-4 -> x**4 in relator list
            if w.array_form[0][1] < 0:
                rels[i] = w**-1
            one_syllable_rels.add(rels[i])

        # since modifies the array rep., so should be
        # added a list
        nw[i] = list(rels[i].array_form)

    # bound the exponent of relators, making use of the single
    # syllable relators
    for i in range(l_rels):
        k = nw[i]
        rels_i = rels[i]
        for gen in one_syllable_rels:
            n = gen.array_form[0][1]
            gen_arr0 = gen.array_form[0][0]
            j = len(k) - 1
            while j >= 0:
                if gen_arr0 == k[j][0] and gen is not rels_i:
                    t = Mod(k[j][1], n)

                    # multiple of one syllable relator
                    if t == 0:
                        del k[j]
                        zero_mul_simp(k, j - 1)
                        j = len(k)

                    # power should be bounded by (-n/2, n/2]
                    elif t <= n/2:
                        k[j] = k[j][0], Mod(k[j][1], n)
                    elif t > n/2:
                        k[j] = k[j][0], Mod(k[j][1], n) - n
                j -= 1

    return nw

def reidemeister_presentation(fp_grp, H, elm_rounds=2, simp_rounds=2):
    """
    fp_group: A finitely presented group, an instance of FpGroup
    H: A subgroup whose presentation is to be found, given as a list
    of words in generators of `fp_grp`

    Examples
    ========

    >>> from sympy.combinatorics.free_groups import free_group
    >>> from sympy.combinatorics.fp_groups import FpGroup, reidemeister_presentation
    >>> F, x, y = free_group("x, y")

    Example 5.6 Pg. 177 from [1]
    >>> f = FpGroup(F, [x**3, y**5, (x*y)**2])
    >>> H = [x*y, x**-1*y**-1*x*y*x]
    >>> reidemeister_presentation(f, H)
    ((y_1, y_2), (y_1**2, y_2**3, y_2*y_1*y_2*y_1*y_2*y_1))

    Example 5.8 Pg. 183 from [1]
    >>> f = FpGroup(F, [x**3, y**3, (x*y)**3])
    >>> H = [x*y, x*y**-1]
    >>> reidemeister_presentation(f, H)
    ((x_0, y_0), (x_0**3, y_0**3, x_0*y_0*x_0*y_0*x_0*y_0))

    Exercises Q2. Pg 187 from [1]
    >>> f = FpGroup(F, [x**2*y**2, y**-1*x*y*x**-3])
    >>> H = [x]
    >>> reidemeister_presentation(f, H)
    ((x_0,), (x_0**4,))

    Example 5.9 Pg. 183 from [1]
    >>> f = FpGroup(F, [x**3*y**-3, (x*y)**3, (x*y**-1)**2])
    >>> H = [x]
    >>> reidemeister_presentation(f, H)
    ((x_0,), (x_0**6,))

    """
    C = coset_enumeration_r(fp_grp, H)
    C.compress(); C.standardize()
    define_schreier_generators(C)
    reidemeister_relators(C)
    for i in range(20):
        elimination_technique_1(C)
        simplify_presentation(C)
    C.schreier_generators = tuple(C._schreier_generators)
    C.reidemeister_relators = tuple(C._reidemeister_relators)
    return C.schreier_generators, C.reidemeister_relators


FpGroupElement = FreeGroupElement

from sympy.ntheory import isprime, sieve
from sympy.combinatorics.permutations import _new_from_array_form

############################################
###
### Utilities for computational group theory
###
############################################

def _check_cycles_alt_sym(perm):
    """
    Checks for cycles of prime length p with n/2 < p < n-2.

    Here `n` is the degree of the permutation. This is a helper function for
    the function is_alt_sym from sympy.combinatorics.perm_groups.

    """
    n = perm.size
    af = perm.array_form
    current_len = 0
    total_len = 0
    used = set()
    for i in xrange(n//2):
        if not i in used and i < n//2 - total_len:
            current_len = 1
            used.add(i)
            j = i
            while(af[j] != i):
                current_len += 1
                j = af[j]
                used.add(j)
            total_len += current_len
            if current_len > n//2 and current_len < n-2 and isprime(current_len):
                return True
    return False

def _strip(g, base, orbs, transversals):
    """
    Attempt to decompose a group element using a subgroup chain with orbits
    ``orbs``.

    This is done by treating the subgroup chain as a chain of
    stabilizers with respect to the sequence ``base`` (even though the groups
    might only be subgroups of the respective stabilizers),
    and treat ``orbs`` as basic orbits.
    This process is called "sifting". A sift is unsuccessful when a certain
    orbit element is not found or when after the sift the decomposition
    doesn't end with the identity element.
    The argument ``transversals`` is a list of dictionaries that provides
    transversal elements for the orbits ``orbs``.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import SymmetricGroup
    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.util import _strip
    >>> S = SymmetricGroup(5)
    >>> S.schreier_sims()
    >>> g = Permutation([0, 2, 3, 1, 4])
    >>> _strip(g, S.base, S.basic_orbits, S.basic_transversals)
    (Permutation([0, 1, 2, 3, 4]), 5)

    Notes
    =====

    The algorithm is described in [1],pp.89-90. The reason for returning
    both the current state of the element being decomposed and the level
    at which the sifting ends is that they provide important information for
    the randomized version of the Schreier-Sims algorithm.

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of computational group theory"

    See Also
    ========

    sympy.combinatorics.perm_groups.PermutationGroup.schreier_sims

    sympy.combinatorics.perm_groups.PermutationGroup.schreier_sims_random

    """
    h = g
    base_len = len(base)
    for i in range(base_len):
        beta = h(base[i])
        if beta == base[i]:
            continue
        if beta not in orbs[i]:
            return h, i + 1
        u = transversals[i][beta]
        h = ~u*h
    return h, base_len + 1

def _distribute_gens_by_base(base, gens):
    """
    Distribute the group elements ``gens`` in basic stabilizers.

    Here, ``base`` is a sequence of points in `\{0, 1, ..., n-1\}`, and
    ``gens`` is a list of elements of a permutation group of degree `n`.
    Notice that for a base `(b_1, b_2, ..., b_k)`, the basic stabilizers
    are defined as `G^{(i)} = G_{b_1, ..., b_{i-1}}` for
    `i \in\{1, 2, ..., k\}`. The result is a list of length `k`, where `k` is
    the length of ``base``. The `i`-th entry contains those elements in
    ``gens`` which fix the first `i` elements of ``base`` (so that the
    `0`-th entry is equal to ``gens`` itself). If no element fixes the first
    `i` elements of ``base``, the `i`-th element is set to a list containing
    the identity element.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import DihedralGroup
    >>> from sympy.combinatorics.util import _distribute_gens_by_base
    >>> D = DihedralGroup(3)
    >>> D.schreier_sims()
    >>> D.strong_gens
    [Permutation([1, 2, 0]), Permutation([2, 1, 0]), Permutation([0, 2, 1])]
    >>> D.base
    [0, 1]
    >>> _distribute_gens_by_base(D.base, D.strong_gens)
    [[Permutation([1, 2, 0]), Permutation([2, 1, 0]), Permutation([0, 2, 1])],\
    [Permutation([0, 2, 1])]]

    See Also
    ========

    _strong_gens_from_distr, _orbits_transversals_from_bsgs,
    _handle_precomputed_bsgs

    """
    base_len = len(base)
    stabs = []
    degree = gens[0].size
    for i in xrange(base_len):
        stabs.append([])
    num_gens = len(gens)
    max_stab_index = 0
    for i in xrange(num_gens):
        j = 0
        while j < base_len - 1 and gens[i](base[j]) == base[j]:
            j += 1
        if j > max_stab_index:
            max_stab_index = j
        for k in xrange(j + 1):
            stabs[k].append(gens[i])
    for i in range(max_stab_index + 1, base_len):
        stabs[i].append(_new_from_array_form(range(degree)))
    return stabs

def _strong_gens_from_distr(distr_gens):
    """
    Retrieve strong generating set from generators of basic stabilizers.

    This is just the union of the generators of the first and second basic
    stabilizers.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import SymmetricGroup
    >>> from sympy.combinatorics.util import (_strong_gens_from_distr,
    ... _distribute_gens_by_base)
    >>> S = SymmetricGroup(3)
    >>> S.schreier_sims()
    >>> S.strong_gens
    [Permutation([1, 2, 0]), Permutation([1, 0, 2]), Permutation([0, 2, 1])]
    >>> distr_gens = _distribute_gens_by_base(S.base, S.strong_gens)
    >>> _strong_gens_from_distr(distr_gens)
    [Permutation([1, 2, 0]), Permutation([1, 0, 2]), Permutation([0, 2, 1])]

    See Also
    ========

    _distribute_gens_by_base

    """
    if len(distr_gens) == 1:
        return distr_gens[0][:]
    else:
        result = distr_gens[0]
        for gen in distr_gens[1]:
            if gen not in result:
                result.append(gen)
        return result

def _orbits_transversals_from_bsgs(base, distr_gens,\
                                   transversals_only=False):
    """
    Compute basic orbits and transversals from a base and strong generating set.

    The generators are provided as distributed across the basic stabilizers.
    If the optional argument ``transversals_only`` is set to True, only the
    transversals are returned.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import SymmetricGroup
    >>> from sympy.combinatorics.util import _orbits_transversals_from_bsgs
    >>> from sympy.combinatorics.util import (_orbits_transversals_from_bsgs,
    ... _distribute_gens_by_base)
    >>> S = SymmetricGroup(3)
    >>> S.schreier_sims()
    >>> distr_gens = _distribute_gens_by_base(S.base, S.strong_gens)
    >>> _orbits_transversals_from_bsgs(S.base, distr_gens)
    ([[0, 1, 2], [1, 2]], [{0: Permutation([0, 1, 2]),\
    1: Permutation([1, 2, 0]), 2: Permutation([2, 0, 1])},\
    {1: Permutation([0, 1, 2]), 2: Permutation([0, 2, 1])}])

    See Also
    ========

    _distribute_gens_by_base, _handle_precomputed_bsgs

    """
    from sympy.combinatorics.perm_groups import PermutationGroup
    base_len = len(base)
    transversals = [None]*base_len
    if transversals_only is False:
        basic_orbits = [None]*base_len
    for i in xrange(base_len):
        group = PermutationGroup(distr_gens[i])
        transversals[i] = dict(group.orbit_transversal(base[i], pairs=True))
        if transversals_only is False:
            basic_orbits[i] = transversals[i].keys()
    if transversals_only:
        return transversals
    else:
        return basic_orbits, transversals

def _handle_precomputed_bsgs(base, strong_gens, transversals=None,\
                             basic_orbits=None, distr_gens=None):
    """
    Calculate BSGS-related structures from whatever are present.

    The base and strong generating set must be provided; if any of the
    transversals, basic orbits or distributed strong generators are not
    provided, they will be calculated from the base and strong generating set.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import DihedralGroup
    >>> from sympy.combinatorics.util import _handle_precomputed_bsgs
    >>> D = DihedralGroup(3)
    >>> D.schreier_sims()
    >>> _handle_precomputed_bsgs(D.base, D.strong_gens,
    ... basic_orbits=D.basic_orbits)
    ([{0: Permutation([0, 1, 2]), 1: Permutation([1, 2, 0]),\
    2: Permutation([2, 1, 0])}, {1: Permutation([0, 1, 2]),\
    2: Permutation([0, 2, 1])}], [[0, 1, 2], [1, 2]], [[Permutation([1, 2, 0]),\
    Permutation([2, 1, 0]), Permutation([0, 2, 1])], [Permutation([0, 2, 1])]])

    See Also
    ========

    _orbits_transversals_from_bsgs, distribute_gens_by_base

    """
    if distr_gens is None:
        distr_gens = _distribute_gens_by_base(base, strong_gens)
    if transversals is None:
        if basic_orbits is None:
            basic_orbits, transversals =\
            _orbits_transversals_from_bsgs(base, distr_gens)
        else:
            transversals =\
            _orbits_transversals_from_bsgs(base, distr_gens,
                                           transversals_only=True)
    else:
        if basic_orbits is None:
            base_len = len(base)
            basic_orbits = [None]*base_len
            for i in xrange(base_len):
                basic_orbits[i] = transversals[i].keys()
    return transversals, basic_orbits, distr_gens

def _base_ordering(base, degree):
    r"""
    Order `\{0, 1, ..., n-1\}` so that base points come first and in order.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import SymmetricGroup
    >>> from sympy.combinatorics.util import _base_ordering
    >>> S = SymmetricGroup(4)
    >>> S.schreier_sims()
    >>> _base_ordering(S.base, S.degree)
    [0, 1, 2, 3]

    Notes
    =====

    This is used in backtrack searches, when we define a relation `<<` on
    the underlying set for a permutation group of degree `n`,
    `\{0, 1, ..., n-1\}`, so that if `(b_1, b_2, ..., b_k)` is a base we
    have `b_i << b_j` whenever `i<j` and `b_i << a` for all
    `i\in\{1,2, ..., k\}` and `a` is not in the base. The idea is developed
    and applied to backtracking algorithms in [1], pp.108-132. The points
    that are not in the base are taken in increasing order.

    See Also
    ========

    sympy.combinatorics.perm_groups.PermutationGroup.list_lex_by_base

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of computational group theory"

    """
    base_len = len(base)
    ordering = [0]*degree
    for i in xrange(base_len):
        ordering[base[i]] = i
    current = base_len
    for i in xrange(degree):
        if i not in base:
            ordering[i] = current
            current += 1
    return ordering

def _verify_bsgs(group, base, gens):
    """
    Verify the correctness of a base and strong generating set.

    This is a naive implementation using the definition of a base and a strong
    generating set relative to it. There are other procedures for
    verifying a base and strong generating set, but this one will
    serve for more robust testing.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import AlternatingGroup
    >>> from sympy.combinatorics.util import _verify_bsgs
    >>> A = AlternatingGroup(4)
    >>> A.schreier_sims()
    >>> _verify_bsgs(A, A.base, A.strong_gens)
    True

    See Also
    ========

    sympy.combinatorics.perm_groups.PermutationGroup.schreier_sims

    """
    from sympy.combinatorics.perm_groups import PermutationGroup
    distr_gens = _distribute_gens_by_base(base, gens)
    base_len = len(base)
    degree = group.degree
    current_stabilizer = group
    for i in range(base_len):
        candidate = PermutationGroup(distr_gens[i])
        if current_stabilizer.order() != candidate.order():
            return False
        current_stabilizer = current_stabilizer.stabilizer(base[i])
    if current_stabilizer.order() != 1:
        return False
    return True

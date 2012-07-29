from sympy.ntheory import isprime, sieve
from sympy.combinatorics.permutations import _new_from_array_form

############################################
###
### Utilities for computational group theory
###
############################################


def _base_ordering(base, degree):
    r"""
    Order `\{0, 1, ..., n-1\}` so that base points come first and in order.

    Parameters
    ==========

    ``base`` - the base
    ``degree`` - the degree of the associated permutation group

    Returns
    =======

    A list ``base_ordering`` such that ``base_ordering[point]`` is the
    number of ``point`` in the ordering.
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

def _check_cycles_alt_sym(perm):
    """
    Checks for cycles of prime length p with n/2 < p < n-2.

    Here `n` is the degree of the permutation. This is a helper function for
    the function is_alt_sym from sympy.combinatorics.perm_groups.

    Examples
    ========

    >>> from sympy.combinatorics.util import _check_cycles_alt_sym
    >>> from sympy.combinatorics.permutations import Permutation
    >>> a = Permutation([[0,1,2,3,4,5,6,7,8,9,10], [11, 12]])
    >>> _check_cycles_alt_sym(a)
    False
    >>> b = Permutation([[0,1,2,3,4,5,6], [7,8,9,10]])
    >>> _check_cycles_alt_sym(b)
    True

    See Also
    ========

    sympy.combinatorics.perm_groups.PermutationGroup.is_alt_sym

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

def _distribute_gens_by_base(base, gens):
    """
    Distribute the group elements ``gens`` by membership in basic stabilizers.

    Notice that for a base `(b_1, b_2, ..., b_k)`, the basic stabilizers
    are defined as `G^{(i)} = G_{b_1, ..., b_{i-1}}` for
    `i \in\{1, 2, ..., k\}`.

    Parameters
    ==========

    ``base`` - a sequence of points in `\{0, 1, ..., n-1\}`
    ``gens`` - a list of elements of a permutation group of degree `n`.

    Returns
    =======

    List of length `k`, where `k` is
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

def _handle_precomputed_bsgs(base, strong_gens, transversals=None,\
                             basic_orbits=None, distr_gens=None):
    """
    Calculate BSGS-related structures from those present.

    The base and strong generating set must be provided; if any of the
    transversals, basic orbits or distributed strong generators are not
    provided, they will be calculated from the base and strong generating set.

    Parameters
    ==========

    ``base`` - the base
    ``strong_gens`` - the strong generators
    ``transversals`` - basic transversals
    ``basic_orbits`` - basic orbits
    ``distr_gens`` - strong generators distributed by membership in basic
    stabilizers

    Returns
    =======

    ``(transversals, basic_orbits, distr_gens)`` where ``transversals`` are the
    basic transversals, ``basic_orbits`` are the basic orbits, and
    ``distr_gens`` are the strong generators distributed by membership in basic
    stabilizers.

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

def _insert_point_in_base(group, base, strong_gens, pos, point, distr_gens=None, basic_orbits=None, transversals=None):
    from sympy.combinatorics.perm_groups import PermutationGroup
    # initialize basic group properties and BSGS structures
    base_len = len(base)
    degree = group.degree
    identity = _new_from_array_form(range(degree))
    transversals, basic_orbits, distr_gens = _handle_precomputed_bsgs(base, strong_gens, randomized=False, transversals=transversals, basic_orbits=basic_orbits, distr_gens=distr_gens)
    # cut the base at position pos and append the new point
    partial_base = base[: pos + 1]
    partial_base.append(point)
    # cut the generators for the stabilizer chain and amend them accordingly
    if pos == base_len - 1:
        partial_distr_gens = distr_gens[: pos + 1]
        partial_distr_gens.append([identity])
    else:
        partial_distr_gens = distr_gens[: pos + 2]
    # cut the basic orbits and transversals and amend them accordingly
    partial_basic_orbits = basic_orbits[: pos + 1]
    partial_transversals = transversals[: pos + 1]
    last_stab = PermutationGroup(partial_distr_gens[pos + 1])
    last_transversal = dict(last_stab.orbit_transversal(point, pairs=True))
    last_orbit = last_transversal.keys()
    partial_basic_orbits.append(last_orbit)
    partial_transversals.append(last_transversal)
    # baseswap with the partial BSGS structures. Notice that we need only
    # the orbit and transversal of the new point under the last stabilizer
    new_base, new_strong_gens = group.baseswap(partial_base, strong_gens, pos, randomized=False, transversals=partial_transversals, basic_orbits=partial_basic_orbits, distr_gens=partial_distr_gens)
    # amend the basic orbits and transversals
    stab_pos = PermutationGroup(distr_gens[pos])
    new_transversal = dict(stab_pos.orbit_transversal(point, pairs=True))
    transversals[pos] = new_transversal
    basic_orbits[pos] = new_transversal.keys()
    # amend the distributed generators if necessary
    if pos != base_len - 1:
        new_stab_gens = []
        for gen in new_strong_gens:
            if [gen(point) for point in new_base[: pos + 1]] == [point for point in new_base[: pos + 1]]:
                new_stab_gens.append(gen)
        distr_gens[pos + 1] = new_stab_gens
    # return the new partial base and partial strong generating set
    new_base.pop()
    new_base = new_base + base[pos + 1 :]
    while len(base) != 0:
        base.pop()
    for point in new_base:
        base.append(point)
    while len(strong_gens) != 0:
        strong_gens.pop()
    for gen in new_strong_gens:
        strong_gens.append(gen)

def _orbits_transversals_from_bsgs(base, distr_gens,\
                                   transversals_only=False):
    """
    Compute basic orbits and transversals from a base and strong generating set.

    The generators are provided as distributed across the basic stabilizers.
    If the optional argument ``transversals_only`` is set to True, only the
    transversals are returned.

    Parameters
    ==========

    ``base`` - the base
    ``distr_gens`` - strong generators distributed by membership in basic
    stabilizers
    ``transversals_only`` - a flag swithing between returning only the
    transversals/ both orbits and transversals

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

def _remove_gens(base, strong_gens, basic_orbits=None, distr_gens=None):
    from sympy.combinatorics.perm_groups import PermutationGroup
    base_len = len(base)
    degree = strong_gens[0].size
    identity = _new_from_array_form(range(degree))
    if distr_gens is None:
        distr_gens = _distribute_gens_by_base(base, strong_gens)
    temp = distr_gens[:]
    if basic_orbits is None:
        basic_orbits = []
        for i in range(base_len):
            stab = PermutationGroup(distr_gens[i])
            basic_orbit = stab.orbit(base[i])
            basic_orbits.append(basic_orbit)
    distr_gens.append([])
    res = strong_gens[:]
    for i in range(base_len - 1, -1, -1):
        gens_copy = distr_gens[i][:]
        for gen in distr_gens[i]:
            if gen not in distr_gens[i + 1]:
                temp_gens = gens_copy[:]
                temp_gens.remove(gen)
                if temp_gens == []:
                    continue
                temp_group = PermutationGroup(temp_gens)
                temp_orbit = temp_group.orbit(base[i])
                if temp_orbit == basic_orbits[i]:
                    gens_copy.remove(gen)
                    res.remove(gen)
    return res

def _strip(g, base, orbs, transversals):
    """
    Attempt to decompose a permutation using a (possibly partial) BSGS
    structure.

    This is done by treating the sequence ``base`` as an actual base, and
    the orbits ``orbs`` and transversals ``transversals`` as basic orbits and
    transversals relative to it.
    This process is called "sifting". A sift is unsuccessful when a certain
    orbit element is not found or when after the sift the decomposition
    doesn't end with the identity element.
    The argument ``transversals`` is a list of dictionaries that provides
    transversal elements for the orbits ``orbs``.

    Parameters
    ==========

    ``g`` - permutation to be decomposed
    ``base`` - sequence of points
    ``orbs`` - a list in which the ``i``-th entry is an orbit of ``base[i]``
    under some subgroup of the pointwise stabilizer of `
    `base[0], base[1], ..., base[i - 1]``. The groups themselves are implicit
    in this function since the only infromation we need is encoded in the orbits
    and transversals
    ``transversals`` - a list of orbit transversals associated with the orbits
    ``orbs``.

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

def _strong_gens_from_distr(distr_gens):
    """
    Retrieve strong generating set from generators of basic stabilizers.

    This is just the union of the generators of the first and second basic
    stabilizers.

    Parameters
    ==========

    ``distr_gens`` - strong generators distributed by membership in basic
    stabilizers

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

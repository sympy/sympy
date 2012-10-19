from sympy.combinatorics.util import _distribute_gens_by_base
from sympy.combinatorics import Permutation

rmul = Permutation.rmul

def _cmp_perm_lists(first, second):
    """
    Compare two lists of permutations as sets.

    This is used for testing purposes. Since the array form of a
    permutation is currently a list, Permutation is not hashable
    and cannot be put into a set.

    Examples
    ========

    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.testutil import _cmp_perm_lists
    >>> a = Permutation([0, 2, 3, 4, 1])
    >>> b = Permutation([1, 2, 0, 4, 3])
    >>> c = Permutation([3, 4, 0, 1, 2])
    >>> ls1 = [a, b, c]
    >>> ls2 = [b, c, a]
    >>> _cmp_perm_lists(ls1, ls2)
    True

    """
    return set([tuple(a) for a in first]) == \
           set([tuple(a) for a in second])

def _naive_list_centralizer(self, other):
    from sympy.combinatorics.perm_groups import PermutationGroup
    """
    Return a list of elements for the centralizer of a subgroup/set/element.

    This is a brute-force implementation that goes over all elements of the
    group and checks for membership in the centralizer. It is used to
    test ``.centralizer()`` from ``sympy.combinatorics.perm_groups``.

    Examples
    ========
    >>> from sympy.combinatorics.testutil import _naive_list_centralizer
    >>> from sympy.combinatorics.named_groups import DihedralGroup
    >>> D = DihedralGroup(4)
    >>> _naive_list_centralizer(D, D)
    [Permutation([0, 1, 2, 3]), Permutation([2, 3, 0, 1])]

    See Also
    ========

    sympy.combinatorics.perm_groups.centralizer

    """
    if hasattr(other, 'generators'):
        elements = list(self.generate_dimino())
        gens = other.generators
        commutes_with_gens = lambda x: [rmul(x, gen) for gen in gens] ==\
                                       [rmul(gen, x) for gen in gens]
        centralizer_list = []
        for element in elements:
            if commutes_with_gens(element):
                centralizer_list.append(element)
        return centralizer_list
    elif hasattr(other, 'getitem'):
        return _naive_list_centralizer(self, PermutationGroup(other))
    elif hasattr(other, 'array_form'):
        return _naive_list_centralizer(self, PermutationGroup([other]))

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
    >>> from sympy.combinatorics.testutil import _verify_bsgs
    >>> A = AlternatingGroup(4)
    >>> A.schreier_sims()
    >>> _verify_bsgs(A, A.base, A.strong_gens)
    True

    See Also
    ========

    sympy.combinatorics.perm_groups.PermutationGroup.schreier_sims

    """
    from sympy.combinatorics.perm_groups import PermutationGroup
    strong_gens_distr = _distribute_gens_by_base(base, gens)
    current_stabilizer = group
    for i in range(len(base)):
        candidate = PermutationGroup(strong_gens_distr[i])
        if current_stabilizer.order() != candidate.order():
            return False
        current_stabilizer = current_stabilizer.stabilizer(base[i])
    if current_stabilizer.order() != 1:
        return False
    return True

def _verify_centralizer(group, arg, centr=None):
    """
    Verify the centralizer of a group/set/element inside another group.

    This is used for testing ``.centralizer()`` from
    ``sympy.combinatorics.perm_groups``

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import (SymmetricGroup,
    ... AlternatingGroup)
    >>> from sympy.combinatorics.perm_groups import PermutationGroup
    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.testutil import _verify_centralizer
    >>> S = SymmetricGroup(5)
    >>> A = AlternatingGroup(5)
    >>> centr = PermutationGroup([Permutation([0, 1, 2, 3, 4])])
    >>> _verify_centralizer(S, A, centr)
    True

    See Also
    ========

    _naive_list_centralizer,\
    sympy.combinatorics.perm_groups.PermutationGroup.centralizer,\
    _cmp_perm_lists

    """
    if centr is None:
        centr = group.centralizer(arg)
    centr_list = list(centr.generate_dimino())
    centr_list_naive = _naive_list_centralizer(group, arg)
    return _cmp_perm_lists(centr_list, centr_list_naive)

def _verify_normal_closure(group, arg, closure=None):
    from sympy.combinatorics.perm_groups import PermutationGroup
    """
    Verify the normal closure of a subgroup/subset/element in a group.

    This is used to test
    sympy.combinatorics.perm_groups.PermutationGroup.normal_closure

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import (SymmetricGroup,\
    ... AlternatingGroup)
    >>> from sympy.combinatorics.testutil import _verify_normal_closure
    >>> S = SymmetricGroup(3)
    >>> A = AlternatingGroup(3)
    >>> _verify_normal_closure(S, A, closure=A)
    True

    See Also
    ========

    sympy.combinatorics.perm_groups.PermutationGroup.normal_closure

    """
    if closure is None:
        closure = group.normal_closure(arg)
    conjugates = []
    group_els = list(group.generate_dimino())
    if hasattr(arg, 'generators'):
        subgr_gens = arg.generators
    elif hasattr(arg, '__getitem__'):
        subgr_gens = arg
    elif hasattr(arg, 'array_form'):
        subgr_gens = [arg]
    for el in group_els:
        for gen in subgr_gens:
            conjugate = rmul(~el, gen, el)
            if conjugate not in conjugates:
                conjugates.append(conjugate)
    naive_closure = PermutationGroup(conjugates)
    return closure.is_subgroup(naive_closure)

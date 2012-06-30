from sympy.combinatorics.perm_groups import PermutationGroup, DirectProduct
from sympy.combinatorics.permutations import Permutation, _new_from_array_form

def SymmetricGroup(n):
    """
    Generates the symmetric group on ``n`` elements as a permutation group.

    The generators taken are the ``n``-cycle
    ``(0 1 2 ... n-1)`` and the transposition ``(0 1)`` (in cycle notation).
    (See [1]). After the group is generated, some of its basic properties
    are set.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import SymmetricGroup
    >>> G = SymmetricGroup(4)
    >>> G.order()
    24
    >>> list(G.generate_schreier_sims(af=True))
    [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2],
    [0, 3, 2, 1], [1, 2, 3, 0], [1, 2, 0, 3], [1, 3, 2, 0],
    [1, 3, 0, 2], [1, 0, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [2, 3, 1, 0],
    [2, 0, 3, 1], [2, 0, 1, 3], [2, 1, 3, 0], [2, 1, 0, 3], [3, 0, 1, 2],
    [3, 0, 2, 1], [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0]]

    See Also
    ========

    CyclicGroup, DihedralGroup, AlternatingGroup

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Symmetric_group#Generators_and_relations

    """
    if n == 1:
        G = PermutationGroup([Permutation([0])])
    elif n == 2:
        G = PermutationGroup([Permutation([1, 0])])
    else:
        a = range(1,n)
        a.append(0)
        gen1 = _new_from_array_form(a)
        a = range(n)
        a[0], a[1] = a[1], a[0]
        gen2 = _new_from_array_form(a)
        G = PermutationGroup([gen1, gen2])

    if n<3:
        G._is_abelian = True
    else:
        G._is_abelian = False
    G._degree = n
    G._is_transitive = True
    G._is_sym = True
    return G

def CyclicGroup(n):
    """
    Generates the cyclic group of order ``n`` as a permutation group.

    The generator taken is the ``n``-cycle ``(0 1 2 ... n-1)``
    (in cycle notation). After the group is generated, some of its basic
    properties are set.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import CyclicGroup
    >>> G = CyclicGroup(6)
    >>> G.order()
    6
    >>> list(G.generate_schreier_sims(af=True))
    [[0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 0], [2, 3, 4, 5, 0, 1],
    [3, 4, 5, 0, 1, 2], [4, 5, 0, 1, 2, 3], [5, 0, 1, 2, 3, 4]]

    See Also
    ========

    SymmetricGroup, DihedralGroup, AlternatingGroup

    """
    a = range(1, n)
    a.append(0)
    gen = _new_from_array_form(a)
    G = PermutationGroup([gen])

    G._is_abelian = True
    G._degree = n
    G._is_transitive = True
    G._order = n
    return G

def DihedralGroup(n):
    r"""
    Generates the dihedral group `D_n` as a permutation group.

    The dihedral group `D_n` is the group of symmetries of the regular
    ``n``-gon. The generators taken are the ``n``-cycle ``a = (0 1 2 ... n-1)``
    (a rotation of the ``n``-gon) and ``b = (0 n-1)(1 n-2)...``
    (a reflection of the ``n``-gon) in cycle rotation. It is easy to see that
    these satisfy ``a**n = b**2 = 1`` and ``bab = ~a`` so they indeed generate
    `D_n` (See [1]). After the group is generated, some of its basic properties
    are set.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import DihedralGroup
    >>> G = DihedralGroup(5)
    >>> a = list(G.generate_dimino())
    >>> [perm.cyclic_form for perm in a]
    [[[4], [3], [2], [1], [0]], [[0, 1, 2, 3, 4]], [[0, 2, 4, 1, 3]],
    [[0, 3, 1, 4, 2]], [[0, 4, 3, 2, 1]], [[2], [1, 3], [0, 4]],
    [[2, 3], [1, 4], [0]], [[3], [2, 4], [0, 1]], [[3, 4], [1], [0, 2]],
    [[4], [1, 2], [0, 3]]]

    See Also
    ========

    SymmetricGroup, CyclicGroup, AlternatingGroup

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Dihedral_group

    """
    # small cases are special
    if n == 1:
        return PermutationGroup([Permutation([1, 0])])
    if n == 2:
        return PermutationGroup([Permutation([1, 0, 3, 2]),
               Permutation([2, 3, 0, 1]), Permutation([3, 2, 1, 0])])

    a = range(1, n)
    a.append(0)
    gen1 = _new_from_array_form(a)
    a = range(n)
    a.reverse()
    gen2 = _new_from_array_form(a)
    G = PermutationGroup([gen1, gen2])

    G._is_abelian = False
    G._degree = n
    G._is_transitive = True
    G._order = 2*n
    return G

def AlternatingGroup(n):
    """
    Generates the alternating group on ``n`` elements as a permutation group.

    For ``n > 2``, the generators taken are ``(0 1 2), (0 1 2 ... n-1)`` for
    ``n`` odd
    and ``(0 1 2), (1 2 ... n-1)`` for ``n`` even (See [1], p.31, ex.6.9.).
    After the group is generated, some of its basic properties are set.
    The cases ``n = 1, 2`` are handled separately.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import AlternatingGroup
    >>> G = AlternatingGroup(4)
    >>> a = list(G.generate_dimino())
    >>> len(a)
    12
    >>> [perm.is_even for perm in a]
    [True, True, True, True, True, True, True, True, True, True, True, True]

    See Also
    ========

    SymmetricGroup, CyclicGroup, DihedralGroup

    References
    ==========

    [1] Armstrong, M. "Groups and Symmetry"

    """
    # small cases are special
    if n == 1 or n == 2:
        return PermutationGroup([Permutation([0])])

    a = range(n)
    a[0], a[1], a[2] = a[1], a[2], a[0]
    gen1 = _new_from_array_form(a)
    if n % 2 == 1:
        a = range(1, n)
        a.append(0)
        gen2 = _new_from_array_form(a)
    else:
        a = range(2, n)
        a.append(1)
        gen2 = _new_from_array_form([0] + a)
    G = PermutationGroup([gen1, gen2])

    if n<4:
        G._is_abelian = True
    else:
        G._is_abelian = False
    G._degree = n
    G._is_transitive = True
    G._is_alt = True
    return G

def AbelianGroup(*cyclic_orders):
    """
    Returns the direct product of cyclic groups with the given orders.

    According to the structure theorem for finite abelian groups ([1]),
    every finite abelian group can be written as the direct product of finitely
    many cyclic groups.
    [1] http://groupprops.subwiki.org/wiki/Structure_theorem_for_finitely
    _generated_abelian_groups


    Examples
    ========


    >>> from sympy.combinatorics.perm_groups import AbelianGroup
    >>> AbelianGroup(3,4)
    PermutationGroup([Permutation([1, 2, 0, 3, 4, 5, 6]),\
    Permutation([0, 1, 2, 4, 5, 6, 3])])

    See Also
    ========
    DirectProduct
    """
    groups = []
    degree = 0
    order = 1
    for size in cyclic_orders:
        degree += size
        order *= size
        groups.append(CyclicGroup(size))
    G = DirectProduct(*groups)
    G._is_abelian = True
    G._degree = degree
    G._order = order

    return G

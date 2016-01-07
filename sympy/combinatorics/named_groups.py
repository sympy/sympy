from __future__ import print_function, division

from sympy.core.compatibility import range
from sympy.core.basic import Basic
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.group_constructs import DirectProduct
from sympy.combinatorics.permutations import Permutation
from sympy import factorial

_af_new = Permutation._af_new


def AbelianGroup(*cyclic_orders):
    """
    Returns the direct product of cyclic groups with the given orders.

    According to the structure theorem for finite abelian groups ([1]),
    every finite abelian group can be written as the direct product of
    finitely many cyclic groups.

    Examples
    ========

    >>> from sympy.combinatorics import Permutation
    >>> Permutation.print_cyclic = True
    >>> from sympy.combinatorics.named_groups import AbelianGroup
    >>> AbelianGroup(3, 4)
    PermutationGroup([
            (6)(0 1 2),
            (3 4 5 6)])
    >>> _.is_group
    True

    See Also
    ========

    DirectProduct

    References
    ==========

    [1] http://groupprops.subwiki.org/wiki/Structure_theorem_for_finitely_generated_abelian_groups

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
    >>> G.is_group
    True
    >>> a = list(G.generate_dimino())
    >>> len(a)
    12
    >>> all(perm.is_even for perm in a)
    True

    See Also
    ========

    SymmetricGroup, CyclicGroup, DihedralGroup

    References
    ==========

    [1] Armstrong, M. "Groups and Symmetry"

    """
    # small cases are special
    if n in (1, 2):
        return PermutationGroup([Permutation([0])])

    a = list(range(n))
    a[0], a[1], a[2] = a[1], a[2], a[0]
    gen1 = a
    if n % 2:
        a = list(range(1, n))
        a.append(0)
        gen2 = a
    else:
        a = list(range(2, n))
        a.append(1)
        a.insert(0, 0)
        gen2 = a
    gens = [gen1, gen2]
    if gen1 == gen2:
        gens = gens[:1]
    G = PermutationGroup([_af_new(a) for a in gens], dups=False)

    if n < 4:
        G._is_abelian = True
        G._is_nilpotent = True
    else:
        G._is_abelian = False
        G._is_nilpotent = False
    if n < 5:
        G._is_solvable = True
    else:
        G._is_solvable = False
    G._degree = n
    G._is_transitive = True
    G._is_alt = True
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
    >>> G.is_group
    True
    >>> G.order()
    6
    >>> list(G.generate_schreier_sims(af=True))
    [[0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 0], [2, 3, 4, 5, 0, 1],
    [3, 4, 5, 0, 1, 2], [4, 5, 0, 1, 2, 3], [5, 0, 1, 2, 3, 4]]

    See Also
    ========

    SymmetricGroup, DihedralGroup, AlternatingGroup

    """
    a = list(range(1, n))
    a.append(0)
    gen = _af_new(a)
    G = PermutationGroup([gen])

    G._is_abelian = True
    G._is_nilpotent = True
    G._is_solvable = True
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
    >>> G.is_group
    True
    >>> a = list(G.generate_dimino())
    >>> [perm.cyclic_form for perm in a]
    [[], [[0, 1, 2, 3, 4]], [[0, 2, 4, 1, 3]],
    [[0, 3, 1, 4, 2]], [[0, 4, 3, 2, 1]], [[0, 4], [1, 3]],
    [[1, 4], [2, 3]], [[0, 1], [2, 4]], [[0, 2], [3, 4]],
    [[0, 3], [1, 2]]]

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

    a = list(range(1, n))
    a.append(0)
    gen1 = _af_new(a)
    a = list(range(n))
    a.reverse()
    gen2 = _af_new(a)
    G = PermutationGroup([gen1, gen2])
    # if n is a power of 2, group is nilpotent
    if n & (n-1) == 0:
        G._is_nilpotent = True
    else:
        G._is_nilpotent = False
    G._is_abelian = False
    G._is_solvable = True
    G._degree = n
    G._is_transitive = True
    G._order = 2*n
    return G


class SymmetricGroup(Basic):
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
    >>> G.is_group
    True
    >>> G.order()
    24
    >>> list(G.generate_schreier_sims(af=True))
    [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 1, 2, 0], [0, 2, 3, 1],
    [1, 3, 0, 2], [2, 0, 1, 3], [3, 2, 0, 1], [0, 3, 1, 2], [1, 0, 2, 3],
    [2, 1, 3, 0], [3, 0, 1, 2], [0, 1, 3, 2], [1, 2, 0, 3], [2, 3, 1, 0],
    [3, 1, 0, 2], [0, 2, 1, 3], [1, 3, 2, 0], [2, 0, 3, 1], [3, 2, 1, 0],
    [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3], [3, 0, 2, 1]]

    See Also
    ========

    CyclicGroup, DihedralGroup, AlternatingGroup

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Symmetric_group#Generators_and_relations

    """
    is_group = True

    def __new__(cls, degree):
        if degree <= 0 or not isinstance(degree, int):
            raise ValueError("Degree of a SymmetricGroup should be a "
                    "positive integer")

        obj = Basic.__new__(cls, degree)

        if degree < 3:
            obj._is_abelian = True
            obj._is_nilpotent = True
        else:
            obj._is_abelian = False
            obj._is_nilpotent = False

        if degree < 5:
            obj._is_solvable = True
        else:
            obj._is_solvable = False

        obj._is_transitive = True
        obj._is_sym = True
        obj._degree = degree

        return obj

    def __getitem__(self, i):
        return self.generators[i]

    def __contains__(self, i):
        if not isinstance(i, Permutation):
            raise TypeError("A SymmetricGroup contains only Permutations as "
                            "elements, not elements of type %s" % type(i))
        return self.contains(i)

    def __len__(self):
        return len(self.generators)

    def __mul__(self, other):
        pass

    def center(self):
        if self.degree == 1 or self.degree == 2:
            return self
        else:
            return SymmetricGroup(1)

    @property
    def generators(self):
        degree = self.degree
        if degree == 1:
            return [Permutation([0])]
        elif degree == 2:
            return [Permutation([1, 0])]
        else:
            a = list(range(1, degree))
            a.append(0)
            gen1 = _af_new(a)
            a = list(range(degree))
            a[0], a[1] = a[1], a[0]
            gen2 = _af_new(a)
            return [gen1, gen2]

    @property
    def degree(self):
        return self._degree

    def is_subgroup(self, other):
        if isinstance(other, SymmetricGroup):
            return self.degree <= other.degree

    def contains(self, g):
        if not isinstance(g, Permutation):
            return False
        return g.size <= self.degree

    @property
    def is_abelian(self):
        return self._is_abelian

    @property
    def is_trivial(self):
        return self == SymmetricGroup(1)

    def order(self):
        return factorial(self.degree)

    def orbit(self, alpha):
        return set(list(range(self.degree)))

    def stabilizer(self, alpha):
        if alpha >= self.degree or alpha < 0:
            raise ValueError("Such element is not there in the finite set on "
                    "which SymmetricGroup is defined")
        gens = _stabilizer(self.degree, alpha)
        if gens == self.generators:
            return SymmetricGroup(self.degree)
        return PermutationGroup(gens)


def _stabilizer(degree, alpha):
    if degree == 1 or degree == 2:
        return [Permutation(degree - 1)]
    gen2_self = Permutation(degree - 1)(0, 1)
    a = list(range(degree))
    a.remove(alpha)
    gen1 = Permutation(*a)
    if alpha^gen2_self == alpha:
        gen2 = gen2_self
    else:
        if alpha == 0:
            gen2 = Permutation(degree - 1)(2, 1)
        elif alpha == 1:
            gen2 = Permutation(degree - 1)(2, 0)
    return [gen1, gen2]


def RubikGroup(n):
    """Return a group of Rubik's cube generators

    >>> from sympy.combinatorics.named_groups import RubikGroup
    >>> RubikGroup(2).is_group
    True
    """
    from sympy.combinatorics.generators import rubik
    if n <= 1:
        raise ValueError("Invalid cube . n has to be greater than 1")
    return PermutationGroup(rubik(n))


# setting an alias for SymmetricGroup
Sym = SymmetricGroup

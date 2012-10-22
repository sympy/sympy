from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.permutations import cyclic as perm_cyclic
from sympy.utilities.iterables import variations, rotate_left

def symmetric(n):
    """
    Generates the symmetric group of order n, Sn.

    Examples
    ========

    >>> from sympy.combinatorics.generators import symmetric
    >>> list(symmetric(3))
    [Permutation([0, 1, 2]), Permutation([0, 2, 1]), Permutation([1, 0, 2]), \
    Permutation([1, 2, 0]), Permutation([2, 0, 1]), Permutation([2, 1, 0])]
    """
    for perm in variations(range(n), n):
        yield Permutation(perm)

def cyclic(n):
    """
    Generates the cyclic group of order n, Cn.

    Examples
    ========

    >>> from sympy.combinatorics.generators import cyclic
    >>> list(cyclic(5))
    [Permutation([0, 1, 2, 3, 4]), Permutation([1, 2, 3, 4, 0]), \
    Permutation([2, 3, 4, 0, 1]), Permutation([3, 4, 0, 1, 2]), \
    Permutation([4, 0, 1, 2, 3])]

    See Also
    ========
    dihedral
    """
    gen = range(n)
    for i in xrange(n):
        yield Permutation(gen)
        gen = rotate_left(gen, 1)

def alternating(n):
    """
    Generates the alternating group of order n, An.

    Examples
    ========

    >>> from sympy.combinatorics.generators import alternating
    >>> list(alternating(4))
    [Permutation([0, 1, 2, 3]), Permutation([0, 2, 3, 1]), \
    Permutation([0, 3, 1, 2]), Permutation([1, 0, 3, 2]), \
    Permutation([1, 2, 0, 3]), Permutation([1, 3, 2, 0]), \
    Permutation([2, 0, 1, 3]), Permutation([2, 1, 3, 0]), \
    Permutation([2, 3, 0, 1]), Permutation([3, 0, 2, 1]), \
    Permutation([3, 1, 0, 2]), Permutation([3, 2, 1, 0])]
    """
    for perm in variations(range(n), n):
        p = Permutation(perm)
        if p.is_even:
            yield p

def dihedral(n):
    """
    Generates the dihedral group of order 2n, Dn.

    The result is given as a subgroup of Sn, except for the special cases n=1
    (the group S2) and n=2 (the Klein 4-group) where that's not possible
    and embeddings in S2 and S4 respectively are given.

    Examples
    ========

    >>> from sympy.combinatorics.generators import dihedral
    >>> list(dihedral(4))
    [Permutation([0, 1, 2, 3]), Permutation([3, 2, 1, 0]), \
    Permutation([1, 2, 3, 0]), Permutation([0, 3, 2, 1]), \
    Permutation([2, 3, 0, 1]), Permutation([1, 0, 3, 2]), \
    Permutation([3, 0, 1, 2]), Permutation([2, 1, 0, 3])]

    See Also
    ========
    cyclic
    """
    if n == 1:
        yield Permutation([0, 1])
        yield Permutation([1, 0])
    elif n == 2:
        yield Permutation([0, 1, 2, 3])
        yield Permutation([1, 0, 3, 2])
        yield Permutation([2, 3, 0, 1])
        yield Permutation([3, 2, 1, 0])
    else:
        gen = range(n)
        for i in xrange(n):
            yield Permutation(gen)
            yield Permutation(gen[::-1])
            gen = rotate_left(gen, 1)

def rubik_cube_generators():
    """return the generators of the Rubik cube, see
    http://www.gap-system.org/Doc/Examples/rubik.html
    """
    a = [[(1,3,8,6),(2,5,7,4),(9,33,25,17),(10,34,26,18),(11,35,27,19)],
      [(9,11,16,14),(10,13,15,12),(1,17,41,40),(4,20,44,37),(6,22,46,35)],
      [(17,19,24,22),(18,21,23,20),(6,25,43,16),(7,28,42,13),(8,30,41,11)],
      [(25,27,32,30),(26,29,31,28),(3,38,43,19),(5,36,45,21),(8,33,48,24)],
      [(33,35,40,38),(34,37,39,36),(3,9,46,32),(2,12,47,29),(1,14,48,27)],
      [(41,43,48,46),(42,45,47,44),(14,22,30,38),(15,23,31,39),(16,24,32,40)]]
    return [Permutation(perm_cyclic(x, 48)) for x in a]

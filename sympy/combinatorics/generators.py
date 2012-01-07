from sympy.combinatorics.permutations import Permutation
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
    Generates the dihedral group of order 2n, D2n.

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
    gen = range(n)
    for i in xrange(n):
        yield Permutation(gen)
        yield Permutation(gen[::-1])
        gen = rotate_left(gen, 1)

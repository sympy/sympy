"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed SymPy-wide and are located
here for easy import.
"""

from collections import defaultdict

# These are in here because telling if something is an iterable just by calling
# hasattr(obj, "__iter__") behaves differently in Python 2 and Python 3.  In
# particular, hasattr(str, "__iter__") is False in Python 2 and True in Python 3.
# I think putting them here also makes it easier to use them in the core.

def iterable(i, exclude=(basestring, dict)):
    """
    Return a boolean indicating whether i is an iterable in the sympy sense.

    When sympy is working with iterables, it is almost always assuming
    that the iterable is not a string or a mapping, so those are excluded
    by default. If you want a pure python definition, make exclude=None. To
    exclude multiple items, pass them as a tuple.

    See also: is_sequence

    Examples
    ========

    >>> from sympy.utilities.iterables import iterable
    >>> from sympy import Tuple
    >>> things = [[1], (1,), set([1]), Tuple(1), (j for j in [1, 2]), {1:2}, '1', 1]
    >>> for i in things:
    ...     print iterable(i), type(i)
    True <... 'list'>
    True <... 'tuple'>
    True <... 'set'>
    True <class 'sympy.core.containers.Tuple'>
    True <... 'generator'>
    False <... 'dict'>
    False <... 'str'>
    False <... 'int'>

    >>> iterable({}, exclude=None)
    True
    >>> iterable({}, exclude=str)
    True
    >>> iterable("no", exclude=str)
    False

    """
    try:
        iter(i)
    except TypeError:
        return False
    if exclude:
        return not isinstance(i, exclude)
    return True

def is_sequence(i, include=None):
    """
    Return a boolean indicating whether i is a sequence in the sympy
    sense. If anything that fails the test below should be included as
    being a sequence for your application, set 'include' to that object's
    type; multiple types should be passed as a tuple of types.

    Note: although generators can generate a sequence, they often need special
    handling to make sure their elements are captured before the generator is
    exhausted, so these are not included by default in the definition of a
    sequence.

    See also: iterable

    Examples
    ========

    >>> from sympy.utilities.iterables import is_sequence
    >>> from types import GeneratorType
    >>> is_sequence([])
    True
    >>> is_sequence(set())
    False
    >>> is_sequence('abc')
    False
    >>> is_sequence('abc', include=str)
    True
    >>> generator = (c for c in 'abc')
    >>> is_sequence(generator)
    False
    >>> is_sequence(generator, include=(str, GeneratorType))
    True

    """
    return (hasattr(i, '__getitem__') and
            iterable(i) or
            bool(include) and
            isinstance(i, include))

"""
Wrapping some imports in try/except statements to allow the same code to
be used in Python 3+ as well.
"""

try:
    callable = callable
except NameError:
    import collections
    def callable(obj):
        return isinstance(obj, collections.Callable)

try:
    from functools import reduce
except ImportError:
    reduce = reduce

def cmp_to_key(mycmp):
    """
    Convert a cmp= function into a key= function

    This code is included in Python 2.7 and 3.2 in functools.
    """
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

try:
    import __builtin__
    cmp = __builtin__.cmp
except AttributeError:
    def cmp(a,b):
        return (a > b) - (a < b)

try:
    from itertools import product
except ImportError: # Python 2.5
    def product(*args, **kwds):
        """
        Cartesian product of input iterables.

        Equivalent to nested for-loops in a generator expression. For example,
        product(A, B) returns the same as ((x,y) for x in A for y in B).

        The nested loops cycle like an odometer with the rightmost element
        advancing on every iteration. This pattern creates a lexicographic
        ordering so that if the input's iterables are sorted, the product
        tuples are emitted in sorted order.

        To compute the product of an iterable with itself, specify the number
        of repetitions with the optional repeat keyword argument. For example,
        product(A, repeat=4) means the same as product(A, A, A, A).

        Examples
        ========

        >>> from sympy.core.compatibility import product
        >>> [''.join(p) for p in list(product('ABC', 'xy'))]
        ['Ax', 'Ay', 'Bx', 'By', 'Cx', 'Cy']
        >>> list(product(range(2), repeat=2))
        [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)

try:
    from itertools import permutations
except ImportError: # Python 2.5
    def permutations(iterable, r=None):
        """
        Return successive r length permutations of elements in the iterable.

        If r is not specified or is None, then r defaults to the length of
        the iterable and all possible full-length permutations are generated.

        Permutations are emitted in lexicographic sort order. So, if the input
        iterable is sorted, the permutation tuples will be produced in sorted
        order.

        Elements are treated as unique based on their position, not on their
        value. So if the input elements are unique, there will be no repeat
        values in each permutation.

        Examples;
        >>> from sympy.core.compatibility import permutations
        >>> [''.join(p) for p in list(permutations('ABC', 2))]
        ['AB', 'AC', 'BA', 'BC', 'CA', 'CB']
        >>> list(permutations(range(3)))
        [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        """

        pool = tuple(iterable)
        n = len(pool)
        r = n if r is None else r
        if r > n:
            return
        indices = range(n)
        cycles = range(n, n-r, -1)
        yield tuple(pool[i] for i in indices[:r])
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i+1:] + indices[i:i+1]
                    cycles[i] = n - i
                else:
                    j = cycles[i]
                    indices[i], indices[-j] = indices[-j], indices[i]
                    yield tuple(pool[i] for i in indices[:r])
                    break
            else:
                return

try:
    from itertools import combinations, combinations_with_replacement
except ImportError: # < python 2.6
    def combinations(iterable, r):
        """
        Return r length subsequences of elements from the input iterable.

        Combinations are emitted in lexicographic sort order. So, if the
        input iterable is sorted, the combination tuples will be produced
        in sorted order.

        Elements are treated as unique based on their position, not on their
        value. So if the input elements are unique, there will be no repeat
        values in each combination.

        See also: combinations_with_replacement

        Examples
        ========

        >>> from sympy.core.compatibility import combinations
        >>> list(combinations('ABC', 2))
        [('A', 'B'), ('A', 'C'), ('B', 'C')]
        >>> list(combinations(range(4), 3))
        [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        """
        pool = tuple(iterable)
        n = len(pool)
        if r > n:
            return
        indices = range(r)
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i+1, r):
                indices[j] = indices[j-1] + 1
            yield tuple(pool[i] for i in indices)

    def combinations_with_replacement(iterable, r):
        """Return r length subsequences of elements from the input iterable
        allowing individual elements to be repeated more than once.

        Combinations are emitted in lexicographic sort order. So, if the
        input iterable is sorted, the combination tuples will be produced
        in sorted order.

        Elements are treated as unique based on their position, not on their
        value. So if the input elements are unique, the generated combinations
        will also be unique.

        See also: combinations

        Examples
        ========

        >>> from sympy.core.compatibility import combinations_with_replacement
        >>> list(combinations_with_replacement('AB', 2))
        [('A', 'A'), ('A', 'B'), ('B', 'B')]
        """
        pool = tuple(iterable)
        n = len(pool)
        if not n and r:
            return
        indices = [0] * r
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != n - 1:
                    break
            else:
                return
            indices[i:] = [indices[i] + 1] * (r - i)
            yield tuple(pool[i] for i in indices)

def set_intersection(*sets):
    """Return the intersection of all the given sets.

    As of Python 2.6 you can write ``set.intersection(*sets)``.

    Examples
    ========

    >>> from sympy.core.compatibility import set_intersection
    >>> set_intersection(set([1, 2]), set([2, 3]))
    set([2])
    >>> set_intersection()
    set()
    """
    if not sets:
        return set()
    rv = sets[0]
    for s in sets:
        rv &= s
    return rv

def set_union(*sets):
    """Return the union of all the given sets.

    As of Python 2.6 you can write ``set.union(*sets)``.

    >>> from sympy.core.compatibility import set_union
    >>> set_union(set([1, 2]), set([2, 3]))
    set([1, 2, 3])
    >>> set_union()
    set()
    """
    rv = set()
    for s in sets:
        rv |= s
    return rv

try:
    bin = bin
except NameError: # Python 2.5
    def bin(x):
        """
        bin(number) -> string

        Stringifies an int or long in base 2.
        """
        if x < 0: return '-' + bin(-x)
        out = []
        if x == 0: out.append('0')
        while x > 0:
            out.append('01'[x & 1])
            x >>= 1
            pass
        return '0b' + ''.join(reversed(out))

try:
    next = next
except NameError: # Python 2.5
    def next(*args):
        """
        next(iterator[, default])

        Return the next item from the iterator. If default is given and the
        iterator is exhausted, it is returned instead of raising StopIteration.
        """
        if len(args) == 1:
            return args[0].next()
        elif len(args) == 2:
            try:
                return args[0].next()
            except StopIteration:
                return args[1]
        else:
            raise TypeError('Expected 1 or 2 arguments, got %s' % len(args))

def as_int(n):
    """
    Convert the argument to a builtin integer.

    The return value is guaranteed to be equal to the input. ValueError is
    raised if the input has a non-integral value.

    Examples
    ========

    >>> from sympy.core.compatibility import as_int
    >>> from sympy import sqrt
    >>> 3.0
    3.0
    >>> as_int(3.0) # convert to int and test for equality
    3
    >>> int(sqrt(10))
    3
    >>> as_int(sqrt(10))
    Traceback (most recent call last):
    ...
    ValueError: ... is not an integer

    """
    result = int(n)
    if result != n:
        raise ValueError('%s is not an integer' % n)
    return result

def quick_sort(seq, quick=True):
    """Sort by hash and break ties with default_sort_key (default)
    or entirely by default_sort_key if ``quick`` is False.

    When sorting for consistency between systems, ``quick`` should be
    False; if sorting is just needed to give consistent orderings during
    a given session ``quick`` can be True.

    >>> from sympy.core.compatibility import quick_sort
    >>> from sympy.abc import x

    For PYTHONHASHSEED=3923375334 the x came first; for
    PYTHONHASHSEED=158315900 the x came last (on a 32-bit system).

    >>> quick_sort([x, 1, 3]) in [(1, 3, x), (x, 1, 3)]
    True
    """
    from sympy.utilities.iterables import default_sort_key

    if not quick:
        seq = list(seq)
        seq.sort(key=default_sort_key)
    else:
        d = defaultdict(list)
        for a in seq:
            d[hash(a)].append(a)
        seq = []
        for k in sorted(d.keys()):
          if len(d[k]) > 1:
              seq.extend(sorted(d[k], key=default_sort_key))
          else:
              seq.extend(d[k])
    return tuple(seq)

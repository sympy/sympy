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
    Return a boolean indicating whether ``i`` is SymPy iterable.

    When SymPy is working with iterables, it is almost always assuming
    that the iterable is not a string or a mapping, so those are excluded
    by default. If you want a pure Python definition, make exclude=None. To
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
    Return a boolean indicating whether ``i`` is a sequence in the SymPy
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
    def cmp(a, b):
        return (a > b) - (a < b)

try:
    from itertools import product
except ImportError:  # Python 2.5
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
            result = [x + [y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)

try:
    from itertools import permutations
except ImportError:  # Python 2.5
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
        cycles = range(n, n - r, -1)
        yield tuple(pool[i] for i in indices[:r])
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i + 1:] + indices[i:i + 1]
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
except ImportError:  # < python 2.6
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
            for j in range(i + 1, r):
                indices[j] = indices[j - 1] + 1
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
except NameError:  # Python 2.5
    def bin(x):
        """
        bin(number) -> string

        Stringifies an int or long in base 2.
        """
        if x < 0:
            return '-' + bin(-x)
        out = []
        if x == 0:
            out.append('0')
        while x > 0:
            out.append('01'[x & 1])
            x >>= 1
            pass
        return '0b' + ''.join(reversed(out))

try:
    next = next
except NameError:  # Python 2.5
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

try:
    from __builtin__ import bin
except ImportError:  # Python 2.5
    _hexDict = {
        '0': '0000', '1': '0001', '2': '0010', '3': '0011', '4': '0100', '5': '0101',
        '6': '0110', '7': '0111', '8': '1000', '9': '1001', 'a': '1010', 'b': '1011',
        'c': '1100', 'd': '1101', 'e': '1110', 'f': '1111', 'L': ''}

    def bin(n):
        """Return the equivalent to Python 2.6's bin function.

        Examples
        ========

        >>> from sympy.core.compatibility import bin
        >>> bin(-123)
        '-0b1111011'
        >>> bin(0) # this is the only time a 0 will be to the right of 'b'
        '0b0'

        See Also
        ========
        sympy.physics.quantum.shor.arr

        Modified from http://code.activestate.com/recipes/576847/
        """
        # =========================================================
        # create hex of int, remove '0x'. now for each hex char,
        # look up binary string, append in list and join at the end.
        # =========================================================
        if n < 0:
            return '-%s' % bin(-n)
        return '0b%s' % (''.join([_hexDict[hstr] for hstr in hex(n)[2:].lower()
            ]).lstrip('0') or '0')


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


def default_sort_key(item, order=None):
    """Return a key that can be used for sorting.

    The key has the structure:

    (class_key, (len(args), args), exponent.sort_key(), coefficient)

    This key is supplied by the sort_key routine of Basic objects when
    ``item`` is a Basic object or an object (other than a string) that
    sympifies to a Basic object. Otherwise, this function produces the
    key.

    The ``order`` argument is passed along to the sort_key routine and is
    used to determe how the terms *within* an expression are ordered.
    (See examples below) ``order`` options are: 'lex', 'grlex', 'grevlex',
    and reversed values of the same (e.g. 'rev-lex'). The default order
    value is None (which translates to 'lex').

    Examples
    ========

    >>> from sympy import Basic, S, I, default_sort_key
    >>> from sympy.core.function import UndefinedFunction
    >>> from sympy.abc import x

    The following are eqivalent ways of getting the key for an object:

    >>> x.sort_key() == default_sort_key(x)
    True

    Here are some examples of the key that is produced:

    >>> default_sort_key(UndefinedFunction('f'))
    ((0, 0, 'UndefinedFunction'), (1, ('f',)), ((1, 0, 'Number'),
        (0, ()), (), 1), 1)
    >>> default_sort_key('1')
    ((0, 0, 'str'), (1, ('1',)), ((1, 0, 'Number'), (0, ()), (), 1), 1)
    >>> default_sort_key(S.One)
    ((1, 0, 'Number'), (0, ()), (), 1)
    >>> default_sort_key(2)
    ((1, 0, 'Number'), (0, ()), (), 2)


    While sort_key is a method only defined for SymPy objects,
    default_sort_key will accept anything as an argument so it is
    more robust as a sorting key. For the following, using key=
    lambda i: i.sort_key() would fail because 2 doesn't have a sort_key
    method; that's why default_sort_key is used. Note, that it also
    handles sympification of non-string items likes ints:

    >>> a = [2, I, -I]
    >>> sorted(a, key=default_sort_key)
    [2, -I, I]

    The returned key can be used anywhere that a key can be specified for
    a function, e.g. sort, min, max, etc...:

    >>> a.sort(key=default_sort_key); a[0]
    2
    >>> min(a, key=default_sort_key)
    2

    Note
    ----

    The key returned is useful for getting items into a canonical order
    that will be the same across platforms. It is not directly useful for
    sorting lists of expressions:

    >>> a, b = x, 1/x

    Since ``a`` has only 1 term, its value of sort_key is unaffected by
    ``order``:

    >>> a.sort_key() == a.sort_key('rev-lex')
    True

    If ``a`` and ``b`` are combined then the key will differ because there
    are terms that can be ordered:

    >>> eq = a + b
    >>> eq.sort_key() == eq.sort_key('rev-lex')
    False
    >>> eq.as_ordered_terms()
    [x, 1/x]
    >>> eq.as_ordered_terms('rev-lex')
    [1/x, x]

    But since the keys for each of these terms are independent of ``order``'s
    value, they don't sort differently when they appear separately in a list:

    >>> sorted(eq.args, key=default_sort_key)
    [1/x, x]
    >>> sorted(eq.args, key=lambda i: default_sort_key(i, order='rev-lex'))
    [1/x, x]

    The order of terms obtained when using these keys is the order that would
    be obtained if those terms were *factors* in a product.

    See Also
    ========

    sympy.core.expr.as_ordered_factors, sympy.core.expr.as_ordered_terms

    """

    from sympy.core import S, Basic
    from sympy.core.sympify import sympify, SympifyError
    from sympy.core.compatibility import iterable

    if isinstance(item, Basic):
        return item.sort_key(order=order)

    if iterable(item, exclude=basestring):
        if isinstance(item, dict):
            args = item.items()
            unordered = True
        elif isinstance(item, set):
            args = item
            unordered = True
        else:
            # e.g. tuple, list
            args = list(item)
            unordered = False

        args = [default_sort_key(arg, order=order) for arg in args]

        if unordered:
            # e.g. dict, set
            args = sorted(args)

        cls_index, args = 10, (len(args), tuple(args))
    else:
        if not isinstance(item, basestring):
            try:
                item = sympify(item)
            except SympifyError:
                # e.g. lambda x: x
                pass
            else:
                if isinstance(item, Basic):
                    # e.g int -> Integer
                    return default_sort_key(item)
                # e.g. UndefinedFunction

        # e.g. str
        cls_index, args = 0, (1, (str(item),))

    return (cls_index, 0, item.__class__.__name__
            ), args, S.One.sort_key(), S.One


# define keys for use with lazyDSU_sort; define as tuples
# so they can't be inadvertently changed

small_first_keys = (lambda _: _.count_ops(), default_sort_key)
big_first_keys = (lambda _: -(_.count_ops()), default_sort_key)


def lazyDSU_sort(seq, keys, warn=False):
    """Return sorted seq, breaking ties by applying keys only when needed.

    If ``warn`` is True then an error will be raised if there were no
    keys remaining to break ties.

    Examples
    ========

    >>> from sympy.utilities.iterables import lazyDSU_sort, default_sort_key
    >>> from sympy.abc import x, y

    The count_ops is not sufficient to break ties in this list and the first
    two items appear in their original order: Python sorting is stable.

    >>> lazyDSU_sort([y + 2, x + 2, x**2 + y + 3],
    ...    [lambda x: x.count_ops()], warn=False)
    ...
    [y + 2, x + 2, x**2 + y + 3]

    The use of default_sort_key allows the tie to be broken (and warn can
    be safely left False).

    >>> lazyDSU_sort([y + 2, x + 2, x**2 + y + 3],
    ...    [lambda x: x.count_ops(),
    ...     default_sort_key])
    ...
    [x + 2, y + 2, x**2 + y + 3]

    Here, sequences are sorted by length, then sum:

    >>> seq, keys = [[[1, 2, 1], [0, 3, 1], [1, 1, 3], [2], [1]], [
    ...    lambda x: len(x),
    ...    lambda x: sum(x)]]
    ...
    >>> lazyDSU_sort(seq, keys, warn=False)
    [[1], [2], [1, 2, 1], [0, 3, 1], [1, 1, 3]]

    If ``warn`` is True, an error will be raised if there were not
    enough keys to break ties:

    >>> lazyDSU_sort(seq, keys, warn=True)
    Traceback (most recent call last):
    ...
    ValueError: not enough keys to break ties


    Notes
    =====

    The decorated sort is one of the fastest ways to sort a sequence for
    which special item comparison is desired: the sequence is decorated,
    sorted on the basis of the decoration (e.g. making all letters lower
    case) and then undecorated. If one wants to break ties for items that
    have the same decorated value, a second key can be used. But if the
    second key is expensive to compute then it is inefficient to decorate
    all items with both keys: only those items having identical first key
    values need to be decorated. This function applies keys successively
    only when needed to break ties.
    """
    d = defaultdict(list)
    keys = list(keys)
    f = keys.pop(0)
    for a in seq:
        d[f(a)].append(a)
    seq = []
    for k in sorted(d.keys()):
        if len(d[k]) > 1:
            if keys:
                d[k] = lazyDSU_sort(d[k], keys, warn)
            elif warn:
                raise ValueError('not enough keys to break ties')
            seq.extend(d[k])
        else:
            seq.append(d[k][0])
    return seq


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

class _Link(object):
    __slots__ = 'right', 'left', 'key', '__weakref__'

class oset(object):
    """
    Return a set with the elements in the same order that they were added.

    Examples
    ========

    >>> from sympy.core.compatibility import oset
    >>> oset([2, 1, 3, 1, 2])
    oset([2, 1, 3])

    Notes
    =====
    Big-O running times for all methods are the same as for regular sets.
    The internal self.__map dictionary maps keys to links in a doubly linked list.
    The circular doubly linked list starts and ends with a sentinel element.
    The sentinel element never gets deleted (this simplifies the algorithm).
    The right/left links are weakref proxies (to prevent circular references).
    Individual links are kept alive by the hard reference in self.__map.
    Those hard references disappear when a key is deleted from an oset.

    References
    ==========
    1. modified from http://code.activestate.com/recipes/576696/ (r5)

    """
    def __init__(self, dat=None):
        self.__root = root = _Link()         # sentinel node for doubly linked list
        root.left = root.right = root
        self.__map = {}                     # key --> link
        if dat is not None:
            self.update(dat)

    def __len__(self):
        """
        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 1]); a
        oset([1, 2])
        >>> len(a)
        2

        """
        return len(self.__map)

    def __contains__(self, key):
        """
        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 1]); a
        oset([1, 2])
        >>> 2 in a
        True

        """
        return key in self.__map

    def add(self, key):
        """
        Store new key in a new link at the end of the linked list.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset()
        >>> a.add(1); a
        oset([1])
        >>> a.add(2); a
        oset([1, 2])

        """
        from weakref import proxy

        if key not in self.__map:
            self.__map[key] = link = _Link()
            root = self.__root
            last = root.left
            link.left, link.right, link.key = last, root, key
            last.right = root.left = proxy(link)

    def discard(self, key):
        """
        Remove an existing item.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3])
        >>> a.discard(2); a
        oset([1, 3])

        """
        # use self.__map to find the link which is then removed
        # by updating the links in the predecessor and successors.
        if key in self.__map:
            link = self.__map.pop(key)
            link.left.right = link.right
            link.right.left = link.left

    def __iter__(self):
        """
        Traverse the linked list in order.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2])
        >>> list(a)
        [1, 2]

        """
        root = self.__root
        curr = root.right
        while curr is not root:
            yield curr.key
            curr = curr.right

    def __reversed__(self):
        """
        Return an iterator that returns the linked list in reversed order.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2])
        >>> list(reversed(a))
        [2, 1]

        """
        root = self.__root
        curr = root.left
        while curr is not root:
            yield curr.key
            curr = curr.left

    def pop(self, last=True):
        """
        Return the last element of the linked list; to return the first
        element, set ``last`` to False.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3])
        >>> a.pop()
        3
        >>> a.pop(0)
        1

        """
        if not self:
            raise KeyError('set is empty')
        key = next(reversed(self)) if last else next(iter(self))
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        """
        Return True if ``other`` is the same as self.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3])
        >>> b = oset([1, 2, 3])
        >>> a == b
        True
        >>> a == list(b)
        False
        >>> b.add(4)
        >>> a == b
        False

        """
        if not isinstance(other, oset):
            return False
        if len(self) != len(other):
            return False
        for s, o in zip(self, other):
            if s != o:
                return False
        return True

    def update(self, other):
        """
        Add all elements of other to self in place; to actually create a
        union of self and other use self.union(other).

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3])
        >>> a.update([1, 4, 7]); a
        oset([1, 2, 3, 4, 7])

        See Also
        ========
        union
        """
        for o in other:
            self.add(o)

    def intersection(self, other):
        """
        Return the items of self that are also in other. s.intersection(o)
        will have the same items as o.intersection(s) but perhaps in different
        order.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3])
        >>> b = oset([3, 1])
        >>> a.intersection(b)
        oset([1, 3])
        >>> b.intersection(a)
        oset([3, 1])

        """
        if not isinstance(other, oset):
            other = oset(other)
        rv = oset()
        for s in self:
            if s in other:
                rv.add(s)
        return rv

    def difference(self, other):
        """
        Return the items of self that are not in other, in the same order
        that they were in self.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3, 4])
        >>> b = oset([3, 1])
        >>> a.difference(b)
        oset([2, 4])
        >>> b.difference(a)
        oset()

        """
        if not isinstance(other, oset):
            other = oset(other)
        rv = oset()
        for s in self:
            if s not in other:
                rv.add(s)
        return rv

    def union(self, other):
        """
        Return the items that are in both self and other. other.union(self)
        will contain the same items but perhaps in a different order.

        Examples
        ========

        >>> from sympy.core.compatibility import oset
        >>> a = oset([1, 2, 3, 4])
        >>> b = oset([3, 1])
        >>> a.union(b) == a
        True
        >>> c = b.union(a)
        >>> c == a
        False
        >>> c
        oset([3, 1, 2, 4])

        See Also
        ========
        update
        """
        rv = oset()
        rv.update(self)
        rv.update(other)
        return rv

"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed SymPy-wide and are located
here for easy import.
"""

from collections import defaultdict
from sympy.external import import_module


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
    def product(*args, **kwargs):
        """
        Cartesian product of input iterables.

        Equivalent to nested for-loops in a generator expression. For example,
        cartes(A, B) returns the same as ((x,y) for x in A for y in B).

        The nested loops cycle like an odometer with the rightmost element
        advancing on every iteration. This pattern creates a lexicographic
        ordering so that if the input's iterables are sorted, the product
        tuples are emitted in sorted order.

        To compute the product of an iterable with itself, specify the number
        of repetitions with the optional repeat keyword argument. For example,
        product(A, repeat=4) means the same as product(A, A, A, A).

        Examples
        ========

        >>> from sympy.utilities.iterables import cartes
        >>> [''.join(p) for p in list(cartes('ABC', 'xy'))]
        ['Ax', 'Ay', 'Bx', 'By', 'Cx', 'Cy']
        >>> list(cartes(range(2), repeat=2))
        [(0, 0), (0, 1), (1, 0), (1, 1)]

        See Also
        ========
        variations
        """
        pools = map(tuple, args) * kwargs.get('repeat', 1)
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
    from collections import namedtuple
except ImportError:
    # code from http://code.activestate.com/recipes/500261-named-tuples/
    # PSF license
    # code is Copyright 2007-2013 Raymond Hettinger
    from operator import itemgetter as _itemgetter
    from keyword import iskeyword as _iskeyword
    import sys as _sys

    # For some reason, doctest will test namedtuple's docstring if we simply
    # have a def namedtuple() here
    def _namedtuple(typename, field_names, verbose=False, rename=False):
        if isinstance(field_names, basestring):
            field_names = field_names.replace(',', ' ').split()
        field_names = tuple(map(str, field_names))
        if rename:
            names = list(field_names)
            seen = set()
            for i, name in enumerate(names):
                if (not min(c.isalnum() or c=='_' for c in name) or _iskeyword(name)
                    or not name or name[0].isdigit() or name.startswith('_')
                    or name in seen):
                        names[i] = '_%d' % i
                seen.add(name)
            field_names = tuple(names)
        for name in (typename,) + field_names:
            if not min(c.isalnum() or c=='_' for c in name):
                raise ValueError('Type names and field names can only contain'
                                 ' alphanumeric characters and underscores: %r' % name)
            if _iskeyword(name):
                raise ValueError('Type names and field names cannot be a'
                                 ' keyword: %r' % name)
            if name[0].isdigit():
                raise ValueError('Type names and field names cannot start'
                                 ' with a number: %r' % name)
        seen_names = set()
        for name in field_names:
            if name.startswith('_') and not rename:
                raise ValueError('Field names cannot start with an'
                                 ' underscore: %r' % name)
            if name in seen_names:
                raise ValueError('Encountered duplicate field name: %r' % name)
            seen_names.add(name)

        # Create and fill-in the class template
        numfields = len(field_names)
        argtxt = repr(field_names).replace("'", "")[1:-1]
        reprtxt = ', '.join('%s=%%r' % name for name in field_names)
        template = '''class %(typename)s(tuple):
    '%(typename)s(%(argtxt)s)' \n
    __slots__ = () \n
    _fields = %(field_names)r \n
    def __new__(_cls, %(argtxt)s):
        return _tuple.__new__(_cls, (%(argtxt)s)) \n
    @classmethod
    def _make(cls, iterable, new=tuple.__new__, len=len):
        'Make a new %(typename)s object from a sequence or iterable'
        result = new(cls, iterable)
        if len(result) != %(numfields)d:
            raise TypeError('Expected %(numfields)d arguments, got %%d' %% len(result))
        return result \n
    def __repr__(self):
        return '%(typename)s(%(reprtxt)s)' %% self \n
    def _asdict(self):
        'Return a new dict which maps field names to their values'
        return dict(zip(self._fields, self)) \n
    def _replace(_self, **kwds):
        'Return a new %(typename)s object replacing specified fields with new values'
        result = _self._make(map(kwds.pop, %(field_names)r, _self))
        if kwds:
            raise ValueError('Got unexpected field names: %%r' %% kwds.keys())
        return result \n
    def __getnewargs__(self):
        return tuple(self) \n\n''' % locals()
        for i, name in enumerate(field_names):
            template += '    %s = _property(_itemgetter(%d))\n' % (name, i)
        if verbose:
            print template

        # Execute the template string in a temporary namespace
        namespace = dict(_itemgetter=_itemgetter, __name__='namedtuple_%s' % typename,
                         _property=property, _tuple=tuple)
        try:
            exec template in namespace
        except SyntaxError, e:
            raise SyntaxError(e.message + ':\n' + template)
        result = namespace[typename]

        # For pickling to work, the __module__ variable needs to be set to the frame
        # where the named tuple is created.  Bypass this step in enviroments where
        # sys._getframe is not defined (Jython for example) or sys._getframe is not
        # defined for arguments greater than 0 (IronPython).
        try:
            result.__module__ = _sys._getframe(1).f_globals.get('__name__', '__main__')
        except (AttributeError, ValueError):
            pass

        return result

    namedtuple = _namedtuple

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
    try:
        result = int(n)
        if result != n:
            raise TypeError
    except TypeError:
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
    used to determine how the terms *within* an expression are ordered.
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


def _nodes(e):
    """
    A helper for ordered() which returns the node count of ``e`` which
    for Basic object is the number of Basic nodes in the expression tree
    but for other object is 1 (unless the object is an iterable or dict
    for which the sum of nodes is returned).
    """
    from basic import Basic

    if isinstance(e, Basic):
        return e.count(Basic)
    elif iterable(e):
        return 1 + sum(_nodes(ei) for ei in e)
    elif isinstance(e, dict):
        return 1 + sum(_nodes(k) + _nodes(v) for k, v in e.iteritems())
    else:
        return 1


def ordered(seq, keys=None, default=True, warn=False):
    """Return an iterator of the seq where keys are used to break ties.
    Two default keys will be applied after and provided unless ``default``
    is False. The two keys are _nodes and default_sort_key which will
    place smaller expressions before larger ones (in terms of Basic nodes)
    and where there are ties, they will be broken by the default_sort_key.

    If ``warn`` is True then an error will be raised if there were no
    keys remaining to break ties. This can be used if it was expected that
    there should be no ties.

    Examples
    ========

    >>> from sympy.utilities.iterables import ordered, default_sort_key
    >>> from sympy import count_ops
    >>> from sympy.abc import x, y

    The count_ops is not sufficient to break ties in this list and the first
    two items appear in their original order (i.e. the sorting is stable):

    >>> list(ordered([y + 2, x + 2, x**2 + y + 3],
    ...    count_ops, default=False, warn=False))
    ...
    [y + 2, x + 2, x**2 + y + 3]

    The default_sort_key allows the tie to be broken:

    >>> list(ordered([y + 2, x + 2, x**2 + y + 3]))
    ...
    [x + 2, y + 2, x**2 + y + 3]

    Here, sequences are sorted by length, then sum:

    >>> seq, keys = [[[1, 2, 1], [0, 3, 1], [1, 1, 3], [2], [1]], [
    ...    lambda x: len(x),
    ...    lambda x: sum(x)]]
    ...
    >>> list(ordered(seq, keys, default=False, warn=False))
    [[1], [2], [1, 2, 1], [0, 3, 1], [1, 1, 3]]

    If ``warn`` is True, an error will be raised if there were not
    enough keys to break ties:

    >>> list(ordered(seq, keys, default=False, warn=True))
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
    only when needed to break ties. By yielding an iterator, use of the
    tie-breaker is delayed as long as possible.

    This function is best used in cases when use of the first key is
    expected to be a good hashing function; if there are no unique hashes
    from application of a key then that key should not have been used. The
    exception, however, is that even if there are many collisions, if the
    first group is small and one does not need to process all items in the
    list then time will not be wasted sorting what one was not interested
    in. For example, if one were looking for the minimum in a list and
    there were several criteria used to define the sort order, then this
    function would be good at returning that quickly if the first group
    of candidates is small relative to the number of items being processed.

    """
    d = defaultdict(list)
    if keys:
        if not isinstance(keys, (list, tuple)):
            keys = [keys]
        keys = list(keys)

        f = keys.pop(0)
        for a in seq:
            d[f(a)].append(a)
    else:
        if not default:
            raise ValueError('if default=False then keys must be provided')
        d[None].extend(seq)

    for k in sorted(d.keys()):
        if len(d[k]) > 1:
            if keys:
                d[k] = ordered(d[k], keys, default, warn)
            elif default:
                d[k] = ordered(d[k], (_nodes, default_sort_key,),
                               default=False, warn=warn)
            elif warn:
                raise ValueError('not enough keys to break ties')
        for v in d[k]:
            yield v
        d.pop(k)

try:
    next = next
except NameError:
    def next(x):
        return x.next()

# If HAS_GMPY is 0, no supported version of gmpy is available. Otherwise,
# HAS_GMPY contains the major version number of gmpy; i.e. 1 for gmpy, and
# 2 for gmpy2.

# Versions of gmpy prior to 1.03 do not work correctly with int(largempz)
# For example, int(gmpy.mpz(2**256)) would raise OverflowError.
# See issue 1881.

# Minimum version of gmpy changed to 1.13 to allow a single code base to also
# work with gmpy2.

def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

GROUND_TYPES = _getenv('SYMPY_GROUND_TYPES', 'auto').lower()

HAS_GMPY = 0

if GROUND_TYPES != 'python':

    # Don't try to import gmpy2 if ground types is set to gmpy1. This is
    # primarily intended for testing.

    if GROUND_TYPES != 'gmpy1':
        gmpy = import_module('gmpy2', min_module_version='2.0.0',
            module_version_attr='version', module_version_attr_call_args=())
        if gmpy:
            HAS_GMPY = 2
    else:
        GROUND_TYPES = 'gmpy'

    if not HAS_GMPY:
        gmpy = import_module('gmpy', min_module_version='1.13',
            module_version_attr='version', module_version_attr_call_args=())
        if gmpy:
            HAS_GMPY = 1

if GROUND_TYPES == 'auto':
    if HAS_GMPY:
        GROUND_TYPES = 'gmpy'
    else:
        GROUND_TYPES = 'python'

if GROUND_TYPES == 'gmpy' and not HAS_GMPY:
    from warnings import warn
    warn("gmpy library is not installed, switching to 'python' ground types")
    GROUND_TYPES = 'python'

# SYMPY_INTS is a tuple containing the base types for valid integer types.

import sys

if sys.version_info[0] == 2:
    SYMPY_INTS = (int, long)
else:
    SYMPY_INTS = (int,)

if GROUND_TYPES == 'gmpy':
    SYMPY_INTS += (type(gmpy.mpz(0)),)

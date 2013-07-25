"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed SymPy-wide and are located
here for easy import.
"""

import operator
from collections import defaultdict
from sympy.external import import_module


"""
Python 2 and Python 3 compatable imports
"""

import sys
PY2 = sys.version_info[0] == 2

if not PY2:
    import collections

    class_types = type,
    integer_types = (int,)
    string_types = (str,)
    long = int

    # String / unicode compatibility
    unicode = str
    def u(x):
        return x

    Iterator = object

    # Moved definitions
    get_function_code = operator.attrgetter("__code__")
    get_function_globals = operator.attrgetter("__globals__")
    get_function_name = operator.attrgetter("__name__")

    import builtins
    def callable(obj):
        return isinstance(obj, collections.Callable)
    def cmp(a, b):
        return (a > b) - (a < b)
    filter = filter
    from functools import reduce
    from io import StringIO
    cStringIO = StringIO

    exec_ = getattr(builtins, "exec")
else:
    import codecs
    import types

    class_types = (type, types.ClassType)
    integer_types = (int, long)
    string_types = (str, unicode)
    long = long

    # String / unicode compatibility
    unicode = unicode
    def u(x):
        return codecs.unicode_escape_decode(x)[0]

    class Iterator(object):
        def next(self):
            return type(self).__next__(self)

    # Moved definitions
    get_function_code = operator.attrgetter("func_code")
    get_function_globals = operator.attrgetter("func_globals")
    get_function_name = operator.attrgetter("func_name")

    import __builtin__ as builtins
    callable  = callable
    cmp = cmp
    from itertools import ifilter as filter
    reduce = reduce
    from StringIO import StringIO
    from cStringIO import StringIO as cStringIO

    def exec_(_code_, _globs_=None, _locs_=None):
        """Execute code in a namespace."""
        if _globs_ is None:
            frame = sys._getframe(1)
            _globs_ = frame.f_globals
            if _locs_ is None:
                _locs_ = frame.f_locals
            del frame
        elif _locs_ is None:
            _locs_ = _globs_
        exec("exec _code_ in _globs_, _locs_")

def with_metaclass(meta, *bases):
    """Create a base class with a metaclass."""
    class metaclass(meta):
        __call__ = type.__call__
        __init__ = type.__init__
        def __new__(cls, name, this_bases, d):
            if this_bases is None:
                return type.__new__(cls, name, (), d)
            return meta(name, bases, d)
    return metaclass("NewBase", None, {})


# These are in here because telling if something is an iterable just by calling
# hasattr(obj, "__iter__") behaves differently in Python 2 and Python 3.  In
# particular, hasattr(str, "__iter__") is False in Python 2 and True in Python 3.
# I think putting them here also makes it easier to use them in the core.


def iterable(i, exclude=(string_types, dict)):
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
    ...     print('%s %s' % (iterable(i), type(i)))
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
    from itertools import combinations_with_replacement
except ImportError:  # <= Python 2.6
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

    >>> from sympy import S, I, default_sort_key
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

    if iterable(item, exclude=string_types):
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
        if not isinstance(item, string_types):
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
    from .basic import Basic

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

    >>> from sympy.utilities.iterables import ordered
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

# check_output() is new in Python 2.7
import os

try:
    from subprocess import CalledProcessError
    try:
        from subprocess import check_output
    except ImportError:
        from subprocess import check_call
        def check_output(*args, **kwargs):
            with open(os.devnull, 'w') as fh:
                kwargs['stdout'] = fh
                try:
                    return check_call(*args, **kwargs)
                except CalledProcessError as e:
                    e.output = ("program output is not available for Python 2.6.x")
                    raise e
except ImportError:
    # running on platform like App Engine, no subprocess at all
    pass

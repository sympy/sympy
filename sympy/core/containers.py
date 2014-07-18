"""Module for SymPy containers

    (SymPy objects that store other SymPy objects)

    The containers implemented in this module are subclassed to Basic.
    They are supposed to work seamlessly within the SymPy framework.
"""

from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.compatibility import as_int, integer_types
from sympy.core.sympify import sympify, converter
from sympy.utilities.iterables import iterable


class Tuple(Basic):
    """
    Wrapper around the builtin tuple object

    The Tuple is a subclass of Basic, so that it works well in the
    SymPy framework.  The wrapped tuple is available as self.args, but
    you can also access elements or slices with [:] syntax.

    >>> from sympy import symbols
    >>> from sympy.core.containers import Tuple
    >>> a, b, c, d = symbols('a b c d')
    >>> Tuple(a, b, c)[1:]
    (b, c)
    >>> Tuple(a, b, c).subs(a, d)
    (d, b, c)

    """

    def __new__(cls, *args, **assumptions):
        args = [ sympify(arg) for arg in args ]
        obj = Basic.__new__(cls, *args, **assumptions)
        return obj

    def __getitem__(self, i):
        if isinstance(i, slice):
            indices = i.indices(len(self))
            return Tuple(*[self.args[j] for j in range(*indices)])
        return self.args[i]

    def __len__(self):
        return len(self.args)

    def __contains__(self, item):
        return item in self.args

    def __iter__(self):
        return iter(self.args)

    def __add__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(self.args + other.args))
        elif isinstance(other, tuple):
            return Tuple(*(self.args + other))
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(other.args + self.args))
        elif isinstance(other, tuple):
            return Tuple(*(other + self.args))
        else:
            return NotImplemented

    def __mul__(self, other):
        try:
            n = as_int(other)
        except ValueError:
            raise TypeError("Can't multiply sequence by non-integer of type '%s'" % type(other))
        return self.func(*(self.args*n))

    __rmul__ = __mul__

    def __eq__(self, other):
        if isinstance(other, Basic):
            return super(Tuple, self).__eq__(other)
        return self.args == other

    def __ne__(self, other):
        if isinstance(other, Basic):
            return super(Tuple, self).__ne__(other)
        return self.args != other

    def __hash__(self):
        return hash(self.args)

    def _to_mpmath(self, prec):
        return tuple([a._to_mpmath(prec) for a in self.args])

    def __lt__(self, other):
        return sympify(self.args < other.args)

    def __le__(self, other):
        return sympify(self.args <= other.args)

    # XXX: Basic defines count() as something different, so we can't
    # redefine it here. Originally this lead to cse() test failure.
    def tuple_count(self, value):
        """T.count(value) -> integer -- return number of occurrences of value"""
        return self.args.count(value)

    def index(self, value, start=None, stop=None):
        """T.index(value, [start, [stop]]) -> integer -- return first index of value.
           Raises ValueError if the value is not present."""
        # XXX: One would expect:
        #
        # return self.args.index(value, start, stop)
        #
        # here. Any trouble with that? Yes:
        #
        # >>> (1,).index(1, None, None)
        # Traceback (most recent call last):
        #   File "<stdin>", line 1, in <module>
        # TypeError: slice indices must be integers or None or have an __index__ method
        #
        # See: http://bugs.python.org/issue13340

        if start is None and stop is None:
            return self.args.index(value)
        elif stop is None:
            return self.args.index(value, start)
        else:
            return self.args.index(value, start, stop)

converter[tuple] = lambda tup: Tuple(*tup)


def tuple_wrapper(method):
    """
    Decorator that converts any tuple in the function arguments into a Tuple.

    The motivation for this is to provide simple user interfaces.  The user can
    call a function with regular tuples in the argument, and the wrapper will
    convert them to Tuples before handing them to the function.

    >>> from sympy.core.containers import tuple_wrapper
    >>> def f(*args):
    ...    return args
    >>> g = tuple_wrapper(f)

    The decorated function g sees only the Tuple argument:

    >>> g(0, (1, 2), 3)
    (0, (1, 2), 3)

    """
    def wrap_tuples(*args, **kw_args):
        newargs = []
        for arg in args:
            if type(arg) is tuple:
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)
    return wrap_tuples


class Dict(Basic):
    """
    Wrapper around the builtin dict object

    The Dict is a subclass of Basic, so that it works well in the
    SymPy framework.  Because it is immutable, it may be included
    in sets, but its values must all be given at instantiation and
    cannot be changed afterwards.  Otherwise it behaves identically
    to the Python dict.

    >>> from sympy.core.containers import Dict

    >>> D = Dict({1: 'one', 2: 'two'})
    >>> for key in D:
    ...    if key == 1:
    ...        print('%s %s' % (key, D[key]))
    1 one

    The args are sympified so the 1 and 2 are Integers and the values
    are Symbols. Queries automatically sympify args so the following work:

    >>> 1 in D
    True
    >>> D.has('one') # searches keys and values
    True
    >>> 'one' in D # not in the keys
    False
    >>> D[1]
    one

    """

    def __new__(cls, *args):
        if len(args) == 1 and ((args[0].__class__ is dict) or
                             (args[0].__class__ is Dict)):
            items = [Tuple(k, v) for k, v in args[0].items()]
        elif iterable(args) and all(len(arg) == 2 for arg in args):
            items = [Tuple(k, v) for k, v in args]
        else:
            raise TypeError('Pass Dict args as Dict((k1, v1), ...) or Dict({k1: v1, ...})')
        elements = frozenset(items)
        obj = Basic.__new__(cls, elements)
        obj.elements = elements
        obj._dict = dict(items)  # In case Tuple decides it wants to sympify
        return obj

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        return self._dict[sympify(key)]

    def __setitem__(self, key, value):
        raise NotImplementedError("SymPy Dicts are Immutable")

    @property
    def args(self):
        return tuple(self.elements)

    def items(self):
        '''D.items() -> list of D's (key, value) pairs, as 2-tuples'''
        return self._dict.items()

    def keys(self):
        '''D.keys() -> list of D's keys'''
        return self._dict.keys()

    def values(self):
        '''D.values() -> list of D's values'''
        return self._dict.values()

    def __iter__(self):
        '''x.__iter__() <==> iter(x)'''
        return iter(self._dict)

    def __len__(self):
        '''x.__len__() <==> len(x)'''
        return self._dict.__len__()

    def get(self, key, default=None):
        '''D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.'''
        return self._dict.get(sympify(key), default)

    def __contains__(self, key):
        '''D.__contains__(k) -> True if D has a key k, else False'''
        return sympify(key) in self._dict

    def __lt__(self, other):
        return sympify(self.args < other.args)

    @property
    def _sorted_args(self):
        from sympy.utilities import default_sort_key
        return sorted(self.args, key=default_sort_key)


from itertools import imap, islice
import operator

class Stream(Basic):
    """
    Stream container

    Examples
    ========

    >>> from sympy.core.containers import Stream
    >>> from sympy.abc import x
    >>> def fib(x, y):
    ...     while True:
    ...         yield x
    ...         x, y = y, x + y
    >>> s = Stream(fib(0, 1))
    >>> list(s[0:5])
    [0, 1, 1, 2, 3]
    """

    __slots__ = ("_last", "_collection", "_origin")

    class _StreamIterator(object):

        __slots__ = ("_stream", "_position")

        def __init__(self, stream):
            self._stream = stream
            self._position = -1 # not started yet

        def __next__(self):
            self._position += 1
            if len(self._stream._collection) > self._position or self._stream._fill_to(self._position):
                return self._stream._collection[self._position]

            raise StopIteration()
        next = __next__

    def __init__(self, origin=[]):
        self._collection = []
        self._last = -1 # not started yet
        self._origin = iter(origin) if origin else []

    def _fill_to(self, index):
        while self._last < index:
            try:
                n = next(self._origin)
            except StopIteration:
                return False

            self._last += 1
            self._collection.append(n)

        return True

    def _eval_subs(self, old, new):
        return self.__class__(imap(lambda t: t.subs(old, new), self))

    def __iter__(self):
        return self._StreamIterator(self)

    def __getitem__(self, index):
        if isinstance(index, integer_types):
            if index < 0:
                raise TypeError("Invalid argument type")
            self._fill_to(index)
            return self._collection[index]
        elif isinstance(index, slice):
            if index.step == 0:
                raise ValueError("Step must not be 0")
            if not index.stop:
                return self.__class__(imap(self.__getitem__, islice(index.start, index.stop, index.step or 1)))
            return self.__class__(imap(self.__getitem__, xrange(index.start, index.stop, index.step or 1)))
        else:
            raise TypeError("Invalid argument type")

    def __add__(self, other):
        return self.__class__(imap(operator.add, self, other))

    def __sub__(self, other):
        from sympy.core import S
        return self + other*S.NegativeOne

    def __mul__(self, other):
        from sympy.core import Number, Symbol, Basic, S
        if not isinstance(other, Stream):
            return self.__class__(imap(lambda x: x*other, self))
        def mul():
            k = 0
            while True:
                p = S.Zero
                for i, a in enumerate(self):
                    p += a*other[k - i]
                    if i >= k:
                        break
                yield p
                k += 1
        return self.__class__(mul())

    def __truediv__(self, other):
        from sympy.core import S
        if not isinstance(other, Stream):
            return self*(S.One/other)
        def invert_unit_series(s):
            r = S.One
            for t in s:
                yield r
                r = S.NegativeOne*r*t
        c = S.One/other[0]
        return self*self.__class__(invert_unit_series(other*c))*c
    __div__ = __truediv__

    def shift(self, n):
        """ Right shift stream by n items """
        return self.__class__(islice(self, n, None))

    @property
    def is_empty(self):
        try:
            if self.__getitem__(0) is not None:
                return False
        except IndexError:
            return True

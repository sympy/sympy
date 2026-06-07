"""Module for SymPy containers

    (SymPy objects that store other SymPy objects)

    The containers implemented in this module are subclassed to Basic.
    They are supposed to work seamlessly within the SymPy framework.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import MutableSet
from typing import Callable, Iterable, Iterator, TypeVar, ParamSpec, overload, Generic, TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.core.expr import Expr
    from sympy.logic.boolalg import Boolean

from .basic import Basic
from .sorting import default_sort_key, ordered
from .sympify import _sympify, sympify, _sympy_converter, SympifyError
from sympy.core.kind import Kind
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import as_int

_P = ParamSpec('_P')
_R = TypeVar('_R')
_T = TypeVar('_T', bound=Basic)
_Texpr = TypeVar('_Texpr', bound='Expr')


class Tuple(Basic, Generic[_T]):
    """
    Wrapper around the builtin tuple object.

    Explanation
    ===========

    The Tuple is a subclass of Basic, so that it works well in the
    SymPy framework.  The wrapped tuple is available as self.args, but
    you can also access elements or slices with [:] syntax.

    Parameters
    ==========

    sympify : bool
        If ``False``, ``sympify`` is not called on ``args``. This
        can be used for speedups for very large tuples where the
        elements are known to already be SymPy objects.

    Examples
    ========

    >>> from sympy import Tuple, symbols
    >>> a, b, c, d = symbols('a b c d')
    >>> Tuple(a, b, c)[1:]
    (b, c)
    >>> Tuple(a, b, c).subs(a, d)
    (d, b, c)

    """

    if TYPE_CHECKING:
        @property
        def args(self) -> tuple[_T, ...]: ...
    else:
        @property
        def args(self):
            return self._args

    def __new__(cls, *args: _T, **kwargs: bool) -> Tuple[_T]:
        if kwargs.get('sympify', True):
            sym_args = tuple(sympify(arg) for arg in args)
            return Basic.__new__(cls, *sym_args)
        return Basic.__new__(cls, *args)

    @overload
    def __getitem__(self, i: int) -> _T: ...

    @overload
    def __getitem__(self, i: slice) -> Tuple[_T]: ...

    def __getitem__(self, i: int | slice) -> _T | Tuple[_T]:
        if isinstance(i, slice):
            indices = i.indices(len(self))
            return Tuple(*(self.args[j] for j in range(*indices)))
        return self.args[i]

    def __len__(self) -> int:
        return len(self.args)

    def __contains__(self, item) -> bool:
        return item in self.args

    def __iter__(self) -> Iterator[_T]:
        return iter(self.args)

    def __add__(self, other: Tuple[_T] | tuple[_T, ...]) -> Tuple[_T]:
        if isinstance(other, Tuple):
            return Tuple(*(self.args + other.args))
        elif isinstance(other, tuple):
            return Tuple(*(self.args + other))
        else:
            return NotImplemented

    def __radd__(self, other: Tuple[_T] | tuple[_T, ...]) -> Tuple[_T]:
        if isinstance(other, Tuple):
            return Tuple(*(other.args + self.args))
        elif isinstance(other, tuple):
            return Tuple(*(other + self.args))
        else:
            return NotImplemented

    def __mul__(self, other: int) -> Tuple[_T]:
        try:
            n = as_int(other)
        except ValueError:
            raise TypeError("Can't multiply sequence by non-integer of type '%s'" % type(other))
        return self.func(*(self.args*n))

    __rmul__ = __mul__

    def __eq__(self, other) -> bool:
        if isinstance(other, Basic):
            return super().__eq__(other)
        return self.args == other

    def __ne__(self, other) -> bool:
        if isinstance(other, Basic):
            return super().__ne__(other)
        return self.args != other

    def __hash__(self) -> int:
        return hash(self.args)

    def _to_mpmath(self: Tuple[_Texpr], prec: int) -> tuple:
        return tuple(a._to_mpmath(prec) for a in self.args)

    def __lt__(self: Tuple[_T], other: Tuple[_T]) -> Boolean:
        return _sympify(self.args < other.args)

    def __le__(self: Tuple[_T], other: Tuple[_T]) -> Boolean:
        return _sympify(self.args <= other.args)

    def tuple_count(self, value) -> int:
        """Return number of occurrences of value."""
        return self.args.count(value)

    def index(self, value, start: int | None = None, stop: int | None = None) -> int:
        """Searches and returns the first index of the value."""
        if start is None:
            return self.args.index(value)
        elif stop is None:
            return self.args.index(value, start)
        else:
            return self.args.index(value, start, stop)

    @property
    def kind(self) -> TupleKind:
        return TupleKind(*(i.kind for i in self.args))

_sympy_converter[tuple] = lambda tup: Tuple(*tup)


def tuple_wrapper(method: Callable[_P, _R]) -> Callable[_P, _R]:
    """
    Decorator that converts any tuple in the function arguments into a Tuple.

    Explanation
    ===========

    The motivation for this is to provide simple user interfaces.  The user can
    call a function with regular tuples in the argument, and the wrapper will
    convert them to Tuples before handing them to the function.

    Explanation
    ===========

    >>> from sympy.core.containers import tuple_wrapper
    >>> def f(*args):
    ...    return args
    >>> g = tuple_wrapper(f)

    The decorated function g sees only the Tuple argument:

    >>> g(0, (1, 2), 3)
    (0, (1, 2), 3)

    """
    def wrap_tuples(*args: _P.args, **kw_args: _P.kwargs) -> _R:
        newargs = []
        for arg in args:
            if isinstance(arg, tuple):
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)
    return wrap_tuples


class Dict(Basic):
    """
    Wrapper around the builtin dict object.

    Explanation
    ===========

    The Dict is a subclass of Basic, so that it works well in the
    SymPy framework.  Because it is immutable, it may be included
    in sets, but its values must all be given at instantiation and
    cannot be changed afterwards.  Otherwise it behaves identically
    to the Python dict.

    Examples
    ========

    >>> from sympy import Dict, Symbol

    >>> D = Dict({1: 'one', 2: 'two'})
    >>> for key in D:
    ...    if key == 1:
    ...        print('%s %s' % (key, D[key]))
    1 one

    The args are sympified so the 1 and 2 are Integers and the values
    are Symbols. Queries automatically sympify args so the following work:

    >>> 1 in D
    True
    >>> D.has(Symbol('one')) # searches keys and values
    True
    >>> 'one' in D # not in the keys
    False
    >>> D[1]
    one

    """

    elements: frozenset[Tuple[Basic]]
    _dict: dict[Basic, Basic]

    def __new__(cls, *args) -> Dict:
        if len(args) == 1 and isinstance(args[0], (dict, Dict)):
            items = [Tuple(k, v) for k, v in args[0].items()]
        elif iterable(args) and all(len(arg) == 2 for arg in args):
            items = [Tuple(k, v) for k, v in args]
        else:
            raise TypeError('Pass Dict args as Dict((k1, v1), ...) or Dict({k1: v1, ...})')
        elements = frozenset(items)
        obj = Basic.__new__(cls, *ordered(items))
        obj.elements = elements
        obj._dict = dict(items)
        return obj

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        try:
            key = _sympify(key)
        except SympifyError:
            raise KeyError(key)

        return self._dict[key]

    def __setitem__(self, key, value) -> None:
        raise NotImplementedError("SymPy Dicts are Immutable")

    def items(self):
        '''Returns a set-like object providing a view on dict's items.
        '''
        return self._dict.items()

    def keys(self):
        '''Returns the list of the dict's keys.'''
        return self._dict.keys()

    def values(self):
        '''Returns the list of the dict's values.'''
        return self._dict.values()

    def __iter__(self):
        '''x.__iter__() <==> iter(x)'''
        return iter(self._dict)

    def __len__(self) -> int:
        '''x.__len__() <==> len(x)'''
        return len(self._dict)

    def get(self, key, default=None):
        '''Returns the value for key if the key is in the dictionary.'''
        try:
            key = _sympify(key)
        except SympifyError:
            return default
        return self._dict.get(key, default)

    def __contains__(self, key) -> bool:
        '''D.__contains__(k) -> True if D has a key k, else False'''
        try:
            key = _sympify(key)
        except SympifyError:
            return False
        return key in self._dict

    def __lt__(self, other: Dict) -> Boolean:
        return _sympify(self.args < other.args)

    @property
    def _sorted_args(self) -> tuple:
        return tuple(sorted(self.args, key=default_sort_key))

    def __eq__(self, other) -> bool:
        if isinstance(other, dict):
            return self == Dict(other)
        return super().__eq__(other)

    def __hash__(self) -> int:
        return Basic.__hash__(self)

# this handles dict, defaultdict, OrderedDict
_sympy_converter[dict] = lambda d: Dict(*d.items())


class OrderedSet(MutableSet[_T]):
    def __init__(self, iterable: Iterable[_T] | None = None) -> None:
        if iterable:
            self.map: OrderedDict[_T, None] = OrderedDict((item, None) for item in iterable)
        else:
            self.map = OrderedDict()

    def __len__(self) -> int:
        return len(self.map)

    def __contains__(self, key) -> bool:
        return key in self.map

    def add(self, value: _T) -> None:
        self.map[value] = None

    def discard(self, value: _T) -> None:
        self.map.pop(value, None)

    def pop(self, last: bool = True) -> _T:
        return self.map.popitem(last=last)[0]

    def __iter__(self) -> Iterator[_T]:
        yield from self.map.keys()

    def __repr__(self) -> str:
        if not self.map:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self.map.keys()))

    def intersection(self, other: Iterable[_T]) -> OrderedSet[_T]:
        return self.__class__([val for val in self if val in other])

    def difference(self, other: Iterable[_T]) -> OrderedSet[_T]:
        return self.__class__([val for val in self if val not in other])

    def update(self, iterable: Iterable[_T]) -> None:
        for val in iterable:
            self.add(val)


class TupleKind(Kind):
    """
    TupleKind is a subclass of Kind, which is used to define Kind of ``Tuple``.

    Parameters of TupleKind will be kinds of all the arguments in Tuples, for
    example

    Parameters
    ==========

    args : tuple(element_kind)
       element_kind is kind of element.
       args is tuple of kinds of element

    Examples
    ========

    >>> from sympy import Tuple, Matrix
    >>> Tuple(1, 2).kind
    TupleKind(NumberKind, NumberKind)
    >>> Tuple(Matrix([1, 2]), 1).kind
    TupleKind(MatrixKind(NumberKind), NumberKind)
    >>> Tuple(1, 2).kind.element_kind
    (NumberKind, NumberKind)

    See Also
    ========

    sympy.matrices.kind.MatrixKind
    sympy.core.kind.NumberKind
    sympy.sets.sets.SetKind
    """
    element_kind: tuple[Kind, ...]
    def __new__(cls, *args: Kind) -> TupleKind:
        obj = super().__new__(cls, *args)
        obj.element_kind = args
        return obj

    def __repr__(self) -> str:
        return "TupleKind{}".format(self.element_kind)
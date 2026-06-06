"""Module for SymPy containers

    (SymPy objects that store other SymPy objects)

    The containers implemented in this module are subclassed to Basic.
    They are supposed to work seamlessly within the SymPy framework.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import MutableSet, KeysView, ValuesView, ItemsView
from typing import Callable, Iterable, Iterator, TypeVar, ParamSpec, overload

from .basic import Basic
from .sorting import default_sort_key, ordered
from .sympify import _sympify, sympify, _sympy_converter, SympifyError
from sympy.core.kind import Kind
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import as_int

P = ParamSpec('P')
R = TypeVar('R')


class Tuple(Basic):
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

    @property
    def args(self) -> tuple[Basic, ...]:
        return self._args

    def __new__(cls, *args: object, **kwargs: bool) -> Tuple:
        if kwargs.get('sympify', True):
            args = tuple(sympify(arg) for arg in args)
        obj = Basic.__new__(cls, *args)
        return obj

    @overload
    def __getitem__(self, i: int) -> Basic: ...

    @overload
    def __getitem__(self, i: slice) -> Tuple: ...

    def __getitem__(self, i: int | slice) -> Basic | Tuple:
        if isinstance(i, slice):
            indices = i.indices(len(self))
            return Tuple(*(self.args[j] for j in range(*indices)))
        return self.args[i]

    def __len__(self) -> int:
        return len(self.args)

    def __contains__(self, item: object) -> bool:
        return item in self.args

    def __iter__(self) -> Iterator[Basic]:
        return iter(self.args)

    def __add__(self, other: Tuple | tuple[Basic, ...]) -> Tuple:
        if isinstance(other, Tuple):
            return Tuple(*(self.args + other.args))
        elif isinstance(other, tuple):
            return Tuple(*(self.args + other))
        else:
            return NotImplemented

    def __radd__(self, other: Tuple | tuple[Basic, ...]) -> Tuple:
        if isinstance(other, Tuple):
            return Tuple(*(other.args + self.args))
        elif isinstance(other, tuple):
            return Tuple(*(other + self.args))
        else:
            return NotImplemented

    def __mul__(self, other: int) -> Tuple:
        try:
            n = as_int(other)
        except ValueError:
            raise TypeError("Can't multiply sequence by non-integer of type '%s'" % type(other))
        return self.func(*(self.args*n))

    __rmul__ = __mul__

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Basic):
            return super().__eq__(other)
        return self.args == other

    def __ne__(self, other: object) -> bool:
        if isinstance(other, Basic):
            return super().__ne__(other)
        return self.args != other

    def __hash__(self) -> int:
        return hash(self.args)

    def _to_mpmath(self, prec: int) -> tuple[object, ...]:
        return tuple(a._to_mpmath(prec) for a in self.args)  # type: ignore[attr-defined]

    def __lt__(self, other: Tuple) -> Basic:
        return _sympify(self.args < other.args)

    def __le__(self, other: Tuple) -> Basic:
        return _sympify(self.args <= other.args)

    # XXX: Basic defines count() as something different, so we can't
    # redefine it here. Originally this lead to cse() test failure.
    def tuple_count(self, value: object) -> int:
        """Return number of occurrences of value."""
        return self.args.count(value)

    def index(self, value: object, start: int | None = None, stop: int | None = None) -> int:
        """Searches and returns the first index of the value."""
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

        if start is None:
            return self.args.index(value)
        elif stop is None:
            return self.args.index(value, start)
        else:
            return self.args.index(value, start, stop)

    @property
    def kind(self) -> TupleKind:  # type: ignore[override]
        """
        The kind of a Tuple instance.

        The kind of a Tuple is always of :class:`TupleKind` but
        parametrised by the number of elements and the kind of each element.

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
        """
        return TupleKind(*(i.kind for i in self.args))  # type: ignore[attr-defined]

_sympy_converter[tuple] = lambda tup: Tuple(*tup)


def tuple_wrapper(method: Callable[P, R]) -> Callable[P, R]:
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
    def wrap_tuples(*args: P.args, **kw_args: P.kwargs) -> R:
        newargs: list[object] = []
        for arg in args:
            if isinstance(arg, tuple):
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)  # type: ignore[arg-type]
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

    elements: frozenset[Tuple]
    _dict: dict[Basic, Basic]

    def __new__(cls, *args: object) -> Dict:
        if len(args) == 1 and isinstance(args[0], (dict, Dict)):
            items = [Tuple(k, v) for k, v in args[0].items()]
        elif iterable(args) and all(len(arg) == 2 for arg in args):  # type: ignore[arg-type]
            items = [Tuple(k, v) for k, v in args]  # type: ignore[misc]
        else:
            raise TypeError('Pass Dict args as Dict((k1, v1), ...) or Dict({k1: v1, ...})')
        elements = frozenset(items)
        obj = Basic.__new__(cls, *ordered(items))
        obj.elements = elements
        obj._dict = dict(items)  # type: ignore[arg-type]
        return obj

    def __getitem__(self, key: object) -> Basic:
        """x.__getitem__(y) <==> x[y]"""
        try:
            key = _sympify(key)
        except SympifyError:
            raise KeyError(key)

        return self._dict[key]  # type: ignore[index]

    def __setitem__(self, key: object, value: object) -> None:
        raise NotImplementedError("SymPy Dicts are Immutable")

    def items(self) -> ItemsView[Basic, Basic]:
        '''Returns a set-like object providing a view on dict's items.
        '''
        return self._dict.items()

    def keys(self) -> KeysView[Basic]:
        '''Returns the list of the dict's keys.'''
        return self._dict.keys()

    def values(self) -> ValuesView[Basic]:
        '''Returns the list of the dict's values.'''
        return self._dict.values()

    def __iter__(self) -> Iterator[Basic]:
        '''x.__iter__() <==> iter(x)'''
        return iter(self._dict)

    def __len__(self) -> int:
        '''x.__len__() <==> len(x)'''
        return self._dict.__len__()

    def get(self, key: object, default: object = None) -> Basic | object:
        '''Returns the value for key if the key is in the dictionary.'''
        try:
            key = _sympify(key)
        except SympifyError:
            return default
        return self._dict.get(key, default)  # type: ignore[call-overload]

    def __contains__(self, key: object) -> bool:
        '''D.__contains__(k) -> True if D has a key k, else False'''
        try:
            key = _sympify(key)
        except SympifyError:
            return False
        return key in self._dict

    def __lt__(self, other: Dict) -> Basic:
        return _sympify(self.args < other.args)

    @property
    def _sorted_args(self) -> tuple[Basic, ...]:
        return tuple(sorted(self.args, key=default_sort_key))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, dict):
            return self == Dict(other)
        return super().__eq__(other)

    def __hash__(self) -> int:
        return Basic.__hash__(self)

# this handles dict, defaultdict, OrderedDict
_sympy_converter[dict] = lambda d: Dict(*d.items())


class OrderedSet(MutableSet[Basic]):
    def __init__(self, iterable: Iterable[Basic] | None = None) -> None:
        if iterable:
            self.map: OrderedDict[Basic, None] = OrderedDict((item, None) for item in iterable)
        else:
            self.map = OrderedDict()

    def __len__(self) -> int:
        return len(self.map)

    def __contains__(self, key: object) -> bool:
        return key in self.map

    def add(self, value: Basic) -> None:
        self.map[value] = None

    def discard(self, value: Basic) -> None:
        self.map.pop(value, None)

    def pop(self, last: bool = True) -> Basic:
        return self.map.popitem(last=last)[0]

    def __iter__(self) -> Iterator[Basic]:
        yield from self.map.keys()

    def __repr__(self) -> str:
        if not self.map:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self.map.keys()))

    def intersection(self, other: Iterable[Basic]) -> OrderedSet:
        return self.__class__([val for val in self if val in other])

    def difference(self, other: Iterable[Basic]) -> OrderedSet:
        return self.__class__([val for val in self if val not in other])

    def update(self, iterable: Iterable[Basic]) -> None:
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
        obj = super().__new__(cls, *args)  # type: ignore[call-arg]
        obj.element_kind = args            # type: ignore[attr-defined]
        return obj                         # type: ignore[return-value]

    def __repr__(self) -> str:
        return "TupleKind{}".format(self.element_kind)
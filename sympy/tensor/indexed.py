"""Module that defines indexed objects

The classes ``IndexedBase``, ``Indexed``, and ``Idx`` represent a
matrix element ``M[i, j]`` as in the following diagram::

   1) The Indexed class represents the entire indexed object.
              |
           ___|___
          '       '
           M[i, j]
          /   \__\______
          |             |
          |             |
          |     2) The Idx class represents indices; each Idx can
          |        optionally contain information about its range.
          |
    3) IndexedBase represents the 'stem' of an Indexed object, here M.

There can be any number of indices on an ``Indexed`` object.  No symmetry
properties are implemented in these objects. Implicit contraction of
repeated indices is supported via the ``EinsteinSum`` function.

Note that the support for complicated (i.e. non-atomic) integer
expressions as indices is limited.  (This should be improved in
future releases.)

Basic Examples
==============

To express the above matrix element example:

>>> from sympy import symbols, IndexedBase, Idx
>>> M = IndexedBase('M')
>>> i, j = symbols('i j', cls=Idx)
>>> M[i, j]
M[i, j]

Repeated indices in a monomial imply summation when enclosed by the
``EinsteinSum`` function. To express a matrix-vector product in terms of
``Indexed`` objects:

>>> from sympy import EinsteinSum
>>> x = IndexedBase('x')
>>> EinsteinSum(M[i, j] * x[j])
EinsteinSum(x[j]*M[i, j])

``EinsteinSum`` expressions can be differentiated:

>>> EinsteinSum(x[j] * x[j]).diff(x[i])
EinsteinSum(2*x[i])

They can also be converted into Python functions accepting NumPy arrays:

>>> from sympy import lambdify
>>> func = lambdify([M, x], EinsteinSum(M[i, j] * x[j]))  # doctest: +SKIP
>>> func([[1, -1], [0, 1]], [1, 2])  # doctest: +SKIP
[-1 2]

"""

from __future__ import print_function, division

from sympy.core import (Expr, Tuple, Symbol, sympify, S)
from sympy.core.compatibility import (is_sequence, string_types, NotIterable,
                                      default_sort_key)


class IndexException(Exception):
    pass


class Indexed(Expr):
    """Represents a mathematical object with indices.

    >>> from sympy import Indexed, IndexedBase, Idx, symbols
    >>> i, j = symbols('i j', cls=Idx)
    >>> Indexed('A', i, j)
    A[i, j]

    It is recommended that ``Indexed`` objects be created via ``IndexedBase``:

    >>> A = IndexedBase('A')
    >>> Indexed('A', i, j) == A[i, j]
    True

    """
    is_commutative = True
    is_Symbol = True

    @property
    def _base_Symbol(self):
        return self.base.label

    def __new__(cls, base, *args, **kw_args):
        from sympy.utilities.misc import filldedent

        if not args:
            raise IndexException("Indexed needs at least one index.")
        if isinstance(base, (string_types, Symbol)):
            base = IndexedBase(base)
        elif not isinstance(base, IndexedBase):
            raise TypeError(filldedent("""
                Indexed expects string, Symbol, or IndexedBase as base."""))
        args = list(map(sympify, args))
        return Expr.__new__(cls, base, *args, **kw_args)

    @property
    def _diff_wrt(self):
        """Allow derivatives with respect to an ``Indexed`` object."""
        return True

    def _eval_derivative(self, wrt):
        if isinstance(wrt, Indexed) and wrt.base == self.base:
            if len(self.indices) != len(wrt.indices):
                msg = "Different # of indices: d({!s})/d({!s})".format(self,
                                                                       wrt)
                raise IndexException(msg)
            result = S.One
            for index1, index2 in zip(self.indices, wrt.indices):
                result *= DeltaIndexedBase()[index1, index2]
            return result
        else:
            return S.Zero

    @property
    def base(self):
        """Returns the ``IndexedBase`` of the ``Indexed`` object.

        Examples
        ========

        >>> from sympy import Indexed, IndexedBase, Idx, symbols
        >>> i, j = symbols('i j', cls=Idx)
        >>> Indexed('A', i, j).base
        A
        >>> B = IndexedBase('B')
        >>> B == B[i, j].base
        True

        """
        return self.args[0]

    @property
    def indices(self):
        """
        Returns the indices of the ``Indexed`` object.

        Examples
        ========

        >>> from sympy import Indexed, Idx, symbols
        >>> i, j = symbols('i j', cls=Idx)
        >>> Indexed('A', i, j).indices
        (i, j)

        """
        return self.args[1:]

    @property
    def rank(self):
        """
        Returns the rank of the ``Indexed`` object.

        Examples
        ========

        >>> from sympy import Indexed, Idx, symbols
        >>> i, j, k, l, m = symbols('i:m', cls=Idx)
        >>> Indexed('A', i, j).rank
        2
        >>> q = Indexed('A', i, j, k, l, m)
        >>> q.rank
        5
        >>> q.rank == len(q.indices)
        True

        """
        return len(self.args) - 1

    @property
    def shape(self):
        """Returns a list with dimensions of each index.

        Dimensions is a property of the array, not of the indices.  Still, if
        the ``IndexedBase`` does not define a shape attribute, it is assumed
        that the ranges of the indices correspond to the shape of the array.

        >>> from sympy import IndexedBase, Idx, symbols
        >>> n, m = symbols('n m', integer=True)
        >>> i = Idx('i', m)
        >>> j = Idx('j', m)
        >>> A = IndexedBase('A', shape=(n, n))
        >>> B = IndexedBase('B')
        >>> A[i, j].shape
        (n, n)
        >>> B[i, j].shape
        (m, m)
        """
        from sympy.utilities.misc import filldedent

        if self.base.shape:
            return self.base.shape
        try:
            return Tuple(*[i.upper - i.lower + 1 for i in self.indices])
        except AttributeError:
            raise IndexException(filldedent("""
                Range is not defined for all indices in: %s""" % self))
        except TypeError:
            raise IndexException(filldedent("""
                Shape cannot be inferred from Idx with
                undefined range: %s""" % self))

    @property
    def ranges(self):
        """Returns a list of tuples with lower and upper range of each index.

        If an index does not define the data members upper and lower, the
        corresponding slot in the list contains ``None`` instead of a tuple.

        Examples
        ========

        >>> from sympy import Indexed,Idx, symbols
        >>> Indexed('A', Idx('i', 2), Idx('j', 4), Idx('k', 8)).ranges
        [(0, 1), (0, 3), (0, 7)]
        >>> Indexed('A', Idx('i', 3), Idx('j', 3), Idx('k', 3)).ranges
        [(0, 2), (0, 2), (0, 2)]
        >>> x, y, z = symbols('x y z', integer=True)
        >>> Indexed('A', x, y, z).ranges
        [None, None, None]

        """
        ranges = []
        for i in self.indices:
            try:
                ranges.append(Tuple(i.lower, i.upper))
            except AttributeError:
                ranges.append(None)
        return ranges

    def _sympystr(self, p):
        indices = list(map(p.doprint, self.indices))
        return "%s[%s]" % (p.doprint(self.base), ", ".join(indices))


class IndexedBase(Expr, NotIterable):
    """Represent the base or stem of an ``Indexed`` object

    The ``IndexedBase`` class represent an array that contains elements. The
    main purpose of this class is to allow the convenient creation of objects of
    the ``Indexed`` class.  The ``__getitem__`` method of ``IndexedBase``
    returns an instance of ``Indexed``.

    >>> from sympy import IndexedBase, Idx, symbols
    >>> A = IndexedBase('A'); A
    A
    >>> type(A)
    <class 'sympy.tensor.indexed.IndexedBase'>

    When an ``IndexedBase`` object receives indices, it returns an array with
    named axes, represented by an ``Indexed`` object:

    >>> i, j = symbols('i j', integer=True)
    >>> A[i, j, 2]
    A[i, j, 2]
    >>> type(A[i, j, 2])
    <class 'sympy.tensor.indexed.Indexed'>

    The ``IndexedBase`` constructor takes an optional shape argument.  If given,
    it overrides any shape information in the indices. (But not the index
    ranges!)

    >>> m, n, o, p = symbols('m n o p', integer=True)
    >>> i = Idx('i', m)
    >>> j = Idx('j', n)
    >>> A[i, j].shape
    (m, n)
    >>> B = IndexedBase('B', shape=(o, p))
    >>> B[i, j].shape
    (o, p)

    """
    is_commutative = True
    is_Symbol = True

    def __new__(cls, label, shape=None, **kw_args):
        if isinstance(label, string_types):
            label = Symbol(label)
        elif isinstance(label, Symbol):
            pass
        else:
            raise TypeError("Base label should be a string or Symbol.")

        obj = Expr.__new__(cls, label, **kw_args)
        if is_sequence(shape):
            obj._shape = Tuple(*shape)
        else:
            obj._shape = sympify(shape)
        return obj

    @property
    def args(self):
        """Returns the arguments used to create this ``IndexedBase`` object.

        Examples
        ========

        >>> from sympy import IndexedBase
        >>> from sympy.abc import x, y
        >>> IndexedBase('A', shape=(x, y)).args
        (A, (x, y))

        """
        if self._shape:
            return self._args + (self._shape,)
        else:
            return self._args

    def _hashable_content(self):
        if self._shape:
            return super(IndexedBase, self)._hashable_content() + (self._shape,)
        else:
            return super(IndexedBase, self)._hashable_content()

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            if self.shape and len(self.shape) != len(indices):
                raise IndexException("Rank mismatch.")
            return Indexed(self, *indices, **kw_args)
        else:
            if self.shape and len(self.shape) != 1:
                raise IndexException("Rank mismatch.")
            return Indexed(self, indices, **kw_args)

    @property
    def shape(self):
        """Returns the shape of the ``IndexedBase`` object.

        Examples
        ========

        >>> from sympy import IndexedBase, Idx, Symbol
        >>> from sympy.abc import x, y
        >>> IndexedBase('A', shape=(x, y)).shape
        (x, y)

        Note: If the shape of the ``IndexedBase`` is specified, it will override
        any shape information given by the indices.

        >>> A = IndexedBase('A', shape=(x, y))
        >>> B = IndexedBase('B')
        >>> i = Idx('i', 2)
        >>> j = Idx('j', 1)
        >>> A[i, j].shape
        (x, y)
        >>> B[i, j].shape
        (2, 1)

        """
        return self._shape

    @property
    def label(self):
        """Returns the label of the ``IndexedBase`` object.

        Examples
        ========

        >>> from sympy import IndexedBase
        >>> from sympy.abc import x, y
        >>> IndexedBase('A', shape=(x, y)).label
        A

        """
        return self.args[0]

    def _sympystr(self, p):
        return p.doprint(self.label)


class DeltaIndexedBase(IndexedBase):
    """Class to represent the ``IndexedBase`` of a `Kronecker delta tensor`_.

    Simplification rules exist that involve ``DeltaIndexedBase`` objects; see
    ``EinsteinSum.simplify_deltas``.

    Examples
    ========

    >>> from sympy import DeltaIndexedBase, Idx, symbols
    >>> i, j = symbols('i j', cls=Idx)
    >>> delta = DeltaIndexedBase()
    >>> delta[i, j]
    delta[i, j]

    A custom label can also be given:

    >>> DeltaIndexedBase('delta1')[i, j]
    delta1[i, j]

    .. _Kronecker delta tensor: http://en.wikipedia.org/wiki/Kronecker_delta

    See Also
    ========

    sympy.functions.special.tensor_functions.KroneckerDelta

    """

    def __new__(cls, label='delta', **kw_args):
        return IndexedBase.__new__(cls, label=label, **kw_args)

    def __getitem__(self, indices, **kw_args):
        if not len(indices) == 2:
            raise IndexException(
                "Kronecker delta tensor should have exactly two indices")
        indices = sorted(indices, key=default_sort_key)
        return super(DeltaIndexedBase, self).__getitem__(indices, **kw_args)


class Idx(Expr):
    """Represents an integer index as an ``Integer`` or integer expression.

    There are a number of ways to create an ``Idx`` object.  The constructor
    takes two arguments:

    ``label``
        An integer or a symbol that labels the index.
    ``range``
        Optionally you can specify a range as either

        * ``Symbol`` or integer: This is interpreted as a dimension. Lower and
          upper bounds are set to ``0`` and ``range - 1``, respectively.
        * ``tuple``: The two elements are interpreted as the lower and upper
          bounds of the range, respectively.

    Note: the ``Idx`` constructor is rather pedantic in that it only accepts
    integer arguments.  The only exception is that you can use ``-oo`` and
    ``oo`` to specify an unbounded range.  For all other cases, both label and
    bounds must be declared as integers, e.g. if ``n`` is given as an argument
    then ``n.is_integer`` must return ``True``.

    For convenience, if the label is given as a string it is automatically
    converted to an integer symbol.  (Note: this conversion is not done for
    range or dimension arguments.)

    Examples
    ========

    >>> from sympy import IndexedBase, Idx, symbols, oo
    >>> n, i, L, U = symbols('n i L U', integer=True)

    If a string is given for the label an integer ``Symbol`` is created and the
    bounds are both ``None``:

    >>> idx = Idx('qwerty'); idx
    qwerty
    >>> idx.lower, idx.upper
    (None, None)

    Both upper and lower bounds can be specified:

    >>> idx = Idx(i, (L, U)); idx
    i
    >>> idx.lower, idx.upper
    (L, U)

    When only a single bound is given it is interpreted as the dimension
    and the lower bound defaults to 0:

    >>> idx = Idx(i, n); idx.lower, idx.upper
    (0, n - 1)
    >>> idx = Idx(i, 4); idx.lower, idx.upper
    (0, 3)
    >>> idx = Idx(i, oo); idx.lower, idx.upper
    (0, oo)

    The label can be a literal integer instead of a string/``Symbol``:

    >>> idx = Idx(2, n); idx.lower, idx.upper
    (0, n - 1)
    >>> idx.label
    2

    """

    is_integer = True

    def __new__(cls, label, range=None, **kw_args):
        from sympy.utilities.misc import filldedent

        if isinstance(label, string_types):
            label = Symbol(label, integer=True)
        label, range = list(map(sympify, (label, range)))

        if not label.is_integer:
            raise TypeError("Idx object requires an integer label.")

        elif is_sequence(range):
            if len(range) != 2:
                msg = "Idx range tuple must have length 2, but got {}"
                msg = msg.format(range)
                raise ValueError(msg)
            for bound in range:
                if not (bound.is_integer or abs(bound) is S.Infinity):
                    raise TypeError("Idx object requires integer bounds.")
            args = label, Tuple(*range)
        elif isinstance(range, Expr):
            if not (range.is_integer or range is S.Infinity):
                raise TypeError("Idx object requires an integer dimension.")
            args = label, Tuple(0, range - 1)
        elif range:
            raise TypeError(filldedent("""
                The range must be an ordered iterable or
                integer SymPy expression."""))
        else:
            args = label,

        obj = Expr.__new__(cls, *args, **kw_args)
        return obj

    @property
    def label(self):
        """Returns the label of the ``Idx`` object.

        Examples
        ========

        >>> from sympy import Idx, Symbol
        >>> Idx(2).label
        2
        >>> j = Symbol('j', integer=True)
        >>> Idx(j).label
        j
        >>> Idx(j + 1).label
        j + 1

        """
        return self.args[0]

    @property
    def lower(self):
        """Returns the lower bound of the ``Idx``.

        Examples
        ========

        >>> from sympy import Idx
        >>> Idx('j', 2).lower
        0
        >>> Idx('j', 5).lower
        0
        >>> Idx('j').lower is None
        True

        """
        try:
            return self.args[1][0]
        except IndexError:
            return

    @property
    def upper(self):
        """Returns the upper bound of the ``Idx``.

        Examples
        ========

        >>> from sympy import Idx
        >>> Idx('j', 2).upper
        1
        >>> Idx('j', 5).upper
        4
        >>> Idx('j').upper is None
        True

        """
        try:
            return self.args[1][1]
        except IndexError:
            return

    def _sympystr(self, p):
        return p.doprint(self.label)

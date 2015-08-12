"""Module that defines indexed objects

    The classes IndexedBase, Indexed and Idx would represent a matrix element
    M[i, j] as in the following graph::

       1) The Indexed class represents the entire indexed object.
                  |
               ___|___
              '       '
               M[i, j]
              /   \__\______
              |             |
              |             |
              |     2) The Idx class represent indices and each Idx can
              |        optionally contain information about its range.
              |
        3) IndexedBase represents the `stem' of an indexed object, here `M'.
           The stem used by itself is usually taken to represent the entire
           array.

    There can be any number of indices on an Indexed object.  No
    transformation properties are implemented in these Base objects, but
    implicit contraction of repeated indices is supported.

    Note that the support for complicated (i.e. non-atomic) integer
    expressions as indices is limited.  (This should be improved in
    future releases.)

    Examples
    ========

    To express the above matrix element example you would write:

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols
    >>> M = IndexedBase('M')
    >>> i, j = symbols('i j', cls=Idx)
    >>> M[i, j]
    M[i, j]

    Repeated indices in a product implies a summation, so to express a
    matrix-vector product in terms of Indexed objects:

    >>> x = IndexedBase('x')
    >>> M[i, j]*x[j]
    x[j]*M[i, j]

    If the indexed objects will be converted to component based arrays, e.g.
    with the code printers or the autowrap framework, you also need to provide
    (symbolic or numerical) dimensions.  This can be done by passing an
    optional shape parameter to IndexedBase upon construction:

    >>> dim1, dim2 = symbols('dim1 dim2', integer=True)
    >>> A = IndexedBase('A', shape=(dim1, 2*dim1, dim2))
    >>> A.shape
    (dim1, 2*dim1, dim2)
    >>> A[i, j, 3].shape
    (dim1, 2*dim1, dim2)

    If an IndexedBase object has no shape information, it is assumed that the
    array is as large as the ranges of its indices:

    >>> n, m = symbols('n m', integer=True)
    >>> i = Idx('i', m)
    >>> j = Idx('j', n)
    >>> M[i, j].shape
    (m, n)
    >>> M[i, j].ranges
    [(0, m - 1), (0, n - 1)]

    The above can be compared with the following:

    >>> A[i, 2, j].shape
    (dim1, 2*dim1, dim2)
    >>> A[i, 2, j].ranges
    [(0, m - 1), None, (0, n - 1)]

    To analyze the structure of indexed expressions, you can use the methods
    get_indices() and get_contraction_structure():

    >>> from sympy.tensor import get_indices, get_contraction_structure
    >>> get_indices(A[i, j, j])
    (set([i]), {})
    >>> get_contraction_structure(A[i, j, j])
    {(j,): set([A[i, j, j]])}

    See the appropriate docstrings for a detailed explanation of the output.

"""

#   TODO:  (some ideas for improvement)
#
#   o test and guarantee numpy compatibility
#      - implement full support for broadcasting
#      - strided arrays
#
#   o more functions to analyze indexed expressions
#      - identify standard constructs, e.g matrix-vector product in a subexpression
#
#   o functions to generate component based arrays (numpy and sympy.Matrix)
#      - generate a single array directly from Indexed
#      - convert simple sub-expressions
#
#   o sophisticated indexing (possibly in subclasses to preserve simplicity)
#      - Idx with range smaller than dimension of Indexed
#      - Idx with stepsize != 1
#      - Idx with step determined by function call

from __future__ import print_function, division

from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable, range


class IndexException(Exception):
    pass


class Indexed(Expr):
    """Represents a mathematical object with indices.

    >>> from sympy.tensor import Indexed, IndexedBase, Idx
    >>> from sympy import symbols
    >>> i, j = symbols('i j', cls=Idx)
    >>> Indexed('A', i, j)
    A[i, j]

    It is recommended that Indexed objects are created via IndexedBase:

    >>> A = IndexedBase('A')
    >>> Indexed('A', i, j) == A[i, j]
    True

    """
    is_commutative = True

    def __new__(cls, base, *args):
        from sympy.utilities.misc import filldedent

        if not args:
            raise IndexException("Indexed needs at least one index.")
        if isinstance(base, (string_types, Symbol)):
            base = IndexedBase(base)
        elif not hasattr(base, '__getitem__') and not isinstance(base, IndexedBase):
            raise TypeError(filldedent("""
                Indexed expects string, Symbol or IndexedBase as base."""))
        args = list(map(sympify, args))
        return Expr.__new__(cls, base, *args)

    @property
    def base(self):
        """Returns the IndexedBase of the Indexed object.

        Examples
        ========

        >>> from sympy.tensor import Indexed, IndexedBase, Idx
        >>> from sympy import symbols
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
        Returns the indices of the Indexed object.

        Examples
        ========

        >>> from sympy.tensor import Indexed, Idx
        >>> from sympy import symbols
        >>> i, j = symbols('i j', cls=Idx)
        >>> Indexed('A', i, j).indices
        (i, j)

        """
        return self.args[1:]

    @property
    def rank(self):
        """
        Returns the rank of the Indexed object.

        Examples
        ========

        >>> from sympy.tensor import Indexed, Idx
        >>> from sympy import symbols
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
        the IndexedBase does not define a shape attribute, it is assumed that
        the ranges of the indices correspond to the shape of the array.

        >>> from sympy.tensor.indexed import IndexedBase, Idx
        >>> from sympy import symbols
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
    """Represent the base or stem of an indexed object

    The IndexedBase class represent an array that contains elements. The main purpose
    of this class is to allow the convenient creation of objects of the Indexed
    class.  The __getitem__ method of IndexedBase returns an instance of
    Indexed.  Alone, without indices, the IndexedBase class can be used as a
    notation for e.g. matrix equations, resembling what you could do with the
    Symbol class.  But, the IndexedBase class adds functionality that is not
    available for Symbol instances:

      -  An IndexedBase object can optionally store shape information.  This can
         be used in to check array conformance and conditions for numpy
         broadcasting.  (TODO)
      -  An IndexedBase object implements syntactic sugar that allows easy symbolic
         representation of array operations, using implicit summation of
         repeated indices.
      -  The IndexedBase object symbolizes a mathematical structure equivalent
         to arrays, and is recognized as such for code generation and automatic
         compilation and wrapping.

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols
    >>> A = IndexedBase('A'); A
    A
    >>> type(A)
    <class 'sympy.tensor.indexed.IndexedBase'>

    When an IndexedBase object receives indices, it returns an array with named
    axes, represented by an Indexed object:

    >>> i, j = symbols('i j', integer=True)
    >>> A[i, j, 2]
    A[i, j, 2]
    >>> type(A[i, j, 2])
    <class 'sympy.tensor.indexed.Indexed'>

    The IndexedBase constructor takes an optional shape argument.  If given,
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
        """Returns the arguments used to create this IndexedBase object.

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
        return Expr._hashable_content(self) + (self._shape,)

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
        """Returns the shape of the IndexedBase object.

        Examples
        ========

        >>> from sympy import IndexedBase, Idx, Symbol
        >>> from sympy.abc import x, y
        >>> IndexedBase('A', shape=(x, y)).shape
        (x, y)

        Note: If the shape of the IndexedBase is specified, it will override
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
        """Returns the label of the IndexedBase object.

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


class Idx(Expr):
    """Represents an integer index as an Integer or integer expression.

    There are a number of ways to create an Idx object.  The constructor
    takes two arguments:

    ``label``
        An integer or a symbol that labels the index.
    ``range``
        Optionally you can specify a range as either

    - Symbol or integer: This is interpreted as a dimension. Lower and
      upper bounds are set to 0 and range - 1, respectively.
    - tuple: The two elements are interpreted as the lower and upper
      bounds of the range, respectively.

    Note: the Idx constructor is rather pedantic in that it only accepts
    integer arguments.  The only exception is that you can use oo and -oo to
    specify an unbounded range.  For all other cases, both label and bounds
    must be declared as integers, e.g. if n is given as an argument then
    n.is_integer must return True.

    For convenience, if the label is given as a string it is automatically
    converted to an integer symbol.  (Note: this conversion is not done for
    range or dimension arguments.)

    Examples
    ========

    >>> from sympy.tensor import Idx
    >>> from sympy import symbols, oo
    >>> n, i, L, U = symbols('n i L U', integer=True)

    If a string is given for the label an integer Symbol is created and the
    bounds are both None:

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

    """

    is_integer = True

    def __new__(cls, label, range=None, **kw_args):
        from sympy.utilities.misc import filldedent

        if isinstance(label, string_types):
            label = Symbol(label, integer=True)
        label, range = list(map(sympify, (label, range)))

        if label.is_Number:
            if not label.is_integer:
                raise TypeError("Index is not an integer number.")
            return label

        if not label.is_integer:
            raise TypeError("Idx object requires an integer label.")

        elif is_sequence(range):
            if len(range) != 2:
                raise ValueError(filldedent("""
                    Idx range tuple must have length 2, but got %s""" % len(range)))
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
        """Returns the label (Integer or integer expression) of the Idx object.

        Examples
        ========

        >>> from sympy import Idx, Symbol
        >>> x = Symbol('x', integer=True)
        >>> Idx(x).label
        x
        >>> j = Symbol('j', integer=True)
        >>> Idx(j).label
        j
        >>> Idx(j + 1).label
        j + 1

        """
        return self.args[0]

    @property
    def lower(self):
        """Returns the lower bound of the Index.

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
        """Returns the upper bound of the Index.

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

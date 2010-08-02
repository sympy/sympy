"""Module that defines indexed objects with arbitrary transformation properties.

    The classes IndexedBase, Indexed and Idx represent a matrix element M(i, j)
    as in the following graph:

       1) The Indexed class represents an element of the indexed object.
                  |
               ___|___
              '       '
               M(i, j)
              /   \__\____
              |           \
              |            \
              |     2) The Idx class represent indices and each Idx can
              |        optionally contain information about its range.
              |
    3) The IndexedBase class represents the `stem' of an indexed object, here `M'.
    The stem used by itself is usually taken to represent the entire array.


    Examples
    ========

    To express the above matrix element example you would write:

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols
    >>> i, j, n, m = symbols('i j n m', integer=True)
    >>> M = IndexedBase('M')
    >>> M[i, j]
    M[i, j]

    If the indexed objects will be converted to component based arrays, e.g.
    numpy arrays, you also need to provide (symbolic) dimensions.  This is done
    with the Idx class:

    >>> ii = Idx(i, m)
    >>> jj = Idx(j, n)
    >>> M[ii, jj]
    M[i, j]
    >>> M[ii, jj].shape
    (m, n)
    >>> M[ii, jj].ranges
    [(0, -1 + m), (0, -1 + n)]

    To express a matrix-vector product in terms of IndexedBase objects:

    >>> x = IndexedBase('x')
    >>> M[ii, jj]*x[jj]     #doctest: +SKIP
    M[i, j]*x[j]


    TODO:  (some ideas for improvement)

    o test and guarantee numpy compatibility

    o functions to operate on the indexed expressions
       - check_conformance()
       - determine_resulting_indices()
       - determine_summation_indices()
       - identify standard constructs, e.g matrix-vector product in a subexpression

    o functions to generate component based arrays (numpy and sympy.Matrix)
       - generate a single array directly from Indexed
       - convert simple sub-expressions

    o sophisticated indexing (possibly in subclasses to preserve simplicity)
       - Idx with range smaller than dimension of Indexed
       - Idx with stepsize != 1
       - Idx with step determined by function call
"""

from sympy import Expr, Basic, Tuple, Symbol, Integer, sympify, S

class IndexException(Exception):
    pass

class IndexedBase(Expr):
    """Represent the stem of an indexed object, e.g. a numpy array.

    The IndexedStem class represent an array that contains elements. An element
    i of an indexed object A is denoted A(i) and is represented by the
    Indexed class.

    The IndexedBase class allows a simple notation for e.g. matrix equations,
    resembling what you could do with the Symbol class.  But, the IndexedBase class
    adds functionality that is not available for Symbol instances:

      -  An IndexedBase object can optionally store shape information.  This can
         be used in functions that check array conformance and conditions for
         numpy broadcasting.
      -  An IndexedBase object implements syntactic sugar that allows easy symbolic
         representation of array elements,  e.g. you can type A(i, j) to
         create a symbol for element i,j of array A.
      -  The IndexedBase object symbolizes a mathematical structure equivalent to
         arrays, and is recognized as such for code generation and
         TODO: conversion to a numpy array and sympy.Matrix objects.

    >>> from sympy.tensor import IndexedBase
    >>> from sympy import symbols
    >>> a = IndexedBase('a'); a
    a
    >>> type(a)
    <class 'sympy.tensor.indexed.IndexedBase'>

    Objects of type Indexed to symbolize elements of the array will by
    created whenever indices are supplied to the stem:

    >>> i, j, k = symbols('i j k', integer=True)
    >>> a[i, j, k]
    a[i, j, k]
    >>> type(a[i, j, k])
    <class 'sympy.tensor.indexed.Indexed'>

    """
    def __new__(cls, label, shape=None, commutative=True, **kw_args):
        if isinstance(label, basestring):
            label = Symbol(label)

        obj = Expr.__new__(cls, label, **kw_args)
        obj._is_commutative = commutative
        if isinstance(shape, (tuple, list)):
            obj._shape = Tuple(*shape)
        else:
            obj._shape = shape
        return obj

    @property
    def args(self):
        if self._shape:
            return self._args + (self._shape,)
        else:
            return self._args

    def _hashable_content(self):
        return Expr._hashable_content(self) + (self._is_commutative, self._shape)

    def __getitem__(self, indices, **kw_args):
        if isinstance(indices, tuple):
            if self.shape and len(self.shape) != len(indices):
                raise IndexException("Rank mismatch")
            return Indexed(self, *indices, **kw_args)
        else:
            if self.shape and len(self.shape) != 1:
                raise IndexException("Rank mismatch")
            return Indexed(self, indices, **kw_args)

    @property
    def is_commutative(self):
        return self._is_commutative

    @property
    def shape(self):
        return self._shape

    @property
    def label(self):
        return self.args[0]

    def _sympystr(self, p):
        return p.doprint(self.label)

#FIXME only needed for 2.4 compatibility
def _ensure_Idx(arg):
    if isinstance(arg, Idx):
        return arg
    else:
        return Idx(arg)

class Indexed(Expr):
    """Represent an indexed element, e.g. a symbolic array element.

    >>> from sympy.tensor import Indexed
    >>> from sympy import symbols
    >>> i, j, k = symbols('i j k', integer=True)
    >>> Indexed('a', i, j)
    a[i, j]

    """

    def __new__(cls, stem, *args, **kw_args):
        if not args: raise IndexException("Indexed needs at least one index")
        if isinstance(stem, (basestring, Symbol)):
            stem = IndexedBase(stem)
        elif not isinstance(stem, IndexedBase):
            raise TypeError("Indexed expects string, Symbol or IndexedBase as stem")
        # FIXME: 2.4 compatibility
        args = map(_ensure_Idx, args)
        # args = tuple([ a if isinstance(a, Idx) else Idx(a) for a in args ])
        return Expr.__new__(cls, stem, *args, **kw_args)

    @property
    def stem(self):
        return self.args[0]

    @property
    def is_commutative(self):
        return self.stem.is_commutative

    @property
    def indices(self):
        return self.args[1:]

    @property
    def rank(self):
        """returns the number of indices"""
        return len(self.args)-1

    @property
    def shape(self):
        """returns a list with dimensions of each index"""
        if self.stem.shape:
            return self.stem.shape
        try:
            return tuple( i.upper - i.lower + 1 for i in self.indices )
        except TypeError:
            # Let's return a more meaningful error
            raise IndexException("Shape is not defined")

    @property
    def ranges(self):
        """returns a list of tuples with lower and upper range of each index"""
        return [ (i.lower, i.upper) for i in self.indices ]

    def _sympystr(self, p):
        indices = map(p.doprint, self.indices)
        return "%s[%s]" % (p.doprint(self.stem), ", ".join(indices))


class Idx(Basic):
    """Represents an index, either symbolic or integer.

    Optionally you can specify a range [default=0]
    .. Symbol, integer  --  interpreted as dimension, lower and upper ranges are
                            set to 0 and range-1
    .. tuple  --  interpreted as lower, upper elements in range.

    Note that the Idx constructor is rather pedantic, and will not accept
    non-integer symbols.  The only exception is that you can use oo and -oo to
    specify an unbounded range.  For all other cases, both label and bounds
    must be declared as integers, in the sense that for a symbol n,
    n.is_integer must return True.

    For convenience, if the label is given as a string, it is automatically
    converted to an integer symbol.  (Note that this conversion is not done for
    range or dimension arguments.)

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols, oo
    >>> n, i, L, U = symbols('n i L U', integer=True)

    0) Construction from a string

    >>> Idx('qwerty')
    qwerty

    1) Both upper and lower bound specified

    >>> idx = Idx(i, (L, U)); idx
    i
    >>> idx.lower, idx.upper
    (L, U)

    2) Only dimension specified, lower bound defaults to 0

    >>> idx = Idx(i, n); idx.lower, idx.upper
    (0, -1 + n)
    >>> idx = Idx(i, 4); idx.lower, idx.upper
    (0, 3)
    >>> idx = Idx(i, oo); idx.lower, idx.upper
    (0, oo)

    3) No bounds given, interpretation of this depends on context.

    >>> idx = Idx(i); idx.lower, idx.upper
    (None, None)

    4) for a literal integer instead of a symbolic label the bounds are still
    there:

    >>> idx = Idx(2, n); idx.lower, idx.upper
    (0, -1 + n)


    """

    is_integer = True

    def __new__(cls, label, range=None, **kw_args):

        if isinstance(label, basestring):
            label = Symbol(label, integer=True)
        label, range = map(sympify, (label, range))

        if not label.is_integer:
            raise TypeError("Idx object requires an integer label")

        elif isinstance(range, (tuple, list, Tuple)):
            assert len(range) == 2, "Idx got range tuple with wrong length"
            for bound in range:
                if not (bound.is_integer or abs(bound) is S.Infinity):
                    raise TypeError("Idx object requires integer bounds")
            args = label, Tuple(*range)
        elif isinstance(range, (Symbol, Integer)) or range is S.Infinity:
            if not (range.is_integer or range is S.Infinity):
                raise TypeError("Idx object requires an integer dimension")
            args = label, Tuple(S.Zero, range-S.One)
        elif range:
            raise TypeError("range must be tuple, symbol or integer")
        else:
            args = label,

        obj = Basic.__new__(cls, *args, **kw_args)
        return obj

    @property
    def label(self):
        return self.args[0]

    @property
    def lower(self):
        try:
            return self.args[1][0]
        except IndexError:
            return


    @property
    def upper(self):
        try:
            return self.args[1][1]
        except IndexError:
            return

    def _sympystr(self, p):
        return p.doprint(self.label)


"""Module that defines indexed objects with arbitrary transformation properties.

    The Idx class represent indices and each Idx can optionally contain
    information about its range.  Transformation properties can be defined in
    subclasses.

    Examples
    ========

    To express a matrix-vector product in terms of Indexed objects:


    >>> from sympy.tensor import Indexed, Idx
    >>> from sympy import symbols, Eq
    >>> i,j,n,m,M,x,y = symbols('i j n m M x y')
    >>> M = Indexed(M, Idx(i, m), Idx(j, n))
    >>> x = Indexed(x, Idx(j, n))
    >>> y = Indexed(y, Idx(i, m))
    >>> Eq(y, M*x)
    y(i) == M(i, j)*x(j)

"""

from sympy import Expr, Basic, SymTuple, Symbol, Integer, sympify, S

class IndexException(Exception):
    pass

class Indexed(Expr):
    """Represent an arbitrary object with indices, e.g. a symbolic array element.

    >>> from sympy.tensor import Indexed
    >>> from sympy import symbols
    >>> i, j, k = symbols('i j k', integer=True)
    >>> a = symbols('a')
    >>> Indexed(a, i, j)
    a(i, j)

    There is also syntactic sugar that allows declaration of the stem independent
    of the indices:

    >>> A = Indexed(a); A
    a
    >>> A(i, j, k)
    a(i, j, k)

    """

    def __call__(self, *indices, **kw_args):
        if self.rank != 0:
            raise IndexException("Indexed with rank > 0 cannot be used as stem")
        return self.func(self.label, *indices, **kw_args)

    @property
    def label(self):
        return self.args[0]

    @property
    def indices(self):
        return self.args[1:]

    @property
    def rank(self):
        """returns the number of indices"""
        return len(args)-1

    @property
    def dimensions(self):
        return [(i.lower, i.upper) for i in self.indices]

    def _sympystr(self, p):
        indices = map(p.doprint, self.indices)
        if indices:
            return "%s(%s)" % (p.doprint(self.label), ", ".join(indices))
        else:
            return "%s" % p.doprint(self.label)


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

    >>> from sympy.tensor import Indexed, Idx
    >>> from sympy import symbols, oo
    >>> n, i, L, U = symbols('n i L U', integer=True)

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

        label, range = map(sympify, (label, range))

        if not label.is_integer:
            raise TypeError("Idx object requires an integer label")

        elif isinstance(range, (tuple, list, SymTuple)):
            assert len(range) == 2, "Idx got range tuple with wrong length"
            for bound in range:
                if not (bound.is_integer or abs(bound) is S.Infinity):
                    raise TypeError("Idx object requires integer bounds")
            args = label, SymTuple(*range)
        elif isinstance(range, (Symbol, Integer)) or range is S.Infinity:
            if not (range.is_integer or range is S.Infinity):
                raise TypeError("Idx object requires an integer dimension")
            args = label, SymTuple(S.Zero, range-S.One)
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


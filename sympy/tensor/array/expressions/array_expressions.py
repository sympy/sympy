import operator
from functools import reduce
import itertools
from sympy import Expr, ImmutableDenseNDimArray, S, Symbol, Integer
from sympy.core.sympify import _sympify


class _ArrayExpr(Expr):
    pass


class ArraySymbol(_ArrayExpr):
    """
    Symbol representing an array expression
    """

    def __new__(cls, symbol, *shape):
        if isinstance(symbol, str):
            symbol = Symbol(symbol)
        # symbol = _sympify(symbol)
        shape = map(_sympify, shape)
        obj = Expr.__new__(cls, symbol, *shape)
        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def shape(self):
        return self._args[1:]

    def __getitem__(self, item):
        return ArrayElement(self, item)

    def as_explicit(self):
        if any(not isinstance(i, (int, Integer)) for i in self.shape):
            raise ValueError("cannot express explicit array with symbolic shape")
        data = [self[i] for i in itertools.product(*[range(j) for j in self.shape])]
        return ImmutableDenseNDimArray(data).reshape(*self.shape)


class ArrayElement(_ArrayExpr):
    """
    An element of an array.
    """
    def __new__(cls, name, indices):
        if isinstance(name, str):
            name = Symbol(name)
        name = _sympify(name)
        indices = _sympify(indices)
        if hasattr(name, "shape"):
            if any([(i >= s) == True for i, s in zip(indices, name.shape)]):
                raise ValueError("shape is out of bounds")
        if any([(i < 0) == True for i in indices]):
            raise ValueError("shape contains negative values")
        obj = Expr.__new__(cls, name, indices)
        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def indices(self):
        return self._args[1]


class ZeroArray(_ArrayExpr):
    """
    Symbolic array of zeros. Equivalent to ``ZeroMatrix`` for matrices.
    """

    def __new__(cls, *shape):
        if len(shape) == 0:
            return S.Zero
        shape = map(_sympify, shape)
        obj = Expr.__new__(cls, *shape)
        return obj

    @property
    def shape(self):
        return self._args

    def as_explicit(self):
        if any(not i.is_Integer for i in self.shape):
            raise ValueError("Cannot return explicit form for symbolic shape.")
        return ImmutableDenseNDimArray.zeros(*self.shape)


class OneArray(_ArrayExpr):
    """
    Symbolic array of ones.
    """

    def __new__(cls, *shape):
        if len(shape) == 0:
            return S.One
        shape = map(_sympify, shape)
        obj = Expr.__new__(cls, *shape)
        return obj

    @property
    def shape(self):
        return self._args

    def as_explicit(self):
        if any(not i.is_Integer for i in self.shape):
            raise ValueError("Cannot return explicit form for symbolic shape.")
        return ImmutableDenseNDimArray([S.One for i in range(reduce(operator.mul, self.shape))]).reshape(*self.shape)

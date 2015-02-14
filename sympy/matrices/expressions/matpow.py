from __future__ import print_function, division

from .matexpr import MatrixExpr, ShapeError, Identity
from sympy.core.sympify import _sympify
from sympy.core.compatibility import range
from sympy.matrices import MatrixBase
from sympy.core import S


class MatPow(MatrixExpr):

    def __new__(cls, base, exp):
        base = _sympify(base)
        if not base.is_Matrix:
            raise TypeError("Function parameter should be a matrix")
        exp = _sympify(exp)
        return super(MatPow, cls).__new__(cls, base, exp)

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    @property
    def shape(self):
        return self.base.shape

    def _entry(self, i, j):
        if self.exp.is_zero:
            if not self.base.is_square:
                raise ShapeError("Power of non-square matrix %s" % self.base)
            T = Identity(self.base.shape[0])
        elif self.exp is S.One:
            T = self.base
        elif isinstance(self.base, MatrixBase):
            # e.g., ImmutableMatrix**S.Half not covered by MatMul cases below
            T = self.base**self.exp
        elif self.exp.is_Integer and self.exp.is_positive:
            # Make an explicit MatMul out of the MatPow
            T = MatMul(*[self.base for k in range(self.exp)])
        #elif self.exp.is_Integer and self.exp.is_negative:
        #    # Note: possible future improvement: in principle we can take
        #    # positive powers of the inverse, but carefully avoid recursion,
        #    # perhaps by adding `_entry` to Inverse (as it is our subclass).
        #    T = self.base.inverse()
        #    T = MatMul(*[T for k in range(-self.exp)])
        else:
            raise NotImplementedError(("(%d, %d) entry" % (int(i), int(j))) +
                    "of matrix power either not defined or not implemented")
        return T._entry(i, j)


    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        base = args[0]
        exp = args[1]
        if isinstance(base, MatrixBase) and exp.is_number:
            if exp is S.One:
                return base
            return base**exp
        # Note: just evaluate cases we know, return unevaluated on others.
        # E.g., MatrixSymbol('x', n, m) to power 0 is not an error.
        if exp.is_zero and base.is_square:
            return Identity(base.shape[0])
        elif exp is S.One:
            return base
        return MatPow(base, exp)


from .matmul import MatMul

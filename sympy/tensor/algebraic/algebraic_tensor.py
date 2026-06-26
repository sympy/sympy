from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.sympify import sympify


class PureTensor(Mul):
    """Pure tensor as an unevaluated (non-commutative) tensor product of factors.

    Extends SymPy's non-commutative Mul so that all existing simplification,
    differentiation, and pattern-matching machinery is available out of the
    box.  The tensor-product operator is non-commutative: order of factors
    is always preserved.

    Each factor must carry ``.shape`` (e.g. any ``MatrixExpr`` or a 1-x1
    wrapper around a non-commutative Symbol).
    """

    __slots__ = ()

    is_PureTensor = True

    _eval_is_commutative = lambda self: False

    @property
    def factors(self):
        """Individual tensor-product factors in left-to-right order."""
        return self.args

    @property
    def num_factors(self):
        return len(self.args)

    @property
    def tensor_shape(self):
        """Outer-product shape: (dim_0-of-first, dim_1-of-last)."""
        shapes = [f.shape for f in self.factors]
        return (shapes[0][0], shapes[-1][1])

    def __str__(self):
        return " \u2297 ".join(str(f) for f in self.factors)

    def __repr__(self):
        return f"PureTensor({', '.join(repr(f) for f in self.factors)})"

    def __new__(cls, *args, evaluate=False):
        if not args:
            raise ValueError("PureTensor requires at least one factor")

        processed = []
        for a in args:
            a = sympify(a)
            if isinstance(a, Number):
                raise TypeError(
                    f"Scalar numbers are not valid tensor factors "
                    f"(use a 1x1 MatrixExpr wrapper): {a}"
                )
            if not hasattr(a, "shape"):
                raise TypeError(
                    f"Tensor factor must have a .shape attribute, "
                    f"got {type(a).__name__}"
                )
            processed.append(a)

        if len(processed) == 1:
            return processed[0]

        obj = Mul.__new__(cls, *processed, evaluate=False)
        return obj


def tensor_product(*args):
    """Convenience constructor for PureTensor.

    >>> from sympy.tensor.algebraic.algebraic_tensor import tensor_product
    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> A = MatrixSymbol("A", 2, 3)
    >>> v = MatrixSymbol("v", 3, 1)
    >>> tensor_product(A, v)
    A \u2297 v
    """
    return PureTensor(*args)

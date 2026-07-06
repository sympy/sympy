from __future__ import annotations

from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify


class ScalarMul(Basic):
    """Scalar multiplication of an AlgebraicPureTensor.

    Represents ``scalar * tensor`` where *scalar* is a commutative SymPy
    expression and *tensor* is an ``AlgebraicPureTensor``.  Used instead of
    ``Mul`` to avoid SymPy's deprecation of ``Mul`` with non-Expr factors
    (AlgebraicPureTensor extends Mul but carries matrix factors).

    Extends ``Basic`` (not ``Mul``) to prevent SymPy's flattening machinery
    from unpacking the tensor factors.
    """

    __slots__ = ()

    is_ScalarMul = True
    is_commutative = False
    is_Mul = False
    is_Add = False

    _op_priority = 12  # Higher than AlgebraicPureTensor so nesting works

    def __new__(cls, scalar, tensor):
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

        scalar = sympify(scalar)
        tensor = sympify(tensor)

        if scalar is S.Zero:
            return AlgebraicZeroTensor(tensor.tensor_shape)

        if scalar is S.One:
            return tensor

        if isinstance(tensor, AlgebraicPureTensor):
            inner_coeff = tensor._get_coeff()
            combined = scalar * inner_coeff
            if combined is S.One:
                factors = tensor.factors
                if len(factors) == 1:
                    return factors[0]
                return AlgebraicPureTensor(*factors)
            if combined is S.Zero:
                return AlgebraicZeroTensor(tensor.tensor_shape)
            if combined is S.NegativeOne:
                factors = tensor.factors
                if len(factors) == 1:
                    return ScalarMul(S.NegativeOne, factors[0])
                return AlgebraicPureTensor(S.NegativeOne, *factors)
            # Try to absorb into PureTensor if combined is a simple Number
            if isinstance(combined, Number):
                factors = tensor.factors
                if len(factors) == 1:
                    return AlgebraicPureTensor(combined, factors[0])
                return AlgebraicPureTensor(combined, *factors)

        # General case: keep as ScalarMul
        obj = Basic.__new__(cls, scalar, tensor)
        return obj

    @property
    def scalar(self):
        """The commutative scalar factor."""
        return self.args[0]

    @property
    def tensor(self):
        """The AlgebraicPureTensor factor."""
        return self.args[1]

    @property
    def factors(self):
        """Tensor-product factors (same as the inner tensor's factors)."""
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        if isinstance(self.tensor, AlgebraicPureTensor):
            return self.tensor.factors
        return (self.tensor,)

    @property
    def tensor_shape(self):
        """Full tensor shape as a tuple of per-factor (rows, cols) pairs."""
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        if isinstance(self.tensor, AlgebraicPureTensor):
            return self.tensor.tensor_shape
        if hasattr(self.tensor, "shape"):
            return (self.tensor.shape,)
        raise AttributeError("Cannot determine tensor_shape")

    @property
    def commutativity_shape(self):
        """Same as the inner tensor's commutativity_shape."""
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        if isinstance(self.tensor, AlgebraicPureTensor):
            return self.tensor.commutativity_shape
        if hasattr(self.tensor, "free_symbols"):
            has_nc = any(
                getattr(s, "is_commutative", True) is False
                for s in getattr(self.tensor, "free_symbols", set())
            )
            return (0 if has_nc else 1,)
        return (0,)

    def _get_coeff(self):
        """Extract the combined scalar coefficient."""
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        scalar_part = self.scalar
        if isinstance(self.tensor, AlgebraicPureTensor):
            return scalar_part * self.tensor._get_coeff()
        return scalar_part

    def __neg__(self):
        """Negate the scalar part."""
        new_scalar = -self.scalar
        if new_scalar is S.One:
            return self.tensor
        if new_scalar is S.Zero:
            from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
            return AlgebraicZeroTensor(self.tensor_shape)
        return ScalarMul(new_scalar, self.tensor)

    def __mul__(self, other):
        """Compose or scale this ScalarMul."""
        from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors

        other = sympify(other)

        if isinstance(other, Number) or (
            hasattr(other, "is_commutative") and other.is_commutative
        ):
            new_scalar = self.scalar * other
            if new_scalar is S.One:
                return self.tensor
            if new_scalar is S.Zero:
                from sympy.tensor.algebraic.algebraic_zero_tensor import (
                    AlgebraicZeroTensor,
                )

                return AlgebraicZeroTensor(self.tensor_shape)
            return ScalarMul(new_scalar, self.tensor)

        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        """Compose or scale this ScalarMul from the left."""
        from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors

        other = sympify(other)

        if isinstance(other, Number) or (
            hasattr(other, "is_commutative") and other.is_commutative
        ):
            new_scalar = other * self.scalar
            if new_scalar is S.One:
                return self.tensor
            if new_scalar is S.Zero:
                from sympy.tensor.algebraic.algebraic_zero_tensor import (
                    AlgebraicZeroTensor,
                )

                return AlgebraicZeroTensor(self.tensor_shape)
            return ScalarMul(new_scalar, self.tensor)

        return compose_algebraic_tensors(other, self)

    def __add__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(other, -self)

    def expand(self, deep=True, **hints):
        """Expand the inner tensor, carrying the scalar through."""
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        if hasattr(self.tensor, 'expand'):
            expanded = self.tensor.expand(deep=deep, **hints)
        else:
            expanded = self.tensor

        if expanded is self.tensor:
            return self

        # If expansion produced an AlgebraicTensor, distribute the scalar
        if isinstance(expanded, AlgebraicTensor):
            return AlgebraicTensor(*(self.scalar * a for a in expanded.args))

        # Single expanded tensor
        return ScalarMul(self.scalar, expanded)

    def simplify(self):
        """Simplify this ScalarMul."""
        from sympy.tensor.algebraic.simplify import tensorsimplify
        return tensorsimplify(self)

    def display(self, mode="latex"):
        """Display this tensor using IPython display or fallback to print."""
        try:
            from IPython.display import display, Latex
            if mode == "latex":
                display(Latex(self._repr_latex_()))
            else:
                display(self, plain=True)
        except ImportError:
            if mode == "latex":
                print(self._repr_latex_())
            else:
                print(self)

    def __str__(self):
        return f"{self.scalar}*{self.tensor}"

    def __repr__(self):
        return f"ScalarMul({repr(self.scalar)}, {repr(self.tensor)})"

from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.zero_tensor import ZeroTensor


def _factor_shapes(factors):
    """Return the full tensor shape as a tuple of per-factor shapes.

    E.g. for factors [A(3,4), C(4,5)] returns ``((3, 4), (4, 5))``.
    """
    return tuple(f.shape for f in factors)


class PureTensor(Mul):
    """Pure tensor as an unevaluated (non-commutative) tensor product of factors.

    Extends SymPy's non-commutative Mul so that all existing simplification,
    differentiation, and pattern-matching machinery is available out of the
    box.  The tensor-product operator is non-commutative: order of factors
    is always preserved.

    Each factor must carry ``.shape`` (e.g. any ``MatrixExpr`` or a 1-x1
    wrapper around a non-commutative Symbol).

    The tensor shape is the full sequence of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for ``PureTensor(A_3x4, C_4x5)``.  No
    contraction is performed.
    """

    __slots__ = ()

    is_PureTensor = True

    _eval_is_commutative = lambda self: False

    @property
    def factors(self):
        """Individual tensor-product factors in left-to-right order.

        Excludes a leading coefficient (Number or commutative symbol) if
        one is stored as the first arg.
        """
        if self.args and (self.args[0].is_Number or
                (hasattr(self.args[0], 'is_commutative') and
                 self.args[0].is_commutative and
                 not isinstance(self.args[0], PureTensor))):
            return self.args[1:]
        return self.args

    @property
    def num_factors(self):
        return len(self.args)

    @property
    def tensor_shape(self):
        """Full tensor shape as a tuple of per-factor (rows, cols) pairs.

        E.g. ``PureTensor(A_3x4, C_4x5).tensor_shape == ((3, 4), (4, 5))``.
        """
        return _factor_shapes(self.factors)

    def __str__(self):
        coeff = self._get_coeff()
        factor_str = " \u2297 ".join(str(f) for f in self.factors)
        if coeff is S.One:
            return factor_str
        if coeff is S.NegativeOne:
            return f"-{factor_str}"
        return f"{coeff}*{factor_str}"

    def __repr__(self):
        coeff = self._get_coeff()
        if coeff is S.One:
            return f"PureTensor({', '.join(repr(f) for f in self.factors)})"
        return f"PureTensor({repr(coeff)}, {', '.join(repr(f) for f in self.factors)})"

    def __new__(cls, *args, evaluate=False):
        if not args:
            raise ValueError("PureTensor requires at least one factor")

        # Separate leading coefficient from tensor factors
        args_list = list(args)
        coeff = S.One

        if args_list and isinstance(args_list[0], Number):
            coeff = args_list[0]
            args_list = args_list[1:]
        elif args_list and len(args_list) > 1 and hasattr(args_list[0], 'is_commutative') \
                and args_list[0].is_commutative and not isinstance(args_list[0], PureTensor):
            coeff = sympify(args_list[0])
            args_list = args_list[1:]

        processed = []
        for a in args_list:
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

        if not processed:
            if coeff is S.Zero:
                return ZeroTensor(())
            raise ValueError("PureTensor requires at least one tensor factor")

        if coeff is S.Zero:
            return ZeroTensor(_factor_shapes(processed))

        if len(processed) == 1 and coeff is S.One:
            return processed[0]

        if coeff is S.One:
            obj = Mul.__new__(cls, *processed, evaluate=False)
        else:
            obj = Mul.__new__(cls, coeff, *processed, evaluate=False)
        return obj

    def __neg__(self):
        """Return a PureTensor with negated coefficient."""
        coeff = self._get_coeff()
        factors = self.factors
        new_coeff = coeff * S.NegativeOne
        if new_coeff is S.One:
            if len(factors) == 1:
                return factors[0]
            return PureTensor(*factors)
        if new_coeff is S.NegativeOne:
            if len(factors) == 1:
                return Mul(S.NegativeOne, factors[0], evaluate=False)
            return PureTensor(S.NegativeOne, *factors)
        return PureTensor(new_coeff, *factors)

    def __mul__(self, other):
        """PureTensor * scalar/symbol -> PureTensor with absorbed coefficient."""
        from sympy.core.singleton import S
        other = sympify(other)
        if isinstance(other, Number):
            if other is S.One:
                return self
            if other is S.Zero:
                return ZeroTensor(self.tensor_shape)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, PureTensor)):
            coeff, factors = self._get_coeff(), self.factors
            new_coeff = coeff * other
            if new_coeff is S.One:
                if len(factors) == 1:
                    return factors[0]
                return PureTensor(*factors)
            if isinstance(new_coeff, Number) and new_coeff is S.Zero:
                return ZeroTensor(_factor_shapes(factors))
            if len(factors) == 0:
                return new_coeff
            return PureTensor(new_coeff, *factors)
        if isinstance(other, PureTensor):
            c1 = self._get_coeff()
            c2 = other._get_coeff()
            combined = c1 * c2
            all_factors = self.factors + other.factors
            if combined is S.One:
                if len(all_factors) == 1:
                    return all_factors[0]
                return PureTensor(*all_factors)
            return PureTensor(combined, *all_factors)
        return Mul(self, other)

    def _get_coeff(self):
        """Extract the leading coefficient from args, defaulting to S.One."""
        from sympy.core.singleton import S
        if self.args and (self.args[0].is_Number or
                (hasattr(self.args[0], 'is_commutative') and
                 self.args[0].is_commutative and
                 not isinstance(self.args[0], PureTensor))):
            return self.args[0]
        return S.One

    def __rmul__(self, other):
        """scalar/symbol * PureTensor -> PureTensor with absorbed coefficient."""
        from sympy.core.singleton import S
        if other == 0:
            return ZeroTensor(self.tensor_shape)
        if other == 1:
            return self
        other = sympify(other)
        if isinstance(other, Number):
            if other is S.One:
                return self
            if other is S.Zero:
                return ZeroTensor(self.tensor_shape)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, PureTensor)):
            coeff = self._get_coeff()
            factors = self.factors
            new_coeff = other * coeff
            if new_coeff is S.One:
                if len(factors) == 1:
                    return factors[0]
                return PureTensor(*factors)
            if isinstance(new_coeff, Number) and new_coeff is S.Zero:
                return ZeroTensor(_factor_shapes(factors))
            if len(factors) == 0:
                return new_coeff
            return PureTensor(new_coeff, *factors)
        if isinstance(other, PureTensor):
            c1 = other._get_coeff()
            c2 = self._get_coeff()
            combined = c1 * c2
            all_factors = other.factors + self.factors
            if combined is S.One:
                if len(all_factors) == 1:
                    return all_factors[0]
                return PureTensor(*all_factors)
            return PureTensor(combined, *all_factors)
        return Mul(other, self)

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

    def simplify(self, **kwargs):
        """Simplify this PureTensor.

        Simplifies the leading coefficient and each tensor factor using
        SymPy's general :func:`simplify`.  See
        :func:`sympy.tensor.algebraic.simplify.tensorsimplify` for details.
        """
        from sympy.tensor.algebraic.simplify import tensorsimplify
        return tensorsimplify(self, **kwargs)


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

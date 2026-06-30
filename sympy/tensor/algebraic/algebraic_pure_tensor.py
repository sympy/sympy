from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor


def _factor_shapes(factors):
    """Return the full tensor shape as a tuple of per-factor shapes.

    E.g. for factors [A(3,4), C(4,5)] returns ``((3, 4), (4, 5))``.
    """
    return tuple(f.shape for f in factors)


class AlgebraicPureTensor(Mul):
    """Pure tensor as an unevaluated (non-commutative) tensor product of factors.

    Extends SymPy's non-commutative Mul so that all existing simplification,
    differentiation, and pattern-matching machinery is available out of the
    box.  The tensor-product operator is non-commutative: order of factors
    is always preserved.

    Each factor must carry ``.shape`` (e.g. any ``MatrixExpr`` or a 1-x1
    wrapper around a non-commutative Symbol).

    The tensor shape is the full sequence of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for ``AlgebraicPureTensor(A_3x4, C_4x5)``.  No
    contraction is performed.
    """

    __slots__ = ()

    is_AlgebraicPureTensor = True
    is_Mul = False  # Prevent Mul.flatten from unpacking AlgebraicPureTensor

    _op_priority = 11  # Higher than Symbol/Expr so x * pt delegates to pt.__rmul__(x)

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
                 not isinstance(self.args[0], AlgebraicPureTensor))):
            return self.args[1:]
        return self.args

    @property
    def num_factors(self):
        return len(self.args)

    @property
    def tensor_shape(self):
        """Full tensor shape as a tuple of per-factor (rows, cols) pairs.

        E.g. ``AlgebraicPureTensor(A_3x4, C_4x5).tensor_shape == ((3, 4), (4, 5))``.
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
            return f"AlgebraicPureTensor({', '.join(repr(f) for f in self.factors)})"
        return f"AlgebraicPureTensor({repr(coeff)}, {', '.join(repr(f) for f in self.factors)})"

    def __new__(cls, *args, evaluate=False):
        if not args:
            raise ValueError("AlgebraicPureTensor requires at least one factor")

        # Separate leading coefficient from tensor factors
        args_list = list(args)
        coeff = S.One

        if args_list:
            first = args_list[0]
            if isinstance(first, Number):
                coeff = first
                args_list = args_list[1:]
            elif len(args_list) > 1:
                first_s = sympify(first)
                if isinstance(first_s, Number) or (
                    hasattr(first_s, 'is_commutative') and
                    first_s.is_commutative and
                    not isinstance(first_s, AlgebraicPureTensor)
                ):
                    coeff = first_s
                    args_list = args_list[1:]
            elif len(args_list) > 1:
                first_s = sympify(first)
                if isinstance(first_s, Number) or (
                    hasattr(first_s, 'is_commutative') and
                    first_s.is_commutative and
                    not isinstance(first_s, AlgebraicPureTensor)
                ):
                    coeff = first_s
                    args_list = args_list[1:]
            elif len(args_list) > 1:
                first_s = sympify(first)
                if isinstance(first_s, Number) or (
                    hasattr(first_s, 'is_commutative') and
                    first_s.is_commutative and
                    not isinstance(first_s, AlgebraicPureTensor)
                ):
                    coeff = first_s
                    args_list = args_list[1:]
            elif len(args_list) > 1:
                first_s = sympify(first)
                if isinstance(first_s, Number) or (
                    hasattr(first_s, 'is_commutative') and
                    first_s.is_commutative and
                    not isinstance(first_s, AlgebraicPureTensor)
                ):
                    coeff = first_s
                    args_list = args_list[1:]
            elif len(args_list) > 1:
                first_s = sympify(first)
                if isinstance(first_s, Number) or (
                    hasattr(first_s, 'is_commutative') and
                    first_s.is_commutative and
                    not isinstance(first_s, AlgebraicPureTensor)
                ):
                    coeff = first_s
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
                return AlgebraicZeroTensor(())
            raise ValueError("AlgebraicPureTensor requires at least one tensor factor")

        if coeff is S.Zero:
            return AlgebraicZeroTensor(_factor_shapes(processed))

        if len(processed) == 1 and coeff is S.One:
            return processed[0]

        if coeff is S.One:
            obj = Mul.__new__(cls, *processed, evaluate=False)
        else:
            obj = Mul.__new__(cls, coeff, *processed, evaluate=False)
        return obj

    def __neg__(self):
        """Return an AlgebraicPureTensor with negated coefficient."""
        coeff = self._get_coeff()
        factors = self.factors
        new_coeff = coeff * S.NegativeOne
        if new_coeff is S.One:
            if len(factors) == 1:
                return factors[0]
            return AlgebraicPureTensor(*factors)
        if new_coeff is S.NegativeOne:
            if len(factors) == 1:
                return Mul(S.NegativeOne, factors[0], evaluate=False)
            return AlgebraicPureTensor(S.NegativeOne, *factors)
        return AlgebraicPureTensor(new_coeff, *factors)

    def __mul__(self, other):
        """Compose or scale this AlgebraicPureTensor.

        For commutative scalars/symbols the scalar is absorbed as a
        coefficient.  For AlgebraicPureTensor, AlgebraicTensor, or bare
        matrices the result is the tensor composition (factor-wise matrix
        multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The right operand.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor, or Mul
            The scaled/composed result.
        """
        from sympy.core.singleton import S
        other = sympify(other)
        if isinstance(other, Number):
            if other is S.One:
                return self
            if other is S.Zero:
                return AlgebraicZeroTensor(self.tensor_shape)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, AlgebraicPureTensor)):
            coeff, factors = self._get_coeff(), self.factors
            new_coeff = coeff * other
            if new_coeff is S.One:
                if len(factors) == 1:
                    return factors[0]
                return AlgebraicPureTensor(*factors)
            if isinstance(new_coeff, Number) and new_coeff is S.Zero:
                return AlgebraicZeroTensor(_factor_shapes(factors))
            if len(factors) == 0:
                return new_coeff
            return AlgebraicPureTensor(new_coeff, *factors)
        # Non-commutative operand: use tensor composition.
        from sympy.tensor.algebraic.algebraic_tensor import (
            compose_algebraic_tensors,
        )
        return compose_algebraic_tensors(self, other)

    def _get_coeff(self):
        """Extract the leading coefficient from args, defaulting to S.One."""
        from sympy.core.singleton import S
        if self.args and (self.args[0].is_Number or
                (hasattr(self.args[0], 'is_commutative') and
                 self.args[0].is_commutative and
                 not isinstance(self.args[0], AlgebraicPureTensor))):
            return self.args[0]
        return S.One

    def __rmul__(self, other):
        """Compose or scale this AlgebraicPureTensor from the left.

        For commutative scalars/symbols the scalar is absorbed as a
        coefficient.  For AlgebraicPureTensor, AlgebraicTensor, or bare
        matrices the result is the tensor composition (factor-wise matrix
        multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The left operand.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor, or Mul
            The scaled/composed result.
        """
        from sympy.core.singleton import S
        if other == 0:
            return AlgebraicZeroTensor(self.tensor_shape)
        if other == 1:
            return self
        other = sympify(other)
        if isinstance(other, Number):
            if other is S.One:
                return self
            if other is S.Zero:
                return AlgebraicZeroTensor(self.tensor_shape)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, AlgebraicPureTensor)):
            coeff = self._get_coeff()
            factors = self.factors
            new_coeff = other * coeff
            if new_coeff is S.One:
                if len(factors) == 1:
                    return factors[0]
                return AlgebraicPureTensor(*factors)
            if isinstance(new_coeff, Number) and new_coeff is S.Zero:
                return AlgebraicZeroTensor(_factor_shapes(factors))
            if len(factors) == 0:
                return new_coeff
            return AlgebraicPureTensor(new_coeff, *factors)
        # Non-commutative operand: use tensor composition.
        from sympy.tensor.algebraic.algebraic_tensor import (
            compose_algebraic_tensors,
        )
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

    def simplify(self):
        from sympy.tensor.algebraic.simplify import _simplify_algebraic_pure_tensor
        return _simplify_algebraic_pure_tensor(self)

    def _eval_expand_mul(self, **hints):
        from sympy.core.add import Add
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        from itertools import product

        deep = hints.pop('deep', True)
        coeff = self._get_coeff()
        factors = self.factors

        expanded_factors = []
        for f in factors:
            if hasattr(f, 'expand'):
                expanded_factors.append(f.expand(deep=deep, **hints))
            else:
                expanded_factors.append(f)

        add_slots = [i for i, f in enumerate(expanded_factors)
                     if isinstance(f, Add)]

        if not add_slots:
            if expanded_factors == list(factors):
                return self
            if coeff is S.One:
                if len(expanded_factors) == 1:
                    return expanded_factors[0]
                return AlgebraicPureTensor(*expanded_factors)
            if not expanded_factors:
                return coeff
            return AlgebraicPureTensor(coeff, *expanded_factors)

        choices = []
        for f in expanded_factors:
            if isinstance(f, Add):
                choices.append(f.args)
            else:
                choices.append((f,))

        terms = []
        for combination in product(*choices):
            if coeff is S.One:
                if len(combination) == 1:
                    terms.append(combination[0])
                else:
                    terms.append(AlgebraicPureTensor(*combination))
            else:
                terms.append(AlgebraicPureTensor(coeff, *combination))

        if len(terms) == 1:
            return terms[0]
        return AlgebraicTensor(*terms)


def algebraic_tensor_product(*args):
    """Convenience constructor for AlgebraicPureTensor.

    >>> from sympy.tensor.algebraic.algebraic_pure_tensor import algebraic_tensor_product
    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> A = MatrixSymbol("A", 2, 3)
    >>> v = MatrixSymbol("v", 3, 1)
    >>> algebraic_tensor_product(A, v)
    A ⊗ v
    """
    return AlgebraicPureTensor(*args)


def compose_algebraic_pure_tensors(left, right):
    """Compose two AlgebraicPureTensors by matrix-multiplying corresponding factors.

    Given ``left = F1 ⊗ F2 ⊗ … ⊗ Fn`` and ``right = G1 ⊗ G2 ⊗ … ⊗ Gn``,
    returns ``H1 ⊗ H2 ⊗ … ⊗ Hn`` where each ``Hj = Fj * Gj`` (matrix product).

    For each factor pair ``(Fj, Gj)`` the column dimension of ``Fj`` must equal
    the row dimension of ``Gj``.  If ``Fj`` has shape ``(a_j, b_j)`` then ``Gj``
    must have shape ``(b_j, c_j)`` so that ``Hj`` has shape ``(a_j, c_j)``.

    Parameters
    ----------
    left : AlgebraicPureTensor
        The left operand.  A bare matrix-like object (single-factor tensor)
        is also accepted.
    right : AlgebraicPureTensor
        The right operand.  A bare matrix-like object is also accepted.

    Returns
    -------
    AlgebraicPureTensor or matrix-like object
        The composed tensor.  If the result has a single factor and
        coefficient 1, the bare factor is returned.

    Raises
    ------
    TypeError
        If either argument is not an ``AlgebraicPureTensor`` or a bare
        matrix-like object with a ``.shape`` attribute.
    ValueError
        If the two tensors have a different number of factors, or if any
        corresponding factor pair has incompatible inner dimensions for
        matrix multiplication.
    """
    # Handle unwrapped single-factor tensors
    if isinstance(left, AlgebraicPureTensor):
        left_factors = left.factors
        left_coeff = left._get_coeff()
    elif hasattr(left, "shape"):
        left_factors = (left,)
        left_coeff = S.One
    else:
        raise TypeError(
            f"Expected AlgebraicPureTensor or a matrix-like object on the left, "
            f"got {type(left).__name__}"
        )

    if isinstance(right, AlgebraicPureTensor):
        right_factors = right.factors
        right_coeff = right._get_coeff()
    elif hasattr(right, "shape"):
        right_factors = (right,)
        right_coeff = S.One
    else:
        raise TypeError(
            f"Expected AlgebraicPureTensor or a matrix-like object on the right, "
            f"got {type(right).__name__}"
        )

    if len(left_factors) != len(right_factors):
        raise ValueError(
            f"Cannot compose AlgebraicPureTensors with different numbers of "
            f"factors: {len(left_factors)} vs {len(right_factors)}"
        )

    # Combine coefficients from both sides
    combined_coeff = left_coeff * right_coeff

    composed_factors = []
    for j, (lf, rf) in enumerate(zip(left_factors, right_factors)):
        lshape = lf.shape
        rshape = rf.shape
        if lshape[1] != rshape[0]:
            raise ValueError(
                f"Cannot compose factor {j}: left factor has shape {lshape} "
                f"but right factor has shape {rshape}; inner dimensions "
                f"{lshape[1]} and {rshape[0]} do not match"
            )
        # Use MatMul for the matrix product of corresponding factors
        from sympy.matrices.expressions.matexpr import MatMul
        composed_factors.append(MatMul(lf, rf, evaluate=False))

    if combined_coeff is S.One:
        if len(composed_factors) == 1:
            return composed_factors[0]
        return AlgebraicPureTensor(*composed_factors)
    if combined_coeff is S.Zero:
        return AlgebraicZeroTensor(_factor_shapes(composed_factors))
    if len(composed_factors) == 1:
        return AlgebraicPureTensor(combined_coeff, composed_factors[0])
    return AlgebraicPureTensor(combined_coeff, *composed_factors)

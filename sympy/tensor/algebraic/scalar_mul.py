from __future__ import annotations

from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify


"""Scalar multiplication of an AlgebraicPureTensor by a symbolic factor.

This module defines :class:`ScalarMul`, which represents ``scalar * tensor``
where *scalar* is a commutative SymPy expression (typically a Symbol or
symbolic expression, **not** a plain :class:`~sympy.core.numbers.Number`)
and *tensor* is an
:class:`~sympy.tensor.algebraic.algebraic_pure_tensor.AlgebraicPureTensor`.

It exists to avoid SymPy's deprecation of ``Mul`` with non-``Expr`` factors:
:class:`~sympy.tensor.algebraic.algebraic_pure_tensor.AlgebraicPureTensor`
extends ``Mul`` but carries matrix factors that are not plain ``Expr``
objects in all contexts.

Examples
========

Create a scalar multiplication with a symbolic coefficient:

>>> from sympy.abc import x
>>> from sympy.matrices.expressions import MatrixSymbol
>>> from sympy.tensor.algebraic import AlgebraicPureTensor
>>> A = MatrixSymbol("A", 3, 4)
>>> B = MatrixSymbol("B", 4, 5)
>>> T = AlgebraicPureTensor(A, B)
>>> S = ScalarMul(x, T)
>>> S
x*A ⊗ B
>>> S.scalar
x
>>> S.tensor
A ⊗ B

Numeric coefficients are absorbed into the PureTensor:

>>> ScalarMul(2, T)
2*A ⊗ B

Zero scalar produces a zero tensor:

>>> ScalarMul(0, T)
0_((3, 4), (4, 5))
"""


class ScalarMul(Basic):
    """Scalar multiplication of an AlgebraicPureTensor.

    Represents ``scalar * tensor`` where *scalar* is a commutative SymPy
    expression and *tensor* is an ``AlgebraicPureTensor``.  Used instead of
    ``Mul`` to avoid SymPy's deprecation of ``Mul`` with non-Expr factors
    (AlgebraicPureTensor extends Mul but carries matrix factors).

    Extends ``Basic`` (not ``Mul``) to prevent SymPy's flattening machinery
    from unpacking the tensor factors.

    Examples
    ========

    Create a scalar multiplication with a symbolic coefficient:

    >>> from sympy.abc import x
    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> ScalarMul(x, T)
    x*A ⊗ B

    Numeric coefficients are absorbed into the PureTensor:

    >>> ScalarMul(2, T)
    2*A ⊗ B

    Zero scalar produces a zero tensor:

    >>> ScalarMul(0, T)
    0_((3, 4), (4, 5))

    Unit scalar unwraps to the bare tensor:

    >>> ScalarMul(1, T)
    A ⊗ B
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

        if isinstance(tensor, AlgebraicZeroTensor):
            return tensor

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
        """The commutative scalar factor.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T).scalar
        x
        """
        return self.args[0]

    @property
    def tensor(self):
        """The AlgebraicPureTensor factor.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T).tensor
        A ⊗ B
        """
        return self.args[1]

    @property
    def factors(self):
        """Tensor-product factors (same as the inner tensor's factors).

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T).factors
        (A, B)
        """
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        if isinstance(self.tensor, AlgebraicPureTensor):
            return self.tensor.factors
        return (self.tensor,)

    @property
    def tensor_shape(self):
        """Full tensor shape as a tuple of per-factor (rows, cols) pairs.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T).tensor_shape
        ((3, 4), (4, 5))
        """
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        if isinstance(self.tensor, AlgebraicPureTensor):
            return self.tensor.tensor_shape
        if hasattr(self.tensor, "shape"):
            return (self.tensor.shape,)
        raise AttributeError("Cannot determine tensor_shape")

    @property
    def commutativity_shape(self):
        """Same as the inner tensor's commutativity_shape.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T).commutativity_shape
        (0, 0)
        """
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
        """Negate the scalar part.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> -ScalarMul(x, T)
        -x*A ⊗ B
        """
        new_scalar = -self.scalar
        if new_scalar is S.One:
            return self.tensor
        if new_scalar is S.Zero:
            from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
            return AlgebraicZeroTensor(self.tensor_shape)
        return ScalarMul(new_scalar, self.tensor)

    def __mul__(self, other):
        """Compose or scale this ScalarMul.

        Examples
        ========

        Scalar multiplication:

        >>> from sympy.abc import x, y
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> S = ScalarMul(x, T)
        >>> S * y
        x*y*A ⊗ B
        >>> S * 0
        0_((3, 4), (4, 5))
        """
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
        """Compose or scale this ScalarMul from the left.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> S = ScalarMul(x, T)
        >>> y * S
        x*y*A ⊗ B
        """
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
        """Add another tensor expression to this ScalarMul.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T) + AlgebraicPureTensor(C, D)
        x*A ⊗ B + C ⊗ D
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        """Right-add: ``other + self``.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> AlgebraicPureTensor(C, D) + ScalarMul(x, T)
        C ⊗ D + x*A ⊗ B
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        """Subtract another tensor expression from this ScalarMul.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x, T) - AlgebraicPureTensor(C, D)
        x*A ⊗ B - C ⊗ D
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        """Right-subtract: ``other - self``.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> AlgebraicPureTensor(C, D) - ScalarMul(x, T)
        C ⊗ D - x*A ⊗ B
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        return AlgebraicTensor(other, -self)

    def expand(self, deep=True, **hints):
        """Expand ScalarMul by distributing over Add in scalar and tensor factors.

        1. Expand the symbolic prefactor as a SymPy Expr.
        2. Expand every tensor factor.
        3. If the scalar prefactor is an Add, split into an Add of ScalarMuls.
        4. For each ScalarMul, if any tensor factor is an Add, distribute the
           tensor product linearly into multiple AlgebraicPureTensors/ScalarMuls.

        Examples
        ========

        Distribute over an Add in a tensor factor:

        >>> from sympy.abc import x
        >>> from sympy.matrices.expressions import (
        ...     MatrixSymbol, MatAdd)
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 4, 5)
        >>> T = AlgebraicPureTensor(A, MatAdd(B, C))
        >>> ScalarMul(x, T).expand()
        x*A ⊗ B + x*A ⊗ C

        Distribute over an Add in the scalar:

        >>> from sympy.abc import y
        >>> T2 = AlgebraicPureTensor(A, B)
        >>> ScalarMul(x + y, T2).expand()
        x*A ⊗ B + y*A ⊗ B
        """
        from sympy.core.add import Add
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        from itertools import product

        # Step 1: Expand the scalar prefactor
        scalar = self.scalar
        if hasattr(scalar, 'expand'):
            scalar = scalar.expand(deep=deep, **hints)

        # Step 2: Expand every tensor factor
        factors = list(self.factors)
        expanded_factors = []
        for f in factors:
            if hasattr(f, 'expand'):
                expanded_factors.append(f.expand(deep=deep, **hints))
            else:
                expanded_factors.append(f)

        # Check if anything actually changed and there's nothing to distribute
        if scalar is self.scalar and expanded_factors == factors:
            if not isinstance(scalar, Add) and not any(
                isinstance(f, Add) for f in expanded_factors
            ):
                return self

        # Step 3: If scalar is an Add, split into multiple ScalarMuls
        if isinstance(scalar, Add):
            terms = []
            for s in scalar.args:
                sm = ScalarMul(s, AlgebraicPureTensor(*expanded_factors))
                # Recursively expand each resulting ScalarMul in case tensor
                # factors also contain Adds
                if hasattr(sm, 'expand'):
                    sm = sm.expand(deep=deep, **hints)
                terms.append(sm)
            if len(terms) == 1:
                return terms[0]
            return AlgebraicTensor(*terms)

        # Step 4: Distribute tensor product over Add in factors
        add_slots = [i for i, f in enumerate(expanded_factors)
                     if isinstance(f, Add)]

        if not add_slots:
            # No Add in factors; build single term
            inner = AlgebraicPureTensor(*expanded_factors)
            if scalar is S.One:
                return inner
            return ScalarMul(scalar, inner)

        # Build cartesian product of all factor choices
        choices = []
        for f in expanded_factors:
            if isinstance(f, Add):
                choices.append(f.args)
            else:
                choices.append((f,))

        terms = []
        for combination in product(*choices):
            inner = AlgebraicPureTensor(*combination)
            if scalar is S.One:
                terms.append(inner)
            else:
                terms.append(ScalarMul(scalar, inner))

        if len(terms) == 1:
            return terms[0]
        return AlgebraicTensor(*terms)

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

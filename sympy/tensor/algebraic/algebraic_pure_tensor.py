from __future__ import annotations

from sympy.core.basic import Basic
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor


"""Pure tensor as a non-commutative tensor product of matrix factors.

This module defines :class:`AlgebraicPureTensor`, the building block for
algebraic tensor expressions.  A pure tensor represents a single term in
a tensor expression -- the non-commutative tensor product of matrix-like
factors, optionally scaled by a commutative coefficient.

Each factor must carry ``.shape`` (e.g. any
:class:`~sympy.matrices.expressions.MatrixSymbol`).  The tensor shape is
the full sequence of per-factor shapes.
"""


class ShapeMismatchError(TypeError):
    """Raised when tensor operations encounter incompatible shapes.

    Used for both addition (shape mismatch between terms) and composition
    (inner dimension mismatch between factor pairs).
    """
    pass


def _factor_shapes(factors):
    """Return the full tensor shape as a tuple of per-factor shapes.

    E.g. for factors [A(3,4), C(4,5)] returns ``((3, 4), (4, 5))``.
    """
    return tuple(f.shape for f in factors)


def _factor_has_noncommutative(factor):
    """Return True if *factor* contains any noncommutative symbol.

    Checks ``factor.free_symbols`` for any symbol with
    ``is_commutative is False`` (e.g., MatrixSymbol or Symbol
    with commutative=False).  Concrete numeric
    matrices have empty free_symbols and are considered commutative.
    """
    for sym in getattr(factor, 'free_symbols', set()):
        if getattr(sym, 'is_commutative', True) is False:
            return True
    return False


def _is_zero_like(expr):
    """Return True if *expr* is effectively zero.

    First checks for ``S.Zero`` directly and the cached ``is_zero``
    assumption (avoiding the descriptor to prevent triggering SymPy's
    assumption inference on complex expressions), then falls back to
    the test ``expr == 0 * expr`` which works for ``ZeroMatrix``,
    zero ``Matrix`` instances, and other zero-like objects.
    """
    if expr is S.Zero:
        return True
    # Check cached assumption directly to avoid triggering _ask() on
    # complex Mul/Add expressions which can cause shape errors.
    if getattr(expr, '_assumptions', None) and expr._assumptions.get('zero') is True:
        return True
    return expr == S.Zero * expr


class AlgebraicPureTensor(Basic):
    """Pure tensor as an unevaluated (non-commutative)
    tensor product of factors.

    Extends ``Basic`` to provide a clean tensor-specific implementation
    without inheriting Mul's flatten() transformations.  The tensor-product
    operator is non-commutative: order of factors is always preserved.

    Each factor must carry ``.shape`` (e.g. any ``MatrixExpr`` or a
    1x1 Matrix wrapper around a non-commutative Symbol).

    The tensor shape is the full sequence of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for ``AlgebraicPureTensor(A_3x4, C_4x5)``.  No
    contraction is performed.

    Coefficients are stored internally as the first argument when present.

    Examples
    ========

    Create a pure tensor from two matrix symbols:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> print(T)
    A ⊗ B
    >>> T.shape
    ((3, 4), (4, 5))

    Numeric and symbolic coefficients are stored internally:

    >>> print(AlgebraicPureTensor(2, A, B))
    2*A ⊗ B
    >>> from sympy.abc import x
    >>> print(AlgebraicPureTensor(x, A, B))
    x*A ⊗ B

    Zero coefficient produces a zero tensor:

    >>> print(AlgebraicPureTensor(0, A, B))
    0_{(3x4), (4x5)}

    Single factor with coefficient 1 unwraps to the bare factor:

    >>> print(AlgebraicPureTensor(A))
    A
    """

    __slots__ = ()

    is_AlgebraicPureTensor = True

    _op_priority = 11  # Higher than Symbol/Expr so x * pt delegates to pt.__rmul__(x)

    _eval_is_commutative = lambda self: False

    @property
    def coeff(self):
        """The commutative coefficient of this tensor.

        Returns S.One if there is no explicit coefficient.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> AlgebraicPureTensor(A, B).coeff
        1
        >>> AlgebraicPureTensor(2, A, B).coeff
        2
        >>> from sympy.abc import x
        >>> AlgebraicPureTensor(x, A, B).coeff
        x
        """
        if len(self.args) > 0:
            first = self.args[0]
            if isinstance(first, Number) or (
                hasattr(first, 'is_commutative') and first.is_commutative
            ):
                return first
        return S.One

    @property
    def factors(self):
        """Individual tensor-product factors in left-to-right order.

        The coefficient (if present) is not included in the factors.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> AlgebraicPureTensor(A, B).factors
        (A, B)
        >>> AlgebraicPureTensor(2, A, B).factors
        (A, B)
        """
        if len(self.args) > 0:
            first = self.args[0]
            if isinstance(first, Number) or (
                hasattr(first, 'is_commutative') and first.is_commutative
            ):
                return self.args[1:]
        return self.args

    @property
    def num_factors(self):
        return len(self.factors)

    @property
    def shape(self):
        """Full tensor shape as a tuple of per-factor (rows, cols) pairs.

        E.g. ``AlgebraicPureTensor(A_3x4, C_4x5).shape == ((3, 4), (4, 5))``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> AlgebraicPureTensor(A, B).shape
        ((3, 4), (4, 5))
        >>> AlgebraicPureTensor(2, A, B).shape
        ((3, 4), (4, 5))
        """
        return _factor_shapes(self.factors)

    @property
    def commutativity_pattern(self):
        """Tuple of binary entries indicating per-factor commutativity.

        Entry i is 1 if the i-th tensor factor contains no
        noncommutative symbols, 0 otherwise.  Same length as ``shape``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> AlgebraicPureTensor(A, B).commutativity_pattern
        (0, 0)

        A numeric matrix factor is commutative in its slot:

        >>> from sympy.matrices import ImmutableDenseMatrix
        >>> M = ImmutableDenseMatrix([[1, 2, 3, 4],
        ...                           [5, 6, 7, 8],
        ...                           [9, 10, 11, 12]])
        >>> AlgebraicPureTensor(M, B).commutativity_pattern
        (1, 0)
        """
        return tuple(
            0 if _factor_has_noncommutative(f) else 1
            for f in self.factors
        )

    def __new__(cls, *args, evaluate=False):
        """Construct an AlgebraicPureTensor from factors.

        The first argument may be a commutative coefficient
        (Number or symbolic).

        Parameters
        ----------
        *args : commutative coefficient (optional), then matrix-like factors
            Optional leading commutative coefficient followed by one or more
            matrix-like factors, each carrying a ``.shape`` attribute.
        evaluate : bool, default False
            Not used; retained for compatibility.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicZeroTensor, or bare factor
            The result depends on the arguments.  A single factor with
            coefficient 1 unwraps to the bare factor.  Zero coefficient
            produces an AlgebraicZeroTensor.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> print(AlgebraicPureTensor(A, B))
        A ⊗ B
        >>> print(AlgebraicPureTensor(2, A, B))
        2*A ⊗ B
        >>> print(AlgebraicPureTensor(0, A, B))
        0_{(3x4), (4x5)}
        >>> print(AlgebraicPureTensor(A))
        A
        """
        if not args:
            raise ValueError("AlgebraicPureTensor requires at least one factor")

        if len(args) == 1:
            return args[0]

        # Separate leading commutative coefficient from tensor factors
        args_list = list(args)
        coeff = S.One

        if args_list:
            first = args_list[0]
            if isinstance(first, Number):
                coeff = first
                args_list = args_list[1:]
            else:
                first_s = sympify(first)
                if isinstance(first_s, Number):
                    coeff = first_s
                    args_list = args_list[1:]
                elif (hasattr(first_s, 'is_commutative') and
                        first_s.is_commutative and
                        not isinstance(first_s, AlgebraicPureTensor)):
                    # Commutative coefficient: extract it
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
            raise ValueError(
                "AlgebraicPureTensor requires at least one tensor factor "
                "(a zero coefficient alone cannot determine a shape)"
            )

        if _is_zero_like(coeff) or any(_is_zero_like(f) for f in processed):
            return AlgebraicZeroTensor(_factor_shapes(processed))

        # Two args with commutative first arg: return coeff * factor (plain Mul/MatMul)
        if len(args) == 2 and len(processed) == 1 and coeff is not S.One:
            return coeff * processed[0]

        if len(processed) == 1 and coeff is S.One:
            return processed[0]

        # Build the AlgebraicPureTensor with coefficient as first arg
        if coeff is not S.One:
            obj = Basic.__new__(cls, coeff, *processed)
        else:
            obj = Basic.__new__(cls, *processed)
        return obj

    def _get_coeff(self):
        """Extract the coefficient: same as coeff property.

        Provided for compatibility with code that uses _get_coeff().
        """
        return self.coeff

    def __neg__(self):
        """Return a negated tensor.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> print(-AlgebraicPureTensor(A, B))
        -1*A ⊗ B
        """
        if self.coeff is S.One:
            return AlgebraicPureTensor(S.NegativeOne, *self.factors)
        return AlgebraicPureTensor(-self.coeff, *self.factors)

    def __mul__(self, other):
        """Compose or scale this AlgebraicPureTensor.

        For commutative scalars/symbols the scalar is absorbed into the
        coefficient.  For AlgebraicPureTensor, AlgebraicTensor, or bare
        matrices the result is the tensor composition (factor-wise matrix
        multiplication) or matrix multiplication.

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The right operand.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor
            The scaled/composed result.

        Examples
        ========

        Scalar multiplication:

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> print(T * 3)
        3*A ⊗ B
        >>> print(T * 0)
        0_{(3x4), (4x5)}

        Composition with another tensor:

        >>> C = MatrixSymbol("C", 4, 2)
        >>> D = MatrixSymbol("D", 5, 3)
        >>> T2 = AlgebraicPureTensor(C, D)
        >>> print(T * T2)  # doctest: +NORMALIZE_WHITESPACE
        A*C ⊗ B*D
        """
        other = sympify(other)
        if isinstance(other, AlgebraicZeroTensor):
            return other
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, AlgebraicPureTensor)):
            if other is S.One:
                return self
            if other is S.Zero:
                return AlgebraicZeroTensor(self.shape)
            # Multiply into the coefficient
            return AlgebraicPureTensor(self.coeff * other, *self.factors)
        # Non-commutative operand: use tensor composition.
        from sympy.tensor.algebraic.algebraic_tensor import (
            compose_algebraic_tensors,
        )
        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        """Compose or scale this AlgebraicPureTensor from the left.

        For commutative scalars/symbols the scalar is absorbed into the
        coefficient.  For AlgebraicPureTensor, AlgebraicTensor, or bare
        matrices the result is the tensor composition (factor-wise matrix
        multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The left operand.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor
            The scaled/composed result.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> print(3 * T)
        3*A ⊗ B
        >>> print(0 * T)
        0_{(3x4), (4x5)}
        """
        if other == 0:
            return AlgebraicZeroTensor(self.shape)
        if other == 1:
            return self
        other = sympify(other)
        if isinstance(other, AlgebraicZeroTensor):
            return other
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, AlgebraicPureTensor)):
            if other is S.One:
                return self
            if other is S.Zero:
                return AlgebraicZeroTensor(self.shape)
            # Multiply into the coefficient
            return AlgebraicPureTensor(other * self.coeff, *self.factors)
        # Non-commutative operand: use tensor composition.
        from sympy.tensor.algebraic.algebraic_tensor import (
            compose_algebraic_tensors,
        )
        return compose_algebraic_tensors(other, self)

    def __add__(self, other):
        """Add another tensor or scalar to this pure tensor.

        Returns an :class:`~sympy.tensor.algebraic.AlgebraicTensor`
        containing both terms.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(A, B)
        >>> T2 = AlgebraicPureTensor(C, D)
        >>> print(T1 + T2)
        C ⊗ D + A ⊗ B
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        if isinstance(other, AlgebraicZeroTensor):
            if other.shape == self.shape:
                return self
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        """Right-add: ``other + self``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(A, B)
        >>> T2 = AlgebraicPureTensor(C, D)
        >>> print(T2 + T1)
        A ⊗ B + C ⊗ D
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        if isinstance(other, AlgebraicZeroTensor):
            if other.shape == self.shape:
                return self
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        """Subtract another tensor or scalar from this pure tensor.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(A, B)
        >>> T2 = AlgebraicPureTensor(C, D)
        >>> print(T1 - T2)
        -1*C ⊗ D + A ⊗ B
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        """Right-subtract: ``other - self``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(A, B)
        >>> T2 = AlgebraicPureTensor(C, D)
        >>> print(T2 - T1)
        -1*A ⊗ B + C ⊗ D
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, -self)

    def has_zero_term(self):
        return False

    @property
    def T(self):
        """Transpose of this pure tensor.

        Applies ``.T`` to every tensor factor, keeping the coefficient
        unchanged.  The tensor shape is updated so each ``(rows, cols)``
        becomes ``(cols, rows)``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> T.shape
        ((3, 4), (4, 5))
        >>> TT = T.T
        >>> TT.shape
        ((4, 3), (5, 4))
        >>> TT.coeff
        1

        Coefficient is preserved through transpose:

        >>> T2 = AlgebraicPureTensor(2, A, B)
        >>> T2T = T2.T
        >>> T2T.coeff
        2
        >>> T2T.shape
        ((4, 3), (5, 4))
        """
        transposed_factors = tuple(getattr(f, 'T', f) for f in self.factors)
        if self.coeff is S.One:
            return AlgebraicPureTensor(*transposed_factors)
        return AlgebraicPureTensor(self.coeff, *transposed_factors)

    def conjugate(self):
        """Return the complex conjugate of this pure tensor.

        Applies ``.conjugate()`` to the commutative coefficient and to
        every tensor-product factor, then reassembles the result.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
            The conjugated tensor.

        Examples
        ========

        Conjugate a tensor with a complex coefficient:

        >>> from sympy import I
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(1 + I, A, B)
        >>> TC = T.conjugate()
        >>> TC.coeff
        1 - I

        MatrixSymbol factors are conjugated via their .conjugate() method
        (which returns Adjoint(Transpose(self)) for symbolic matrices):

        >>> TC.factors[0] == A.conjugate()
        True
        >>> TC.factors[1] == B.conjugate()
        True

        Conjugate a tensor with a numeric matrix factor:

        >>> from sympy.matrices import ImmutableDenseMatrix
        >>> M = ImmutableDenseMatrix([[1 + I, 2], [3, 4 - I]])
        >>> N = MatrixSymbol("N", 2, 3)
        >>> T2 = AlgebraicPureTensor(M, N)
        >>> T2C = T2.conjugate()
        >>> T2C.factors[0][0, 0]
        1 - I
        """
        result = self._eval_conjugate()
        if result is not None:
            return result
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self)

    def _eval_conjugate(self):
        """Return the complex conjugate of this pure tensor.

        Applies ``.conjugate()`` to the commutative coefficient and to
        every tensor-product factor, then reassembles the result.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
            The conjugated tensor.

        Examples
        ========

        Conjugate a tensor with a complex coefficient:

        >>> from sympy import I
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(1 + I, A, B)
        >>> TC = T.conjugate()
        >>> TC.coeff
        1 - I

        MatrixSymbol factors are conjugated via their .conjugate() method
        (which returns Adjoint(Transpose(self)) for symbolic matrices):

        >>> TC.factors[0] == A.conjugate()
        True
        >>> TC.factors[1] == B.conjugate()
        True

        Conjugate a tensor with a numeric matrix factor:

        >>> from sympy.matrices import ImmutableDenseMatrix
        >>> M = ImmutableDenseMatrix([[1 + I, 2], [3, 4 - I]])
        >>> N = MatrixSymbol("N", 2, 3)
        >>> T2 = AlgebraicPureTensor(M, N)
        >>> T2C = T2.conjugate()
        >>> T2C.factors[0][0, 0]
        1 - I
        """
        coeff = self.coeff
        if hasattr(coeff, 'conjugate'):
            coeff = coeff.conjugate()

        conjugated_factors = []
        for f in self.factors:
            if hasattr(f, 'conjugate'):
                conjugated_factors.append(f.conjugate())
            else:
                conjugated_factors.append(f)

        if coeff is S.One:
            return AlgebraicPureTensor(*conjugated_factors)
        return AlgebraicPureTensor(coeff, *conjugated_factors)

    def simplify(self):
        from sympy.tensor.algebraic.simplify import _simplify_algebraic_pure_tensor
        return _simplify_algebraic_pure_tensor(self)

    def doit(self, **hints):
        """Evaluate the coefficient and each tensor factor.

        Applies ``.doit()`` to the commutative coefficient and to every
        tensor-product factor, then reassembles the result.  If any
        factor evaluates to zero, the result is an ``AlgebraicZeroTensor``.

        Parameters
        ----------
        **hints : dict
            Passed through to the ``doit()`` calls on factors.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicZeroTensor, or bare factor
            The reassembled tensor with evaluated components.

        Examples
        ========

        A tensor with no unevaluated sub-expressions returns an equal result:

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> T.doit() == T
        True

        With a symbolic coefficient:

        >>> from sympy.abc import x, y
        >>> T2 = AlgebraicPureTensor(x + y, A, B)
        >>> T2.doit().coeff
        x + y

        With deep=False, sub-expressions are not evaluated:

        >>> from sympy.matrices.expressions import MatAdd
        >>> C = MatrixSymbol("C", 4, 5)
        >>> T3 = AlgebraicPureTensor(A, MatAdd(B, C))
        >>> T3.doit(deep=False) is T3
        True
        """
        deep = hints.get('deep', True)
        
        # Check if any arg needs doit evaluation
        changed = False
        new_args = []
        for arg in self.args:
            if hasattr(arg, 'doit') and deep:
                da = arg.doit(**hints)
                new_args.append(da)
                if da is not arg:
                    changed = True
            else:
                new_args.append(arg)

        if not changed:
            return self

        # Rebuild from new args
        if len(new_args) == 0:
            return self
        if len(new_args) == 1:
            return new_args[0]
        # Check if first arg is a coefficient
        first = new_args[0]
        if isinstance(first, Number) or (hasattr(first, 'is_commutative') and first.is_commutative):
            if first is S.One:
                return AlgebraicPureTensor(*new_args[1:])
            return AlgebraicPureTensor(first, *new_args[1:])
        return AlgebraicPureTensor(*new_args)

    def diff(self, *symbols, **assumptions):
        """Differentiate this pure tensor with respect to *symbols*.

        Applies the Leibniz rule: the derivative acts on the coefficient
        and each tensor factor individually, summing the results.  For
        ``T = c * f1 ⊗ f2 ⊗ … ⊗ fn``:

            dT/dx = (dc/dx)*f1⊗…⊗fn
                   + c*(df1/dx)⊗f2⊗…⊗fn
                   + …
                   + c*f1⊗…⊗(dfn/dx)

        Parameters
        ----------
        *symbols : Symbol or str
            Symbol(s) to differentiate with respect to.
        **assumptions : dict
            Passed through to the underlying ``.diff()`` calls.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor, or matrix
            The differentiated tensor.  If only one nonzero term survives
            the Leibniz expansion, a single term is returned.

        Examples
        ========

        Differentiate a tensor with symbolic coefficient:

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> from sympy.abc import x
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(x**2, A, B)
        >>> print(T.diff(x))
        2*x*A ⊗ B

        Differentiate a tensor factor that depends on the symbol:

        >>> from sympy.matrices import ImmutableDenseMatrix
        >>> M = ImmutableDenseMatrix([[x, x**2], [1, x]])
        >>> N = MatrixSymbol("N", 2, 3)
        >>> T2 = AlgebraicPureTensor(M, N)
        >>> T2.diff(x)  # doctest: +SKIP

        Leibniz rule -- symbol appears in both coefficient and a factor:

        >>> M2 = ImmutableDenseMatrix([[x, 1], [0, x]])
        >>> T3 = AlgebraicPureTensor(x, M2, N)
        >>> T3_diff = T3.diff(x)
        >>> len(T3_diff.args)
        2
        """
        from sympy.core.singleton import S
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor

        coeff = self.coeff
        factors = self.factors

        terms = []

        # Differentiate the coefficient
        dc = coeff.diff(*symbols, **assumptions)
        if dc is not S.Zero:
            if dc is S.One:
                terms.append(AlgebraicPureTensor(*factors))
            else:
                terms.append(AlgebraicPureTensor(dc, *factors))

        # Differentiate each factor
        for i, f in enumerate(factors):
            df = f.diff(*symbols, **assumptions)
            if df is not S.Zero and not _is_zero_like(df):
                new_factors = list(factors)
                new_factors[i] = df
                if coeff is S.One:
                    terms.append(AlgebraicPureTensor(*new_factors))
                else:
                    terms.append(AlgebraicPureTensor(coeff, *new_factors))

        if not terms:
            return AlgebraicZeroTensor(self.shape)

        if len(terms) == 1:
            return terms[0]

        return AlgebraicTensor(*terms)

    def display(self, mode="latex"):
        """Display this tensor using IPython display or fallback to print.

        Parameters
        ----------
        mode : str, default 'latex'
            'latex' for LaTeX rendering, 'text' for plain text.
        """
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

    def expand(self, deep=True, **hints):
        """Expand this pure tensor by distributing over ``Add``
        in factors and coefficient.

        Parameters
        ----------
        deep : bool, default True
            Whether to expand subexpressions.
        **hints : dict
            Additional expansion hints.

        Returns
        -------
        AlgebraicPureTensor, AlgebraicTensor, or bare factor
            The expanded expression.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol, MatAdd
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 4, 5)
        >>> T = AlgebraicPureTensor(A, MatAdd(B, C))
        >>> print(T.expand())
        A ⊗ C + A ⊗ B
        """
        return self._eval_expand_mul(deep=deep, **hints)

    def _eval_expand_mul(self, **hints):
        """Expand this pure tensor by distributing over ``Add``
        in factors and coefficient.

        Expands the coefficient and each factor individually, then distributes
        the tensor product linearly across all addend combinations when either
        the coefficient or any factor is a SymPy ``Add``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol, MatAdd
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 4, 5)
        >>> T = AlgebraicPureTensor(A, MatAdd(B, C))
        >>> print(T.expand())
        A ⊗ C + A ⊗ B
        """
        from sympy.core.add import Add
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        from itertools import product

        deep = hints.pop('deep', True)
        factors = self.factors
        coeff = self.coeff

        expanded_coeff = coeff
        if hasattr(coeff, 'expand'):
            expanded_coeff = coeff.expand(deep=deep, **hints)

        expanded_factors = []
        for f in factors:
            if hasattr(f, 'expand'):
                expanded_factors.append(f.expand(deep=deep, **hints))
            else:
                expanded_factors.append(f)

        add_slots = [i for i, f in enumerate(expanded_factors)
                     if isinstance(f, Add)]

        has_add_coeff = isinstance(expanded_coeff, Add)

        if not add_slots and not has_add_coeff:
            if expanded_factors == list(factors) and expanded_coeff is coeff:
                return self
            if len(expanded_factors) == 1:
                if expanded_coeff is S.One:
                    return expanded_factors[0]
                return AlgebraicPureTensor(expanded_coeff, expanded_factors[0])
            if expanded_coeff is S.One:
                return AlgebraicPureTensor(*expanded_factors)
            return AlgebraicPureTensor(expanded_coeff, *expanded_factors)

        choices = []
        for f in expanded_factors:
            if isinstance(f, Add):
                choices.append(f.args)
            else:
                choices.append((f,))

        if has_add_coeff:
            coeff_choices = expanded_coeff.args
        else:
            coeff_choices = (expanded_coeff,)

        terms = []
        for coeff_combination in coeff_choices:
            for combination in product(*choices):
                if len(combination) == 1:
                    if coeff_combination is S.One:
                        terms.append(combination[0])
                    else:
                        terms.append(AlgebraicPureTensor(coeff_combination, combination[0]))
                else:
                    if coeff_combination is S.One:
                        terms.append(AlgebraicPureTensor(*combination))
                    else:
                        terms.append(AlgebraicPureTensor(coeff_combination, *combination))

        if len(terms) == 1:
            return terms[0]
        return AlgebraicTensor(*terms)


def algebraic_tensor_product(*args):
    """Convenience constructor for AlgebraicPureTensor.
    
    Examples
    ========

    >>> from sympy.tensor.algebraic import algebraic_tensor_product
    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> A = MatrixSymbol("A", 2, 3)
    >>> v = MatrixSymbol("v", 3, 1)
    >>> print(algebraic_tensor_product(A, v))
    A ⊗ v
    >>> print(algebraic_tensor_product(2, A, v))
    2*A ⊗ v
    """
    from sympy.tensor.algebraic.algebraic_tensor import (
        algebraic_tensor_product as _atp,
    )
    return _atp(*args)


def compose_algebraic_pure_tensors(left, right):
    """Compose two AlgebraicPureTensors by matrix-multiplying
    corresponding factors.

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

    Examples
    ========

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import (
    ...     AlgebraicPureTensor, compose_algebraic_pure_tensors)
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> C = MatrixSymbol("C", 4, 2)
    >>> D = MatrixSymbol("D", 5, 3)
    >>> T1 = AlgebraicPureTensor(A, B)
    >>> T2 = AlgebraicPureTensor(C, D)
    >>> print(compose_algebraic_pure_tensors(T1, T2))  # doctest: +NORMALIZE_WHITESPACE
    A*C ⊗ B*D
    """
    # --- AlgebraicZeroTensor shortcuts ---
    if isinstance(left, AlgebraicZeroTensor) and isinstance(right, AlgebraicZeroTensor):
        # Both zero: compute composed shape from both
        if len(left.shape) != len(right.shape):
            raise ShapeMismatchError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left.shape)} vs {len(right.shape)}"
            )
        for i, (l, r) in enumerate(zip(left.shape, right.shape)):
            if l[1] != r[0]:
                raise ShapeMismatchError(
                    f"Cannot compose factor {i}: left shape {l} "
                    f"but right shape {r}; inner dimensions "
                    f"{l[1]} and {r[0]} do not match"
                )
        composed_shape = tuple(
            (l[0], r[1]) for l, r in zip(left.shape, right.shape)
        )
        return AlgebraicZeroTensor(composed_shape)

    if isinstance(left, AlgebraicZeroTensor):
        # Extract factors from right to compute composed shape
        if isinstance(right, AlgebraicPureTensor):
            right_factors = right.factors
        elif hasattr(right, "shape"):
            right_factors = (right,)
        else:
            raise TypeError(
                f"Expected AlgebraicPureTensor or a matrix-like object on the right, "
                f"got {type(right).__name__}"
            )
        if len(left.shape) != len(right_factors):
            raise ShapeMismatchError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left.shape)} vs {len(right_factors)}"
            )
        for i, (l, rf) in enumerate(zip(left.shape, right_factors)):
            if l[1] != rf.shape[0]:
                raise ShapeMismatchError(
                    f"Cannot compose factor {i}: left shape {l} "
                    f"but right factor shape {rf.shape}; inner dimensions "
                    f"{l[1]} and {rf.shape[0]} do not match"
                )
        composed_shape = tuple(
            (l[0], rf.shape[1]) for l, rf in zip(left.shape, right_factors)
        )
        return AlgebraicZeroTensor(composed_shape)

    if isinstance(right, AlgebraicZeroTensor):
        # Extract factors from left to compute composed shape
        if isinstance(left, AlgebraicPureTensor):
            left_factors = left.factors
        elif hasattr(left, "shape"):
            left_factors = (left,)
        else:
            raise TypeError(
                f"Expected AlgebraicPureTensor or a matrix-like object on the left, "
                f"got {type(left).__name__}"
            )
        if len(left_factors) != len(right.shape):
            raise ShapeMismatchError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left_factors)} vs {len(right.shape)}"
            )
        for i, (lf, r) in enumerate(zip(left_factors, right.shape)):
            if lf.shape[1] != r[0]:
                raise ShapeMismatchError(
                    f"Cannot compose factor {i}: left factor shape {lf.shape} "
                    f"but right shape {r}; inner dimensions "
                    f"{lf.shape[1]} and {r[0]} do not match"
                )
        composed_shape = tuple(
            (lf.shape[0], r[1]) for lf, r in zip(left_factors, right.shape)
        )
        return AlgebraicZeroTensor(composed_shape)

    # Handle unwrapped single-factor tensors
    if isinstance(left, AlgebraicPureTensor):
        left_factors = left.factors
        left_coeff = left.coeff
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
        right_coeff = right.coeff
    elif hasattr(right, "shape"):
        right_factors = (right,)
        right_coeff = S.One
    else:
        raise TypeError(
            f"Expected AlgebraicPureTensor or a matrix-like object on the right, "
            f"got {type(right).__name__}"
        )

    if len(left_factors) != len(right_factors):
        raise ShapeMismatchError(
            f"Cannot compose tensors with different numbers of "
            f"factors: {len(left_factors)} vs {len(right_factors)}"
        )

    # Combine coefficients from both sides
    combined_coeff = left_coeff * right_coeff

    composed_factors = []
    for j, (lf, rf) in enumerate(zip(left_factors, right_factors)):
        lshape = lf.shape
        rshape = rf.shape
        if lshape[1] != rshape[0]:
            raise ShapeMismatchError(
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

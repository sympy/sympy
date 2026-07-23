from __future__ import annotations

from sympy.core.basic import Basic
from sympy.core.numbers import Number
from sympy.core.sympify import sympify


class AlgebraicZeroTensor(Basic):
    r"""Zero tensor carrying a specific tensor shape.

    An AlgebraicZeroTensor of shape ``((m0, n0), (m1, n1), ...)``
    acts as the additive identity for all AlgebraicTensors and
    AlgebraicPureTensors with the same sequence of factor shapes.
    AlgebraicZeroTensors with different shape sequences belong to
    different tensor spaces and are not summable.

    The shape is a tuple of ``(rows, cols)`` pairs, one per tensor-product
    factor. For a single-matrix tensor of shape ``(m, n)``, the constructor
    delegates to :class:`~sympy.matrices.expressions.special.ZeroMatrix`.

    Examples
    ========

    >>> from sympy.tensor.algebraic import AlgebraicZeroTensor
    >>> Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    >>> Z.shape
    ((3, 4), (4, 5))
    >>> print(-Z)
    0_{(3x4), (4x5)}

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> print(T + Z)
    A ⊗ B
    """

    __slots__ = ('_shape',)

    is_AlgebraicZeroTensor = True
    is_zero = True
    is_commutative = True

    _op_priority = 11

    def __new__(cls, shape):
        shape = tuple(shape)
        if len(shape) == 2 and not isinstance(shape[0], (tuple, list)):
            shape = (shape,)
        shape = tuple(tuple(s) for s in shape)

        if len(shape) == 1:
            rows, cols = shape[0]
            from sympy.matrices.expressions import ZeroMatrix
            return ZeroMatrix(rows, cols)

        obj = Basic.__new__(cls)
        obj._shape = shape
        return obj

    @property
    def args(self):
        return ()

    def _hashable_content(self):
        return (self._shape,)

    def __getnewargs__(self):
        return (self._shape,)

    def conjugate(self):
        return self

    @property
    def free_symbols(self):
        return set()

    @property
    def shape(self):
        """The tensor shape stored in _shape."""
        return self._shape

    @property
    def commutativity_pattern(self):
        """All-1s tuple -- a zero tensor is commutative in every slot."""
        return tuple(1 for _ in self._shape)

    def __neg__(self):
        return self

    def __add__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import (
            _validate_addition_shape,
        )
        _validate_addition_shape(self._shape, other)
        return other

    def __radd__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import (
            _validate_addition_shape,
        )
        _validate_addition_shape(self._shape, other)
        return other

    def __sub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import (
            _validate_addition_shape,
        )
        _validate_addition_shape(self._shape, other)
        return -other

    def __rsub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import (
            _validate_addition_shape,
        )
        _validate_addition_shape(self._shape, other)
        return other

    def __mul__(self, other):
        other = sympify(other)
        from sympy.matrices import Matrix, ImmutableDenseMatrix
        from sympy.matrices.expressions import MatrixSymbol
        from sympy.matrices.expressions.special import ZeroMatrix
        if isinstance(other, (Matrix, ImmutableDenseMatrix, MatrixSymbol, ZeroMatrix)):
            raise TypeError(
                f"Cannot multiply AlgebraicZeroTensor with {type(other).__name__}"
            )
        if isinstance(other, Number) or (
            hasattr(other, 'is_commutative') and not (hasattr(other, 'is_AlgebraicZeroTensor') or \
            hasattr(other, 'is_AlgebraicPureTensor') or hasattr(other, 'is_AlgebraicTensor'))
        ):
            return self
        from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors
        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        other = sympify(other)
        from sympy.matrices import Matrix, ImmutableDenseMatrix
        from sympy.matrices.expressions import MatrixSymbol
        from sympy.matrices.expressions.special import ZeroMatrix
        if isinstance(other, (Matrix, ImmutableDenseMatrix, MatrixSymbol, ZeroMatrix)):
            raise TypeError(
                f"Cannot multiply {type(other).__name__} with AlgebraicZeroTensor"
            )
        if isinstance(other, Number) or (
            hasattr(other, 'is_commutative') and not (hasattr(other, 'is_AlgebraicZeroTensor') or \
            hasattr(other, 'is_AlgebraicPureTensor') or hasattr(other, 'is_AlgebraicTensor'))
        ):
            return self
        from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors
        return compose_algebraic_tensors(other, self)

    def __bool__(self):
        return False

    @property
    def T(self):
        """Transpose of this zero tensor.

        Returns a new AlgebraicZeroTensor with every factor shape
        ``(rows, cols)`` reversed to ``(cols, rows)``.

        Examples
        ========

        >>> from sympy.tensor.algebraic import AlgebraicZeroTensor
        >>> Z = AlgebraicZeroTensor(((1, 2), (3, 4)))
        >>> Z.T.shape
        ((2, 1), (4, 3))
        """
        transposed_shape = tuple((c, r) for r, c in self._shape)
        return AlgebraicZeroTensor(transposed_shape)

    def _eval_conjugate(self):
        return self

    def expand(self, **kwargs):
        return self

    def doit(self, **hints):
        return self

    def diff(self, *symbols, **assumptions):
        return self

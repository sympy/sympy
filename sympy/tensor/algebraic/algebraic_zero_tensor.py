from __future__ import annotations

from sympy.core.expr import AtomicExpr
from sympy.core.numbers import Number
from sympy.core.sympify import sympify


"""Zero tensor carrying a specific tensor shape.

This module defines :class:`AlgebraicZeroTensor`, the additive identity for
tensors of a given shape.  Unlike SymPy's ``S.Zero``, an
:class:`AlgebraicZeroTensor` carries shape information, so zero tensors of
different shapes are distinct objects.

Examples
========

Create a zero tensor for a single-factor shape:

>>> from sympy.tensor.algebraic import AlgebraicZeroTensor
>>> Z = AlgebraicZeroTensor((3, 4))
>>> print(Z)
0_((3, 4),)
>>> Z.shape
((3, 4),)

Create a zero tensor for a multi-factor shape:

>>> Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
>>> print(Z2)
0_((3, 4), (4, 5))

Zero tensor acts as additive identity:

>>> from sympy.matrices.expressions import MatrixSymbol
>>> from sympy.tensor.algebraic import AlgebraicPureTensor
>>> A = MatrixSymbol("A", 3, 4)
>>> B = MatrixSymbol("B", 4, 5)
>>> T = AlgebraicPureTensor(A, B)
>>> Z3 = AlgebraicZeroTensor(((3, 4), (4, 5)))
>>> print(T + Z3)
A ⊗ B
"""


class AlgebraicZeroTensor(AtomicExpr):
    """Zero tensor carrying a specific tensor shape.

    An AlgebraicZeroTensor of shape ``((m0, n0), (m1, n1), ...)`` acts as the additive
    identity for all AlgebraicTensors and AlgebraicPureTensors with the same sequence
    of factor shapes.  AlgebraicZeroTensors with different shape sequences belong to
    different tensor spaces and are not summable.

    The shape is a tuple of (rows, cols) pairs, one per tensor-product factor.
    For a single-matrix tensor of shape ``(m, n)`` the canonical form is
    ``((m, n),)`` (a one-element tuple).

    Extends ``AtomicExpr`` (``Atom`` + ``Expr``) to integrate with SymPy's
    expression system: ``sympify()``, tree traversal (``atoms()``, ``has()``,
    ``replace()``), the assumptions system (``is_commutative``),
    generic operations (``subs()``, ``xreplace()``, ``doit()``), and the
    ``Expr`` machinery (``as_base_exp()``, ``as_coeff_Mul()``, ``as_coeff_Add()``).

    Examples
    ========

    Create a zero tensor for a single-factor shape:

    >>> Z = AlgebraicZeroTensor((3, 4))
    >>> print(Z)
    0_((3, 4),)
    >>> Z.shape
    ((3, 4),)

    Create a zero tensor for a multi-factor shape:

    >>> Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    >>> print(Z2)
    0_((3, 4), (4, 5))

    Zero tensor is its own negation:

    >>> print(-Z)
    0_((3, 4),)

    Zero tensor is the additive identity:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> Z3 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    >>> print(T + Z3)
    A ⊗ B
    """

    __slots__ = ('_shape',)

    is_AlgebraicZeroTensor = True
    is_zero = None
    is_commutative = True

    _op_priority = 11  # Higher than Expr (10) so x * zt dispatches to zt.__rmul__(x)

    def __new__(cls, shape):
        # Normalise to a tuple of (rows, cols) tuples.
        # Accept: ((3,4), (4,5)), [(3,4)], plain (3,4), [3,4], etc.
        shape = tuple(shape)
        if len(shape) == 2 and not isinstance(shape[0], tuple):
            # Bare (m, n) or [m, n] -> wrap as ((m, n),)
            shape = (shape,)
        shape = tuple(tuple(s) for s in shape)
        obj = AtomicExpr.__new__(cls)
        obj._shape = shape
        return obj

    def _hashable_content(self):
        """Return content for hashing and equality comparison."""
        return (self._shape,)

    def __getnewargs__(self):
        """Return args for pickle reconstruction."""
        return (self._shape,)

    @property
    def free_symbols(self):
        """A zero tensor has no free symbols."""
        return set()

    @property
    def shape(self):
        """The tensor shape stored in _shape.

        Examples
        ========

        >>> Z = AlgebraicZeroTensor((3, 4))
        >>> Z.shape
        ((3, 4),)
        >>> AlgebraicZeroTensor(((3, 4), (4, 5))).shape
        ((3, 4), (4, 5))
        """
        return self._shape

    @property
    def commutativity_pattern(self):
        """All-1s tuple -- a zero tensor is commutative in every slot.

        Examples
        ========

        >>> AlgebraicZeroTensor((3, 4)).commutativity_pattern
        (1,)
        >>> AlgebraicZeroTensor(((3, 4), (4, 5))).commutativity_pattern
        (1, 1)
        """
        return tuple(1 for _ in self._shape)

    def __neg__(self):
        """Negation of zero is zero.

        Examples
        ========

        >>> Z = AlgebraicZeroTensor((3, 4))
        >>> print(-Z)
        0_((3, 4),)
        """
        return self

    def __add__(self, other):
        """Additive identity -- return the other operand unchanged.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
        >>> print(Z + T)
        A ⊗ B
        """
        return other

    def __radd__(self, other):
        """Additive identity -- return the other operand unchanged.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
        >>> print(T + Z)
        A ⊗ B
        """
        return other

    def __sub__(self, other):
        """Subtract *other* from the zero tensor.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
        >>> print(Z - T)
        -1*A ⊗ B + 0_((3, 4), (4, 5))
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        """Right-subtract: ``other - self`` (returns *other* unchanged).

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> T = AlgebraicPureTensor(A, B)
        >>> Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
        >>> print(T - Z)
        A ⊗ B + 0_((3, 4), (4, 5))
        """
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, -self)

    def __mul__(self, other):
        """Scaling by a commutative scalar returns self.

        For non-commutative operands, delegates to composition.

        Examples
        ========

        >>> Z = AlgebraicZeroTensor((3, 4))
        >>> print(5 * Z)
        0_((3, 4),)
        >>> print(Z * 5)
        0_((3, 4),)
        """
        other = sympify(other)
        if isinstance(other, Number) or (
            hasattr(other, 'is_commutative') and other.is_commutative
        ):
            return self
        from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors
        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        """Scaling by a commutative scalar returns self.

        For non-commutative operands, delegates to composition.

        Examples
        ========

        >>> Z = AlgebraicZeroTensor((3, 4))
        >>> print(5 * Z)
        0_((3, 4),)
        """
        other = sympify(other)
        if isinstance(other, Number) or (
            hasattr(other, 'is_commutative') and other.is_commutative
        ):
            return self
        from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors
        return compose_algebraic_tensors(other, self)

    def __bool__(self):
        return False

    def copy(self):
        return self

    def expand(self, **kwargs):
        """Return self unchanged. A zero tensor is already in expanded form."""
        return self

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


def algebraic_zero_tensor(shape):
    """Convenience constructor for AlgebraicZeroTensor.

    Parameters
    ----------
    shape : tuple of tuples, or a single (rows, cols) pair
        The full tensor shape, e.g. ``((3, 4), (4, 5))`` for a product
        of a 3x4 matrix and a 4x5 matrix, or ``((3, 4),)`` for a single
        3x4 factor.  A bare pair ``(3, 4)`` is accepted and wrapped as
        ``((3, 4),)``.

    Returns
    -------
    AlgebraicZeroTensor

    Examples
    ========

    >>> from sympy.tensor.algebraic import algebraic_zero_tensor
    >>> print(algebraic_zero_tensor((3, 4)))
    0_((3, 4),)
    >>> print(algebraic_zero_tensor(((3, 4), (4, 5))))
    0_((3, 4), (4, 5))
    """
    return AlgebraicZeroTensor(shape)

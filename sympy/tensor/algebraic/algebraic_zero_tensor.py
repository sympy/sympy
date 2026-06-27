from __future__ import annotations


class AlgebraicZeroTensor:
    """Zero tensor carrying a specific tensor shape.

    An AlgebraicZeroTensor of shape ``((m0, n0), (m1, n1), ...)`` acts as the additive
    identity for all AlgebraicTensors and AlgebraicPureTensors with the same sequence
    of factor shapes.  AlgebraicZeroTensors with different shape sequences belong to
    different tensor spaces and are not summable.

    The shape is a tuple of (rows, cols) pairs, one per tensor-product factor.
    For a single-matrix tensor of shape ``(m, n)`` the canonical form is
    ``((m, n),)`` (a one-element tuple).
    """

    __slots__ = ("_shape",)
    is_AlgebraicZeroTensor = True

    def __init__(self, shape):
        # Normalise to a tuple of (rows, cols) tuples.
        # Accept: ((3,4), (4,5)), [(3,4)], plain (3,4), [3,4], etc.
        shape = tuple(shape)
        if len(shape) == 2 and not isinstance(shape[0], tuple):
            # Bare (m, n) or [m, n] → wrap as ((m, n),)
            shape = (shape,)
        self._shape = tuple(tuple(s) for s in shape)

    @property
    def shape(self):
        return self._shape

    @property
    def tensor_shape(self):
        return self._shape

    def __neg__(self):
        return self

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

    def __eq__(self, other):
        if isinstance(other, AlgebraicZeroTensor):
            return self._shape == other._shape
        return NotImplemented

    def __hash__(self):
        return hash(("AlgebraicZeroTensor", self._shape))

    def __repr__(self):
        return f"AlgebraicZeroTensor{self._shape}"

    def __str__(self):
        return f"0_{self._shape}"

    def __bool__(self):
        return False

    def simplify(self, **kwargs):
        """An AlgebraicZeroTensor is already in simplest form; return self."""
        return self


def algebraic_zero_tensor(shape):
    """Convenience constructor for AlgebraicZeroTensor.

    Parameters
    ----------
    shape : tuple of tuples, or a single (rows, cols) pair
        The full tensor shape, e.g. ``((3, 4), (4, 5))`` for a product
        of a 3×4 matrix and a 4×5 matrix, or ``((3, 4),)`` for a single
        3×4 factor.  A bare pair ``(3, 4)`` is accepted and wrapped as
        ``((3, 4),)``.

    Returns
    -------
    AlgebraicZeroTensor
    """
    return AlgebraicZeroTensor(shape)

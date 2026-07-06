from __future__ import annotations

from sympy.core.basic import Atom


class AlgebraicZeroTensor(Atom):
    """Zero tensor carrying a specific tensor shape.

    An AlgebraicZeroTensor of shape ``((m0, n0), (m1, n1), ...)`` acts as the additive
    identity for all AlgebraicTensors and AlgebraicPureTensors with the same sequence
    of factor shapes.  AlgebraicZeroTensors with different shape sequences belong to
    different tensor spaces and are not summable.

    The shape is a tuple of (rows, cols) pairs, one per tensor-product factor.
    For a single-matrix tensor of shape ``(m, n)`` the canonical form is
    ``((m, n),)`` (a one-element tuple).

    Extends ``Atom`` (a leaf ``Basic`` subclass) to integrate with SymPy's
    expression system: ``sympify()``, tree traversal (``atoms()``, ``has()``,
    ``replace()``), the assumptions system (``is_zero``, ``is_commutative``),
    and generic operations (``subs()``, ``xreplace()``, ``doit()``).
    """

    __slots__ = ('_shape',)

    is_AlgebraicZeroTensor = True
    is_zero = True
    is_commutative = True

    def __new__(cls, shape):
        # Normalise to a tuple of (rows, cols) tuples.
        # Accept: ((3,4), (4,5)), [(3,4)], plain (3,4), [3,4], etc.
        shape = tuple(shape)
        if len(shape) == 2 and not isinstance(shape[0], tuple):
            # Bare (m, n) or [m, n] -> wrap as ((m, n),)
            shape = (shape,)
        shape = tuple(tuple(s) for s in shape)
        obj = Atom.__new__(cls)
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
        """The tensor shape stored in _shape."""
        return self._shape

    @property
    def tensor_shape(self):
        """Alias for ``shape`` -- the full tensor shape tuple."""
        return self._shape

    @property
    def commutativity_shape(self):
        """All-1s tuple -- a zero tensor is commutative in every slot."""
        return tuple(1 for _ in self._shape)

    def __neg__(self):
        """Negation of zero is zero."""
        return self

    def __add__(self, other):
        """Additive identity -- return the other operand unchanged."""
        return other

    def __radd__(self, other):
        """Additive identity -- return the other operand unchanged."""
        return other

    def __sub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, -self)

    def __mul__(self, other):
        """Composition of a zero tensor returns a zero tensor of the same shape."""
        return self

    def __rmul__(self, other):
        """Composition of a zero tensor from the left returns the zero tensor."""
        return self

    def __repr__(self):
        return f"AlgebraicZeroTensor{self._shape}"

    def __str__(self):
        return f"0_{self._shape}"

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
    """
    return AlgebraicZeroTensor(shape)

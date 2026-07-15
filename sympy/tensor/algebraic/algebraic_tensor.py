from __future__ import annotations

from sympy.core.add import add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor, _factor_has_noncommutative
)
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor


"""Linear combinations of same-shape algebraic pure tensors.

This module defines :class:`AlgebraicTensor`, a sum of
:class:`~sympy.tensor.algebraic.algebraic_pure_tensor.AlgebraicPureTensor`
terms that all share the same tensor shape.  It also provides
:func:`compose_algebraic_tensors`, the linear composition operator that
extends factor-wise matrix multiplication to sums of tensors.

Examples
========

Create a sum of two pure tensors with the same shape:

>>> from sympy.matrices.expressions import MatrixSymbol
>>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor, AlgebraicTensor
>>> A = MatrixSymbol("A", 3, 4)
>>> B = MatrixSymbol("B", 4, 5)
>>> C = MatrixSymbol("C", 3, 4)
>>> D = MatrixSymbol("D", 4, 5)
>>> T1 = AlgebraicPureTensor(A, B)
>>> T2 = AlgebraicPureTensor(C, D)
>>> S = AlgebraicTensor(T1, T2)
>>> print(S)
C ⊗ D + A ⊗ B
>>> S.shape
((3, 4), (4, 5))

Addition of pure tensors is routed to AlgebraicTensor:

>>> print(T1 + T2)
C ⊗ D + A ⊗ B
"""


# ---------------------------------------------------------------------------
# Shape helpers
# ---------------------------------------------------------------------------

def _shape_of(expr):
    """Return the full tensor shape of *expr* as a tuple of factor shapes.

    Handles AlgebraicPureTensor, AlgebraicZeroTensor,
    Mul(coeff, AlgebraicPureTensor), and bare matrix-like objects
    (whose single-factor shape is wrapped).
    """
    if isinstance(expr, AlgebraicPureTensor):
        return expr.shape
    if isinstance(expr, AlgebraicZeroTensor):
        return expr.shape
    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                return f.shape
            if isinstance(f, AlgebraicZeroTensor):
                return f.shape
    if hasattr(expr, "shape"):
        return (expr.shape,)
    return None


def _commutativity_pattern_of(expr):
    """Return the commutativity_pattern of *expr*, or None if not determinable.

    Handles AlgebraicPureTensor, AlgebraicZeroTensor,
    AlgebraicTensor, Mul(coeff, AlgebraicPureTensor),
    and bare matrix-like objects.
    """
    if isinstance(expr, AlgebraicPureTensor):
        return expr.commutativity_pattern
    if isinstance(expr, AlgebraicZeroTensor):
        return tuple(1 for _ in expr.shape)
    if isinstance(expr, AlgebraicTensor):
        return expr.commutativity_pattern
    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                return f.commutativity_pattern
            if isinstance(f, AlgebraicZeroTensor):
                return tuple(1 for _ in f.shape)
    if hasattr(expr, "shape"):
        if _factor_has_noncommutative(expr):
            return (0,)
        return (1,)
    return None


# ---------------------------------------------------------------------------

class ShapeMismatchError(TypeError):
    """Raised when attempting to add tensors of incompatible shapes.

    Examples
    ========

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
    >>> from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> C = MatrixSymbol("C", 2, 3)
    >>> D = MatrixSymbol("D", 3, 4)
    >>> try:
    ...     AlgebraicTensor(AlgebraicPureTensor(A, B),
    ...                     AlgebraicPureTensor(C, D))
    ... except ShapeMismatchError:
    ...     print("shape mismatch")
    shape mismatch
    """
    pass


def _compose_reassemble(results, shape):
    """Reassemble a list of composition results into a single expression.

    Strips AlgebraicZeroTensor anchors, handles cancellation, and wraps
    back into an AlgebraicTensor if needed.
    """
    from sympy.core.singleton import S
    from sympy.tensor.algebraic.algebraic_zero_tensor import (
        AlgebraicZeroTensor as _ZT,
    )

    real = []
    zero_shape = None
    for r in results:
        if isinstance(r, _ZT):
            zero_shape = r.shape
        elif r is not S.Zero:
            real.append(r)

    if not real:
        if zero_shape:
            return _ZT(zero_shape)
        return _ZT(shape)

    if len(real) == 1:
        return real[0]

    return AlgebraicTensor(*real, _sympify=False)


def compose_algebraic_tensors(left, right):
    """Compose two algebraic-tensor expressions by linearity.

    Composes *left* and *right* by factor-wise matrix multiplication,
    extending ``compose_algebraic_pure_tensors`` by linearity
    to sums of tensors.

    Parameters
    ----------
    left : AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
        The left operand.
    right : AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
        The right operand.

    Returns
    -------
    AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
        The composition result.

    Raises
    ------
    TypeError
        If either argument is not a recognised tensor or matrix type.
    ValueError
        If the factor structures are incompatible for composition.

    Examples
    ========

    Compose two pure tensors:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import (
    ...     AlgebraicPureTensor, compose_algebraic_tensors)
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> C = MatrixSymbol("C", 4, 2)
    >>> D = MatrixSymbol("D", 5, 3)
    >>> T1 = AlgebraicPureTensor(A, B)
    >>> T2 = AlgebraicPureTensor(C, D)
    >>> print(compose_algebraic_tensors(T1, T2))  # doctest: +NORMALIZE_WHITESPACE
    A*C ⊗ B*D

    Compose two sums:

    >>> G = MatrixSymbol("G", 3, 4)
    >>> H = MatrixSymbol("H", 4, 5)
    >>> E = MatrixSymbol("E", 4, 2)
    >>> F = MatrixSymbol("F", 5, 3)
    >>> S1 = T1 + AlgebraicPureTensor(G, H)
    >>> S2 = T2 + AlgebraicPureTensor(E, F)
    >>> print(compose_algebraic_tensors(S1, S2))  # doctest: +NORMALIZE_WHITESPACE
    A*C ⊗ B*D + A*E ⊗ B*F + G*C ⊗ H*D + G*E ⊗ H*F
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import (
        AlgebraicPureTensor as _PT,
        compose_algebraic_pure_tensors,
    )
    from sympy.tensor.algebraic.algebraic_zero_tensor import (
        AlgebraicZeroTensor as _ZT,
    )

    # --- AlgebraicZeroTensor shortcuts ---
    if isinstance(left, _ZT):
        return _ZT(left.shape)
    if isinstance(right, _ZT):
        return _ZT(right.shape)

    # --- AlgebraicTensor × AlgebraicTensor (check first before single-side) ---
    if isinstance(left, AlgebraicTensor) and isinstance(right, AlgebraicTensor):
        results = []
        for la in left.args:
            if isinstance(la, _ZT):
                results.append(la)
                continue
            for ra in right.args:
                if isinstance(ra, _ZT):
                    results.append(ra)
                else:
                    comp = compose_algebraic_tensors(la, ra)
                    results.append(comp)
        return _compose_reassemble(results, left.shape)

    # --- PureTensor / bare matrix × AlgebraicTensor ---
    if isinstance(right, AlgebraicTensor):
        if isinstance(left, _PT) or hasattr(left, "shape"):
            results = []
            for a in right.args:
                if isinstance(a, _ZT):
                    results.append(a)
                else:
                    results.append(
                        compose_algebraic_tensors(left, a)
                    )
            return _compose_reassemble(results, right.shape)
        raise TypeError(
            f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
            f"the left, got {type(left).__name__}"
        )

    # --- AlgebraicTensor × PureTensor / bare matrix ---
    if isinstance(left, AlgebraicTensor):
        if isinstance(right, _PT) or hasattr(right, "shape"):
            return left._compose_with_term(right)
        raise TypeError(
            f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
            f"the right, got {type(right).__name__}"
        )

    # --- PureTensor / bare matrix × PureTensor / bare matrix ---
    if isinstance(left, _PT) and isinstance(right, _PT):
        return compose_algebraic_pure_tensors(left, right)
    if isinstance(left, _PT) and hasattr(right, "shape"):
        return compose_algebraic_pure_tensors(left, right)
    if hasattr(left, "shape") and isinstance(right, _PT):
        return compose_algebraic_pure_tensors(left, right)
    if hasattr(left, "shape") and hasattr(right, "shape"):
        return compose_algebraic_pure_tensors(left, right)

    raise TypeError(
        f"Cannot compose {type(left).__name__} with {type(right).__name__}"
    )


class AlgebraicTensor(Basic):
    """Sum of AlgebraicPureTensors (and/or AlgebraicZeroTensors)
    sharing the same tensor shape.

    Built on top of ``Basic`` (not ``Add``) to avoid the is_commutative
    descriptor conflict in ``AssocOp._from_args``, while still exposing
    ``is_Add = True`` so that the wider SymPyecosystem recognises this
    as an additive expression.

    Shape enforcement and AlgebraicZeroTensor handling are done in ``__new__``.

    Tensor shapes are full sequences of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for the product ``A_3x4 ⊗ C_4x5``.

    Examples
    ========

    Create a sum of two pure tensors with the same shape:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> C = MatrixSymbol("C", 3, 4)
    >>> D = MatrixSymbol("D", 4, 5)
    >>> T1 = AlgebraicPureTensor(A, B)
    >>> T2 = AlgebraicPureTensor(C, D)
    >>> print(AlgebraicTensor(T1, T2))
    C ⊗ D + A ⊗ B

    Addition is routed to AlgebraicTensor:

    >>> print(T1 + T2)
    C ⊗ D + A ⊗ B

    Single term unwraps to the bare term:

    >>> print(AlgebraicTensor(T1))
    A ⊗ B

    Cancelling terms produce a zero tensor:

    >>> print(AlgebraicTensor(T1, -T1))
    0_{(3x4), (4x5)}
    """

    __slots__ = ()

    is_AlgebraicTensor = True
    is_Add = True

    identity = None  # no single identity; use AlgebraicZeroTensor(shape) instead

    def __new__(cls, *args, evaluate=True, _sympify=True):
        if not args:
            raise ValueError("AlgebraicTensor requires at least one argument")

        if _sympify:
            args = tuple(
                a if isinstance(a, AlgebraicZeroTensor) else sympify(a)
                for a in args
            )

        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (AlgebraicPureTensor, AlgebraicZeroTensor)):
                return arg
            if isinstance(arg, AlgebraicTensor):
                return arg
            if isinstance(arg, (list, tuple)) and len(arg) == 1:
                return arg[0]

        # Flatten nested AlgebraicTensors and collect all leaf terms
        flat, shape, zero_term, comm_cs = cls._flatten_args(args)

        # Remember whether zero_term was user-provided (from _flatten_args)
        # vs. created later by coefficient cancellation.
        zero_term_was_user_provided = zero_term is not None

        # Collect coefficients of PureTensor terms with identical factors
        flat, zero_term = cls._collect_coefficients(flat, shape, zero_term)

        if not flat:
            if zero_term is not None:
                return zero_term
            raise ValueError("AlgebraicTensor resulted in zero terms")

        # Remove identity Numbers (S.Zero) from flat list
        flat = [a for a in flat if a is not S.Zero]

        if not flat:
            if zero_term is not None:
                return zero_term
            raise ValueError("AlgebraicTensor resulted in zero terms")

        # If an AlgebraicZeroTensor was explicitly provided by the user,
        # always keep it as an anchor.  Only unwrap to a single term when
        # there is exactly one non-zero term AND no AlgebraicZeroTensor was
        # explicitly given.
        if len(flat) == 1 and not zero_term_was_user_provided and not isinstance(flat[0], AlgebraicZeroTensor):
            return flat[0]

        # Append user-provided AlgebraicZeroTensor (for shape anchoring).
        # Cancellation-created zero_term is NOT appended when there are
        # surviving non-zero terms -- it was only needed for the all-cancelled
        # case handled above.
        if zero_term_was_user_provided:
            flat.append(zero_term)

        obj = Basic.__new__(cls, *flat)
        return obj

    @classmethod
    def _collect_coefficients(cls, terms, shape, zero_term):
        """Combine coefficients of AlgebraicPureTensor terms
        with identical factors.

        Uses a dictionary keyed by the factor-only AlgebraicPureTensor to
        accumulate coefficients in O(N) time.  AlgebraicZeroTensor entries
        that somehow ended up in the terms list are also removed.

        Parameters
        ----------
        terms : list
            Flattened term list from ``_flatten_args``.
        shape : tuple | None
            The common tensor shape.
        zero_term : AlgebraicZeroTensor | None
            Existing zero-tensor anchor (if any).

        Returns
        -------
        (new_terms : list, new_zero_term : AlgebraicZeroTensor | None)
        """
        result = []
        coeff_map = {}

        # Separate AlgebraicPureTensor terms from others
        for t in terms:
            if isinstance(t, AlgebraicZeroTensor):
                if zero_term is None:
                    zero_term = t
                continue

            if not isinstance(t, AlgebraicPureTensor):
                result.append(t)
                continue

            # Build key_pure_tensor from factors (without coefficient)
            factors = t.factors
            if len(factors) == 1:
                key_pure_tensor = factors[0]
            else:
                key_pure_tensor = AlgebraicPureTensor(*factors)

            c = t.coeff

            if key_pure_tensor in coeff_map:
                coeff_map[key_pure_tensor] = coeff_map[key_pure_tensor] + c
            else:
                coeff_map[key_pure_tensor] = c

        # Rebuild terms from the dictionary
        for key_pure_tensor, combined_coeff in coeff_map.items():
            if combined_coeff == 0:
                if zero_term is None and shape is not None:
                    zero_term = AlgebraicZeroTensor(shape)
            elif combined_coeff is S.One:
                result.append(key_pure_tensor)
            else:
                result.append(AlgebraicPureTensor(combined_coeff, *key_pure_tensor.factors))

        return result, zero_term

    @classmethod
    def _flatten_args(cls, args):
        """Flatten nested AlgebraicTensors, validate shapes,
        collect AlgebraicZeroTensors.

        Returns
        -------
        (terms : list, shape : tuple | None,
         zero_term : AlgebraicZeroTensor | None,
         commutativity_pattern : tuple | None)

        *shape* is a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.
        *commutativity_pattern* is a tuple of binary entries,
        same length as shape, representing the component-wise AND
        of all term commutativity_patterns.
        """
        shape = None
        zero_term = None
        comm_cs = None
        terms = []

        work = list(args)
        while work:
            o = work.pop()

            if isinstance(o, AlgebraicZeroTensor):
                candidate = o.shape
                if shape is None:
                    shape = candidate
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                zero_term = o
                continue

            if isinstance(o, AlgebraicTensor):
                if comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, o.commutativity_pattern))
                work.extend(o.args)
                continue

            if isinstance(o, AlgebraicPureTensor):
                candidate = o.shape
                if shape is None:
                    shape = candidate
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                if comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, o.commutativity_pattern))
                terms.append(o)
                continue

            if isinstance(o, Mul) and not isinstance(o, AlgebraicPureTensor):
                candidate = _shape_of(o)
                if candidate is not None:
                    if shape is None:
                        shape = candidate
                        comm_cs = tuple(1 for _ in shape)
                    elif candidate != shape:
                        raise ShapeMismatchError(
                            f"Cannot add tensors of different shapes: "
                            f"{shape} vs {candidate}"
                        )
                cs = _commutativity_pattern_of(o)
                if cs is not None and comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, cs))
                terms.append(o)
                continue

            if o.is_Number:
                terms.append(o)
                continue

            # Catch-all: bare matrix-like objects or anything else
            candidate = _shape_of(o)
            if candidate is not None:
                if shape is None:
                    shape = candidate
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
            cs = _commutativity_pattern_of(o)
            if cs is not None and comm_cs is not None:
                comm_cs = tuple(r & c for r, c in zip(comm_cs, cs))
            terms.append(o)

        return terms, shape, zero_term, comm_cs

    @property
    def shape(self):
        """Shape shared by every term in this sum.

        Returns a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> S.shape
        ((3, 4), (4, 5))
        """
        for arg in self.args:
            ts = _shape_of(arg)
            if ts is not None:
                return ts
        raise AttributeError("Cannot determine shape")

    @property
    def commutativity_pattern(self):
        """Component-wise AND of commutativity_pattern
        over all terms in this sum.

        For AlgebraicTensor operands, uses their stored commutativity_pattern
        rather than recomputing from individual PureTensor factors.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> S.commutativity_pattern
        (0, 0)
        """
        ts = self.shape
        if ts is None:
            raise AttributeError("Cannot determine commutativity_pattern without shape")
        result = tuple(1 for _ in ts)
        for arg in self.args:
            cs = _commutativity_pattern_of(arg)
            if cs is not None:
                result = tuple(r & c for r, c in zip(result, cs))
        return result

    @property
    def terms(self):
        """Non-zero, non-coefficient terms in this sum.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> len(S.terms)
        2
        """
        return tuple(
            a for a in self.args
            if not isinstance(a, Number) and not isinstance(a, AlgebraicZeroTensor)
        )

    def has_zero_term(self):
        """Return True if an AlgebraicZeroTensor anchors this sum."""
        return any(isinstance(a, AlgebraicZeroTensor) for a in self.args)

    @property
    def T(self):
        """Transpose of this algebraic tensor.

        Transposes every term in the sum.  For each ``AlgebraicPureTensor``
        term, ``.T`` is applied factor-wise.  ``AlgebraicZeroTensor``
        anchors are transposed by reversing each factor shape.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> S.shape
        ((3, 4), (4, 5))
        >>> ST = S.T
        >>> ST.shape
        ((4, 3), (5, 4))
        """
        transposed_args = []
        for a in self.args:
            if isinstance(a, AlgebraicZeroTensor):
                transposed_args.append(a.T)
            elif hasattr(a, 'T'):
                transposed_args.append(a.T)
            else:
                transposed_args.append(a)
        return AlgebraicTensor(*transposed_args)

    # ---- Arithmetic delegation ----

    def __neg__(self):
        """Negate each term in the sum.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(-S)
        -1*A ⊗ B - 1*C ⊗ D
        """
        return AlgebraicTensor(*(-a for a in self.args))

    def _compose_with_term(self, other):
        """Compose this AlgebraicTensor with a single term
        (AlgebraicPureTensor or bare matrix).

        Composes each term of self with *other* by factor-wise matrix
        multiplication.

        Parameters
        ----------
        other : AlgebraicPureTensor or bare matrix-like object
            The right operand.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, or AlgebraicZeroTensor
            The composition result.
        """
        from sympy.core.singleton import S
        from sympy.tensor.algebraic.algebraic_pure_tensor import (
            AlgebraicPureTensor as _PT,
            compose_algebraic_pure_tensors,
        )
        from sympy.tensor.algebraic.algebraic_zero_tensor import (
            AlgebraicZeroTensor as _ZT,
        )

        results = []
        for a in self.args:
            if isinstance(a, _ZT):
                results.append(a)
                continue

            if isinstance(a, _PT):
                comp = compose_algebraic_pure_tensors(a, other)
                results.append(comp)
            elif isinstance(a, Mul) and not isinstance(a, _PT):
                # Mul(coeff, AlgebraicPureTensor) — compose the PureTensor part
                for f in a.args:
                    if isinstance(f, _PT):
                        coeff = Mul(*[x for x in a.args if x is not f])
                        comp = compose_algebraic_pure_tensors(f, other)
                        # Re-wrap with coefficient
                        if coeff is S.One:
                            results.append(comp)
                        else:
                            results.append(AlgebraicPureTensor(coeff, *comp.factors if isinstance(comp, _PT) else (comp,)))
                        break
                else:
                    results.append(a)
            elif hasattr(a, "shape"):
                comp = compose_algebraic_pure_tensors(a, other)
                results.append(comp)
            else:
                results.append(a)

        return _compose_reassemble(results, self.shape)

    def __add__(self, other):
        """Add another tensor expression to this sum.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> E = MatrixSymbol("E", 3, 4)
        >>> F = MatrixSymbol("F", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(S + AlgebraicPureTensor(E, F))
        E ⊗ F + A ⊗ B + C ⊗ D
        """
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        """Right-add: ``other + self``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> E = MatrixSymbol("E", 3, 4)
        >>> F = MatrixSymbol("F", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(AlgebraicPureTensor(E, F) + S)
        A ⊗ B + C ⊗ D + E ⊗ F
        """
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        """Subtract another tensor expression from this sum.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(S - AlgebraicPureTensor(A, B))
        C ⊗ D
        """
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        """Right-subtract: ``other - self``.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(AlgebraicPureTensor(A, B) - S)
        -1*C ⊗ D
        """
        return AlgebraicTensor(other, -self)

    def __mul__(self, other):
        """Compose or scale this AlgebraicTensor.

        For commutative scalars/symbols the scalar is absorbed as a
        coefficient into each term.  For AlgebraicPureTensor,
        AlgebraicTensor, or bare matrices the result is the tensor
        composition (factor-wise matrix multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The right operand.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or Mul
            The scaled/composed result.

        Examples
        ========

        Scalar multiplication:

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(S * 2)
        2*A ⊗ B + 2*C ⊗ D
        """
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(a * other for a in self.args))
        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        """Compose or scale this AlgebraicTensor from the left.

        For commutative scalars/symbols the scalar is absorbed as a
        coefficient into each term.  For AlgebraicPureTensor,
        AlgebraicTensor, or bare matrices the result is the tensor
        composition (factor-wise matrix multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The left operand.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or Mul
            The scaled/composed result.

        Examples
        ========

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> S = AlgebraicTensor(AlgebraicPureTensor(A, B),
        ...                     AlgebraicPureTensor(C, D))
        >>> print(2 * S)
        2*A ⊗ B + 2*C ⊗ D
        """
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(other * a for a in self.args))
        return compose_algebraic_tensors(other, self)

    def simplify(self):
        """Simplify this AlgebraicTensor using the tensor
        simplification pipeline.

        Applies proportionality factoring and per-term simplification.
        """
        from sympy.tensor.algebraic.simplify import _simplify_algebraic_tensor
        return _simplify_algebraic_tensor(self)

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
        """Expand this AlgebraicTensor by expanding each term.

        Each term is expanded individually, and the results are reassembled
        into a single AlgebraicTensor (flattening any nested sums).

        Examples
        ========

        >>> from sympy.matrices.expressions import (
        ...     MatrixSymbol, MatAdd)
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 4, 5)
        >>> T = AlgebraicPureTensor(A, MatAdd(B, C))
        >>> S = AlgebraicTensor(T)
        >>> print(S.expand())
        A ⊗ C + A ⊗ B
        """
        from sympy.core.add import Add

        expanded_args = []
        for a in self.args:
            if isinstance(a, AlgebraicZeroTensor):
                expanded_args.append(a)
            elif hasattr(a, 'expand'):
                expanded_term = a.expand(deep=deep, **hints)
                if isinstance(expanded_term, Add):
                    expanded_args.extend(expanded_term.args)
                else:
                    expanded_args.append(expanded_term)
            else:
                expanded_args.append(a)

        if expanded_args == list(self.args):
            return self

        return AlgebraicTensor(*expanded_args)


# Register AlgebraicTensor with the SymPy add dispatcher so that
# expressions like (A + B) where A/B involve PureTensor or AlgebraicTensor
# are routed through AlgebraicTensor.__new__ instead of Add.__new__.
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)

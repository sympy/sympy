from __future__ import annotations

from sympy.core.add import add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor, _factor_has_noncommutative, _is_zero_like
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

    Handles AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor,
    Mul(coeff, AlgebraicPureTensor), and bare matrix-like objects
    (whose single-factor shape is wrapped).
    """
    if isinstance(expr, AlgebraicPureTensor):
        return expr.shape
    if isinstance(expr, AlgebraicZeroTensor):
        return expr.shape
    if isinstance(expr, AlgebraicTensor):
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
    if isinstance(left, _ZT) and isinstance(right, _ZT):
        if len(left.shape) != len(right.shape):
            raise ValueError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left.shape)} vs {len(right.shape)}"
            )
        composed_shape = tuple(
            (l[0], r[1]) for l, r in zip(left.shape, right.shape)
        )
        return _ZT(composed_shape)

    if isinstance(left, _ZT):
        # Extract shape from right operand
        if isinstance(right, AlgebraicTensor):
            right_shape = right.shape
        elif isinstance(right, _PT):
            right_shape = right.shape
        elif hasattr(right, "shape"):
            right_shape = (right.shape,)
        else:
            raise TypeError(
                f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
                f"the right, got {type(right).__name__}"
            )
        if len(left.shape) != len(right_shape):
            raise ValueError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left.shape)} vs {len(right_shape)}"
            )
        composed_shape = tuple(
            (l[0], r[1]) for l, r in zip(left.shape, right_shape)
        )
        return _ZT(composed_shape)

    if isinstance(right, _ZT):
        # Extract shape from left operand
        if isinstance(left, AlgebraicTensor):
            left_shape = left.shape
        elif isinstance(left, _PT):
            left_shape = left.shape
        elif hasattr(left, "shape"):
            left_shape = (left.shape,)
        else:
            raise TypeError(
                f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
                f"the left, got {type(left).__name__}"
            )
        if len(left_shape) != len(right.shape):
            raise ValueError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left_shape)} vs {len(right.shape)}"
            )
        composed_shape = tuple(
            (l[0], r[1]) for l, r in zip(left_shape, right.shape)
        )
        return _ZT(composed_shape)

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

    __slots__ = ('_coeff_map',)

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

 

        # Collect coefficients of PureTensor terms with identical factors
        flat, zero_term, coeff_map = cls._collect_coefficients(flat, shape, zero_term)

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

        # Unwrap single term (AlgebraicZeroTensor already handled above).
        if len(flat) == 1 and not isinstance(flat[0], AlgebraicZeroTensor):
            return flat[0]

        # Drop AlgebraicZeroTensor anchor when non-zero terms survive;
        # it is only returned when flat is empty (handled above).

        obj = Basic.__new__(cls, *flat)
        obj._coeff_map = coeff_map
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
        (new_terms : list, new_zero_term : AlgebraicZeroTensor | None,
         coeff_map : dict)
            The coeff_map keys are always AlgebraicPureTensor (unit coefficient)
            representing the factor structure, and values are the combined
            commutative coefficients.
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

            # Build key_pure_tensor from factors (without coefficient).
            # Always wrap in AlgebraicPureTensor for consistent keys.
            factors = t.factors
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

        return result, zero_term, coeff_map

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

    def __getstate__(self):
        """Exclude _coeff_map from pickled state (recomputable from args)."""
        state = super().__getstate__()
        if state is not None and '_coeff_map' in state:
            del state['_coeff_map']
        return state

    def _has_simple_terms(self):
        """Check if all args are AlgebraicPureTensor or AlgebraicZeroTensor."""
        return all(
            isinstance(a, (AlgebraicPureTensor, AlgebraicZeroTensor))
            for a in self.args
        )

    @classmethod
    def _combine_coeff_maps(cls, left, right, negate_right=False):
        """Combine coefficient maps of two AlgebraicTensors.

        Shared helper for _merge and _subtract. Merges the coefficient
        maps by iterating over right's args and adding (or subtracting,
        when negate_right is True) their coefficients into left's map.

        Both operands must have only AlgebraicPureTensor and AlgebraicZeroTensor
        terms (checked by _has_simple_terms before calling).

        Parameters
        ----------
        left : AlgebraicTensor
            The base tensor whose coeff_map is copied and accumulated into.
        right : AlgebraicTensor
            The tensor whose terms are merged into left's map.
        negate_right : bool, default False
            If True, negate each coefficient from right before merging
            (used for subtraction).

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
            The combined result, with the same unwrapping semantics as __new__.
        """
        combined = left._coeff_map.copy()
        zero_term = None
        zero_term_was_user_provided = left.has_zero_term()

        for a in right.args:
            if isinstance(a, AlgebraicTensor):
                nested = a._coeff_map
                for key, c in nested.items():
                    if negate_right:
                        c = -c
                    if key in combined:
                        combined[key] = combined[key] + c
                    else:
                        combined[key] = c
                if a.has_zero_term():
                    zero_term_was_user_provided = True
            elif isinstance(a, AlgebraicPureTensor):
                factors = a.factors
                key = AlgebraicPureTensor(*factors)
                c = a.coeff
                if negate_right:
                    c = -c
                if key in combined:
                    combined[key] = combined[key] + c
                else:
                    combined[key] = c
            elif isinstance(a, AlgebraicZeroTensor):
                zero_term = a
                zero_term_was_user_provided = True

        real = []
        for key, c in combined.items():
            if c == 0:
                if zero_term is None:
                    zero_term = AlgebraicZeroTensor(left.shape)
            elif c is S.One:
                real.append(key)
            else:
                real.append(AlgebraicPureTensor(c, *key.factors))

        if not real:
            if zero_term is not None:
                return zero_term
            return AlgebraicZeroTensor(left.shape)

        if len(real) == 1 and not zero_term_was_user_provided and not isinstance(real[0], AlgebraicZeroTensor):
            return real[0]

        if zero_term_was_user_provided:
            real.append(zero_term)

        obj = Basic.__new__(cls, *real)
        obj._coeff_map = combined
        return obj

    @classmethod
    def _merge(cls, left, right):
        """Merge two AlgebraicTensors by combining their coefficient maps.

        Both operands must have only AlgebraicPureTensor and AlgebraicZeroTensor
        terms (checked by _has_simple_terms before calling).
        """
        return cls._combine_coeff_maps(left, right, negate_right=False)

    @classmethod
    def _subtract(cls, left, right):
        """Subtract right from left by merging with negated coefficients.

        Both operands must have only AlgebraicPureTensor and AlgebraicZeroTensor
        terms (checked by _has_simple_terms before calling).
        """
        return cls._combine_coeff_maps(left, right, negate_right=True)

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

    def conjugate(self):
        """Return the complex conjugate of this algebraic tensor.

        Applies ``.conjugate()`` to each term in the sum by linearity.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, or AlgebraicZeroTensor
            The conjugated tensor.  If the result collapses to a single
            term, that term is returned directly.

        Examples
        ========

        Conjugate a sum of two pure tensors:

        >>> from sympy import I
        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(1 + I, A, B)
        >>> T2 = AlgebraicPureTensor(2 - I, C, D)
        >>> S = AlgebraicTensor(T1, T2)
        >>> SC = S.conjugate()
        >>> isinstance(SC, AlgebraicTensor)
        True

        Coefficients are conjugated in each term:

        >>> coeff_vals = set()
        >>> for arg in SC.args:
        ...     if isinstance(arg, AlgebraicPureTensor):
        ...         coeff_vals.add(arg.coeff)
        >>> (1 - I) in coeff_vals
        True
        >>> (2 + I) in coeff_vals
        True

        Double conjugate returns the original tensor:

        >>> S.conjugate().conjugate() == S
        True
        """
        result = self._eval_conjugate()
        if result is not None:
            return result
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self)

    def _eval_conjugate(self):
        """Apply conjugate to each term in the sum by linearity.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, or AlgebraicZeroTensor
            The conjugated tensor.
        """
        conjugated_args = []
        for a in self.args:
            if hasattr(a, 'conjugate'):
                conjugated_args.append(a.conjugate())
            else:
                conjugated_args.append(a)

        if conjugated_args == list(self.args):
            return self

        return AlgebraicTensor(*conjugated_args)

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
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._merge(self, other)
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
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._merge(other, self)
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
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._subtract(self, other)
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
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._subtract(other, self)
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

    def doit(self, **hints):
        """Evaluate each term in the sum by linearity.

        Applies ``.doit()`` to every term (``AlgebraicPureTensor`` or
        ``AlgebraicZeroTensor``) in the sum, then reassembles the
        results into a single ``AlgebraicTensor``.

        Parameters
        ----------
        **hints : dict
            Passed through to the ``doit()`` calls on individual terms.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, or AlgebraicZeroTensor
            The reassembled tensor with evaluated terms.  If the result
            collapses to a single term, that term is returned directly.

        Examples
        ========

        A sum with no unevaluated sub-expressions returns itself:

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(A, B)
        >>> T2 = AlgebraicPureTensor(C, D)
        >>> S = AlgebraicTensor(T1, T2)
        >>> S.doit() is S
        True

        With symbolic coefficients, doit evaluates the coefficients:

        >>> from sympy.abc import x, y
        >>> T3 = AlgebraicPureTensor(x, A, B)
        >>> T4 = AlgebraicPureTensor(y, A, B)
        >>> S2 = AlgebraicTensor(T3, T4)
        >>> S2.doit().coeff
        x + y

        With deep=False, sub-expressions are not evaluated:

        >>> S.doit(deep=False) is S
        True
        """
        deep = hints.get('deep', True)
        evaluated_args = []
        for a in self.args:
            if hasattr(a, 'doit') and deep:
                evaluated_args.append(a.doit(**hints))
            else:
                evaluated_args.append(a)

        if evaluated_args == list(self.args):
            return self

        return AlgebraicTensor(*evaluated_args)

    def diff(self, *symbols, **assumptions):
        """Differentiate this tensor sum with respect to *symbols*.

        Applies ``.diff()`` to each term in the sum by linearity and
        reassembles the results into a single ``AlgebraicTensor``.

        Parameters
        ----------
        *symbols : Symbol or str
            Symbol(s) to differentiate with respect to.
        **assumptions : dict
            Passed through to the underlying ``.diff()`` calls.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
            The differentiated tensor sum.  If the result collapses to
            a single term, that term is returned directly.

        Examples
        ========

        Differentiate a sum of two pure tensors:

        >>> from sympy.matrices.expressions import MatrixSymbol
        >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
        >>> from sympy.abc import x
        >>> A = MatrixSymbol("A", 3, 4)
        >>> B = MatrixSymbol("B", 4, 5)
        >>> C = MatrixSymbol("C", 3, 4)
        >>> D = MatrixSymbol("D", 4, 5)
        >>> T1 = AlgebraicPureTensor(x**2, A, B)
        >>> T2 = AlgebraicPureTensor(x, C, D)
        >>> S = AlgebraicTensor(T1, T2)
        >>> S_diff = S.diff(x)
        >>> len(S_diff.args)
        2

        Differentiating a sum with respect to a symbol not present:

        >>> from sympy.abc import y
        >>> S.diff(y)
        0_{(3x4), (4x5)}
        """
        from sympy.core.singleton import S

        differentiated = []
        for a in self.args:
            if hasattr(a, 'diff'):
                differentiated.append(a.diff(*symbols, **assumptions))
            else:
                differentiated.append(S.Zero)

        return AlgebraicTensor(*differentiated)

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


def algebraic_tensor_product(*args):
    """Tensor product of matrix factors, pure tensors, and tensor sums.

    Generalizes the tensor-product constructor to accept
    ``AlgebraicPureTensor`` and ``AlgebraicTensor`` arguments.

    For each argument the function extracts (coefficient, factors)::

        scalar / Number          -> (value, ())
        AlgebraicPureTensor      -> (coeff, factors)
        AlgebraicTensor          -> list of (coeff, factors) per term
        AlgebraicZeroTensor      -> zero result with combined shape
        bare matrix with .shape  -> (1, (matrix,))

    The result is the sum of all tensor-product combinations, with
    coefficients multiplied together.

    Parameters
    ----------
    *args
        Any combination of scalars, matrix-like objects,
        ``AlgebraicPureTensor``, ``AlgebraicTensor``, or
        ``AlgebraicZeroTensor``.

    Returns
    -------
    AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor, or matrix
        The tensor product.  A single-factor result with coefficient 1
        unwraps to the bare factor.  If any argument is a zero tensor,
        the result is a zero tensor of the combined shape.

    Raises
    ------
    TypeError
        If an argument is not a recognised type.

    Examples
    ========

    Basic tensor product of matrix factors:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import algebraic_tensor_product
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> print(algebraic_tensor_product(A, B))
    A ⊗ B

    With a scalar coefficient:

    >>> print(algebraic_tensor_product(2, A, B))
    2*A ⊗ B

    PureTensor argument -- coefficient is extracted and factors are
    flattened into the product:

    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> from sympy.abc import x, y
    >>> pt = AlgebraicPureTensor(y, B, A.T)
    >>> print(algebraic_tensor_product(x, A, pt))
    x*y*A ⊗ B ⊗ A.T

    AlgebraicTensor argument -- the product distributes over the sum:

    >>> from sympy.tensor.algebraic import AlgebraicTensor
    >>> C = MatrixSymbol("C", 3, 4)
    >>> D = MatrixSymbol("D", 4, 5)
    >>> at = AlgebraicTensor(AlgebraicPureTensor(A, B),
    ...                      AlgebraicPureTensor(C, D))
    >>> E = MatrixSymbol("E", 5, 3)
    >>> print(algebraic_tensor_product(at, E))
    A ⊗ B ⊗ E + C ⊗ D ⊗ E

    Zero tensor argument produces a zero tensor of the combined shape:

    >>> from sympy.tensor.algebraic import AlgebraicZeroTensor
    >>> zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    >>> result = algebraic_tensor_product(A, zt)
    >>> result
    0_{(3x4), (3x4), (4x5)}
    """
    from itertools import product as iterproduct

    # ---- helpers ----

    def _terms_of_arg(arg):
        """Return list of (coeff, factors_tuple) for *arg*.

        Returns [S.Zero] as a sentinel when the argument is an
        AlgebraicZeroTensor or an AlgebraicTensor with no real terms.
        """
        if isinstance(arg, AlgebraicPureTensor):
            return [(arg.coeff, arg.factors)]

        if isinstance(arg, AlgebraicTensor):
            terms = []
            for t in arg.args:
                if isinstance(t, AlgebraicZeroTensor):
                    continue
                if isinstance(t, AlgebraicPureTensor):
                    terms.append((t.coeff, t.factors))
                else:
                    terms.append((S.One, (t,)))
            if not terms:
                return [S.Zero]
            return terms

        if isinstance(arg, AlgebraicZeroTensor):
            return [S.Zero]

        if isinstance(arg, Number):
            return [(arg, ())]

        arg_s = sympify(arg)
        if isinstance(arg_s, Number):
            return [(arg_s, ())]

        # Commutative symbols and expressions are treated as coefficients
        if (hasattr(arg_s, 'is_commutative') and arg_s.is_commutative
                and not isinstance(arg_s, AlgebraicPureTensor)):
            return [(arg_s, ())]

        if hasattr(arg_s, "shape"):
            return [(S.One, (arg_s,))]

        raise TypeError(
            f"algebraic_tensor_product does not accept {type(arg).__name__}"
        )

    def _combined_shape(args):
        """Concatenate the shape tuples of all arguments."""
        parts = []
        for a in args:
            s = _shape_of(a)
            if s is None:
                continue
            parts.append(s)
        return tuple(__i for __p in parts for __i in __p)

    # ---- main logic ----

    if not args:
        raise ValueError(
            "algebraic_tensor_product requires at least one argument"
        )

    if len(args) == 1:
        return args[0]

    if len(args) == 2:
        first = args[0]
        if isinstance(first, Number):
            return first * sympify(args[1])
        first_s = sympify(first)
        if isinstance(first_s, Number) or (
            hasattr(first_s, 'is_commutative') and first_s.is_commutative
            and not isinstance(first_s, AlgebraicZeroTensor)
        ):
            return first_s * sympify(args[1])

    arg_lists = [_terms_of_arg(a) for a in args]

    # If any argument is zero-like (0 * a == a), the whole product is zero.
    for a in args:
        if _is_zero_like(a):
            return AlgebraicZeroTensor(_combined_shape(args))

    # Check for zero-tensor sentinel
    for al in arg_lists:
        if len(al) == 1 and al[0] is S.Zero:
            return AlgebraicZeroTensor(_combined_shape(args))

    combined_shape = _combined_shape(args)

    # Build all combinations
    result_terms = []
    for combination in iterproduct(*arg_lists):
        coeff = S.One
        all_factors = ()
        for c, fs in combination:
            coeff = coeff * c
            all_factors = all_factors + fs
        result_terms.append((coeff, all_factors))

    if not result_terms:
        return AlgebraicZeroTensor(combined_shape)

    # Build individual AlgebraicPureTensor terms
    pt_terms = []
    for coeff, factors in result_terms:
        if coeff is S.Zero:
            pt_terms.append(AlgebraicZeroTensor(combined_shape))
        elif len(factors) == 0:
            continue
        elif len(factors) == 1 and coeff is S.One:
            pt_terms.append(factors[0])
        else:
            pt_terms.append(
                AlgebraicPureTensor(coeff, *factors)
            )

    if not pt_terms:
        if combined_shape:
            return AlgebraicZeroTensor(combined_shape)
        raise ValueError(
            "algebraic_tensor_product requires at least one "
            "tensor factor (matrix-like object)"
        )

    if len(pt_terms) == 1:
        return pt_terms[0]

    return AlgebraicTensor(*pt_terms)


# Register AlgebraicTensor with the SymPy add dispatcher so that
# expressions like (A + B) where A/B involve PureTensor or AlgebraicTensor
# are routed through AlgebraicTensor.__new__ instead of Add.__new__.
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)

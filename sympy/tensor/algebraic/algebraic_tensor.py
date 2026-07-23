from __future__ import annotations

from sympy.core.add import add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor, ShapeMismatchError,
    _factor_has_noncommutative, _is_zero_like
)
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor


"""Linear combinations of same-shape algebraic pure tensors.

This module defines :class:`AlgebraicTensor`, a sum of
:class:`~sympy.tensor.algebraic.algebraic_pure_tensor.AlgebraicPureTensor`
terms that all share the same tensor shape.  It also provides
:func:`compose_algebraic_tensors`, the linear composition operator that
extends factor-wise matrix multiplication to sums of tensors.
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


def _validate_addition_shape(self_shape, other):
    """Validate that *other* can be added to a tensor of *self_shape*.

    Raises ``ShapeMismatchError`` if *other* lacks a ``.shape`` attribute,
    or if both have shapes that differ.
    """
    if not hasattr(other, 'shape'):
        raise ShapeMismatchError(
            "Addition of tensors and non-tensors is not defined"
        )
    if other.shape != self_shape:
        raise ShapeMismatchError(
            f"Cannot add tensors of different shapes: "
            f"{self_shape} vs {other.shape}"
        )


# ---------------------------------------------------------------------------

def _compose_reassemble(results, shape):
    """Reassemble composition results into a single expression."""
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
        return _ZT(composed_shape)

    if isinstance(left, _ZT):
        # Extract shape from right operand
        if isinstance(right, AlgebraicTensor):
            right_shape = right.shape
        elif isinstance(right, _PT):
            right_shape = right.shape
        else:
            raise TypeError(
                f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
                f"the right, got {type(right).__name__}"
            )
        if len(left.shape) != len(right_shape):
            raise ShapeMismatchError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left.shape)} vs {len(right_shape)}"
            )
        for i, (l, r) in enumerate(zip(left.shape, right_shape)):
            if l[1] != r[0]:
                raise ShapeMismatchError(
                    f"Cannot compose factor {i}: left shape {l} "
                    f"but right shape {r}; inner dimensions "
                    f"{l[1]} and {r[0]} do not match"
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
        else:
            raise TypeError(
                f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
                f"the left, got {type(left).__name__}"
            )
        if len(left_shape) != len(right.shape):
            raise ShapeMismatchError(
                f"Cannot compose tensors with different numbers of "
                f"factors: {len(left_shape)} vs {len(right.shape)}"
            )
        for i, (l, r) in enumerate(zip(left_shape, right.shape)):
            if l[1] != r[0]:
                raise ShapeMismatchError(
                    f"Cannot compose factor {i}: left shape {l} "
                    f"but right shape {r}; inner dimensions "
                    f"{l[1]} and {r[0]} do not match"
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

    # --- PureTensor × AlgebraicTensor ---
    if isinstance(right, AlgebraicTensor):
        if isinstance(left, _PT):
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

    # --- AlgebraicTensor × PureTensor ---
    if isinstance(left, AlgebraicTensor):
        if isinstance(right, _PT):
            return left._compose_with_term(right)
        raise TypeError(
            f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
            f"the right, got {type(right).__name__}"
        )

    # --- PureTensor × PureTensor ---
    if isinstance(left, _PT) and isinstance(right, _PT):
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
        flat, shape, zero_term, _ = cls._flatten_args(args)

        # Collect coefficients of PureTensor terms with identical factors
        flat, zero_term, coeff_map = cls._collect_coefficients(flat, shape, zero_term)

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
        over all terms in this sum."""
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
        """Non-zero, non-coefficient terms in this sum."""
        return tuple(
            a for a in self.args
            if not isinstance(a, Number) and not isinstance(a, AlgebraicZeroTensor)
        )

    def has_zero_term(self):
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

        Shared helper for _merge and _subtract. Both operands must
        have only AlgebraicPureTensor and AlgebraicZeroTensor terms.
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

        # Intentionally bypass AlgebraicTensor.__new__ (and its _sympify
        # parameter) because all terms are already properly constructed
        # AlgebraicPureTensor/AlgebraicZeroTensor objects. Sympification
        # would be redundant and could trigger unwanted re-normalization.
        obj = Basic.__new__(cls, *real)
        obj._coeff_map = combined
        return obj

    @classmethod
    def _merge(cls, left, right):
        return cls._combine_coeff_maps(left, right, negate_right=False)

    @classmethod
    def _subtract(cls, left, right):
        return cls._combine_coeff_maps(left, right, negate_right=True)

    @property
    def T(self):
        """Transpose of this algebraic tensor.

        Transposes every term in the sum.

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
        >>> S.T.shape
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
        return AlgebraicTensor(*(-a for a in self.args))

    def _compose_with_term(self, other):
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
                        elif isinstance(comp, _PT):
                            results.append(AlgebraicPureTensor(coeff, *comp.factors))
                        else:
                            # Bare matrix (single-factor unwrapped result)
                            results.append(AlgebraicPureTensor(coeff, comp))
                        break
                else:
                    results.append(a)
            else:
                results.append(a)

        return _compose_reassemble(results, self.shape)

    def __add__(self, other):
        _validate_addition_shape(self.shape, other)
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._merge(self, other)
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        _validate_addition_shape(self.shape, other)
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._merge(other, self)
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        _validate_addition_shape(self.shape, other)
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._subtract(self, other)
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        _validate_addition_shape(self.shape, other)
        if isinstance(other, AlgebraicTensor):
            if self._has_simple_terms() and other._has_simple_terms():
                return self._subtract(other, self)
        return AlgebraicTensor(other, -self)

    def __mul__(self, other):
        other = sympify(other)
        from sympy.matrices import Matrix, ImmutableDenseMatrix
        from sympy.matrices.expressions import MatrixSymbol
        from sympy.matrices.expressions.special import ZeroMatrix
        if isinstance(other, (Matrix, ImmutableDenseMatrix, MatrixSymbol, ZeroMatrix)):
            raise TypeError(
                f"Cannot multiply AlgebraicTensor with {type(other).__name__}"
            )
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor))):
            return AlgebraicTensor(*(a * other for a in self.args))
        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        other = sympify(other)
        from sympy.matrices import Matrix, ImmutableDenseMatrix
        from sympy.matrices.expressions import MatrixSymbol
        from sympy.matrices.expressions.special import ZeroMatrix
        if isinstance(other, (Matrix, ImmutableDenseMatrix, MatrixSymbol, ZeroMatrix)):
            raise TypeError(
                f"Cannot multiply {type(other).__name__} with AlgebraicTensor"
            )
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor))):
            return AlgebraicTensor(*(other * a for a in self.args))
        return compose_algebraic_tensors(other, self)

    def simplify(self):
        """Simplify this AlgebraicTensor.

        Applies proportionality factoring and per-term simplification.
        """
        from sympy.tensor.algebraic.simplify import _simplify_algebraic_tensor
        return _simplify_algebraic_tensor(self)

    def doit(self, **hints):
        """Evaluate each term in the sum by linearity.

        Applies ``.doit()`` to every term in the sum, then reassembles
        the results into a single ``AlgebraicTensor``.

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
    ``AlgebraicPureTensor`` and ``AlgebraicTensor`` arguments.  The result
    is the sum of all tensor-product combinations, with coefficients
    multiplied together.

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

    # If any argument is zero-like or is an AlgebraicZeroTensor (detected
    # by the [S.Zero] sentinel from _terms_of_arg), the whole product is zero.
    for a, al in zip(args, arg_lists):
        if _is_zero_like(a) or (len(al) == 1 and al[0] is S.Zero):
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

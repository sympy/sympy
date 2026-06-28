from __future__ import annotations

from sympy.core.add import Add, add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor, _factor_shapes
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor, algebraic_zero_tensor


# ---------------------------------------------------------------------------
# Shape helpers
# ---------------------------------------------------------------------------

def _normalize_shape(s):
    """Normalize any shape-like input to a tuple of (rows, cols) tuples.

    ``((3,4), (4,5))``  -> ``((3, 4), (4, 5))``  (unchanged)
    ``(3, 4)``          -> ``((3, 4),)``          (bare pair wrapped)
    ``[(3,4)]``         -> ``((3, 4),)``
    """
    if not s:
        return ()
    if len(s) == 2 and not isinstance(s[0], (tuple, list)):
        # Bare (m, n) -> ((m, n),)
        return (tuple(s),)
    return tuple(tuple(x) for x in s)


def _tensor_shape_of(expr):
    """Return the full tensor shape of *expr* as a tuple of factor shapes.

    Handles AlgebraicPureTensor, AlgebraicZeroTensor, Mul(coeff, AlgebraicPureTensor), and bare
    matrix-like objects (whose single-factor shape is wrapped).
    """
    if isinstance(expr, AlgebraicPureTensor):
        return expr.tensor_shape
    if isinstance(expr, AlgebraicZeroTensor):
        return expr.shape
    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                return f.tensor_shape
            if isinstance(f, AlgebraicZeroTensor):
                return f.shape
    if hasattr(expr, "shape"):
        return (expr.shape,)
    return None


# ---------------------------------------------------------------------------

class ShapeMismatchError(TypeError):
    """Raised when attempting to add tensors of incompatible shapes."""
    pass


def _compose_reassemble(results, tensor_shape):
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
        return _ZT(tensor_shape)

    if len(real) == 1:
        return real[0]

    return AlgebraicTensor(*real, _sympify=False)


def compose_algebraic_tensors(left, right):
    """Compose two algebraic-tensor expressions by linearity.

    Composes *left* and *right* by factor-wise matrix multiplication,
    extending ``compose_algebraic_pure_tensors`` by linearity to sums.

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
                    comp = compose_algebraic_pure_tensors(la, ra)
                    results.append(comp)
        return _compose_reassemble(results, left.tensor_shape)

    # --- PureTensor / bare matrix × AlgebraicTensor ---
    if isinstance(right, AlgebraicTensor):
        if isinstance(left, _PT) or hasattr(left, "shape"):
            results = []
            for a in right.args:
                if isinstance(a, _ZT):
                    results.append(a)
                else:
                    results.append(
                        compose_algebraic_pure_tensors(left, a)
                    )
            return _compose_reassemble(results, right.tensor_shape)
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
    """Sum of AlgebraicPureTensors (and/or AlgebraicZeroTensors) sharing the same tensor shape.

    Built on top of ``Basic`` (not ``Add``) to avoid the is_commutative
    descriptor conflict in ``AssocOp._from_args``, while still exposing
    ``is_Add = True`` so that the wider SymPy ecosystem recognises this
    as an additive expression.

    Shape enforcement and AlgebraicZeroTensor handling are done in ``__new__``.
    The flattening logic is kept intentionally simple so that a future
    ``.simplify`` extension can layer noncommutative-polynomial factoring
    on top without fighting against ``Add.flatten`` internals.

    Tensor shapes are full sequences of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for the product ``A_3x4 ⊗ C_4x5``.
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
        flat, shape, zero_term = cls._flatten_args(args)

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

        # If an AlgebraicZeroTensor was provided, always keep it as an anchor.
        # Only unwrap to a single term when there is exactly one non-zero
        # term AND no AlgebraicZeroTensor was explicitly given.
        if len(flat) == 1 and zero_term is None and not isinstance(flat[0], AlgebraicZeroTensor):
            return flat[0]

        # Append AlgebraicZeroTensor if present (for shape anchoring)
        if zero_term is not None:
            flat.append(zero_term)

        obj = Basic.__new__(cls, *flat)
        return obj

    @classmethod
    def _flatten_args(cls, args):
        """Flatten nested AlgebraicTensors, validate shapes, collect AlgebraicZeroTensors.

        Returns
        -------
        (terms : list, shape : tuple | None, zero_term : AlgebraicZeroTensor | None)

        *shape* is a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.
        """
        shape = None
        zero_term = None
        terms = []

        work = list(args)
        while work:
            o = work.pop()

            if isinstance(o, AlgebraicZeroTensor):
                candidate = o.shape
                if shape is None:
                    shape = candidate
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                zero_term = o
                continue

            if isinstance(o, AlgebraicTensor):
                work.extend(o.args)
                continue

            if isinstance(o, AlgebraicPureTensor):
                candidate = o.tensor_shape
                if shape is None:
                    shape = candidate
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                terms.append(o)
                continue

            if isinstance(o, Mul) and not isinstance(o, AlgebraicPureTensor):
                candidate = _tensor_shape_of(o)
                if candidate is not None:
                    if shape is None:
                        shape = candidate
                    elif candidate != shape:
                        raise ShapeMismatchError(
                            f"Cannot add tensors of different shapes: "
                            f"{shape} vs {candidate}"
                        )
                terms.append(o)
                continue

            if o.is_Number:
                terms.append(o)
                continue

            # Catch-all: bare matrix-like objects or anything else
            candidate = _tensor_shape_of(o)
            if candidate is not None:
                if shape is None:
                    shape = candidate
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
            terms.append(o)

        return terms, shape, zero_term

    @property
    def tensor_shape(self):
        """Shape shared by every term in this sum.

        Returns a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.
        """
        for arg in self.args:
            ts = _tensor_shape_of(arg)
            if ts is not None:
                return ts
        raise AttributeError("Cannot determine tensor_shape")

    @property
    def terms(self):
        """Non-zero, non-coefficient terms in this sum."""
        return tuple(
            a for a in self.args
            if not isinstance(a, Number) and not isinstance(a, AlgebraicZeroTensor)
        )

    def has_zero_term(self):
        """Return True if an AlgebraicZeroTensor anchors this sum."""
        return any(isinstance(a, AlgebraicZeroTensor) for a in self.args)

    # ---- Arithmetic delegation ----

    def __neg__(self):
        """Negate each term in the sum."""
        return AlgebraicTensor(*(-a for a in self.args))

    def _compose_with_term(self, other):
        """Compose this AlgebraicTensor with a single term (AlgebraicPureTensor or bare matrix).

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
                        results.append(coeff * comp)
                        break
                else:
                    results.append(a)
            elif hasattr(a, "shape"):
                comp = compose_algebraic_pure_tensors(a, other)
                results.append(comp)
            else:
                results.append(a)

        return _compose_reassemble(results, self.tensor_shape)

    def __add__(self, other):
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
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
        """
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(other * a for a in self.args))
        return compose_algebraic_tensors(other, self)

    # ---- Factorization helpers (infrastructure for future .simplify) ----

    def as_common_left(self):
        """Extract common left factors from all AlgebraicPureTensor terms.

        Handles factors wrapped in MatMul with scalar coefficients (e.g. 2*A
        and 3*A are recognized as sharing the same matrix base A).

        Returns
        -------
        (left_factors, rest, right_factors)
            *left_factors* is a tuple of shared leading matrix bases (with
            scalar coefficients stripped out).
            *rest* is an AlgebraicTensor of the middle parts (or a single
            PureTensor / sum of PureTensors when there is only one term).
            *right_factors* is always () for left extraction.

        If no common left factors exist, returns ((), self, ()).
        """
        from sympy.tensor.algebraic.simplify import _proportionality_ratio, \
            _matrix_base

        pure_terms = [
            a for a in self.args
            if isinstance(a, AlgebraicPureTensor) or (isinstance(a, Mul) and
                any(isinstance(f, AlgebraicPureTensor) for f in a.args))
        ]
        if not pure_terms:
            return ((), self, ())

        def _pure_part(expr):
            if isinstance(expr, AlgebraicPureTensor):
                return expr.factors
            if isinstance(expr, Mul):
                for f in expr.args:
                    if isinstance(f, AlgebraicPureTensor):
                        return f.factors
            return ()

        factor_lists = [_pure_part(t) for t in pure_terms]
        if not factor_lists[0]:
            return ((), self, ())

        common_len = 0
        for pos in range(min(len(fl) for fl in factor_lists)):
            refs = factor_lists[0][pos]
            if all(_proportionality_ratio(refs, fl[pos]) is not None
                   for fl in factor_lists[1:]):
                common_len += 1
            else:
                break

        if common_len == 0:
            return ((), self, ())

        # Extract matrix base (scalar-stripped) for the common left factors
        left_factors = tuple(_matrix_base(f) for f in factor_lists[0][:common_len])

        def _scalar_from_factor(f):
            """Extract scalar coefficient from a factor like MatMul(k, M) -> k."""
            from sympy.core.mul import Mul as _Mul
            if not isinstance(f, _Mul):
                return S.One
            commutative_parts = [a for a in f.args
                                 if hasattr(a, 'is_commutative') and a.is_commutative]
            if not commutative_parts:
                return S.One
            return Mul(*commutative_parts, evaluate=True)

        # Collect proportionality ratios for each term's common factors
        middles = []
        for t, fl in zip(pure_terms, factor_lists):
            # Accumulate coefficient from scalar parts of common left factors
            term_coeff = S.One
            for i in range(common_len):
                if t is pure_terms[0]:
                    # Pivot term: extract scalar from its own factors
                    term_coeff = term_coeff * _scalar_from_factor(fl[i])
                else:
                    # Non-pivot: ratio accounts for scalar difference vs pivot base
                    r = _proportionality_ratio(factor_lists[0][i], fl[i])
                    if r is not None:
                        term_coeff = term_coeff / r

            mid = fl[common_len:]
            # Also extract scalars from non-common factors into the coefficient
            for f in mid:
                term_coeff = term_coeff * _scalar_from_factor(f)
            # Replace mid factors with their base matrices
            mid = [_matrix_base(f) for f in mid]

            if isinstance(t, Mul) and not isinstance(t, AlgebraicPureTensor):
                coeff = S.One
                for f in t.args:
                    if isinstance(f, AlgebraicPureTensor):
                        coeff = t / f
                        break
                term_coeff = term_coeff * coeff

            if mid:
                if term_coeff is S.One:
                    middles.append(AlgebraicPureTensor(*mid))
                else:
                    middles.append(term_coeff * AlgebraicPureTensor(*mid))
            else:
                middles.append(term_coeff)

        if len(middles) == 1:
            rest = middles[0]
        else:
            rest = AlgebraicTensor(*middles)

        return (left_factors, rest, ())

    def as_common_right(self):
        """Extract common right factors from all AlgebraicPureTensor terms.

        Symmetric to :meth:`as_common_left` but for trailing factors.
        Handles factors wrapped in MatMul with scalar coefficients.

        Returns
        -------
        (left_factors, rest, right_factors)
        """
        from sympy.tensor.algebraic.simplify import _proportionality_ratio, \
            _matrix_base

        pure_terms = [
            a for a in self.args
            if isinstance(a, AlgebraicPureTensor) or (isinstance(a, Mul) and
                any(isinstance(f, AlgebraicPureTensor) for f in a.args))
        ]
        if not pure_terms:
            return ((), self, ())

        def _pure_part(expr):
            if isinstance(expr, AlgebraicPureTensor):
                return expr.factors
            if isinstance(expr, Mul):
                for f in expr.args:
                    if isinstance(f, AlgebraicPureTensor):
                        return f.factors
            return ()

        factor_lists = [_pure_part(t) for t in pure_terms]
        if not factor_lists[0]:
            return ((), self, ())

        common_len = 0
        for pos in range(1, min(len(fl) for fl in factor_lists) + 1):
            refs = factor_lists[0][-pos]
            if all(_proportionality_ratio(refs, fl[-pos]) is not None
                   for fl in factor_lists[1:]):
                common_len += 1
            else:
                break

        if common_len == 0:
            return ((), self, ())

        right_factors = tuple(
            _matrix_base(f) for f in factor_lists[0][-common_len:])

        def _scalar_from_factor(f):
            """Extract scalar coefficient from a factor like MatMul(k, M) -> k."""
            from sympy.core.mul import Mul as _Mul
            if not isinstance(f, _Mul):
                return S.One
            commutative_parts = [a for a in f.args
                                 if hasattr(a, 'is_commutative') and a.is_commutative]
            if not commutative_parts:
                return S.One
            return Mul(*commutative_parts, evaluate=True)

        middles = []
        for t, fl in zip(pure_terms, factor_lists):
            term_coeff = S.One
            for i in range(common_len):
                if t is pure_terms[0]:
                    term_coeff = term_coeff * _scalar_from_factor(fl[-common_len + i])
                else:
                    r = _proportionality_ratio(factor_lists[0][-common_len + i],
                                               fl[-common_len + i])
                    if r is not None:
                        term_coeff = term_coeff / r

            mid = fl[:-common_len] if common_len < len(fl) else ()
            # Also extract scalars from non-common factors into the coefficient
            for f in mid:
                term_coeff = term_coeff * _scalar_from_factor(f)
            mid = [_matrix_base(f) for f in mid]

            if isinstance(t, Mul) and not isinstance(t, AlgebraicPureTensor):
                coeff = S.One
                for f in t.args:
                    if isinstance(f, AlgebraicPureTensor):
                        coeff = t / f
                        break
                term_coeff = term_coeff * coeff

            if mid:
                if term_coeff is S.One:
                    middles.append(AlgebraicPureTensor(*mid))
                else:
                    middles.append(term_coeff * AlgebraicPureTensor(*mid))
            else:
                middles.append(term_coeff)

        if len(middles) == 1:
            rest = middles[0]
        else:
            rest = AlgebraicTensor(*middles)

        return ((), rest, right_factors)

    def as_common_factors(self):
        """Extract both common left and right factors simultaneously.

        Returns
        -------
        (left_factors, middle_expr, right_factors)
        """
        left_factors, rest, _ = self.as_common_left()
        if isinstance(rest, AlgebraicTensor):
            _, rest2, right_factors = rest.as_common_right()
        else:
            rest2, right_factors = rest, ()
        return (left_factors, rest2, right_factors)

    def __str__(self):
        if not self.args:
            return ""
        parts = []
        first = True
        for a in self.args:
            s = str(a)
            if first:
                parts.append(s)
                first = False
            elif s.startswith("-"):
                parts.append(f"- {s[1:]}")
            else:
                parts.append(f"+ {s}")
        return " ".join(parts)

    def simplify(self, **kwargs):
        """Simplify this AlgebraicTensor.

        Combines like terms, extracts common left/right factors, and uses
        SymPy's simplification machinery on coefficients and middle parts.
        See :func:`sympy.tensor.algebraic.simplify.tensorsimplify` for details.
        """
        from sympy.tensor.algebraic.simplify import tensorsimplify
        return tensorsimplify(self, **kwargs)

    def __repr__(self):
        return f"AlgebraicTensor({', '.join(repr(a) for a in self.args)})"


# Register AlgebraicTensor with the SymPy add dispatcher so that
# expressions like (A + B) where A/B involve PureTensor or AlgebraicTensor
# are routed through AlgebraicTensor.__new__ instead of Add.__new__.
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)

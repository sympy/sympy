from __future__ import annotations

from collections import defaultdict

from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify


def _get_sympy_simplify():
    """Lazy import of SymPy's top-level simplify to avoid circular imports."""
    from sympy.simplify.simplify import simplify as _s
    return _s


def _tensor_key(expr):
    """Return a hashable key identifying the tensor factor structure of *expr*.

    For a ``PureTensor``, the key is the tuple of ``id(f)`` for each factor,
    so that two PureTensors with the exact same factor objects share a key.

    For ``Mul(coeff, PureTensor)`` the PureTensor part is extracted first.
    For anything else (plain matrix, symbol, number) the key is just the
    object identity.
    """
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    if isinstance(expr, PureTensor):
        return ("pt", tuple(id(f) for f in expr.factors))
    if isinstance(expr, Mul) and not isinstance(expr, PureTensor):
        for f in expr.args:
            if isinstance(f, PureTensor):
                return ("pt", tuple(id(f2) for f2 in f.factors))
    return ("other", id(expr))


def _extract_coeff_and_pt(expr):
    """Return (coeff, unit_puretensor) from *expr*.

    Returns the fully combined coefficient and a **unit** PureTensor
    (coefficient stripped to S.One) so the caller can re-apply a new
    coefficient without doubling up.

    For non-tensor expressions return (expr, None).
    """
    from sympy.tensor.algebraic.pure_tensor import PureTensor

    if isinstance(expr, PureTensor):
        coeff = expr._get_coeff()
        factors = expr.factors
        # Build a unit PureTensor (coeff = 1).
        if len(factors) == 1:
            unit_pt = factors[0]
        else:
            unit_pt = PureTensor(*factors)
        return (coeff, unit_pt)

    if isinstance(expr, Mul) and not isinstance(expr, PureTensor):
        for f in expr.args:
            if isinstance(f, PureTensor):
                outer_coeff = Mul(*[a for a in expr.args if a is not f])
                inner_coeff = f._get_coeff()
                combined_coeff = outer_coeff * inner_coeff
                factors = f.factors
                if len(factors) == 1:
                    unit_pt = factors[0]
                else:
                    unit_pt = PureTensor(*factors)
                return (combined_coeff, unit_pt)

    return (expr, None)


def _make_scaled_pt(coeff, unit_pt):
    """Construct ``coeff * unit_pt`` as a PureTensor without going through
    SymPy's general dispatch (which would produce MatMul).

    *unit_pt* is either a bare matrix-like object (single-factor PureTensor
    unwrapped) or a multi-factor PureTensor with coefficient S.One.
    """
    from sympy.tensor.algebraic.pure_tensor import PureTensor

    if coeff is S.One:
        return unit_pt
    if coeff is S.Zero:
        if isinstance(unit_pt, PureTensor):
            from sympy.tensor.algebraic.zero_tensor import ZeroTensor
            return ZeroTensor(unit_pt.tensor_shape)
        if hasattr(unit_pt, "shape"):
            from sympy.tensor.algebraic.zero_tensor import ZeroTensor
            return ZeroTensor((unit_pt.shape,))
        return coeff

    # Determine the factor list.
    if isinstance(unit_pt, PureTensor):
        factors = unit_pt.factors
    elif hasattr(unit_pt, "shape"):
        factors = (unit_pt,)
    else:
        return coeff * unit_pt

    if len(factors) == 1:
        return PureTensor(coeff, factors[0])
    return PureTensor(coeff, *factors)


def _tensor_multiply_left(factor, expr):
    """Return ``factor ⊗ expr`` where *expr* is PureTensor, AlgebraicTensor,
    or a plain matrix-like object."""
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.pure_tensor import PureTensor

    if isinstance(expr, AlgebraicTensor):
        return AlgebraicTensor(*(factor * a for a in expr.args))
    if isinstance(expr, PureTensor):
        return PureTensor(factor, expr)
    # Plain matrix-like: create a two-factor PureTensor.
    return PureTensor(factor, expr)


def _tensor_multiply_right(expr, factor):
    """Return ``expr ⊗ factor`` where *expr* is PureTensor, AlgebraicTensor,
    or a plain matrix-like object."""
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.pure_tensor import PureTensor

    if isinstance(expr, AlgebraicTensor):
        return AlgebraicTensor(*(a * factor for a in expr.args))
    if isinstance(expr, PureTensor):
        return PureTensor(expr, factor)
    # Plain matrix-like: create a two-factor PureTensor.
    return PureTensor(expr, factor)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def tensorsimplify(expr, **kwargs):
    """Simplify an algebraic-tensor expression.

    Dispatches to the appropriate handler for ``PureTensor``,
    ``AlgebraicTensor``, ``ZeroTensor``, or falls back to SymPy's
    general :func:`simplify` for everything else.

    Parameters
    ----------
    expr : Any
        The expression to simplify.
    **kwargs
        Extra keyword arguments forwarded to :func:`sympy.simplify.simplify`
        when simplifying individual coefficients.

    Returns
    -------
    simplified expression (same type as *expr* when possible)
    """
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor

    if isinstance(expr, ZeroTensor):
        return expr
    if isinstance(expr, AlgebraicTensor):
        return _simplify_algebraic_tensor(expr, **kwargs)
    if isinstance(expr, PureTensor):
        return _simplify_pure_tensor(expr, **kwargs)

    # Non-tensor expression — use SymPy's general simplify.
    _s = _get_sympy_simplify()
    return _s(expr, **kwargs)


# ---------------------------------------------------------------------------
# PureTensor simplification
# ---------------------------------------------------------------------------

def _simplify_pure_tensor(pt, **kwargs):
    """Simplify a single PureTensor.

    Strategy
    --------
    1. Simplify the leading coefficient with SymPy's ``simplify``.
    2. Attempt to simplify each factor individually (useful when factors
       are MatrixExpr containing simplifiable sub-expressions).
    3. Reconstruct the PureTensor with the (possibly changed) pieces.
    """
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor

    _s = _get_sympy_simplify()

    coeff = pt._get_coeff()
    factors = pt.factors

    # Simplify coefficient
    new_coeff = _s(coeff, **kwargs)
    if new_coeff is S.Zero:
        return ZeroTensor(pt.tensor_shape)

    # Simplify individual factors
    new_factors = []
    changed = new_coeff != coeff
    for f in factors:
        sf = _s(f, **kwargs)
        if sf != f:
            changed = True
        new_factors.append(sf)

    if not changed:
        return pt

    if len(new_factors) == 0:
        return new_coeff

    if new_coeff is S.One:
        if len(new_factors) == 1:
            return new_factors[0]
        return PureTensor(*new_factors)
    if len(new_factors) == 1:
        return PureTensor(new_coeff, new_factors[0])
    return PureTensor(new_coeff, *new_factors)


# ---------------------------------------------------------------------------
# AlgebraicTensor simplification
# ---------------------------------------------------------------------------

def _simplify_algebraic_tensor(at, **kwargs):
    """Simplify an AlgebraicTensor (sum of same-shape tensor terms).

    The simplification pipeline has three phases that each leverage SymPy's
    existing simplification machinery:

    1. **Combine like terms** -- group terms that share the same tensor
       factor sequence and add their coefficients together (using SymPy's
       ``Add`` / ``simplify`` on the accumulated coefficient).

    2. **Extract common left / right factors** -- use the existing
       ``as_common_factors`` infrastructure so that the middle sub-sum can
       be simplified independently.

    3. **Recurse** on the middle sub-sum and reassemble.

    ZeroTensor anchors and S.Zero terms are handled naturally.
    """
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor, zero_tensor

    _s = _get_sympy_simplify()

    args = list(at.args)
    original_shape = at.tensor_shape

    # Strip S.Zero entries but remember ZeroTensor anchors.
    has_zero = False
    real_args = []
    for a in args:
        if isinstance(a, ZeroTensor):
            has_zero = True
        elif a is not S.Zero:
            real_args.append(a)

    if not real_args:
        if has_zero:
            z = next(a for a in args if isinstance(a, ZeroTensor))
            return z
        # All terms cancelled and no ZeroTensor anchor — create one.
        return ZeroTensor(original_shape)

    # ------------------------------------------------------------------
    # Phase 1: combine like terms
    # ------------------------------------------------------------------
    grouped = defaultdict(list)
    non_tensor = []

    for a in real_args:
        key = _tensor_key(a)
        if key[0] == "pt":
            grouped[key].append(a)
        else:
            non_tensor.append(a)

    combined = []
    for key, group in grouped.items():
        if len(group) == 1:
            # Single term — still simplify coefficient / factors
            combined.append(_simplify_pure_tensor_in_context(group[0], **kwargs))
        else:
            # Multiple terms with same factor structure: sum coefficients.
            combined.extend(_combine_like_terms(group, **kwargs))

    # Simplify non-tensor entries (plain matrices, numbers, …)
    for nt in non_tensor:
        snt = _s(nt, **kwargs)
        if snt is not S.Zero:
            combined.append(snt)

    if has_zero:
        combined.append(next(a for a in args if isinstance(a, ZeroTensor)))

    # Handle empty combined list (all terms cancelled, no ZeroTensor anchor)
    if not combined:
        return ZeroTensor(original_shape)

    # Build intermediate AlgebraicTensor for the next phases.
    if len(combined) == 1 and not isinstance(combined[0], ZeroTensor):
        result = combined[0]
    else:
        result = AlgebraicTensor(*combined, _sympify=False)

    # ------------------------------------------------------------------
    # Phase 2 + 3: extract common factors and recurse on the middle
    # ------------------------------------------------------------------
    if isinstance(result, AlgebraicTensor):
        result = _simplify_with_common_factors(result, **kwargs)

    return result


def _simplify_pure_tensor_in_context(term, **kwargs):
    """Simplify a term that carries a PureTensor (with or without coeff)."""
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor

    _s = _get_sympy_simplify()

    if isinstance(term, PureTensor):
        return _simplify_pure_tensor(term, **kwargs)

    if isinstance(term, Mul) and not isinstance(term, PureTensor):
        coeff, unit_pt = _extract_coeff_and_pt(term)
        if unit_pt is not None:
            # unit_pt is a unit PureTensor (coeff stripped) or bare matrix.
            if isinstance(unit_pt, PureTensor):
                factors = unit_pt.factors
            elif hasattr(unit_pt, "shape"):
                factors = (unit_pt,)
            else:
                return _s(term, **kwargs)

            simp_factors = tuple(_s(f, **kwargs) for f in factors)
            simp_coeff = _s(coeff, **kwargs)
            if simp_coeff is S.Zero:
                shape = tuple(f.shape for f in factors)
                return ZeroTensor(shape)
            if simp_coeff is S.One:
                if len(simp_factors) == 1:
                    return simp_factors[0]
                return PureTensor(*simp_factors)
            return _make_scaled_pt(simp_coeff,
                                   PureTensor(*simp_factors) if len(simp_factors) > 1 else simp_factors[0])

    # Fallback: generic simplify
    return _s(term, **kwargs)


def _combine_like_terms(group, **kwargs):
    """Combine a list of terms sharing the same factor sequence.

    Each element of *group* is either a PureTensor or Mul(coeff, PureTensor).
    We extract the (fully combined) coefficient from each, sum all coefficients
    into a single SymPy expression, simplify that sum, and return the product
    of the simplified sum with the unit PureTensor.

    Returns a list (length 0 when the summed coefficient is zero, or length 1
    with the combined term).
    """
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor

    _s = _get_sympy_simplify()

    # All terms in the group share the same factor sequence; pick the
    # representative UNIT PureTensor from the first entry.
    _, representative_unit_pt = _extract_coeff_and_pt(group[0])
    assert representative_unit_pt is not None

    coeffs = []
    for term in group:
        c, _ = _extract_coeff_and_pt(term)
        coeffs.append(c)

    # Sum coefficients using SymPy's Add, then simplify.
    total_coeff = _s(Add(*coeffs), **kwargs)

    if total_coeff is S.Zero:
        return []  # terms cancel — caller will handle emptiness

    if total_coeff is S.One:
        return [representative_unit_pt]

    # Use _make_scaled_pt to avoid SymPy dispatch to MatMul.
    return [_make_scaled_pt(total_coeff, representative_unit_pt)]


def _simplify_with_common_factors(at, **kwargs):
    """Extract common left/right factors, simplify the middle, reassemble."""
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor

    _s = _get_sympy_simplify()

    left_factors, middle, right_factors = at.as_common_factors()

    if not left_factors and not right_factors:
        # Nothing to extract — return as-is (already combined).
        return at

    if isinstance(middle, AlgebraicTensor):
        middle = _simplify_algebraic_tensor(middle, **kwargs)
    elif isinstance(middle, PureTensor):
        middle = _simplify_pure_tensor(middle, **kwargs)
    else:
        middle = _s(middle, **kwargs)

    # Reassemble: left_factors ⊗ middle ⊗ right_factors
    # Use _tensor_multiply_left/right so AlgebraicTensor middle parts
    # are handled by distributing the factor across all terms.
    result = middle
    if left_factors:
        for f in left_factors:
            result = _tensor_multiply_left(f, result)
    if right_factors:
        for f in right_factors:
            result = _tensor_multiply_right(result, f)

    return result

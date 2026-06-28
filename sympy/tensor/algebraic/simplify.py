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

    For an ``AlgebraicPureTensor``, the key is the tuple of ``id(f)`` for each factor,
    so that two AlgebraicPureTensors with the exact same factor objects share a key.

    For ``Mul(coeff, AlgebraicPureTensor)`` the AlgebraicPureTensor part is extracted first.
    For anything else (plain matrix, symbol, number) the key is just the
    object identity.
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    if isinstance(expr, AlgebraicPureTensor):
        return ("pt", tuple(id(f) for f in expr.factors))
    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                return ("pt", tuple(id(f2) for f2 in f.factors))
    return ("other", id(expr))


def _extract_coeff_and_pt(expr):
    """Return (coeff, unit_puretensor) from *expr*.

    Returns the fully combined coefficient and a **unit** AlgebraicPureTensor
    (coefficient stripped to S.One) so the caller can re-apply a new
    coefficient without doubling up.

    For non-tensor expressions return (expr, None).
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if isinstance(expr, AlgebraicPureTensor):
        coeff = expr._get_coeff()
        factors = expr.factors
        # Build a unit AlgebraicPureTensor (coeff = 1).
        if len(factors) == 1:
            unit_pt = factors[0]
        else:
            unit_pt = AlgebraicPureTensor(*factors)
        return (coeff, unit_pt)

    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                outer_coeff = Mul(*[a for a in expr.args if a is not f])
                inner_coeff = f._get_coeff()
                combined_coeff = outer_coeff * inner_coeff
                factors = f.factors
                if len(factors) == 1:
                    unit_pt = factors[0]
                else:
                    unit_pt = AlgebraicPureTensor(*factors)
                return (combined_coeff, unit_pt)

    return (expr, None)


def _make_scaled_pt(coeff, unit_pt):
    """Construct ``coeff * unit_pt`` as an AlgebraicPureTensor without going through
    SymPy's general dispatch (which would produce MatMul).

    *unit_pt* is either a bare matrix-like object (single-factor AlgebraicPureTensor
    unwrapped) or a multi-factor AlgebraicPureTensor with coefficient S.One.
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if coeff is S.One:
        return unit_pt
    if coeff is S.Zero:
        if isinstance(unit_pt, AlgebraicPureTensor):
            from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
            return AlgebraicZeroTensor(unit_pt.tensor_shape)
        if hasattr(unit_pt, "shape"):
            from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
            return AlgebraicZeroTensor((unit_pt.shape,))
        return coeff

    # Determine the factor list.
    if isinstance(unit_pt, AlgebraicPureTensor):
        factors = unit_pt.factors
    elif hasattr(unit_pt, "shape"):
        factors = (unit_pt,)
    else:
        return coeff * unit_pt

    if len(factors) == 1:
        return AlgebraicPureTensor(coeff, factors[0])
    return AlgebraicPureTensor(coeff, *factors)


def _tensor_multiply_left(factor, expr):
    """Return ``factor ⊗ expr`` where *expr* is AlgebraicPureTensor, AlgebraicTensor,
    a plain matrix-like object, or a number."""
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if isinstance(expr, AlgebraicTensor):
        return AlgebraicTensor(*(factor * a for a in expr.args))
    if isinstance(expr, AlgebraicPureTensor):
        return AlgebraicPureTensor(factor, expr)
    if isinstance(expr, Number):
        if expr is S.One:
            return factor
        if expr is S.Zero:
            from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
            return AlgebraicZeroTensor((factor.shape,))
        return AlgebraicPureTensor(expr, factor)
    # Plain matrix-like: create a two-factor AlgebraicPureTensor.
    return AlgebraicPureTensor(factor, expr)


def _tensor_multiply_right(expr, factor):
    """Return ``expr ⊗ factor`` where *expr* is AlgebraicPureTensor, AlgebraicTensor,
    a plain matrix-like object, or a number."""
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if isinstance(expr, AlgebraicTensor):
        return AlgebraicTensor(*(a * factor for a in expr.args))
    if isinstance(expr, AlgebraicPureTensor):
        return AlgebraicPureTensor(expr, factor)
    if isinstance(expr, Number):
        if expr is S.One:
            return factor
        if expr is S.Zero:
            from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
            return AlgebraicZeroTensor((factor.shape,))
        return AlgebraicPureTensor(expr, factor)
    # Plain matrix-like: create a two-factor AlgebraicPureTensor.
    return AlgebraicPureTensor(expr, factor)


# ---------------------------------------------------------------------------
# Proportionality checking helpers
# ---------------------------------------------------------------------------

def _is_true_matrix(expr):
    """Check if expr is a matrix-like object with shape != (1, 1)."""
    shape = getattr(expr, 'shape', None)
    if shape is None:
        return False
    return shape != (1, 1)


def _is_1x1_matrix(expr):
    """Check if expr is a 1x1 matrix-like object."""
    shape = getattr(expr, 'shape', None)
    if shape is None:
        return False
    return shape == (1, 1)


def _get_1x1_element(expr):
    """Extract the single element from a 1x1 matrix expression.

    For MatrixSymbol of shape (1,1), returns the indexed element expr[0,0].
    For dense Matrix, returns the single entry.
    Falls back to expr itself if extraction fails.
    """
    if hasattr(expr, 'shape') and expr.shape == (1, 1):
        if hasattr(expr, '__getitem__'):
            try:
                elem = expr[0, 0]
                return elem
            except (TypeError, IndexError, NotImplementedError):
                pass
    return expr


def _matrix_proportionality_ratio(m1, m2):
    """Check element-wise proportionality for true matrices (shape != (1,1)).

    Returns the commutative ratio k such that m1 = k * m2, or None.
    """
    if m1.shape != m2.shape:
        return None

    ratio = None
    _s = _get_sympy_simplify()

    for i in range(m1.shape[0]):
        for j in range(m1.shape[1]):
            e1 = m1[i, j]
            e2 = m2[i, j]

            e1_is_zero = (e1 == S.Zero) or (e1.is_zero is True)
            e2_is_zero = (e2 == S.Zero) or (e2.is_zero is True)

            if e1_is_zero and e2_is_zero:
                continue
            if e1_is_zero or e2_is_zero:
                return None

            r = _s(e1 / e2)
            if ratio is None:
                ratio = r
            elif _s(ratio - r) != S.Zero:
                return None

    if ratio is None:
        return None
    return ratio


def _nc_symbol_proportionality_ratio(f1, f2):
    """Check proportionality for 1x1 matrices wrapping noncommutative symbols.

    Returns a commutative ratio k such that f1 = k * f2, or None.
    """
    elem1 = _get_1x1_element(f1)
    elem2 = _get_1x1_element(f2)

    if elem1 == elem2:
        return S.One

    _s = _get_sympy_simplify()

    # Try direct division: elem1 / elem2
    try:
        ratio = _s(elem1 / elem2)
        if hasattr(ratio, 'is_commutative') and ratio.is_commutative:
            return ratio
        if isinstance(ratio, Number):
            return ratio
    except (ZeroDivisionError, TypeError, AttributeError):
        pass

    # Try the reverse direction
    try:
        ratio_inv = _s(elem2 / elem1)
        if hasattr(ratio_inv, 'is_commutative') and ratio_inv.is_commutative:
            if ratio_inv != S.Zero:
                return _s(S.One / ratio_inv)
        if isinstance(ratio_inv, Number) and ratio_inv != S.Zero:
            return _s(S.One / ratio_inv)
    except (ZeroDivisionError, TypeError, AttributeError):
        pass

    # Try factoring: check if elem1 = k * elem2 by extracting commutative parts
    try:
        from sympy.core.mul import Mul as _Mul
        if isinstance(elem1, _Mul):
            commutative_part = _Mul(*[a for a in elem1.args if a.is_commutative], evaluate=True)
            noncomm_part = _Mul(*[a for a in elem1.args if not a.is_commutative], evaluate=True)
            if noncomm_part == elem2 or _s(noncomm_part - elem2) == S.Zero:
                return commutative_part
        if isinstance(elem2, _Mul):
            commutative_part = _Mul(*[a for a in elem2.args if a.is_commutative], evaluate=True)
            noncomm_part = _Mul(*[a for a in elem2.args if not a.is_commutative], evaluate=True)
            if noncomm_part == elem1 or _s(noncomm_part - elem1) == S.Zero:
                return _s(S.One / commutative_part)
    except (TypeError, AttributeError):
        pass

    return None


def _matrix_base(factor):
    """Extract the pure matrix base from a factor, stripping scalar coefficients.

    If *factor* is a Mul/MatMul with a commutative (scalar) coefficient
    and a single matrix-like part, returns the matrix part.  Otherwise
    returns *factor* unchanged.

    Examples
    --------
    >>> _matrix_base(A)                # A is a MatrixSymbol
    A
    >>> _matrix_base(MatMul(2, A))    # MatMul with scalar coeff
    A
    """
    from sympy.core.mul import Mul as _Mul
    if not isinstance(factor, _Mul):
        return factor

    commutative_parts = []
    noncommutative_parts = []
    for a in factor.args:
        if hasattr(a, 'is_commutative') and a.is_commutative:
            commutative_parts.append(a)
        else:
            noncommutative_parts.append(a)

    if len(noncommutative_parts) == 1 and commutative_parts:
        return noncommutative_parts[0]
    return factor


def _proportionality_ratio(factor1, factor2):
    """Check if factor1 and factor2 are proportional.

    Returns the commutative ratio k such that factor1 = k * factor2,
    or None if they are not proportional.

    - For true matrices (shape != (1,1)): checks all nonzero elements share
      the same ratio. Avoids division by zero.
    - For 1x1 matrices: checks if there is a commutative proportionality
      constant between the wrapped noncommutative symbols.
    """
    if factor1 == factor2:
        return S.One

    if _is_1x1_matrix(factor1) and _is_1x1_matrix(factor2):
        return _nc_symbol_proportionality_ratio(factor1, factor2)

    if _is_true_matrix(factor1) and _is_true_matrix(factor2):
        return _matrix_proportionality_ratio(factor1, factor2)

    return None


# ---------------------------------------------------------------------------
# Proportionality factoring
# ---------------------------------------------------------------------------

def _extract_pt_and_coeff(term):
    """Extract (unit_pt, coefficient) from a term.

    Returns (unit_pt, coeff) where unit_pt is either a bare matrix-like object,
    an AlgebraicPureTensor, or the term itself (when no extraction is possible).
    coeff is the extracted commutative coefficient (S.One if none).

    Handles:
    - AlgebraicPureTensor directly
    - Mul(coeff, AlgebraicPureTensor)
    - MatMul(coeff, *matrix_factors) with direct matrix factors
    - Anything else -> (term, S.One)
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if isinstance(term, AlgebraicPureTensor):
        return (term, term._get_coeff())

    if isinstance(term, Mul) and not isinstance(term, AlgebraicPureTensor):
        # Check for AlgebraicPureTensor in args
        for f in term.args:
            if isinstance(f, AlgebraicPureTensor):
                outer_coeff = Mul(*[a for a in term.args if a is not f])
                inner_coeff = f._get_coeff()
                return (f, outer_coeff * inner_coeff)

        # Handle MatMul with direct matrix factors: extract commutative coeff
        # and matrix-like factors (objects with .shape attribute)
        matrix_factors = []
        commutative_parts = []
        for a in term.args:
            if isinstance(a, AlgebraicPureTensor):
                continue  # Already handled above
            if hasattr(a, 'shape') and a.shape is not None:
                matrix_factors.append(a)
            elif hasattr(a, 'is_commutative') and a.is_commutative:
                commutative_parts.append(a)

        if matrix_factors:
            coeff = Mul(*commutative_parts, evaluate=True) if commutative_parts else S.One
            coeff = _get_sympy_simplify()(coeff)
            return (None, coeff, matrix_factors)

    return (term, S.One)


def _build_pt(coeff, factors):
    """Build an AlgebraicPureTensor from a coefficient and a list of factors.

    Handles special cases: coeff=1 drops coefficient, single factor with
    coeff=1 unwraps to bare factor, coeff=0 produces AlgebraicZeroTensor.
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    if coeff is S.Zero:
        if len(factors) > 0 and hasattr(factors[0], 'shape'):
            return AlgebraicZeroTensor(tuple(f.shape for f in factors))
        return coeff

    if not factors:
        return coeff

    if coeff is S.One:
        if len(factors) == 1:
            return factors[0]
        return AlgebraicPureTensor(*factors)

    return AlgebraicPureTensor(coeff, *factors)
 

def proportionality_factoring(at):
    """Simplify an AlgebraicTensor by merging proportional AlgebraicPureTensor terms.

    Algorithm
    ---------
    1. Pick a pivot AlgebraicPureTensor from the sum.
    2. Compare the pivot with every other AlgebraicPureTensor in the sum.
    3. For each pair, walk through their tensor factors slot by slot:

       - If ALL factor slots are proportional (each factor of the selected
         term is a commutative multiple of the corresponding factor of the
         pivot), extract all ratios, multiply them into the selected term's
         prefactor, and add the two terms by summing their prefactors.

       - If EXACTLY ONE factor slot is not proportional, gather the ratios
         from all proportional slots, multiply them into the selected term's
         prefactor, and create a new combined AlgebraicPureTensor where the
         non-proportional slot holds the linear combination of the two factors
         weighted by their respective (accumulated) prefactors.

    4. After a successful merge, reset the pivot to the beginning of the list.
    5. If the pivot reaches the end without any merge, increment the pivot.

    Parameters
    ----------
    at : AlgebraicTensor
        The tensor sum to simplify.

    Returns
    -------
    AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or other
        The simplified expression.
    """
    from sympy.core.add import Add as _Add
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    args = list(at.args)

    # Separate AlgebraicPureTensor (and Mul containing one) from other terms
    tensor_terms = []
    non_tensor = []
    has_zero = False

    for a in args:
        if isinstance(a, AlgebraicZeroTensor):
            has_zero = True
            non_tensor.append(a)
        elif isinstance(a, AlgebraicPureTensor):
            tensor_terms.append(a)
        elif isinstance(a, Mul) and not isinstance(a, AlgebraicPureTensor):
            has_pt = any(isinstance(f, AlgebraicPureTensor) for f in a.args)
            has_matrix = any(hasattr(f, 'shape') for f in a.args)
            if has_pt or has_matrix:
                tensor_terms.append(a)
            else:
                non_tensor.append(a)
        elif a is not S.Zero:
            non_tensor.append(a)

    if len(tensor_terms) < 2:
        return at

    # Convert to entries with (coeff, factors) representation
    entries = []
    for t in tensor_terms:
        extracted = _extract_pt_and_coeff(t)
        if len(extracted) == 3:
            # MatMul with direct matrix factors: (None, coeff, [factors])
            _, coeff, factors = extracted
            entries.append({
                'coeff': coeff,
                'factors': factors,
            })
        else:
            pt, coeff = extracted
            if isinstance(pt, AlgebraicPureTensor):
                entries.append({
                    'coeff': coeff,
                    'factors': list(pt.factors),
                })
            elif hasattr(pt, 'shape'):
                entries.append({
                    'coeff': coeff,
                    'factors': [pt],
                })
            else:
                entries.append({
                    'coeff': coeff,
                    'factors': [pt],
                })

    changed = True
    while changed and len(entries) >= 2:
        changed = False
        pivot_idx = 0

        while pivot_idx < len(entries):
            merged = False
            pivot = entries[pivot_idx]

            for sel_idx in range(pivot_idx + 1, len(entries)):
                selected = entries[sel_idx]

                if len(pivot['factors']) != len(selected['factors']):
                    continue

                # Compare factors slot by slot
                ratios = []
                diff_slot = None

                for slot in range(len(pivot['factors'])):
                    r = _proportionality_ratio(
                        selected['factors'][slot],
                        pivot['factors'][slot]
                    )
                    if r is not None:
                        ratios.append(r)
                    else:
                        if diff_slot is not None:
                            ratios = None
                            break
                        diff_slot = slot
                        ratios.append(None)

                if ratios is None:
                    continue

                # Product of all proportional ratios
                combined_ratio = S.One
                for r in ratios:
                    if r is not None:
                        combined_ratio = combined_ratio * r

                add_coeff = _s(combined_ratio * selected['coeff'])

                if diff_slot is None:
                    # All factors proportional -> merge by adding coefficients
                    new_coeff = _s(pivot['coeff'] + add_coeff)

                    if new_coeff is S.Zero:
                        entries.pop(sel_idx)
                        entries.pop(pivot_idx)
                        pivot_idx = 0
                        changed = True
                        merged = True
                        break

                    pivot['coeff'] = new_coeff
                    entries.pop(sel_idx)
                    pivot_idx = 0
                    changed = True
                    merged = True
                    break

                else:
                    # Exactly one non-proportional slot -> linear combination
                    new_coeff_pivot = _s(pivot['coeff'])
                    new_coeff_sel = add_coeff

                    pivot_factor = pivot['factors'][diff_slot]
                    sel_factor = selected['factors'][diff_slot]

                    combined_factor = _s(
                        _Add(new_coeff_pivot * pivot_factor,
                             new_coeff_sel * sel_factor)
                    )

                    if combined_factor is S.Zero:
                        entries.pop(sel_idx)
                        entries.pop(pivot_idx)
                        pivot_idx = 0
                        changed = True
                        merged = True
                        break

                    new_factors = list(pivot['factors'])
                    new_factors[diff_slot] = combined_factor

                    result = _build_pt(S.One, new_factors)

                    if isinstance(result, AlgebraicZeroTensor):
                        entries.pop(sel_idx)
                        entries.pop(pivot_idx)
                        pivot_idx = 0
                        changed = True
                        merged = True
                        break

                    new_extracted = _extract_pt_and_coeff(result)
                    if len(new_extracted) == 3:
                        _, new_c, new_factors_list = new_extracted
                    else:
                        new_pt, new_c = new_extracted
                        if isinstance(new_pt, AlgebraicPureTensor):
                            new_factors_list = list(new_pt.factors)
                        elif hasattr(new_pt, 'shape'):
                            new_factors_list = [new_pt]
                        else:
                            new_factors_list = [new_pt]

                    entries[pivot_idx] = {
                        'coeff': new_c,
                        'factors': new_factors_list,
                    }

                    entries.pop(sel_idx)
                    pivot_idx = 0
                    changed = True
                    merged = True
                    break

            if not merged:
                pivot_idx += 1

    # Rebuild terms from entries
    rebuilt_terms = []
    for entry in entries:
        coeff = entry['coeff']
        factors = entry['factors']
        term = _build_pt(coeff, factors)
        rebuilt_terms.append(term)

    all_terms = rebuilt_terms + non_tensor

    if not all_terms:
        from sympy.tensor.algebraic.algebraic_zero_tensor import algebraic_zero_tensor
        return algebraic_zero_tensor(at.tensor_shape)

    if len(all_terms) == 1:
        return all_terms[0]

    return AlgebraicTensor(*all_terms, _sympify=False)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def tensorsimplify(expr, **kwargs):
    """Simplify an algebraic-tensor expression.

    Dispatches to the appropriate handler for ``AlgebraicPureTensor``,
    ``AlgebraicTensor``, ``AlgebraicZeroTensor``, or falls back to SymPy's
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
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    if isinstance(expr, AlgebraicZeroTensor):
        return expr
    if isinstance(expr, AlgebraicTensor):
        return _simplify_algebraic_tensor(expr, **kwargs)
    if isinstance(expr, AlgebraicPureTensor):
        return _simplify_algebraic_pure_tensor(expr, **kwargs)

    # Non-tensor expression — use SymPy's general simplify.
    _s = _get_sympy_simplify()
    return _s(expr, **kwargs)


# ---------------------------------------------------------------------------
# PureTensor simplification
# ---------------------------------------------------------------------------

def _simplify_algebraic_pure_tensor(pt, **kwargs):
    """Simplify a single AlgebraicPureTensor.

    Strategy
    --------
    1. Simplify the leading coefficient with SymPy's ``simplify``.
    2. Attempt to simplify each factor individually (useful when factors
       are MatrixExpr containing simplifiable sub-expressions).
    3. Reconstruct the PureTensor with the (possibly changed) pieces.
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    coeff = pt._get_coeff()
    factors = pt.factors

    # Simplify coefficient
    new_coeff = _s(coeff, **kwargs)
    if new_coeff is S.Zero:
        return AlgebraicZeroTensor(pt.tensor_shape)

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
        return AlgebraicPureTensor(*new_factors)
    if len(new_factors) == 1:
        return AlgebraicPureTensor(new_coeff, new_factors[0])
    return AlgebraicPureTensor(new_coeff, *new_factors)


# ---------------------------------------------------------------------------
# AlgebraicTensor simplification
# ---------------------------------------------------------------------------

def _simplify_algebraic_tensor(at, **kwargs):
    """Simplify an AlgebraicTensor (sum of same-shape tensor terms).

    The simplification pipeline:

    0. **Proportionality factoring** (primary tool) -- merges proportional
       AlgebraicPureTensor terms by combining coefficients or creating linear
       combinations at non-proportional factor slots.

    1. **Combine like terms** -- group terms that share the same tensor
       factor sequence and add their coefficients together (using SymPy's
       ``Add`` / ``simplify`` on the accumulated coefficient).

    2. **Extract common left / right factors** -- use the existing
       ``as_common_factors`` infrastructure so that the middle sub-sum can
       be simplified independently.

    3. **Recurse** on the middle sub-sum and reassemble.

    AlgebraicZeroTensor anchors and S.Zero terms are handled naturally.
    """
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor, algebraic_zero_tensor

    _s = _get_sympy_simplify()

    # ------------------------------------------------------------------
    # Phase 0: proportionality factoring (primary simplification tool)
    # ------------------------------------------------------------------
    result = proportionality_factoring(at)

    # If proportionality_factoring returned a non-AlgebraicTensor, simplify it
    if not isinstance(result, AlgebraicTensor):
        if isinstance(result, AlgebraicPureTensor):
            return _simplify_algebraic_pure_tensor(result, **kwargs)
        return result

    # If the result is still an AlgebraicTensor, run the remaining phases
    args = list(result.args)
    original_shape = result.tensor_shape

    # Strip S.Zero entries but remember AlgebraicZeroTensor anchors.
    has_zero = False
    real_args = []
    for a in args:
        if isinstance(a, AlgebraicZeroTensor):
            has_zero = True
        elif a is not S.Zero:
            real_args.append(a)

    if not real_args:
        if has_zero:
            z = next(a for a in args if isinstance(a, AlgebraicZeroTensor))
            return z
        return algebraic_zero_tensor(original_shape)

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
            combined.append(_simplify_pure_tensor_in_context(group[0], **kwargs))
        else:
            combined.extend(_combine_like_terms(group, **kwargs))

    # Simplify non-tensor entries (plain matrices, numbers, …)
    for nt in non_tensor:
        snt = _s(nt, **kwargs)
        if snt is not S.Zero:
            combined.append(snt)

    if has_zero:
        combined.append(next(a for a in args if isinstance(a, AlgebraicZeroTensor)))

    if not combined:
        return algebraic_zero_tensor(original_shape)

    if len(combined) == 1 and not isinstance(combined[0], AlgebraicZeroTensor):
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
    """Simplify a term that carries an AlgebraicPureTensor (with or without coeff)."""
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    if isinstance(term, AlgebraicPureTensor):
        return _simplify_algebraic_pure_tensor(term, **kwargs)

    if isinstance(term, Mul) and not isinstance(term, AlgebraicPureTensor):
        coeff, unit_pt = _extract_coeff_and_pt(term)
        if unit_pt is not None:
            # unit_pt is a unit AlgebraicPureTensor (coeff stripped) or bare matrix.
            if isinstance(unit_pt, AlgebraicPureTensor):
                factors = unit_pt.factors
            elif hasattr(unit_pt, "shape"):
                factors = (unit_pt,)
            else:
                return _s(term, **kwargs)

            simp_factors = tuple(_s(f, **kwargs) for f in factors)
            simp_coeff = _s(coeff, **kwargs)
            if simp_coeff is S.Zero:
                shape = tuple(f.shape for f in factors)
                return AlgebraicZeroTensor(shape)
            if simp_coeff is S.One:
                if len(simp_factors) == 1:
                    return simp_factors[0]
                return AlgebraicPureTensor(*simp_factors)
            return _make_scaled_pt(simp_coeff,
                                   AlgebraicPureTensor(*simp_factors) if len(simp_factors) > 1 else simp_factors[0])

    # Fallback: generic simplify
    return _s(term, **kwargs)


def _combine_like_terms(group, **kwargs):
    """Combine a list of terms sharing the same factor sequence.

    Each element of *group* is either an AlgebraicPureTensor or Mul(coeff, AlgebraicPureTensor).
    We extract the (fully combined) coefficient from each, sum all coefficients
    into a single SymPy expression, simplify that sum, and return the product
    of the simplified sum with the unit AlgebraicPureTensor.

    Returns a list (length 0 when the summed coefficient is zero, or length 1
    with the combined term).
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

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
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    left_factors, middle, right_factors = at.as_common_factors()

    if not left_factors and not right_factors:
        # Nothing to extract — return as-is (already combined).
        return at

    if isinstance(middle, AlgebraicTensor):
        middle = _simplify_algebraic_tensor(middle, **kwargs)
    elif isinstance(middle, AlgebraicPureTensor):
        middle = _simplify_algebraic_pure_tensor(middle, **kwargs)
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

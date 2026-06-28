from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.symbol import symbols


def _get_sympy_simplify():
    """Lazy import of SymPy's top-level simplify to avoid circular imports."""
    from sympy.simplify.simplify import simplify as _s
    return _s


# ---------------------------------------------------------------------------
# Proportionality checking helpers
# ---------------------------------------------------------------------------

def _matrix_proportionality_ratio(m1, m2):
    """Check element-wise proportionality for true matrices (shape != (1,1)).

    Returns the commutative ratio k such that m1 = k * m2, or None.
    """
    if m1.shape != m2.shape:
        return None

    ratio = None
    _s = _get_sympy_simplify()
    from sympy.solvers.solvers import solve
    r_symbol = symbols("ratio_symbol")

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

            r = solve(e1 - r_symbol * e2, r_symbol)
            if r is None or len(r) == 0:
                return None
            r = r[0]
            r = _s(r)
            if ratio is None:
                ratio = r
            elif _s(ratio - r) != S.Zero:
                return None

    if ratio is None:
        return None
    return ratio

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

    return _matrix_proportionality_ratio(factor1, factor2)

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
        The tensor sum to simplify.s

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

    **Proportionality factoring** -- merges proportional AlgebraicPureTensor
    terms by combining coefficients or creating linear combinations at
    non-proportional factor slots.

    AlgebraicZeroTensor anchors and S.Zero terms are handled naturally.
    """
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    # Phase 0: proportionality factoring (only simplification tool)
    result = proportionality_factoring(at)

    # If proportionality_factoring returned a non-AlgebraicTensor, simplify it
    if not isinstance(result, AlgebraicTensor):
        if isinstance(result, AlgebraicPureTensor):
            return _simplify_algebraic_pure_tensor(result, **kwargs)
        return result

    return result





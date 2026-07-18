from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.core.symbol import symbols


"""Simplification routines for algebraic tensor expressions.

This module provides :func:`tensorsimplify`, the public entry point for
simplifying algebraic tensor expressions.  The simplification pipeline
consists of two phases:

    1. **Proportionality factoring** -- merges proportional
       :class:`~sympy.tensor.algebraic.AlgebraicPureTensor`
       terms within a sum by combining coefficients or creating linear
       combinations at non-proportional factor slots.

    2. **Commutativity-based simplification** -- decomposes terms by their
       commutativity shape, groups by commutative patterns, applies
       proportionality factoring to non-commutative components, and
       reconstructs the full expression.

Examples
========

Simplify a sum of proportional pure tensors:

>>> from sympy.matrices.expressions import MatrixSymbol
>>> from sympy.tensor.algebraic import AlgebraicPureTensor, tensorsimplify
>>> A = MatrixSymbol("A", 3, 4)
>>> B = MatrixSymbol("B", 4, 5)
>>> C = MatrixSymbol("C", 3, 4)
>>> D = MatrixSymbol("D", 4, 5)
>>> T1 = AlgebraicPureTensor(A, B)
>>> T2 = AlgebraicPureTensor(2, A, B)
>>> print(tensorsimplify(T1 + T2))
3*A ⊗ B
"""


def _get_sympy_simplify():
    """Lazy import of SymPy's top-level simplify to avoid circular imports."""
    from sympy.simplify.simplify import simplify as _s
    return _s


# ---------------------------------------------------------------------------
# Proportionality checking helpers
# ---------------------------------------------------------------------------

def _matrix_proportionality_ratio(m1, m2):
    """Check element-wise proportionality for true matrices (shape != (1,1)).

    Determines the candidate ratio from the first nonzero element pair using
    ``solve`` (handles noncommutative elements where division is undefined),
    then verifies ``simplify(e1 - r * e2) == 0`` against all remaining
    elements.  Calling ``solve`` only once is significantly faster than
    solving per-element.

    Returns the commutative ratio k such that m1 = k * m2, or None.

    Examples
    ========

    >>> from sympy.matrices import ImmutableDenseMatrix
    >>> from sympy.tensor.algebraic.simplify import _matrix_proportionality_ratio
    >>> M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    >>> M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    >>> _matrix_proportionality_ratio(M1, M2)
    1/2
    """
    if m1.shape != m2.shape:
        return None

    _s = _get_sympy_simplify()
    from sympy.solvers.solvers import solve
    r_symbol = symbols("ratio_symbol")

    # Find the first nonzero element pair to determine the candidate ratio.
    ratio = None
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

            # Both nonzero: solve for the candidate ratio once.
            r = solve(e1 - r_symbol * e2, r_symbol)
            if r is None or len(r) == 0:
                return None
            ratio = _s(r[0])
            break
        if ratio is not None:
            break

    if ratio is None:
        # Every element is zero in both matrices — no determinable ratio.
        return None

    # Verify the candidate ratio against every element.
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

            if _s(e1 - ratio * e2) != S.Zero:
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

    Examples
    ========

    >>> from sympy.matrices import ImmutableDenseMatrix
    >>> from sympy.tensor.algebraic.simplify import _proportionality_ratio
    >>> M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    >>> M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    >>> _proportionality_ratio(M1, M2)
    1/2
    >>> _proportionality_ratio(M2, M1)
    2
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
    - AlgebraicPureTensor directly (coeff extracted from term.coeff)
    - Mul with direct matrix factors
    - Anything else -> (term, S.One)

    Examples
    ========

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> from sympy.tensor.algebraic.simplify import _extract_pt_and_coeff
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(2, A, B)
    >>> pt, coeff = _extract_pt_and_coeff(T)
    >>> coeff
    2
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if isinstance(term, AlgebraicPureTensor):
        # Plain AlgebraicPureTensor - extract coefficient
        return (term, term.coeff)

    if isinstance(term, Mul):
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

    Examples
    ========

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic.simplify import _build_pt
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> print(_build_pt(2, [A, B]))
    2*A ⊗ B
    >>> print(_build_pt(1, [A]))
    A
    >>> print(_build_pt(0, [A, B]))
    0_{(3x4), (4x5)}
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

def _try_merge_for_slot(entries, allowed_diff_slot, _s, equality_only=False):
    """Try to merge entries allowing only *allowed_diff_slot* as the
    non-proportional slot (or fully proportional, diff_slot is None).

    Parameters
    ----------
    entries : list of dict
        Mutable list of entry dicts with 'coeff' and 'factors'.
    allowed_diff_slot : int
        The single slot index allowed to be non-proportional.
    _s : callable
        SymPy simplify function.
    equality_only : bool
        If True, only merge when combined_ratio is +1 or -1 and use
        direct add/subtract at the diff slot instead of coefficient-weighted
        linear combination.

    Returns
    -------
    bool
        True if any merge was performed.
    """
    from sympy.core.add import Add as _Add
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    if len(entries) < 2:
        return False

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

                if diff_slot is not None and diff_slot != allowed_diff_slot:
                    continue

                combined_ratio = S.One
                for r in ratios:
                    if r is not None:
                        combined_ratio = combined_ratio * r

                if equality_only:
                    total_ratio = _s(combined_ratio * selected['coeff'] / pivot['coeff'])
                    if total_ratio not in (S.One, S.NegativeOne):
                        continue

                if equality_only and diff_slot is not None:
                    pivot_factor = pivot['factors'][diff_slot]
                    sel_factor = selected['factors'][diff_slot]

                    if total_ratio is S.One:
                        combined_factor = _s(
                            _Add(pivot_factor, sel_factor)
                        )
                    else:
                        combined_factor = _s(
                            _Add(pivot_factor, -sel_factor)
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

                    result = _build_pt(pivot['coeff'], new_factors)

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

                add_coeff = _s(combined_ratio * selected['coeff'])

                if diff_slot is None:
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

    return changed


def _equality_merge_dict(entries, allowed_slot, _s, _factor):
    """Dictionary-based merge for equality factoring at a single allowed slot.

    Groups entries by their canonical key (all factors except the allowed slot,
    with coefficient) and merges entries with matching keys or negated keys.
    Replaces the O(N^2) pivot-based approach with O(N) dictionary lookups.

    Parameters
    ----------
    entries : list of dict
        Entry dicts with 'coeff' and 'factors'.
    allowed_slot : int
        The single slot index allowed to differ between merged entries.
    _s : callable
        SymPy simplify function.
    _factor : callable
        SymPy factor function.

    Returns
    -------
    list of dict
        New entries list with merged entries.
    """
    from sympy.core.mul import Mul
    from sympy.core.add import Add as _Add
    from sympy.matrices.matrixbase import MatrixBase as _MatrixBase

    if len(entries) < 2:
        return entries

    dict_map = {}

    for entry in entries:
        coeff = entry['coeff']
        factors = entry['factors']
        n = len(factors)

        key_factors = [factors[i] for i in range(n) if i != allowed_slot]
        value_factor = factors[allowed_slot]

        # Canonicalize key factors: for each concrete Matrix factor, ensure
        # the first nonzero element does not have -1 as a factor after factoring.
        canonical_coeff = coeff
        canonical_key_factors = list(key_factors)

        for i in range(len(canonical_key_factors)):
            f = canonical_key_factors[i]
            if isinstance(f, _MatrixBase):
                rows, cols = f.shape
                first_nonzero = None
                for r in range(rows):
                    for c in range(cols):
                        elem = f[r, c]
                        if elem != S.Zero and (
                            not hasattr(elem, 'is_zero') or elem.is_zero is not True
                        ):
                            first_nonzero = elem
                            break
                    if first_nonzero is not None:
                        break

                if first_nonzero is not None:
                    ff = _factor(first_nonzero)
                    if isinstance(ff, Mul) and ff.args and ff.args[0] == S.NegativeOne:
                        canonical_key_factors[i] = (-f).applyfunc(_factor)
                        canonical_coeff = S.NegativeOne * canonical_coeff

        # Finalize canonical coefficient
        canonical_coeff = _factor(canonical_coeff)
        canonical_key = _build_pt(canonical_coeff, canonical_key_factors)

        # Dictionary lookup: check canonical_key first, then -1*canonical_key
        found = False
        if canonical_key in dict_map:
            dict_map[canonical_key]['value_sum'] = _s(
                _Add(dict_map[canonical_key]['value_sum'], value_factor)
            )
            found = True

        if not found:
            neg_coeff = _factor(S.NegativeOne * canonical_coeff)
            neg_key = _build_pt(neg_coeff, canonical_key_factors)
            if neg_key in dict_map:
                dict_map[neg_key]['value_sum'] = _s(
                    _Add(dict_map[neg_key]['value_sum'],
                         S.NegativeOne * value_factor)
                )
                found = True

        if not found:
            dict_map[canonical_key] = {
                'coeff': canonical_coeff,
                'factors': canonical_key_factors,
                'value_sum': value_factor,
            }

    # Reconstruct entries from dictionary
    new_entries = []
    for key, data in dict_map.items():
        value_sum = data['value_sum']
        if value_sum == S.Zero or (
            hasattr(value_sum, 'is_zero') and value_sum.is_zero is True
        ):
            continue

        full_factors = list(data['factors'])
        full_factors.insert(allowed_slot, value_sum)

        new_entries.append({
            'coeff': data['coeff'],
            'factors': full_factors,
        })

    return new_entries


def _equality_factoring(at):
    """Simplify an AlgebraicTensor by merging equal or negated
    AlgebraicPureTensor terms (proportionality constant is +1 or -1 only).

    Algorithm
    ---------
    1. Determine the commutativity_pattern of the tensor.
    2. Build an ordered list of allowed diff slots: commutative indices
       (pattern == 1) first, then non-commutative indices (pattern == 0).
    3. For each allowed diff slot in order, run the pivot-based merging:

       a. Pick a pivot entry from the list.
       b. Compare the pivot with every other entry.
       c. Two entries merge only if the combined proportionality ratio
          of all proportional slots is exactly +1 or -1.
       d. For a fully proportional pair (all slots proportional), merge by
          adding coefficients.
       e. For a pair differing in exactly one slot (the allowed diff slot),
          add or subtract the factors directly (no coefficient weighting),
          keeping the pivot coefficient unchanged.
       f. After a merge, reset the pivot to the beginning.

    4. Process commutative slots first, then non-commutative slots.

    Parameters
    ----------
    at : AlgebraicTensor
        The tensor sum to simplify.

    Returns
    -------
    AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or other
        The simplified expression.

    Examples
    ========

    Merge equal terms:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> from sympy.tensor.algebraic.simplify import _equality_factoring
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> C = MatrixSymbol("C", 4, 5)
    >>> T1 = AlgebraicPureTensor(A, B)
    >>> T2 = AlgebraicPureTensor(A, C)
    >>> print(_equality_factoring(T1 + T2))
    A ⊗ (B + C)
    """
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    args = list(at.args)

    tensor_terms = []
    non_tensor = []

    for a in args:
        if isinstance(a, AlgebraicZeroTensor):
            non_tensor.append(a)
        elif isinstance(a, AlgebraicPureTensor):
            tensor_terms.append(a)
        elif isinstance(a, Mul):
            has_matrix = any(hasattr(f, 'shape') for f in a.args)
            if has_matrix:
                tensor_terms.append(a)
            else:
                non_tensor.append(a)
        elif a is not S.Zero:
            non_tensor.append(a)

    if len(tensor_terms) < 2:
        return at

    entries = []
    for t in tensor_terms:
        extracted = _extract_pt_and_coeff(t)
        if len(extracted) == 3:
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

    comm_cs = at.commutativity_pattern
    commutative_indices = [i for i, v in enumerate(comm_cs) if v == 1]
    non_commutative_indices = [i for i, v in enumerate(comm_cs) if v == 0]
    allowed_diff_slots = commutative_indices + non_commutative_indices

    from sympy import factor as _efactor
    for allowed_slot in allowed_diff_slots:
        entries = _equality_merge_dict(entries, allowed_slot, _s, _efactor)

    rebuilt_terms = []
    for entry in entries:
        coeff = entry['coeff']
        factors = entry['factors']
        term = _build_pt(coeff, factors)
        rebuilt_terms.append(term)

    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor as _APT

    for idx, term in enumerate(rebuilt_terms):
        if isinstance(term, _APT):
            new_factors = list(term.factors)
            extra_coeff = S.One
            for fi in range(len(new_factors)):
                ec, ef = _extract_commutative_from_factor(new_factors[fi])
                if ec is not S.One:
                    extra_coeff *= ec
                new_factors[fi] = ef
            if extra_coeff is not S.One:
                new_coeff = _s(extra_coeff * term.coeff)
                rebuilt_terms[idx] = _build_pt(new_coeff, new_factors)

    all_terms = rebuilt_terms + non_tensor

    if not all_terms:
        from sympy.tensor.algebraic.algebraic_zero_tensor import algebraic_zero_tensor
        return algebraic_zero_tensor(at.shape)

    if len(all_terms) == 1:
        return all_terms[0]

    return AlgebraicTensor(*all_terms, _sympify=False)


def _proportionality_factoring(at):
    """Simplify an AlgebraicTensor by merging proportional
    AlgebraicPureTensor terms.

    Algorithm
    ---------
    1. Determine the commutativity_pattern of the tensor.
    2. Build an ordered list of allowed diff slots: commutative indices
       (pattern == 1) first, then non-commutative indices (pattern == 0).
    3. For each allowed diff slot in order, run the pivot-based merging:

       a. Pick a pivot entry from the list.
       b. Compare the pivot with every other entry.
       c. Two entries merge if they are fully proportional (all slots
          proportional) or if their single non-proportimal slot matches
          the current *allowed_diff_slot*.
       d. Merging is either coefficient addition (fully proportional) or
          linear combination at the diff slot.
       e. After a merge, reset the pivot to the beginning.

    4. Process commutative slots first, then non-commutative slots.
       This prioritizes factoring in commutative slots before
       non-commutative ones.

    Parameters
    ----------
    at : AlgebraicTensor
        The tensor sum to simplify.

    Returns
    -------
    AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or other
        The simplified expression.

    Examples
    ========

    Merge proportional terms:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> from sympy.tensor.algebraic.simplify import _proportionality_factoring
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> C = MatrixSymbol("C", 4, 5)
    >>> T1 = AlgebraicPureTensor(A, B)
    >>> T2 = AlgebraicPureTensor(2, A, C)
    >>> print(_proportionality_factoring(T1 + T2))
    A ⊗ (B + 2*C)
    """
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    args = list(at.args)

    tensor_terms = []
    non_tensor = []

    for a in args:
        if isinstance(a, AlgebraicZeroTensor):
            non_tensor.append(a)
        elif isinstance(a, AlgebraicPureTensor):
            tensor_terms.append(a)
        elif isinstance(a, Mul):
            has_matrix = any(hasattr(f, 'shape') for f in a.args)
            if has_matrix:
                tensor_terms.append(a)
            else:
                non_tensor.append(a)
        elif a is not S.Zero:
            non_tensor.append(a)

    if len(tensor_terms) < 2:
        return at

    entries = []
    for t in tensor_terms:
        extracted = _extract_pt_and_coeff(t)
        if len(extracted) == 3:
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

    comm_cs = at.commutativity_pattern
    commutative_indices = [i for i, v in enumerate(comm_cs) if v == 1]
    non_commutative_indices = [i for i, v in enumerate(comm_cs) if v == 0]
    allowed_diff_slots = commutative_indices + non_commutative_indices

    for allowed_slot in allowed_diff_slots:
        _try_merge_for_slot(entries, allowed_slot, _s)

    rebuilt_terms = []
    for entry in entries:
        coeff = entry['coeff']
        factors = entry['factors']
        term = _build_pt(coeff, factors)
        rebuilt_terms.append(term)

    # Post-processing: extract commutative prefactors from every tensor factor
    # of every AlgebraicPureTensor so that commutative scalars appear as
    # the coefficient rather than embedded inside factors.
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor as _APT

    for idx, term in enumerate(rebuilt_terms):
        if isinstance(term, _APT):
            new_factors = list(term.factors)
            extra_coeff = S.One
            for fi in range(len(new_factors)):
                ec, ef = _extract_commutative_from_factor(new_factors[fi])
                if ec is not S.One:
                    extra_coeff *= ec
                new_factors[fi] = ef
            if extra_coeff is not S.One:
                new_coeff = _s(extra_coeff * term.coeff)
                rebuilt_terms[idx] = _build_pt(new_coeff, new_factors)

    all_terms = rebuilt_terms + non_tensor

    if not all_terms:
        from sympy.tensor.algebraic.algebraic_zero_tensor import algebraic_zero_tensor
        return algebraic_zero_tensor(at.shape)

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

    Examples
    ========

    Simplify a sum of proportional pure tensors:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> from sympy.tensor.algebraic.simplify import tensorsimplify
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T1 = AlgebraicPureTensor(A, B)
    >>> T2 = AlgebraicPureTensor(2, A, B)
    >>> print(tensorsimplify(T1 + T2))
    3*A ⊗ B

    Simplify a pure tensor:

    >>> print(tensorsimplify(AlgebraicPureTensor(2, A, B)))
    2*A ⊗ B

    Zero tensor passes through unchanged:

    >>> from sympy.tensor.algebraic import AlgebraicZeroTensor
    >>> Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    >>> print(tensorsimplify(Z))
    0_{(3x4), (4x5)}
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
# Commutativity-based simplification helpers
# ---------------------------------------------------------------------------

def _decompose_commutative_factors(commutative_factors, term_coeff):
    """Decompose commutative matrix factors into sums of basis matrices.

    Iteratively processes each commutative slot, expanding matrix entries
    into basis matrices (0/1 only) while accumulating symbolic coefficients
    into prefactors.

    Parameters
    ----------
    commutative_factors : list of matrix-like objects
        The commutative matrix factors to decompose.
    term_coeff : SymPy expression
        The initial commutative coefficient for the term.

    Returns
    -------
    list of (prefactor, basis_factor_list) tuples
        Each tuple contains a commutative prefactor and a list of basis
        matrices (0/1 entries only), one per commutative slot.

    Examples
    ========

    >>> from sympy.matrices import ImmutableDenseMatrix
    >>> from sympy.tensor.algebraic.simplify import _decompose_commutative_factors
    >>> M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    >>> result = _decompose_commutative_factors([M], 1)
    >>> len(result)
    4
    """
    from sympy.matrices.immutable import ImmutableDenseMatrix

    if not commutative_factors:
        return [(term_coeff, [])]

    terms = [(term_coeff, list(commutative_factors))]

    for s in range(len(commutative_factors)):
        next_terms = []
        for (prefactor, factor_list) in terms:
            M = factor_list[s]
            rows, cols = M.shape[0], M.shape[1]
            for r in range(rows):
                for c in range(cols):
                    entry = M[r, c]
                    if entry == S.Zero or (hasattr(entry, 'is_zero') and entry.is_zero is True):
                        continue
                    data = [[S.Zero for _ in range(cols)] for _ in range(rows)]
                    data[r][c] = S.One
                    E_rc = ImmutableDenseMatrix(data)
                    new_coeff = prefactor * entry
                    new_fl = list(factor_list)
                    new_fl[s] = E_rc
                    next_terms.append((new_coeff, new_fl))
        terms = next_terms
        if not terms:
            break

    return terms


def _reconstruct_term(key, non_commutative_pt, coeff, comm_cs,
                      commutative_indices, non_commutative_indices):
    """Reconstruct a full tensor from commutative and
    non-commutative subtensors.

    Interleaves commutative basis matrices (from *key*) and non-commutative
    factors (from *non_commutative_pt*) back into the original factor order
    specified by *comm_cs*, then builds the result via ``_build_pt``.

    Parameters
    ----------
    key : tuple of matrices or None
        Commutative basis matrices. None when there are no commutative slots.
    non_commutative_pt : AlgebraicPureTensor or None
        Non-commutative subtensor. None when all slots are commutative.
    coeff : SymPy expression
        Combined commutative coefficient.
    comm_cs : tuple of 0/1
        The commutativity_pattern of the original tensor.
    commutative_indices : list of int
        Indices of commutative factor slots.
    non_commutative_indices : list of int
        Indices of non-commutative factor slots.

    Returns
    -------
    AlgebraicPureTensor, bare matrix, or AlgebraicZeroTensor
        The reconstructed tensor term.

    Examples
    ========

    >>> from sympy.matrices import ImmutableDenseMatrix
    >>> from sympy.tensor.algebraic.simplify import _reconstruct_term
    >>> E00 = ImmutableDenseMatrix([[1, 0], [0, 0]])
    >>> print(_reconstruct_term((E00,), None, 1, (1,), [0], []))
    Matrix([[1, 0], [0, 0]])
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    if coeff is S.Zero:
        return None

    comm_iter = iter(key) if key else iter(())
    if non_commutative_pt is not None:
        if isinstance(non_commutative_pt, AlgebraicPureTensor):
            nc_factors = list(non_commutative_pt.factors)
        else:
            nc_factors = [non_commutative_pt]
    else:
        nc_factors = []
    nc_iter = iter(nc_factors)

    merged = []
    for i in range(len(comm_cs)):
        if comm_cs[i] == 1:
            merged.append(next(comm_iter))
        else:
            merged.append(next(nc_iter))

    # Extract commutative prefactors from every merged factor and fold
    # them into the coefficient so that the reconstructed AlgebraicPureTensor
    # carries the fully factored-out coefficient.
    for idx in range(len(merged)):
        ec, ef = _extract_commutative_from_factor(merged[idx])
        if ec is not S.One:
            coeff = coeff * ec
        merged[idx] = ef

    return _build_pt(coeff, merged)


# ---------------------------------------------------------------------------
# Commutative prefactor extraction from non-commutative factors
# ---------------------------------------------------------------------------

def _has_negative_power_or_fraction(expr):
    """Check if *expr* contains Pow with negative exponent or a fraction.

    Returns True if any atom in *expr* is a Pow with a negative exponent,
    or if *expr* as a rational expression has a nontrivial denominator.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.tensor.algebraic.simplify import _has_negative_power_or_fraction
    >>> _has_negative_power_or_fraction(x**2 * y)
    False
    >>> _has_negative_power_or_fraction(x**(-1))
    True
    >>> _has_negative_power_or_fraction(1 / x)
    True
    >>> _has_negative_power_or_fraction(x / y)
    True
    >>> _has_negative_power_or_fraction(2 * x)
    False
    """
    from sympy.core.power import Pow

    for atom in expr.atoms(Pow):
        if hasattr(atom.exp, 'is_negative') and atom.exp.is_negative:
            return True
    num, den = expr.as_numer_denom()
    if den != S.One:
        return True
    return False


def _normalize_factor_sign(f):
    """Normalize the sign of a polynomial factor so -(w-z) and (w-z) match.

    Makes polynomial factors monic in sign by flipping the overall sign when
    the coefficient of the first symbol (sorted by name for determinism) is
    negative.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.tensor.algebraic.simplify import _normalize_factor_sign
    >>> _normalize_factor_sign(x - y)
    x - y
    >>> _normalize_factor_sign(-x + y)
    x - y
    """
    from sympy.core.add import Add as _Add
    from sympy.core.mul import Mul as _Mul

    if isinstance(f, _Mul):
        if f.args[0].is_number and f.args[0].is_negative:
            nf = -f
            if nf != f:
                return nf
    elif isinstance(f, _Add):
        # Sort symbols by string representation for deterministic ordering.
        # Using f.atoms() directly is non-deterministic because atoms()
        # returns a set whose iteration order depends on the Python hash seed.
        syms = sorted(
            (a for a in f.atoms()
             if getattr(a, 'is_Symbol', False)
             and not getattr(a, 'is_number', True)),
            key=lambda s: str(s),
        )
        if syms:
            lc = f.coeff(syms[0])
            if hasattr(lc, 'is_negative') and lc.is_negative:
                nf = -f
                if nf != f:
                    return nf
    return f


def _is_exactly_divisible(entry, candidate):
    """Check if *entry* is exactly divisible by *candidate* with no remainder.

    A zero entry is considered trivially divisible.  For nonzero entries we
    verify that ``candidate`` does not appear in the denominator of the
    simplified quotient ``entry / candidate``.

    Examples
    ========

    >>> from sympy.tensor.algebraic.simplify import _is_exactly_divisible
    >>> _is_exactly_divisible(6, 2)
    True
    >>> _is_exactly_divisible(5, 2)
    False
    >>> _is_exactly_divisible(0, 3)
    True
    """
    from sympy import cancel as _cancel
    from sympy.core.sympify import sympify

    if entry == S.Zero or (hasattr(entry, 'is_zero') and entry.is_zero is True):
        return True
    # Use sympify to avoid Python float division for numeric operands
    q = _cancel(sympify(entry) / sympify(candidate))
    _slocal = _get_sympy_simplify()
    q = _slocal(q)
    if q.has(sympify(1) / sympify(candidate)):
        return False
    # For numeric entries, check if the result is an integer
    if q.is_number:
        denom = q.as_numer_denom()[1]
        return denom == S.One
    # For symbolic entries, check if candidate appears in denominator
    denom = q.as_numer_denom()[1]
    if denom != S.One and denom.has(candidate):
        return False
    return True


def _deduplicate_proportional(factors):
    """Remove factors that are constant multiples of an already-kept factor.

    Given a list of irreducible commutative factors, returns a filtered list
    that contains at most one representative per proportionality class
    (i.e., factors that differ only by a constant scalar such as -1).

    This is needed because _normalize_factor_sign may not catch every sign
    convention, and both ``a`` and ``-a`` could otherwise end up as
    separate survivors, producing ``-a**2`` as the divisor.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.tensor.algebraic.simplify import _deduplicate_proportional
    >>> _deduplicate_proportional([x, -x, x + 1])
    [x, x + 1]
    """
    if not factors:
        return []
    kept = []
    for f in factors:
        proportional = False
        for g in kept:
            # Check if f / g is a nonzero constant
            if g == S.Zero:
                continue
            q = _get_sympy_simplify()(f / g)
            if (hasattr(q, 'is_number') and q.is_number and q != S.Zero):
                proportional = True
                break
            # Also check g / f (in case f == 0)
            if f != S.Zero:
                q2 = _get_sympy_simplify()(g / f)
                if (hasattr(q2, 'is_number') and q2.is_number and q2 != S.Zero):
                    proportional = True
                    break
        if not proportional:
            kept.append(f)
    return kept


def _extract_commutative_from_factor(factor):
    """Extract commutative prefactor from a single tensor factor.

    Handles three cases:
    1. **Add** -- factors the expression with ``factor()`` to pull out a
       common commutative divisor (e.g., ``(x+2)*P + (x+2)*Q -> (x+2)*(P+Q)``),
       then delegates to the Mul handler.
    2. **Mul** -- separates commutative and non-commutative args.
    3. **Matrix** -- factors each entry, finds common commutative divisors
       across all nonzero entries, divides the matrix, and returns the
       divisor as prefactor.

    Returns
    -------
    (commutative_coeff, new_factor)
        *commutative_coeff* is S.One when nothing was extracted.

    Examples
    ========

    >>> from sympy.matrices import ImmutableDenseMatrix
    >>> from sympy.abc import x
    >>> from sympy.tensor.algebraic.simplify import _extract_commutative_from_factor
    >>> M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    >>> coeff, new_M = _extract_commutative_from_factor(M)
    >>> coeff
    x
    """
    from sympy.core.add import Add as _Add
    from sympy.core.mul import Mul
    from sympy.matrices.immutable import ImmutableDenseMatrix
    from sympy import cancel
    from sympy import factor as _factor

    # --- Add: factor to pull out common commutative divisor ---
    if isinstance(factor, _Add):
        # Only attempt factor() if the Add contains non-commutative parts
        # that could share a commutative prefactor.  factor() crashes on
        # expressions mixing MatrixExpr and scalars, so we guard with try/except.
        try:
            factored = _factor(factor)
        except (TypeError, NotImplementedError):
            factored = factor
        if isinstance(factored, Mul) and factored != factor:
            comm_parts = []
            noncomm_parts = []
            for a in factored.args:
                if hasattr(a, 'is_commutative') and a.is_commutative:
                    comm_parts.append(a)
                else:
                    noncomm_parts.append(a)

            if comm_parts and noncomm_parts:
                comm_coeff = Mul(*comm_parts, evaluate=True)
                if _has_negative_power_or_fraction(comm_coeff):
                    pass
                else:
                    if len(noncomm_parts) == 1:
                        new_factor = noncomm_parts[0]
                    else:
                        new_factor = Mul(*noncomm_parts, evaluate=False)
                    return (comm_coeff, new_factor)

    # --- MatAdd / MatrixExpr sum: factor to pull out common commutative divisor ---
    # MatAdd (and other additive MatrixExpr) are NOT instances of sympy Add,
    # so we handle them separately.  SymPy's factor() crashes on MatAdd with
    # non-commutative matrix symbols, so we manually walk the args and find
    # the common commutative divisor.
    from sympy.matrices.expressions.matexpr import MatrixExpr as _MatrixExpr
    from sympy.matrices.expressions.matmul import MatMul as _MatMul
    from sympy.matrices.expressions.matadd import MatAdd as _MatAdd

    if isinstance(factor, _MatrixExpr) and hasattr(factor, 'is_Add') and factor.is_Add:
        # Try SymPy's factor() first (works for some cases)
        try:
            factored = _factor(factor)
            if isinstance(factored, _MatMul) and factored != factor:
                c, rest = factored.as_coeff_Mul()
                if c != S.One and rest != factor:
                    if _has_negative_power_or_fraction(c):
                        pass
                    else:
                        return (c, rest)
        except (TypeError, NotImplementedError):
            pass

        # Manual extraction: walk each term, extract commutative coeff from
        # MatMul.args (not as_coeff_Mul which only extracts numeric coeffs),
        # find the common divisor across all terms.
        terms = list(factor.args)
        coeffs = []
        rests = []
        all_ok = True
        for term in terms:
            if isinstance(term, _MatMul):
                # Separate commutative and non-commutative args from MatMul
                comm_args = []
                noncomm_args = []
                for a in term.args:
                    if hasattr(a, 'is_commutative') and a.is_commutative:
                        comm_args.append(a)
                    else:
                        noncomm_args.append(a)
                if comm_args:
                    cf = Mul(*comm_args, evaluate=True)
                    if len(noncomm_args) == 1:
                        rf = noncomm_args[0]
                    elif len(noncomm_args) > 1:
                        rf = _MatMul(*noncomm_args, evaluate=False)
                    else:
                        rf = S.One
                    coeffs.append(cf)
                    rests.append(rf)
                else:
                    all_ok = False
                    break
            else:
                # Bare matrix symbol or other non-MatMul term
                coeffs.append(S.One)
                rests.append(term)
        if all_ok and len(coeffs) > 1:
            # Find common divisor of all coefficients
            common = coeffs[0]
            for cf in coeffs[1:]:
                if _is_exactly_divisible(cf, common):
                    pass  # common still divides all
                elif _is_exactly_divisible(common, cf):
                    common = cf
                else:
                    # Try to find partial GCD by factoring both
                    from sympy import factor as _fc
                    fa = _fc(common).as_ordered_factors() if isinstance(_fc(common), Mul) else [_fc(common)]
                    fb = _fc(cf).as_ordered_factors() if isinstance(_fc(cf), Mul) else [_fc(cf)]
                    new_common = S.One
                    for a in fa:
                        for b in fb:
                            if cancel(a / b) == S.One:
                                new_common *= a
                                break
                    if new_common != S.One:
                        common = new_common
                    else:
                        common = S.One
                        break
            if common != S.One:
                # Reconstruct: sum of (coeff_i/common)*rest_i
                new_terms = []
                for cf, rf in zip(coeffs, rests):
                    new_cf = cancel(cf / common)
                    if new_cf == S.One:
                        new_terms.append(rf)
                    else:
                        new_terms.append(_MatMul(new_cf, rf))
                if len(new_terms) == 1:
                    new_factor = new_terms[0]
                else:
                    new_factor = _MatAdd(*new_terms)
                if _has_negative_power_or_fraction(common):
                    pass
                else:
                    return (common, new_factor)

    # --- Mul: separate commutative / non-commutative args ---
    if isinstance(factor, Mul):
        comm_parts = []
        noncomm_parts = []
        for a in factor.args:
            if hasattr(a, 'is_commutative') and a.is_commutative:
                comm_parts.append(a)
            else:
                noncomm_parts.append(a)

        if comm_parts and noncomm_parts:
            comm_coeff = Mul(*comm_parts, evaluate=True)
            if _has_negative_power_or_fraction(comm_coeff):
                pass
            else:
                if len(noncomm_parts) == 1:
                    new_factor = noncomm_parts[0]
                else:
                    new_factor = Mul(*noncomm_parts, evaluate=False)
                return (comm_coeff, new_factor)

    # --- Matrix: common commutative divisor across nonzero entries ---
    # Only process concrete MatrixBase instances (ImmutableDenseMatrix, etc.),
    # NOT symbolic MatrixExpr (MatrixSymbol, MatAdd, MatMul).  MatrixBase
    # already excludes symbolic expressions, so no extra check is needed.
    from sympy.matrices.matrixbase import MatrixBase as _MatrixBase

    if isinstance(factor, _MatrixBase):
        rows, cols = factor.shape[0], factor.shape[1]

        nonzero_entries = []
        nonzero_positions = []
        for r in range(rows):
            for c in range(cols):
                entry = factor[r, c]
                if entry != S.Zero and (
                    not hasattr(entry, 'is_zero') or entry.is_zero is not True
                ):
                    nonzero_entries.append(entry)
                    nonzero_positions.append((r, c))

        if not nonzero_entries:
            return (S.One, factor)

        # Collect candidate commutative factors from every nonzero entry.
        # Use factor() to decompose into irreducible factors, then normalize
        # the sign so that -(w-z) and (w-z) are recognized as the same factor.
        from sympy import factor as _factor
        candidates = set()
        for entry in nonzero_entries:
            factored = _factor(entry)
            if isinstance(factored, Mul):
                fac_args = factored.args
            else:
                fac_args = (factored,)
            for a in fac_args:
                if hasattr(a, 'is_number') and a.is_number:
                    continue
                if (hasattr(a, 'is_commutative') and a.is_commutative
                        and a != S.One):
                    candidates.add(_normalize_factor_sign(a))

        # Keep only candidates that divide ALL nonzero entries
        surviving = []
        for cand in candidates:
            ok = True
            for entry in nonzero_entries:
                if not _is_exactly_divisible(entry, cand):
                    ok = False
                    break
            if ok:
                surviving.append(cand)

        # For purely numeric entries, compute the GCD as a candidate
        from sympy.core.numbers import Integer
        from sympy.core.intfunc import igcd
        all_numeric = all(
            isinstance(e, (int, Integer)) or
            (hasattr(e, 'is_integer') and e.is_integer is True)
            for e in nonzero_entries
        )
        if all_numeric and len(nonzero_entries) > 0:
            gcd_val = nonzero_entries[0]
            for e in nonzero_entries[1:]:
                gcd_val = igcd(gcd_val, e)
            if gcd_val != 1 and gcd_val != -1:
                surviving.append(abs(int(gcd_val)))

        if not surviving:
            return (S.One, factor)

        # Remove proportional (constant-scaled) duplicates.  Because
        # _normalize_factor_sign relies on a sign convention that may not
        # catch every case, both ``a`` and ``-a`` could end up in
        # ``surviving``.  Their product ``-a**2`` does not divide the
        # entries, so we must keep at most one representative per
        # proportionality class.
        surviving = _deduplicate_proportional(surviving)

        if not surviving:
            return (S.One, factor)

        # Divide every nonzero entry by the product of survivors
        divisor = Mul(*surviving, evaluate=True)

        # Verify that dividing by the divisor does not introduce any
        # denominators in the matrix entries.  If even one entry would
        # contain a denominator involving the divisor (or any of its
        # factors), bail out and return the original factor unchanged.
        from sympy.core.mul import Mul as _MulCheck
        divisor_factors = set()
        if isinstance(divisor, _MulCheck):
            for df in divisor.args:
                if not (hasattr(df, 'is_number') and df.is_number):
                    divisor_factors.add(df)
        else:
            if not (hasattr(divisor, 'is_number') and divisor.is_number):
                divisor_factors.add(divisor)

        if divisor_factors:
            for r in range(rows):
                for c in range(cols):
                    entry = factor[r, c]
                    if entry == S.Zero or (
                        hasattr(entry, 'is_zero') and entry.is_zero is True
                    ):
                        continue
                    quotient = cancel(entry / divisor)
                    q_denom = quotient.as_numer_denom()[1]
                    if q_denom != S.One:
                        for df in divisor_factors:
                            if q_denom.has(df):
                                return (S.One, factor)

        new_data = []
        for r in range(rows):
            row = []
            for c in range(cols):
                entry = factor[r, c]
                if entry == S.Zero or (
                    hasattr(entry, 'is_zero') and entry.is_zero is True
                ):
                    row.append(S.Zero)
                else:
                    row.append(cancel(entry / divisor))
            new_data.append(row)

        new_factor = ImmutableDenseMatrix(new_data)
        return (divisor, new_factor)

    # Pass through
    return (S.One, factor)


def _extract_commutative_prefactors(pt):
    """Extract commutative prefactors from every factor of a PureTensor.

    Walks each tensor factor, pulls out commutative multiplicative parts,
    and accumulates them into a single prefactor.

    Returns
    -------
    (extracted_coeff, new_factors_list)

    Examples
    ========

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> from sympy.tensor.algebraic.simplify import _extract_commutative_prefactors
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> coeff, factors = _extract_commutative_prefactors(T)
    >>> coeff
    1
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    extracted = S.One

    if isinstance(pt, AlgebraicPureTensor):
        factors = list(pt.factors)
    else:
        factors = [pt]

    new_factors = []

    for factor in factors:
        comm_coeff, new_factor = _extract_commutative_from_factor(factor)
        extracted *= comm_coeff
        new_factors.append(new_factor)

    return (extracted, new_factors)


# ---------------------------------------------------------------------------
# Commutativity-based simplification
# ---------------------------------------------------------------------------

def _commutativity_simplify(at, **kwargs):
    """Simplify an AlgebraicTensor using commutativity-aware decomposition.

    Decomposes each term into commutative and non-commutative subtensors,
    groups by commutative pattern, applies ``_proportionality_factoring`` to
    non-commutative components, and reconstructs the full expression.

    Parameters
    ----------
    at : AlgebraicTensor
        The tensor to simplify.
    **kwargs
        Extra keyword arguments forwarded to SymPy's simplify.

    Returns
    -------
    Simplified tensor expression.
    """
    from sympy.matrices.expressions.special import ZeroMatrix as _ZeroMatrix
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    comm_cs = at.commutativity_pattern
    commutative_indices = [i for i, v in enumerate(comm_cs) if v == 1]
    non_commutative_indices = [i for i, v in enumerate(comm_cs) if v == 0]

    comm_dict = {}

    for term in at.args:
        if isinstance(term, AlgebraicZeroTensor):
            continue

        extracted = _extract_pt_and_coeff(term)

        if len(extracted) == 3:
            _, term_coeff, factors = extracted
        else:
            pt, term_coeff = extracted
            if isinstance(pt, AlgebraicPureTensor):
                factors = list(pt.factors)
            else:
                factors = [pt]
                term_coeff = S.One

        if term_coeff is S.Zero:
            continue

        nc_factors = [factors[i] for i in non_commutative_indices]
        comm_factors = [factors[i] for i in commutative_indices]

        # Extract commutative coefficients from non-commutative factors so
        # they appear as the coefficient of the non-commutative subtensor
        # rather than embedded inside the factors.
        nc_coeff = S.One
        new_nc_factors = []
        for nf in nc_factors:
            ec, ef = _extract_commutative_from_factor(nf)
            nc_coeff *= ec
            new_nc_factors.append(ef)

        if non_commutative_indices:
            non_commutative_pt = _build_pt(nc_coeff, new_nc_factors)
        else:
            non_commutative_pt = None

        decomposed = _decompose_commutative_factors(comm_factors, term_coeff)

        for (prefactor, basis_factor_list) in decomposed:
            if prefactor is S.Zero:
                continue

            if commutative_indices:
                commutative_key = tuple(basis_factor_list)
            else:
                commutative_key = None

            if commutative_key not in comm_dict:
                if non_commutative_indices:
                    shape = tuple(f.shape for f in new_nc_factors)
                    comm_dict[commutative_key] = AlgebraicZeroTensor(shape)
                else:
                    comm_dict[commutative_key] = S.Zero

            if non_commutative_indices:
                scaled = prefactor * non_commutative_pt
                value = comm_dict[commutative_key]
                comm_dict[commutative_key] = value + scaled
            else:
                comm_dict[commutative_key] = comm_dict[commutative_key] + prefactor

    # Normalize ZeroMatrix results (from matrix expression cancellation like -B+B)
    # to AlgebraicZeroTensor so the reconstruction phase handles them correctly.
    nc_shape = tuple(at.shape[i] for i in non_commutative_indices) if non_commutative_indices else ()
    for key in comm_dict:
        value = comm_dict[key]
        if isinstance(value, _ZeroMatrix):
            if non_commutative_indices:
                comm_dict[key] = AlgebraicZeroTensor(nc_shape)
            else:
                comm_dict[key] = S.Zero

    for key in comm_dict:
        value = comm_dict[key]
        if isinstance(value, AlgebraicTensor):
            # Cancel terms like X + (-X) that arise when single-factor
            # non-commutative tensors unwrap to bare matrices and
            # AlgebraicTensor.__new__ doesn't simplify B + (-B).
            cancel_map = {}
            for arg in value.args:
                if isinstance(arg, AlgebraicZeroTensor):
                    continue
                if isinstance(arg, AlgebraicPureTensor):
                    nc_key = arg.factors
                    c = arg.coeff
                else:
                    # Bare non-commutative matrix (unwrapped single-factor tensor)
                    if hasattr(arg, 'is_commutative') and arg.is_commutative:
                        nc_key = ()
                        c = arg
                    else:
                        nc_key = (arg,)
                        c = S.One
                cancel_map[nc_key] = cancel_map.get(nc_key, S.Zero) + c
            surviving = []
            for nc_key, combined_c in cancel_map.items():
                combined_c = _get_sympy_simplify()(combined_c)
                if combined_c == S.Zero:
                    continue
                if nc_key:
                    if len(nc_key) == 1 and combined_c is S.One:
                        surviving.append(nc_key[0])
                    elif len(nc_key) == 1:
                        surviving.append(AlgebraicPureTensor(combined_c, nc_key[0]))
                    elif combined_c is S.One:
                        surviving.append(AlgebraicPureTensor(*nc_key))
                    else:
                        surviving.append(AlgebraicPureTensor(combined_c, *nc_key))
                else:
                    surviving.append(combined_c)
            if not surviving:
                comm_dict[key] = AlgebraicZeroTensor(value.shape)
            elif len(surviving) == 1:
                comm_dict[key] = surviving[0]
            else:
                comm_dict[key] = AlgebraicTensor(*surviving, _sympify=False)
        if isinstance(comm_dict[key], AlgebraicTensor):
            comm_dict[key] = _equality_factoring(comm_dict[key])
    all_reconstructed = []

    for key, value in comm_dict.items():
        if non_commutative_indices:
            if isinstance(value, AlgebraicTensor):
                bodies = list(value.args)
            elif isinstance(value, (AlgebraicPureTensor, AlgebraicZeroTensor)):
                bodies = [value]
            else:
                bodies = [value]

            for body in bodies:
                if isinstance(body, AlgebraicZeroTensor):
                    all_reconstructed.append(AlgebraicZeroTensor(at.shape))
                    continue
                if isinstance(body, _ZeroMatrix):
                    all_reconstructed.append(AlgebraicZeroTensor(at.shape))
                    continue

                if isinstance(body, AlgebraicPureTensor):
                    coeff = body.coeff
                    pt = body
                else:
                    coeff, pt = S.One, body

                reconstructed = _reconstruct_term(key, pt, coeff,
                    comm_cs, commutative_indices, non_commutative_indices)
                if reconstructed is not None:
                    all_reconstructed.append(reconstructed)
        else:
            if value is S.Zero:
                all_reconstructed.append(AlgebraicZeroTensor(at.shape))
                continue
            reconstructed = _reconstruct_term(key, None, value,
                comm_cs, commutative_indices, non_commutative_indices)
            if reconstructed is not None:
                all_reconstructed.append(reconstructed)
    real_terms = [t for t in all_reconstructed if t is not S.Zero]

    if not real_terms:
        final_result = AlgebraicZeroTensor(at.shape)
    elif len(real_terms) == 1:
        final_result = real_terms[0]
    else:
        final_result = AlgebraicTensor(*real_terms, _sympify=False)

    if isinstance(final_result, AlgebraicTensor):
        final_result = _equality_factoring(final_result)

    if isinstance(final_result, AlgebraicTensor):
        new_args = []
        for arg in final_result.args:
            new_args.append(tensorsimplify(arg, **kwargs))
        final_result = AlgebraicTensor(*new_args, _sympify=False)
    elif isinstance(final_result, AlgebraicPureTensor):
        final_result = _simplify_algebraic_pure_tensor(final_result, **kwargs)

    return final_result


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
    3. Extract commutative prefactors from each factor (common divisors
       across matrix entries) and fold them into the coefficient.
    4. Reconstruct the PureTensor with the (possibly changed) pieces.
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    _s = _get_sympy_simplify()

    coeff = pt.coeff
    factors = pt.factors

    # Simplify coefficient
    new_coeff = _s(coeff, **kwargs)
    if new_coeff is S.Zero:
        return AlgebraicZeroTensor(pt.shape)

      # Simplify individual factors and extract commutative prefactors
    from sympy import factor as _factor
    from sympy.matrices.matrixbase import MatrixBase as _MatrixBase

    new_factors = []
    changed = new_coeff != coeff
    for f in factors:
        sf = _s(f, **kwargs)
        if sf != f:
            changed = True
        # Factorize matrix entries so commutative prefactor extraction works
        # better. Only apply to concrete matrices (MatrixBase), not symbolic
        # MatrixExpr (MatrixSymbol, MatAdd, MatMul) which would produce
        # ElementwiseApplyFunction(Lambda(...), ...) wrappers.
        if isinstance(sf, _MatrixBase):
            fsf = sf.applyfunc(_factor)
            if fsf != sf:
                changed = True
            sf = fsf
        # Extract commutative prefactors from this factor
        fc, nf = _extract_commutative_from_factor(sf)
        if fc is not S.One:
            new_coeff = _s(new_coeff * fc)
            changed = True
        if nf != sf:
            sf = nf
        new_factors.append(sf)

    if not changed:
        return pt

    if new_coeff is S.Zero:
        return AlgebraicZeroTensor(pt.shape)

    if len(new_factors) == 0:
        return new_coeff

    if new_coeff is S.One:
        if len(new_factors) == 1:
            return new_factors[0]
        return AlgebraicPureTensor(*new_factors)
    return AlgebraicPureTensor(new_coeff, *new_factors)


# ---------------------------------------------------------------------------
# AlgebraicTensor simplification
# ---------------------------------------------------------------------------

def _simplify_algebraic_tensor(at, **kwargs):
    """Simplify an AlgebraicTensor (sum of same-shape tensor terms).

    The simplification pipeline:

    **Commutativity-based simplification** -- decomposes terms by their
    commutativity shape, groups by commutative patterns, applies
    ``_proportionality_factoring`` to non-commutative components, and
    reconstructs the full expression with simplified factors.

    **Proportionality factoring** -- merges proportional AlgebraicPureTensor
    terms by combining coefficients or creating linear combinations at
    non-proportional factor slots.

    AlgebraicZeroTensor anchors and S.Zero terms are handled naturally.
    """
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor

    # Phase 0: commutativity-based simplification
    result = _commutativity_simplify(at, **kwargs)

    # If commutativity simplification collapsed to a single term, simplify it
    if not isinstance(result, AlgebraicTensor):
        if isinstance(result, AlgebraicPureTensor):
            return _simplify_algebraic_pure_tensor(result, **kwargs)
        return result

    # Phase 1: proportionality factoring on the commutativity-simplified result
    result = _equality_factoring(result)

    if not isinstance(result, AlgebraicTensor):
        if isinstance(result, AlgebraicPureTensor):
            return _simplify_algebraic_pure_tensor(result, **kwargs)
        return result

    return result

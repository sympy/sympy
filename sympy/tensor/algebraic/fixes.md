# Fixes for `sympy.tensor.algebraic`

## CRITICAL BUGS

- [x] **1. `AlgebraicZeroTensor.__add__` / `__radd__` — no shape check**

  **File:** `algebraic_zero_tensor.py:180-238`

  Both methods unconditionally return `other`. A zero tensor of shape `((3,4),)` can be "added" to any object of any shape, silently accepted. Compare with `AlgebraicPureTensor.__add__` (`algebraic_pure_tensor.py:502-506`), which checks `other.shape == self.shape`.

  ```python
  def __add__(self, other):
      return other  # No shape validation!
  ```

  This means `AlgebraicZeroTensor((3,4)) + AlgebraicZeroTensor((5,6))` returns the second zero tensor instead of raising an error, violating the invariant that tensors of different shapes are not summable (stated explicitly in the docstring at line 53).

  **Fix:** Add a shape validation guard. If `other` has a `.shape` attribute and it differs from `self.shape`, raise `ShapeMismatchError`. For `AlgebraicZeroTensor` operands specifically, compare shapes. For other tensor types, delegate to the existing shape-aware logic (e.g., check via `_shape_of(other)`).

  **[FIXED]**

---

- [x] **2. `compose_algebraic_tensors` — wrong shape for zero-tensor shortcuts**

  **File:** `algebraic_tensor.py:221-268`

  ```python
  if isinstance(left, _ZT):
      return _ZT(left.shape)
  if isinstance(right, _ZT):
      return _ZT(right.shape)
  ```

  Composition (`*`) is factor-wise matrix multiplication. The result shape should be `(left_rows, right_cols)` per factor slot, not the input shape. The `compose_algebraic_pure_tensors` function handles this correctly (lines 1006-1058), computing the composed shape. But `compose_algebraic_tensors` preserves the wrong shape, so `Z_((3,4)) * M_((4,5))` returns `Z_((3,4))` instead of `Z_((3,5))`.

  **Fix:** Compute the composed shape from both operands, following the same logic used in `compose_algebraic_pure_tensors` (lines 1006-1058): `composed_shape = tuple((l[0], r[1]) for l, r in zip(left.shape, right.shape))`.

  **[FIXED]**

---

- [x] **3. `_extract_pt_and_coeff` — inconsistent return type**

  **File:** `simplify.py:159-200`

  Returns `(term, term.coeff)` (2-tuple) for `AlgebraicPureTensor`, `(term, S.One)` (2-tuple) for pass-through, but `(None, coeff, matrix_factors)` (3-tuple) for `Mul` with matrix factors. Every caller must do:

  ```python
  extracted = _extract_pt_and_coeff(t)
  if len(extracted) == 3:
      _, coeff, factors = extracted
  else:
      pt, coeff = extracted
      factors = list(pt.factors)  # or [pt]
  ```

  This pattern appears in `_equality_factoring` (lines 365-491), `_commutativity_simplify` (lines 1292-1410), and elsewhere.

  **Fix:** Return a consistent 2-tuple `(coeff, factors)` in all cases. For `AlgebraicPureTensor`, return `(term.coeff, list(term.factors))`. For `Mul`, return `(coeff, matrix_factors)`. For pass-through, return `(S.One, [term])`.

  **[FIXED]**

---

## SIGNIFICANT ISSUES

- [x] **4. `_proportionality_factoring` is dead code**

  **File:** `simplify.py:690-837` (removed)

  The function is documented but never called. `_simplify_algebraic_tensor` (line 1827) describes a pipeline of "Commutativity-based simplification" and "Proportionality factoring," but Phase 1 calls only `_equality_factoring` (line 1856). The proportionality factoring step is entirely missing from the actual pipeline.

  **Fix:** Remove the function and update the docstring.

  **[FIXED]** Removed `_proportionality_factoring`, `_try_merge_for_slot`, related tests, and updated all docstring references.

---

- [x] **5. `AlgebraicTensor._merge` / `_subtract` — near-duplicate methods**

  **File:** `algebraic_tensor.py:703-824`

  These two methods are ~120 lines each and differ only in the sign of `c` being merged. A single helper would eliminate ~110 lines of duplicated logic.

  **Fix:** Extract a shared `_combine_coeff_maps(left, right, negate_right=False)` method. `_merge` calls it with `negate_right=False`, `_subtract` with `negate_right=True`.

  **[FIXED]** Extracted `_combine_coeff_maps` classmethod; `_merge` and `_subtract` now delegate to it.

---

- [x] **6. `_equality_factoring` / `_proportionality_factoring` — duplicated structure**

  **File:** `simplify.py:561-709` and `simplify.py:690-837`

  Both functions share identical entry extraction, commutativity pattern analysis, and term rebuilding logic. Only the inner merge call differs. This is ~150 lines of duplicated boilerplate.

  **Fix:** Extract the common scaffolding (term separation, entry building, slot ordering, rebuilding) into a single function that accepts a `merge_fn` callback for the slot-specific logic.

  **[RESOLVED BY #4]** With `_proportionality_factoring` removed, only `_equality_factoring` remains. No duplication to refactor.

---

- [x] **7. `_extract_commutative_from_factor` — 280+ lines, deep nesting**

  **File:** `simplify.py:861-1176` (refactored)

  This single function handles Add, MatAdd, Mul, and Matrix cases with 5-6 levels of nesting.

  **Fix:** Decompose into `_extract_from_add`, `_extract_from_matadd`, `_extract_from_mul`, `_extract_from_matrix`, each called from a dispatcher at the top of `_extract_commutative_from_factor`. Write docstrings for each new helper function.

  **[FIXED]** Decomposed into dispatcher + 4 specialized helpers (`_extract_from_add`, `_extract_from_matadd`, `_extract_from_mul`, `_extract_from_matrix`) plus supporting helpers (`_find_common_divisor`, `_collect_matrix_factor_candidates`, `_compute_numeric_gcd`, `_verify_safe_division`).

---

- [] **8. `AlgebraicTensor._compose_with_term` — unsafe `.factors` access**

  **File:** `algebraic_tensor.py:974`

  ```python
  results.append(AlgebraicPureTensor(coeff, *comp.factors if isinstance(comp, _PT) else (comp,)))
  ```

  When `comp` is a bare matrix (single-factor result unwrapped by `compose_algebraic_pure_tensors`), the ternary handles it but the expression is parsed ambiguously and is hard to read.

  **Fix:** Split into an explicit `if isinstance(comp, _PT)` / `else` block with clear variable names.

---

## MODERATE ISSUES

- [] **9. `AlgebraicZeroTensor.is_zero = None`**

  **File:** `algebraic_zero_tensor.py:106`

  A zero tensor *is* zero. Setting `is_zero = None` tells SymPy's assumption system the answer is unknown. This affects `simplify`, `expand`, and boolean context checks downstream.

  **Fix:** Change to `is_zero = True`.

---

- [] **10. `AlgebraicPureTensor` extends `Mul` — fragile inheritance**

  **File:** `algebraic_pure_tensor.py:76`

  `Mul` has deeply internal logic (`flatten`, `_from_args`, `dual`, etc.). Setting `is_Mul = False` suppresses `Mul.flatten` unpacking, but other `Mul` internals may assume `is_Mul is True`. This works now but is a maintenance risk.

  **Fix:** Evaluate whether `AlgebraicPureTensor` should extend `Basic` or `Expr` directly (as `AlgebraicTensor` does), storing `(_coeff, factors...)` as `args`. This is a larger refactor; mark as "investigate" rather than "fix now."

---

- [] **11. Duplicate import in module docstring**

  **File:** `algebraic_tensor.py:30`

  ```python
  >>> from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor, AlgebraicTensor
  ```

  `AlgebraicTensor` imported twice.

  **Fix:** Remove the duplicate.

---

- [] **12. `AlgebraicTensor.__new__` unused variable `comm_cs`**

  **File:** `algebraic_tensor.py:408`

  `comm_cs` is returned by `_flatten_args` but never used in `__new__`.

  **Fix:** Assign to `_` or remove from the return tuple if it's not needed elsewhere in `__new__`.

---

- [] **13. `AlgebraicTensor.__new__` redundant S.Zero filtering**

  **File:** `algebraic_tensor.py:415-426`

  `S.Zero` is filtered from `flat` twice with identical logic and the same fallback to `zero_term`. The second block can never execute differently from the first.

  **Fix:** Remove the duplicate block.

---

- [] **14. `_commutativity_simplify` — `S.One` prefactor edge case**

  **File:** `simplify.py:1367`

  ```python
  scaled = prefactor * non_commutative_pt
  value = comm_dict[commutative_key]
  comm_dict[commutative_key] = value + scaled
  ```

  If `value` starts as an `AlgebraicZeroTensor` and `scaled` has a different shape, the shape mismatch is silently ignored (see bug #1).

  **Fix:** Once bug #1 is fixed, this will be caught. As a belt-and-suspenders measure, add an explicit shape assertion here.

---

- [] **15. `AlgebraicZeroTensor.__sub__` / `__rsub__` — unnecessary wrapping**

  **File:** `algebraic_zero_tensor.py:244-280`

  `Z - T` returns `AlgebraicTensor(Z, -T)` which immediately simplifies to `-T`. These should directly return `-other` and `other` respectively.

  **Fix:** `__sub__` should return `-other`. `__rsub__` should return `other` (or equivalently `AlgebraicTensor(other)` which unwraps to `other`).

---

- [] **16. `_is_zero_like` — fragile test**

  **File:** `algebraic_pure_tensor.py:67-73`

  ```python
  def _is_zero_like(expr):
      return expr == S.Zero * expr
  ```

  Relies on `S.Zero * expr` evaluating to a canonical zero. For a custom type with broken `__mul__`, this could silently return `True`.

  **Fix:** Add an explicit check for `expr is S.Zero` or `hasattr(expr, 'is_zero') and expr.is_zero is True` before the multiplication test.

---

- [] **17. `AlgebraicPureTensor.__new__` — zero coeff with no factors**

  **File:** `algebraic_pure_tensor.py:327-329`

  ```python
  if not processed:
      if coeff is S.Zero:
          return AlgebraicZeroTensor(())
  ```

  Empty shape `()` is inconsistent — all other zero tensors carry their actual shape. This path is hard to reach but produces a degenerate object.

  **Fix:** Raise `ValueError` since a zero coefficient with no tensor factors is not a valid zero tensor (there is no shape to anchor it to).

---

## MINOR / STYLE ISSUES

- [] **18. Excessive doctest duplication**

  Module-level docstrings, class docstrings, and method docstrings all repeat identical examples. `AlgebraicZeroTensor` class docstring (lines 68-99) is nearly identical to the module docstring (lines 18-43).

  **Fix:** Keep examples at the class level only. Remove duplicate examples from module-level docstrings, or vice versa.

---

- [] **19. Over-documented trivial methods**

  `copy()` (`algebraic_zero_tensor.py:349-350`), `__bool__` (line 324-325), `has_zero_term` on `AlgebraicPureTensor` (`algebraic_pure_tensor.py:571-573`) are one-liners with full docstrings.

  **Fix:** Study some other modules in Sympy and see how such simple methods are documented. If they have docstrings, keep them also here. If it is not standard practice to include docstrings for trivial methods, remove them from the methods.

---

- [x] **20. `_try_merge_for_slot` — O(N²) algorithm**

  **File:** `simplify.py:248-450` (removed)

  The while-with-pivot loop is O(N²) per slot. `_equality_merge_dict` (line 453) was written as an O(N) replacement for the equality case, but `_proportionality_factoring` still uses the O(N²) version. For tensors with dozens of terms, this is a performance concern.

  **Fix:** Remove any reliance on this function. Remove _proportionality_factoring if it is still present in the codebase.

  **[FIXED]** Removed `_try_merge_for_slot` along with `_proportionality_factoring` (issue #4).

---

- [] **21. Inconsistent use of `_sympify` parameter**

  `AlgebraicTensor.__new__` accepts `_sympify=True` (line 388), but it's only checked at the top level. When `_merge` or `_subtract` construct new `AlgebraicTensor` via `Basic.__new__`, they bypass sympification entirely.

  **Fix:** Add a comment documenting the intentional bypass, or refactor `_merge`/`_subtract` to call `AlgebraicTensor.__new__` with `_sympify=False`.

---

- [] **22. `algebraic_tensor_product` — redundant zero check**

  **File:** `algebraic_tensor.py:1507-1513`

  Zero-like check (line 1507) and zero-tensor sentinel check (line 1511) are nearly identical. Both produce `AlgebraicZeroTensor(_combined_shape(args))`.

  **Fix:** Merge into a single check or extract a shared helper.

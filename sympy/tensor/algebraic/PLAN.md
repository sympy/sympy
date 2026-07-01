# Simplification Pipeline Plan for `_simplify_algebraic_tensor`

## Overview

This document outlines the implementation plan for enhancing `_simplify_algebraic_tensor` in `/tensor/algebraic/simplify.py` with a new multi-step simplification algorithm that leverages the `commutativity_shape` property of `AlgebraicTensor`.

## Current State

Currently, `_simplify_algebraic_tensor` only performs `proportionality_factoring` as its simplification step. This is insufficient for fully simplifying algebraic tensors that contain commutative and non-commutative components that could be factored more aggressively.

## New Algorithm Overview

The new algorithm will decompose an `AlgebraicTensor` based on its `commutativity_shape`, separate commutative and non-commutative components, perform simplifications on each component, and then reconstruct the full tensor expression.

### High-Level Steps

1. **Analyze commutativity structure** - Extract which tensor factor slots are commutative (value 1) vs non-commutative (value 0)
2. **Build commutative subtensor dictionary** - Create a lookup for commutative subtensors (factors with all-1 commutativity_shape)
3. **Process each pure tensor term** - For each term in the input `AlgebraicTensor`:
   - Extract non-commutative subtensor (using commutativity_shape of input `at`)
   - Extract and expand commutative subtensor (fully expand to isolate symbols/numbers)
   - Use expanded terms to populate dictionary with commutative subtensors as keys
4. **Apply proportionality_factoring to dictionary values** - Merge proportional non-commutative subtensors
5. **Reconstruct full tensor expression** - Rebuild terms from dictionary keys and values
6. **Apply final proportionality_factoring** - Final merge pass
7. **Simplify individual factors** - Apply `.simplify()` to each factor in each pure tensor

---

## Detailed Implementation Plan

### Step 0: Read commutativity_shape of input `at`

```python
comm_cs = at.commutativity_shape  # tuple of 0s and 1s, same length as tensor_shape
```

### Step 1: Identify commutative and non-commutative factor positions

Create two lists to track factor indices:
- `commutative_indices = [i for i, val in enumerate(comm_cs) if val == 1]`
- `non_commutative_indices = [i for i, val in enumerate(comm_cs) if val == 0]`

These lists will be used throughout to extract the correct factors from each pure tensor.

### Step 2: Build commutative subtensor dictionary

Create a dictionary where:
- **Keys**: `AlgebraicPureTensor` objects representing commutative subtensors
  - These have `tensor_shape` as a sub-tuple of `at.tensor_shape` for indices where `comm_cs == 1`
  - Each factor in these subtensors is a matrix with exactly one nonzero entry equal to 1
  - Example: `Matrix([[1, 0], [0, 0]])` or `Matrix([[0, 0], [0, 1]])`

- **Values**: `AlgebraicTensor` (or single `AlgebraicPureTensor`/`AlgebraicZeroTensor`) representing accumulated non-commutative subtensors

**Dictionary initialization**:
- For each combination of commutative factor slots, create the canonical basis matrices
- The number of possible keys depends on the number of commutative factor slots and their shapes
- For each commutative slot `i` with shape `(m_i, n_i)`, there are `m_i * n_i` possible basis matrices
- Total possible keys = product of `(m_i * n_i)` over all commutative slots

**Practical approach**:
- Don't pre-populate all possible keys
- Initialize `comm_dict = {}`
- As we process terms, we'll create keys dynamically
- Default value for new keys: `AlgebraicZeroTensor` with shape from non-commutative factor slots

### Step 3: Process each pure tensor term

For each term in `at.args`:

#### Step 3.0: Normalize term to (term_coeff, factors) form

Terms in `at.args` may be `AlgebraicPureTensor`, `Mul(coeff, AlgebraicPureTensor)`, `AlgebraicZeroTensor`, or bare matrices. Normalize before extracting factors:

```python
from sympy.tensor.algebraic.simplify import _extract_pt_and_coeff

# _extract_pt_and_coeff returns (pt_or_none, coeff, optional_factors_list)
pt_or_none, term_coeff, opt_factors = _extract_pt_and_coeff(term)

if isinstance(term, AlgebraicZeroTensor):
    continue  # zero terms contribute nothing to the dictionary

if opt_factors is not None:
    # Direct matrix factors (bypass PureTensor wrapper) — use as-is
    factors = opt_factors
elif pt_or_none is not None:
    factors = list(pt_or_none.factors)
else:
    # Fallback: treat the whole term as a single factor with unit coefficient
    factors = [term]
    term_coeff = S.One
```

This ensures that by Step 3.1, `factors` is a clean list of matrix factors and `term_coeff` is the combined commutative coefficient.

#### Step 3.1: Extract non-commutative factors

Using `non_commutative_indices` from Step 1:
```python
non_commutative_factors = [factors[i] for i in non_commutative_indices]
```

Note: Always use the input `at`'s `commutativity_shape`, not the term's own `commutativity_shape`.

#### Step 3.2: Construct non-commutative subtensor (or handle degenerate cases)

If `non_commutative_indices` is empty (all slots are commutative), there is no non-commutative subtensor — the "value" stored in the dictionary is just a commutative coefficient (a SymPy `Number` or symbolic expression). In this case set a flag and skip construction:
```python
if non_commutative_indices:
    non_commutative_pt = AlgebraicPureTensor(*non_commutative_factors)
    # WARNING: AlgebraicPureTensor.__new__ may unwrap to a bare matrix,
    # AlgebraicZeroTensor, or a raw factor (see DESIGN.md §2 "Constructor behavior").
    # The dictionary values in Step 3.4 must be prepared to accept any of these types.
else:
    non_commutative_pt = None  # values will be plain commutative coefficients
```

#### Step 3.3: Iterative decomposition of commutative factors into basis matrices

Using `commutative_indices` from Step 1:
```python
commutative_factors = [factors[i] for i in commutative_indices]
```

**Goal**: Decompose every commutative factor into a sum of basis matrices (entries are 0 or 1 only), accumulating all symbolic coefficients into a single prefactor. **Do NOT use `.expand()`** — it may destroy pre-existing symbolic factorizations in the matrix entries.

**Data representation for intermediate terms**:

Each intermediate term is a tuple `(prefactor, factor_list)`:
- `prefactor`: commutative SymPy expression (initially the term's coefficient from `_get_coeff()`)
- `factor_list`: list of matrix factors for the commutative slots; entries at positions already processed contain basis matrices (0/1 only), unprocessed positions still hold the original matrix

**Algorithm — iterative slot-by-slot decomposition**:

Start with a single-term seed list:
```python
terms = [(term_coeff, list(commutative_factors))]
```

For each commutative slot `s` (iterate `s = 0, 1, ..., len(commutative_factors)-1`):
1. Take the current `terms` list
2. For each `(prefactor, factor_list)` in `terms`:
   a. Let `M = factor_list[s]` (the matrix at slot `s`)
   b. Iterate over all matrix entries `(r, c)` of `M`:
      - Let `entry = M[r, c]`
      - If `entry == 0`: skip
      - Create basis matrix `E_rc` of same shape as `M`, with 1 at `(r, c)` and 0 elsewhere
      - New prefactor: `new_coeff = prefactor * entry`
      - New factor list: copy of `factor_list` with `factor_list[s]` replaced by `E_rc`
      - Append `(new_coeff, new_factor_list)` to the next round's term list
3. Replace `terms` with the accumulated next round's list

After processing all slots, every `factor_list` contains only basis matrices (0/1 entries), and all symbolic content has been extracted into `prefactor`.

**Example**:
```python
# Input term: AlgebraicPureTensor(3, Matrix([[2, C], [0, 0]]), Matrix([[1, 0], [0, 4]]))
# term_coeff = 3
# commutative_factors = [Matrix([[2, C], [0, 0]]), Matrix([[1, 0], [0, 4]])]

# --- Process slot 0 ---
# M = Matrix([[2, C], [0, 0]])
# Nonzero entries: (0,0)->2, (0,1)->C
# terms after slot 0:
#   (3*2,        [E_00, Matrix([[1,0],[0,4]])])
#   (3*C,        [E_01, Matrix([[1,0],[0,4]])])

# --- Process slot 1 ---
# For (6, [E_00, Matrix([[1,0],[0,4]])]):
#   M = Matrix([[1,0],[0,4]]), nonzero: (0,0)->1, (1,1)->4
#   Produces:
#     (6*1,  [E_00, E_00])
#     (6*4,  [E_00, E_11])
# For (3*C, [E_01, Matrix([[1,0],[0,4]])]):
#   Produces:
#     (3*C*1,  [E_01, E_00])
#     (3*C*4,  [E_01, E_11])

# Final terms:
#   (6,   [E_00, E_00])
#   (24,  [E_00, E_11])
#   (3*C, [E_01, E_00])
#   (12*C,[E_01, E_11])
```

**Important properties**:
- The prefactor always carries the full symbolic coefficient. The factor lists contain only 0/1 basis matrices after the loop completes.
- Zero entries are skipped early, so they never produce terms.
- If a matrix factor has no nonzero entries (all-zero matrix), the slot produces no children and the parent term is effectively eliminated (its contribution is zero).
- The prefactor is multiplied cumulatively. SymPy's own `Mul` simplification handles coefficient arithmetic (e.g., `3*2` → `6`, `3*C*4` → `12*C`).

**Output**: A list of `(prefactor, basis_factor_list)` tuples, one per nonzero combination of basis matrices across all commutative slots.

#### Step 3.4: Populate dictionary

For each `(prefactor, basis_factor_list)` tuple from Step 3.3:

##### Step 3.4.1: Build commutative subtensor key

Construct the dictionary key from the basis matrices. If there are commutative slots, build an `AlgebraicPureTensor` from the basis matrices. If there are no commutative slots (all slots non-commutative), use `None` as the single universal key:
```python
if commutative_indices:
    commutative_key = AlgebraicPureTensor(*basis_factor_list)
    # WARNING: may unwrap to a bare basis matrix if there is only one slot
    # and the coefficient is S.One. Use the raw basis_factor_list tuple as
    # a fallback key if needed, or wrap in a tuple to guarantee hashability:
    commutative_key = tuple(basis_factor_list)  # safe, hashable key
else:
    commutative_key = None  # single bucket — all terms share the same key
```

Using a tuple of matrices as the key avoids all unwrapping issues and is hashable (SymPy matrices implement `__hash__`).

##### Step 3.4.2: Add to dictionary value

**Case A — non-commutative slots exist (`non_commutative_pt` is not None):**

Multiply the non-commutative subtensor by the prefactor and accumulate into the dictionary value. Be careful: `prefactor * non_commutative_pt` dispatches through `__rmul__` / `__mul__` and may return a `Mul(prefactor, AlgebraicPureTensor)` wrapper, a scaled `AlgebraicPureTensor`, or (if prefactor is S.One) the bare `non_commutative_pt`. The addition `value + scaled` similarly may unwrap.

```python
if commutative_key not in comm_dict:
    # Create zero tensor with shape from non_commutative_factors
    shape = tuple(f.shape for f in non_commutative_factors)
    comm_dict[commutative_key] = AlgebraicZeroTensor(shape)

# Scale non_commutative_pt by prefactor. Since non_commutative_pt is
# non-commutative, scalar multiplication goes through __rmul__/__mul__
# which scales the internal coefficient. Result may be a PureTensor,
# Mul(coeff, PureTensor), or AlgebraicZeroTensor (if prefactor is zero).
scaled = prefactor * non_commutative_pt

# Accumulate. The + operator routes through AlgebraicTensor.__new__ via
# the add dispatcher, which handles shape validation and flattening.
value = comm_dict[commutative_key]
comm_dict[commutative_key] = value + scaled
```

If `prefactor` evaluates to `S.Zero`, the term contributes nothing and can be skipped entirely.

**Case B — no non-commutative slots (`non_commutative_pt` is None):**

The value is a plain commutative coefficient. Accumulate arithmetically:
```python
if commutative_key not in comm_dict:
    comm_dict[commutative_key] = S.Zero

comm_dict[commutative_key] = comm_dict[commutative_key] + prefactor
```

### Step 4: Apply proportionality_factoring to dictionary values

For each key in `comm_dict`:
```python
from sympy.tensor.algebraic.simplify import proportionality_factoring

value = comm_dict[key]
if isinstance(value, AlgebraicTensor):
    comm_dict[key] = proportionality_factoring(value)
```

This merges proportional non-commutative subtensors that share the same commutative pattern.

### Step 5: Reconstruct full tensor expression

Collect all reconstructed terms into a flat list, then assemble once at the end.

#### Step 5.1: Iterate dictionary and decompose values into individual terms

For each `(key, value)` pair in `comm_dict`:

```python
all_reconstructed = []

for key, value in comm_dict.items():
    # Normalize value to a list of (coeff, non_commutative_body) pairs.
    # The non_commutative_body is either an AlgebraicPureTensor (Case A)
    # or S.One (Case B, all-commutative slots).

    if non_commutative_indices:
        # Case A: value is a tensor expression (AlgebraicTensor,
        # AlgebraicPureTensor, AlgebraicZeroTensor, or Mul wrapper).
        if isinstance(value, AlgebraicTensor):
            bodies = list(value.args)
        elif isinstance(value, (AlgebraicPureTensor, AlgebraicZeroTensor)):
            bodies = [value]
        elif value.is_Mul:
            # Mul(coeff, AlgebraicPureTensor) — keep as single term
            bodies = [value]
        else:
            bodies = [value]

        for body in bodies:
            if isinstance(body, AlgebraicZeroTensor):
                # Zero non-commutative part: still need to reconstruct a
                # properly-shaped zero tensor with the commutative factors.
                # Build a full-shape AlgebraicZeroTensor and add it.
                all_reconstructed.append(AlgebraicZeroTensor(at.tensor_shape))
                continue

            # Extract coefficient and PureTensor from the body.
            if isinstance(body, AlgebraicPureTensor):
                coeff = body._get_coeff()
                pt = body
            elif body.is_Mul:
                # Mul(outer_coeff, AlgebraicPureTensor(..., inner_coeff, factors))
                # Separate the outer commutative coefficient from the inner PureTensor.
                outer_coeff, inner = body.args[0], body.args[1]
                if isinstance(inner, AlgebraicPureTensor):
                    coeff = outer_coeff * inner._get_coeff()
                    pt = inner
                else:
                    coeff, pt = body, None  # fallback, shouldn't happen
            else:
                coeff, pt = S.One, body

            reconstructed = _reconstruct_term(key, pt, coeff,
                                               comm_cs, commutative_indices,
                                               non_commutative_indices)
            if reconstructed is not None:
                all_reconstructed.append(reconstructed)
    else:
        # Case B: all slots commutative. Value is a plain coefficient.
        # Reconstruct a commutative-only PureTensor from the key.
        if value == S.Zero:
            all_reconstructed.append(AlgebraicZeroTensor(at.tensor_shape))
            continue
        reconstructed = _reconstruct_term(key, None, value,
                                           comm_cs, commutative_indices,
                                           non_commutative_indices)
        if reconstructed is not None:
            all_reconstructed.append(reconstructed)
```

#### Step 5.2: `_reconstruct_term` signature and behavior

```python
def _reconstruct_term(key, non_commutative_pt, coeff, comm_cs,
                      commutative_indices, non_commutative_indices):
```

1. Extract commutative basis matrices from `key` (a tuple of matrices from Step 3.4.1)
2. Extract non-commutative factors from `non_commutative_pt` (or skip if None)
3. Interleave them back into the original factor order:
   - Walk `i = 0, 1, ..., len(comm_cs)-1`
   - If `comm_cs[i] == 1`: take next factor from the commutative key
   - If `comm_cs[i] == 0`: take next factor from `non_commutative_pt.factors`
4. Build the result with `_build_pt(coeff, merged_factors)` (from simplify.py) to handle unwrapping correctly:
   - `coeff == 0` → `AlgebraicZeroTensor`
   - single factor, `coeff == 1` → bare factor
   - otherwise → `AlgebraicPureTensor(coeff, *factors)` or `AlgebraicPureTensor(*factors)`

**Example**:
```python
# Original at.tensor_shape = ((2,2), (2,2), (2,2))
# comm_cs = (1, 0, 1)
# non_commutative_indices = [1]
# commutative_indices = [0, 2]

# key = (Matrix([[1,0],[0,0]]), Matrix([[0,1],[0,0]]))  # tuple, positions 0 and 2
# non_commutative_pt = AlgebraicPureTensor(M)  # position 1
# coeff = 5

# Merged factors: [Matrix([[1,0],[0,0]]), M, Matrix([[0,1],[0,0]])]
# Result: _build_pt(5, merged_factors)
#   → AlgebraicPureTensor(5, Matrix([[1,0],[0,0]]), M, Matrix([[0,1],[0,0]]))
```

**Return**: `AlgebraicPureTensor`, bare matrix, or `AlgebraicZeroTensor` (never returns `None` unless coeff is zero and the caller should skip — use `_build_pt` which handles this).

### Step 6: Assemble final sum from reconstructed terms

After Step 5.1, `all_reconstructed` contains every reconstructed term (possibly including `AlgebraicZeroTensor` shape anchors). Filter and assemble:

```python
# Filter out S.Zero (shouldn't appear, but defensive)
real_terms = [t for t in all_reconstructed if t is not S.Zero]

if not real_terms:
    # Everything cancelled — return the original zero tensor for shape preservation
    final_result = AlgebraicZeroTensor(at.tensor_shape)
elif len(real_terms) == 1:
    final_result = real_terms[0]
else:
    final_result = AlgebraicTensor(*real_terms, _sympify=False)
```

### Step 7: Apply proportionality_factoring

Apply final proportionality factoring to handle any remaining proportional terms across different commutative patterns:
```python
if isinstance(final_result, AlgebraicTensor):
    final_result = proportionality_factoring(final_result)
```

### Step 8: Simplify individual factors

Apply SymPy's `simplify` to the coefficient and `.simplify()` to each factor in each pure tensor. Handle `Mul(coeff, AlgebraicPureTensor)` wrappers and `AlgebraicZeroTensor` pass-through:

```python
from sympy import simplify as sympy_simplify

if isinstance(final_result, AlgebraicTensor):
    new_args = []
    for arg in final_result.args:
        if isinstance(arg, AlgebraicZeroTensor):
            # Zero tensor — nothing to simplify
            new_args.append(arg)
        elif isinstance(arg, AlgebraicPureTensor):
            coeff = sympy_simplify(arg._get_coeff())
            factors = [f.simplify() if hasattr(f, 'simplify') else f
                       for f in arg.factors]
            # Use _build_pt to handle unwrapping correctly
            new_arg = _build_pt(coeff, factors)
            new_args.append(new_arg)
        elif arg.is_Mul:
            # Mul(outer_coeff, AlgebraicPureTensor) — simplify both parts
            outer_coeff, inner = arg.args[0], arg.args[1]
            if isinstance(inner, AlgebraicPureTensor):
                inner_coeff = sympy_simplify(inner._get_coeff())
                factors = [f.simplify() if hasattr(f, 'simplify') else f
                           for f in inner.factors]
                combined_coeff = sympy_simplify(outer_coeff) * inner_coeff
                new_arg = _build_pt(combined_coeff, factors)
                new_args.append(new_arg)
            else:
                new_args.append(arg)  # unexpected structure, pass through
        else:
            new_args.append(arg)  # unexpected type, pass through
    final_result = AlgebraicTensor(*new_args, _sympify=False)
elif isinstance(final_result, AlgebraicPureTensor):
    final_result = _simplify_algebraic_pure_tensor(final_result, **kwargs)
```

Note: `_build_pt` (from simplify.py) handles all unwrapping cases: zero coefficient produces `AlgebraicZeroTensor`, single factor with coeff=1 unwraps to bare matrix, etc. This avoids the pitfall of calling `AlgebraicPureTensor(coeff, *factors)` and getting an unexpected return type.

---

## Helper Functions to Implement

### 1. `_decompose_commutative_factors(commutative_factors, term_coeff)`

Implements the iterative slot-by-slot decomposition from Step 3.3.

**Input**:
- `commutative_factors`: list of commutative matrix factors
- `term_coeff`: initial commutative coefficient for the term

**Output**: List of `(prefactor, basis_factor_list)` tuples

**Implementation**:
- Seed with `[(term_coeff, list(commutative_factors))]`
- For each slot, iterate matrix entries, skip zeros, create basis matrices, accumulate prefactors
- Return final list

### 2. `_reconstruct_term(key, non_commutative_pt, coeff, comm_cs, commutative_indices, non_commutative_indices)`

Reconstructs a full tensor from commutative and non-commutative subtensors.

**Input**:
- `key`: tuple of commutative basis matrices (or `None` if no commutative slots)
- `non_commutative_pt`: `AlgebraicPureTensor` of non-commutative factors (or `None` if all slots commutative)
- `coeff`: combined commutative coefficient
- `comm_cs`: original commutativity_shape tuple
- `commutative_indices`: list of commutative slot indices
- `non_commutative_indices`: list of non-commutative slot indices

**Output**: `AlgebraicPureTensor`, bare matrix, or `AlgebraicZeroTensor` (via `_build_pt`)

**Implementation**:
- Interleave commutative and non-commutative factors back into original order
- Delegate to `_build_pt(coeff, merged_factors)` for safe construction

### 3. `_extract_factors_by_commutativity(pt, comm_cs)`

Helper to extract factors from a pure tensor based on commutativity indices.

**Input**:
- `pt`: `AlgebraicPureTensor` or expression with `.factors` attribute
- `comm_cs`: tuple of 0s and 1s

**Output**: Tuple of (commutative_factors_list, non_commutative_factors_list)

---

## Integration Points

### Modify `_simplify_algebraic_tensor` in `simplify.py`

Replace the current single-step implementation:

```python
def _simplify_algebraic_tensor(at, **kwargs):
    """Simplify an AlgebraicTensor (sum of same-shape tensor terms)."""
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
```

With the new multi-phase pipeline:

```python
def _simplify_algebraic_tensor(at, **kwargs):
    """Simplify an AlgebraicTensor (sum of same-shape tensor terms)."""
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    # Phase 0: proportionality factoring (keep existing behavior)
    result = proportionality_factoring(at)

    # If result is not an AlgebraicTensor, simplify and return
    if not isinstance(result, AlgebraicTensor):
        if isinstance(result, AlgebraicPureTensor):
            return _simplify_algebraic_pure_tensor(result, **kwargs)
        return result

    # Phase 1-8: New commutativity-based simplification
    result = _commutativity_simplify(result, **kwargs)

    return result
```

### New function: `_commutativity_simplify(at, **kwargs)`

This function implements Steps 0-8 of the new algorithm.

---

## Edge Cases to Handle

1. **Empty commutative_indices** (all slots non-commutative)
   - `commutative_factors` is empty; Step 3.3 produces a single `(term_coeff, [])` tuple
   - Dictionary key is `None` (single bucket); Step 3.4 Case A applies normally
   - Reconstruction in Step 5 simply passes through the non-commutative subtensor unchanged

2. **Empty non_commutative_indices** (all slots commutative)
   - `non_commutative_pt` is `None`; Step 3.4 Case B applies (values are plain coefficients)
   - Step 5 reconstructs commutative-only `AlgebraicPureTensor` from dictionary keys
   - No `AlgebraicPureTensor()` with zero factors is ever constructed

3. **Zero prefactors after expansion**
   - If `prefactor` evaluates to `S.Zero`, the term is skipped in Step 3.4.2
   - If an entire commutative factor is all-zero, Step 3.3 produces no children for that parent term, effectively eliminating it

4. **Single-term input `at`**
   - After Phase 0, `proportionality_factoring` may have already unwrapped to `AlgebraicPureTensor`
   - The early return guard (non-AlgebraicTensor check) handles this before reaching `_commutativity_simplify`

5. **Cancellation resulting in zero**
   - Step 6 detects empty `real_terms` and returns `AlgebraicZeroTensor(at.tensor_shape)`
   - Shape information is preserved from the original input

6. **Nested `Mul` terms in `at.args`**
   - Step 3 must extract the coefficient before decomposing. Use `_extract_pt_and_coeff` (from simplify.py) to handle `Mul(coeff, AlgebraicPureTensor)`, bare `AlgebraicPureTensor`, and other forms uniformly
   - Step 5.1 also handles `Mul` wrappers in dictionary values after proportionality_factoring

7. **Term coefficient extraction in Step 3**
   - Before Step 3.1, each `term` from `at.args` must be normalized to `(term_coeff, factors)` form
   - Use the same extraction logic as Step 5.1: handle `AlgebraicPureTensor` directly, unwrap `Mul(outer_coeff, AlgebraicPureTensor)` by multiplying outer and inner coefficients, and skip `AlgebraicZeroTensor` (it contributes nothing to the dictionary)

---

## Testing Strategy

1. **Unit tests for helper functions**:
   - `_decompose_commutative_factors` with various matrix types (numeric, symbolic, mixed, all-zero)
   - `_reconstruct_term` with various index combinations, including degenerate cases (all commutative, all non-commutative)
   - `_extract_factors_by_commutativity` with edge cases
   - `_build_pt` unwrapping behavior (zero coeff, single factor, coefficient S.One)

2. **Integration tests**:
   - Input with all commutative slots (Case B path)
   - Input with all non-commutative slots (Case A, key=None path)
   - Input with mixed commutativity
   - Inputs that should produce cancellation
   - Inputs with zero coefficients
   - Terms with `Mul(coeff, AlgebraicPureTensor)` wrappers in `at.args`

3. **Regression tests**:
   - Verify existing `proportionality_factoring` behavior still works (Phase 0)
   - Verify `.simplify()` on individual factors is applied (Step 8)
   - Verify pre-existing symbolic factorizations in matrix entries are NOT destroyed by the decomposition (regression against the old `.expand()` approach)

---

## Performance Considerations

1. **Dictionary key creation**:
   - May create many keys if many commutative slots with large shapes
   - Consider using sparse representation if memory becomes an issue

2. **Expansion complexity**:
   - Expanding factors with many nonzero entries creates many terms
   - Consider early filtering for zero entries

3. **Proportionality checking**:
   - Applied to dictionary values which may still be large
   - Consider early filtering for obvious non-proportional cases

---

## References

- Current `proportionality_factoring` implementation in `simplify.py:150`
- `commutativity_shape` property in `algebraic_pure_tensor.py:84`
- `tensor_shape` property in `algebraic_pure_tensor.py:76`
- `_simplify_algebraic_pure_tensor` for factor-wise simplification pattern

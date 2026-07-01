# Algebraic Tensor Module — Design Document

## Purpose

This module implements algebraic tensor expressions in SymPy. Tensors are built from
matrix-like factors combined via the **tensor product** (`⊗`) and **addition** (`+`).
A third operation, **composition** (`*`), performs factor-wise matrix multiplication
and serves as the tensor analogue of contraction.

The module provides three core types plus a simplification pipeline:

| Class | Role |
|---|---|
| `AlgebraicPureTensor` | A single tensor-product term (non-commutative Mul of matrix factors, optionally with a commutative coefficient) |
| `AlgebraicZeroTensor` | The additive identity (zero tensor) for a given tensor shape |
| `AlgebraicTensor` | A linear combination (sum) of same-shape `AlgebraicPureTensor` terms |
| `simplify.py` | Proportionality factoring and per-term simplification |

---

## 1. Tensor Shape Model

**Shape is the fundamental invariant.** Every tensor object carries a `tensor_shape`
that is a tuple of per-factor `(rows, cols)` pairs.

```
A ⊗ B    →  tensor_shape = ((m, n), (n, p))    # A is m×n, B is n×p
C        →  tensor_shape = ((m, n),)            # single factor C is m×n
```

- A bare matrix `M` with `M.shape == (m, n)` is treated as a single-factor tensor
  with `tensor_shape = ((m, n),)`.
- `AlgebraicZeroTensor` stores its shape explicitly in `_shape`.
- Two tensors can be **added** only if their `tensor_shape` tuples are identical.
- Two tensors can be **composed** only if they have the same number of factors and
  each corresponding factor pair has matching inner dimensions.

### Shape normalization

`AlgebraicZeroTensor.__new__` normalizes input shapes:
- `(3, 4)` → `((3, 4),)` (bare pair wrapped)
- `[(3, 4)]` → `((3, 4),)` (list converted to tuple of tuples)
- `((3, 4), (4, 5))` → unchanged

**Agent rule:** Always compare full `tensor_shape` tuples for equality. Never
compare only the outer `shape` of individual factors.

---

## 2. AlgebraicPureTensor

### Idea

An `AlgebraicPureTensor` represents a *single term* in a tensor expression — the
non-commutative tensor product of matrix factors, optionally scaled by a
commutative coefficient.

```
A ⊗ B ⊗ C          →  AlgebraicPureTensor(A, B, C)
2 * (A ⊗ B)        →  AlgebraicPureTensor(2, A, B)
-A ⊗ B             →  AlgebraicPureTensor(-1, A, B)
```

### Inheritance: extends `Mul`

`AlgebraicPureTensor` extends `sympy.core.mul.Mul` (not `Basic`) to inherit SymPy's
built-in simplification, differentiation, and pattern-matching machinery.

**Critical:** `is_Mul = False` is set to prevent `Mul.flatten` from unpacking the
tensor product factors. Without this guard, SymPy would try to flatten/commute the
factors, destroying the non-commutative tensor-product ordering.

### Key attributes and properties

| Attribute | Description |
|---|---|
| `is_AlgebraicPureTensor` | Type tag for `isinstance` checks |
| `is_Mul = False` | Prevents `Mul.flatten` from unpacking this object |
| `is_commutative = False` | Tensor product is non-commutative |
| `_op_priority = 11` | Ensures `x * pt` delegates to `pt.__rmul__(x)` for commutative scalars |
| `_eval_is_commutative = lambda self: False` | Descriptor to enforce non-commutativity |

### `commutativity_shape`

The `commutativity_shape` property returns a tuple of binary entries, one per
tensor factor. Entry `i` is `1` if the `i`-th factor contains no non-commutative
symbols (e.g., a concrete numeric matrix), and `0` otherwise (e.g., a
`MatrixSymbol`). Same length as `tensor_shape`.

```
AlgebraicPureTensor(A_3x4, C_4x5).commutativity_shape → (0, 0)   # both symbolic
AlgebraicPureTensor(numeric_3x4, C_4x5).commutativity_shape → (1, 0)
```

This property is used by `AlgebraicTensor._flatten_args` to compute the
component-wise AND of all term commutativity shapes.

### Coefficient storage

A leading commutative coefficient is stored as the **first element of `self.args`**:

```
AlgebraicPureTensor(2, A, B).args → (2, A, B)
AlgebraicPureTensor(A, B).args    → (A, B)          # implicit coefficient S.One
```

The `_get_coeff()` method extracts this coefficient (returns `S.One` if absent).
The `factors` property returns only the tensor-product factors, excluding the
leading coefficient.

**Agent rule:** When iterating factors, always use `.factors` (not `.args`) to
skip the coefficient. When checking for a coefficient, use `._get_coeff()`.

### Constructor behavior (`__new__`)

1. **Empty args** → `ValueError`
2. **All factors are scalar/zero** → returns `AlgebraicZeroTensor`
3. **Single factor with coefficient S.One** → returns the bare factor directly
   (unwraps)
4. **Coefficient is S.Zero** → returns `AlgebraicZeroTensor(factor_shapes)`
5. **Multiple factors** → constructs `Mul` subclass instance with `evaluate=False`

**Agent rule:** Never assume the result of `AlgebraicPureTensor(...)` is an
`AlgebraicPureTensor` instance. It may return a bare matrix, an
`AlgebraicZeroTensor`, or a raw factor. Always check the return type.

### Arithmetic operators

| Operation | Behavior |
|---|---|
| `__neg__` | Negates the coefficient |
| `__mul__` | If `other` is commutative: scales coefficient. If non-commutative: calls `compose_algebraic_tensors(self, other)` |
| `__rmul__` | Same as `__mul__` but with operands reversed |
| `__add__` / `__radd__` | Returns `AlgebraicTensor(self, other)` |
| `__sub__` / `__rsub__` | Returns `AlgebraicTensor(self, -other)` |

**Agent rule:** Multiplication (`*`) on tensors is **composition**, not element-wise
multiplication. For scalars, it behaves as scalar scaling. The dispatch happens
inside `__mul__` / `__rmul__` based on whether the other operand is commutative.

---

## 3. AlgebraicZeroTensor

### Idea

The additive identity for tensors of a specific shape. Unlike SymPy's `S.Zero`,
an `AlgebraicZeroTensor` carries shape information, so `0_((3,4),)` and
`0_((2,2),)` are different objects.

### Inheritance: extends `Atom`

`AlgebraicZeroTensor` extends `sympy.core.basic.Atom` (a leaf `Basic` subclass)
to integrate with SymPy's expression system. This provides:

- **`sympify()` support** — `sympify(z)` returns `z` unchanged.
- **Tree traversal** — `atoms()`, `has()`, `replace()` work correctly.
- **Assumptions system** — `is_zero = True`, `is_commutative = True` are
  recognized by the SymPy assumptions machinery.
- **Generic operations** — `subs()`, `xreplace()`, `doit()` are inherited.

The shape is stored in a dedicated `_shape` slot (not in `_args`). This is
required because `Basic` expects `_args` to contain only `Basic` objects —
storing a plain tuple in `_args` would break `atoms()`, `free_symbols`,
`has()`, and `count_ops()`.

Key implementation details:
- `__slots__ = ('_shape',)` — shape lives in its own slot.
- `Atom.__new__(cls)` is called with **no args**; `_shape` is set afterward.
- `_hashable_content()` returns `(self._shape,)` for hashing and equality.
- `__getnewargs__()` returns `(self._shape,)` for pickle reconstruction.
- `free_symbols` returns `set()` (a zero tensor has no free symbols).
- `copy()` returns `self` (immutable atom).

### `commutativity_shape`

A zero tensor is commutative in every slot. The `commutativity_shape` property
returns an all-1s tuple matching the length of the tensor shape:

```
AlgebraicZeroTensor(((3, 4), (4, 5))).commutativity_shape → (1, 1)
```

### Arithmetic operators

| Operation | Behavior |
|---|---|
| `__neg__` | Returns self (zero is its own negation) |
| `__add__` / `__radd__` | Returns `other` (identity) |
| `__sub__` / `__rsub__` | Delegates to `AlgebraicTensor` constructor |
| `__mul__` / `__rmul__` | Returns self (composition of zero is zero of same shape) |
| `__bool__` | Returns `False` |

**Agent rule:** When `AlgebraicPureTensor` or `AlgebraicTensor` operations produce
a zero result, they return `AlgebraicZeroTensor(shape)`, not `S.Zero`. The shape
information is critical for type safety.

---

## 4. AlgebraicTensor

### Idea

An `AlgebraicTensor` represents a **linear combination** of `AlgebraicPureTensor`
terms that all share the same `tensor_shape`. Conceptually, it is a sum:

```
A ⊗ B + 2*(C ⊗ D)    →  AlgebraicTensor(AlgebraicPureTensor(A, B),
                                          AlgebraicPureTensor(2, C, D))
```

### Inheritance: extends `Basic` (not `Add`)

`AlgebraicTensor` extends `sympy.core.basic.Basic` rather than `Add` to avoid the
`is_commutative` descriptor conflict in `AssocOp._from_args`. Since tensor terms
are non-commutative, using `Add` as a base would cause SymPy's flattening logic
to reorder/merge terms incorrectly.

**Critical:** `is_Add = True` is set so the wider SymPy ecosystem recognizes this
as an additive expression. The `add` dispatcher is registered:

```python
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)
```

This ensures that expressions like `(A⊗B) + (C⊗D)` are routed through
`AlgebraicTensor.__new__` instead of `Add.__new__`.

**Agent rule:** Always use `Basic.__new__` when creating subclasses that need
`is_Add = True` but non-commutative semantics. Never inherit from `Add` for
non-commutative sums.

### Constructor behavior (`__new__`)

The constructor uses `_flatten_args` to:
1. **Flatten** nested `AlgebraicTensor` instances (collect all leaf terms)
2. **Validate** that all terms share the same `tensor_shape`
3. **Collect** `AlgebraicZeroTensor` as a separate anchor term
4. **Filter** out `S.Zero` identity numbers

Return behavior:
- Single term with no zero anchor → returns the bare term (unwraps)
- Single `AlgebraicPureTensor` argument → returns the `AlgebraicPureTensor` directly
- Empty result with zero anchor → returns the `AlgebraicZeroTensor`
- Multiple terms → constructs `Basic` subclass instance

**Agent rule:** `AlgebraicTensor(a, b)` may return an `AlgebraicPureTensor`,
`AlgebraicZeroTensor`, or a bare matrix if there is only one non-zero term.
Always check the return type.

### Shape enforcement

`_flatten_args` raises `ShapeMismatchError` if any two terms have different
`tensor_shape`. The `_tensor_shape_of` helper extracts the shape from various
expression types:
- `AlgebraicPureTensor` → `.tensor_shape`
- `AlgebraicZeroTensor` → `.shape`
- `Mul(coeff, AlgebraicPureTensor)` → inner `AlgebraicPureTensor.tensor_shape`
- Bare matrix with `.shape` → wrapped as `((m, n),)`
- Anything else → `None`

### `commutativity_shape`

The `commutativity_shape` property returns the component-wise AND of the
`commutativity_shape` of every term in the sum. The result is a tuple of
binary entries, same length as `tensor_shape`.

`_flatten_args` computes this incrementally as it walks the arguments, using
the `_commutativity_shape_of` helper which handles `AlgebraicPureTensor`,
`AlgebraicZeroTensor`, `AlgebraicTensor`, `Mul(coeff, AlgebraicPureTensor)`,
and bare matrix-like objects.

```
t = AlgebraicTensor(AlgebraicPureTensor(A, C), AlgebraicPureTensor(B, C))
t.commutativity_shape  →  (0, 0)   # A, B, C are all symbolic

t2 = AlgebraicTensor(AlgebraicPureTensor(numeric, C), AlgebraicPureTensor(B, C))
t2.commutativity_shape  →  (0, 0)  # slot 0: 1 & 0 = 0, slot 1: 0 & 0 = 0
```

### The `_sympify` parameter

`__new__` accepts `_sympify=False` to skip sympification of already-processed
args. This is used internally by `proportionality_factoring` when rebuilding
expressions from already-sympified components.

**Agent rule:** Use `_sympify=False` only when you are certain the arguments are
already valid SymPy objects. Use the default `_sympify=True` for public API calls.

### Composition with a single term: `_compose_with_term`

This method composes each term of the sum with a single right operand (bare matrix
or `AlgebraicPureTensor`). It handles:
- `AlgebraicZeroTensor` terms → pass through unchanged
- `AlgebraicPureTensor` terms → compose via `compose_algebraic_pure_tensors`
- `Mul(coeff, AlgebraicPureTensor)` terms → extract coefficient, compose the
  `PureTensor` part, reattach coefficient
- Bare matrix terms → compose directly

---

## 5. Composition

### `compose_algebraic_pure_tensors(left, right)`

Composes two pure tensors by matrix-multiplying corresponding factors:

```
left  = F₁ ⊗ F₂ ⊗ … ⊗ Fₙ
right = G₁ ⊗ G₂ ⊗ … ⊗ Gₙ
result = (F₁*G₁) ⊗ (F₂*G₂) ⊗ … ⊗ (Fₙ*Gₙ)
```

**Requirements:**
- Same number of factors
- Each factor pair `(Fⱼ, Gⱼ)` must have matching inner dimensions:
  `Fⱼ.shape[1] == Gⱼ.shape[0]`
- Result factor `Hⱼ = MatMul(Fⱼ, Gⱼ, evaluate=False)`

Bare matrix objects (with `.shape`) are accepted and treated as single-factor
tensors. Coefficients from both sides are multiplied together.

### `compose_algebraic_tensors(left, right)`

Extends composition to `AlgebraicTensor` by linearity. Dispatch table:

| left type | right type | Behavior |
|---|---|---|
| `AlgebraicZeroTensor` | anything | Returns zero of left's shape |
| anything | `AlgebraicZeroTensor` | Returns zero of right's shape |
| `AlgebraicTensor` | `AlgebraicTensor` | Pairwise compose all term combinations, reassemble |
| `AlgebraicTensor` | `PureTensor`/matrix | Calls `left._compose_with_term(right)` |
| `PureTensor`/matrix | `AlgebraicTensor` | Compose left with each term of right, reassemble |
| `PureTensor` | `PureTensor` | Calls `compose_algebraic_pure_tensors` |
| matrix | matrix | Calls `compose_algebraic_pure_tensors` |

The `_compose_reassemble` helper:
- Strips `AlgebraicZeroTensor` anchors from results
- Handles cancellation (empty real terms → zero tensor)
- Single real term → returns bare term (unwraps)
- Multiple terms → wraps in `AlgebraicTensor(*real, _sympify=False)`

---

## 6. Simplification

### `tensorsimplify(expr)`

Public entry point. Dispatches based on type:
- `AlgebraicZeroTensor` → return unchanged
- `AlgebraicTensor` → `_simplify_algebraic_tensor`
- `AlgebraicPureTensor` → `_simplify_algebraic_pure_tensor`
- Everything else → falls back to SymPy's `simplify`

### `_simplify_algebraic_pure_tensor(pt)`

1. Simplify the leading coefficient with SymPy's `simplify`
2. Simplify each factor individually (useful when factors are `MatrixExpr` with
   simplifiable sub-expressions)
3. Reconstruct the `PureTensor`
4. If coefficient becomes zero → return `AlgebraicZeroTensor`

### `proportionality_factoring(at)`

The core simplification algorithm. Merges proportional `AlgebraicPureTensor`
terms within an `AlgebraicTensor` sum.

#### Algorithm

1. Separate tensor terms from non-tensor terms in the sum
2. Convert each tensor term to `{coeff, factors}` representation via
   `_extract_pt_and_coeff`
3. Iterate with a pivot index:
   a. Compare pivot with every subsequent term
   b. For each pair, walk factor slots and check proportionality via
      `_proportionality_ratio`
   c. **All slots proportional** → merge by summing coefficients
      (if sum is zero, remove both terms)
   d. **Exactly one slot non-proportional** → replace that slot with a
      linear combination of the two factors, weighted by their accumulated
      prefactors
   e. **More than one slot non-proportional** → skip, do not merge
   f. After any successful merge, reset pivot to 0
4. Rebuild terms from entries via `_build_pt`
5. Combine with non-tensor terms and return

#### `_extract_pt_and_coeff(term)`

Returns `(unit_pt_or_none, coeff, [optional_factors])`:
- `AlgebraicPureTensor` → `(pt, pt._get_coeff())`
- `Mul(coeff, AlgebraicPureTensor)` → `(pt, outer_coeff * inner_coeff)`
- `MatMul` with direct matrix factors → `(None, commutative_coeff, [matrix_factors])`
- Anything else → `(term, S.One)`

**Agent rule:** The 3-tuple return `(None, coeff, factors_list)` indicates direct
matrix factors that bypass the `AlgebraicPureTensor` wrapper. This is used for
terms created by `proportionality_factoring` during the "one non-proportional slot"
merge path.

#### `_proportionality_ratio(factor1, factor2)`

Returns commutative `k` such that `factor1 = k * factor2`, or `None`.

For true matrices (shape ≠ (1,1)):
- Checks all non-zero element pairs share the same ratio
- Uses `solve(e1 - k*e2, k)` per element
- Handles zero/non-zero pairs (returns `None` if one is zero and the other isn't)

For 1×1 matrices: delegates to `_matrix_proportionality_ratio` (same element-wise logic).

#### `_build_pt(coeff, factors)`

Reconstructs a tensor term from a coefficient and factor list:
- `coeff == 0` → `AlgebraicZeroTensor`
- Empty factors → return `coeff`
- `coeff == 1`, single factor → bare factor
- `coeff == 1`, multiple factors → `AlgebraicPureTensor(*factors)`
- Otherwise → `AlgebraicPureTensor(coeff, *factors)`

---

## 7. Expansion

### `.expand()`

All three tensor types support `.expand()`, following the standard SymPy pattern.

#### `AlgebraicPureTensor._eval_expand_mul()`

Implemented via the `_eval_expand_mul` hook (called by SymPy's `Expr.expand`
when `mul=True`, which is the default). The algorithm:

1. **Expand each factor** — call `.expand(deep=True, **hints)` on every tensor
   factor. This expands any `MatAdd`, `MatMul`, or nested `Add`/`Mul` inside
   matrix expressions.
2. **Distribute over `Add`** — if any expanded factor is a SymPy `Add`,
   distribute the tensor product by linearity across all addend combinations
   (cartesian product of addends).

```
AlgebraicPureTensor(A, B + C, D).expand()
→  AlgebraicPureTensor(A, B, D) + AlgebraicPureTensor(A, C, D)
```

The commutative coefficient is carried through to every resulting term.

#### `AlgebraicTensor.expand()`

Expands each term in the sum individually and reassembles the result into a
single flattened `AlgebraicTensor`. `AlgebraicZeroTensor` anchors pass through
unchanged.

#### `AlgebraicZeroTensor.expand()`

Returns self unchanged (a zero tensor is already in expanded form).

---

## 8. Implementation Details for Agents

### 7.1 Type checking

When working with these types, use `isinstance` checks. Do NOT rely on duck
typing with `.is_AlgebraicPureTensor` etc. as flags, since the inheritance from
`Mul`/`Basic` means these flags may not propagate in all contexts.

```python
# Correct:
if isinstance(x, AlgebraicPureTensor): ...
if isinstance(x, AlgebraicZeroTensor): ...
if isinstance(x, AlgebraicTensor): ...

# Also check for bare matrix objects:
if hasattr(x, "shape") and x.shape is not None: ...
```

### 7.2 Commutative vs non-commutative dispatch

The `__mul__` / `__rmul__` methods dispatch based on whether the other operand
is commutative:

```python
isinstance(other, Number) or (
    hasattr(other, 'is_commutative') and
    other.is_commutative and
    not isinstance(other, AlgebraicPureTensor)
)
```

**Agent rule:** When adding new arithmetic operators, always check for
`AlgebraicPureTensor` explicitly in the commutative check, since it has
`is_commutative = False` but might match other heuristics.

### 7.3 The `add` dispatcher

```python
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)
```

This is critical for SymPy to route additions through `AlgebraicTensor.__new__`.
Without it, `PureTensor + PureTensor` would go through `Add.__new__` and lose
shape enforcement.

**Agent rule:** If you add a new tensor-like class, register it with the `add`
dispatcher so that additions involving the new type are handled correctly.

### 7.4 Coefficient handling

Coefficients are always commutative (Numbers or commutative Symbols). They are
stored as the first element of `args` in `AlgebraicPureTensor`. When manipulating
terms, always separate the coefficient from the factors:

```python
coeff = pt._get_coeff()
factors = pt.factors  # excludes coefficient
```

### 7.5 Unwrapping behavior

Both `AlgebraicPureTensor.__new__` and `AlgebraicTensor.__new__` unwrap single-term
results. This means:

```python
AlgebraicPureTensor(A, B)          # might return A⊗B PureTensor, or bare A if B is identity
AlgebraicTensor(A⊗B)               # returns A⊗B directly (bare PureTensor)
AlgebraicTensor(A⊗B, C⊗D)         # returns AlgebraicTensor
AlgebraicTensor(A⊗B, -A⊗B)        # returns AlgebraicZeroTensor
```

**Agent rule:** Always be prepared for constructors to return a different type
than the class you called. Write code that handles all possible return types.

### 7.6 Shape comparison

Always compare full `tensor_shape` tuples. Never compare individual factor shapes
unless you are iterating factor-by-factor (as in composition).

### 7.6.1 `commutativity_shape` convention

All three tensor types expose `commutativity_shape`:

- **`AlgebraicPureTensor`** — per-factor check: `0` if the factor contains any
  non-commutative symbol, `1` otherwise.
- **`AlgebraicZeroTensor`** — all-1s tuple (zero is commutative in every slot).
- **`AlgebraicTensor`** — component-wise AND of all term commutativity shapes.

The `_commutativity_shape_of` helper in `algebraic_tensor.py` extracts the
commutativity shape from any recognized expression type (including `Mul` wrappers
and bare matrices). Use this helper when you need to query commutativity from
an arbitrary tensor expression rather than a known instance.

### 7.7 Import patterns

The module avoids circular imports by using lazy imports inside methods. Follow
this pattern:

```python
def my_method(self, other):
    from sympy.tensor.algebraic.algebraic_tensor import compose_algebraic_tensors
    # ...
```

Do NOT add top-level imports between the algebraic submodules unless there is no
circular dependency risk.

### 7.8 MatMul usage

When composing factors, always use `MatMul(lf, rf, evaluate=False)`. The
`evaluate=False` flag prevents SymPy from eagerly evaluating the matrix product,
which preserves the symbolic structure of the expression.

### 7.9 AlgebraicZeroTensor in sums

`AlgebraicZeroTensor` serves as a shape anchor inside `AlgebraicTensor` sums.
When `AlgebraicTensor.__new__` encounters a zero term, it is kept as an anchor
to preserve shape information even when all non-zero terms cancel out. In
`_flatten_args`, the zero term is separated from regular terms and re-appended
at the end during reconstruction.

**Agent rule:** When iterating `AlgebraicTensor.args`, filter out
`AlgebraicZeroTensor` instances for processing, but preserve them when
rebuilding the sum.

### 7.11 Testing considerations

When writing tests, cover:
1. Constructor unwrapping behavior (single term, zero coefficient, zero tensor factor, empty)
2. Shape mismatch detection in addition and composition
3. Composition with different operand type combinations (the full dispatch table)
4. Scalar multiplication vs composition dispatch (commutative vs non-commutative)
5. Negation behavior
6. `proportionality_factoring` with all/one/more-than-one non-proportional slots
7. `AlgebraicZeroTensor` preservation through operations
8. `tensorsimplify` on each type
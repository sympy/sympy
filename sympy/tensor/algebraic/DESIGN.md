# Algebraic Tensor Library — Design & Idea

## Overview

This library implements symbolic tensor products and linear combinations thereof
on top of SymPy. It is built around three core classes:

| Class | Base | Role |
|---|---|---|
| `AlgebraicZeroTensor` | plain Python class (`__slots__`) | Additive identity for a given tensor shape |
| `AlgebraicPureTensor` | `Mul` (non-commutative) | An unevaluated tensor product `A ⊗ B ⊗ C ⊗ …` |
| `AlgebraicTensor` | `Basic` | A sum of same-shape `AlgebraicPureTensor` terms (and optionally a `AlgebraicZeroTensor` anchor) |

The library lives at `sympy.tensor.algebraic.` and is intended to work with
`MatrixExpr` / `MatrixSymbol` factors and non-commutative `Symbol` (wrapped in
1×1 matrices) as tensor factors.

---

## Core Concept: Tensor Shape

Every tensor-carrying object exposes a `tensor_shape` property.

**Tensor shape** is a *tuple of per-factor `(rows, cols)` pairs*. For example:

| Expression | `tensor_shape` |
|---|---|
| `A` (a single 3×4 matrix) | `((3, 4),)` |
| `A ⊗ B` where A is 3×4, B is 4×5 | `((3, 4), (4, 5))` |
| `A ⊗ B ⊗ C` where C is 5×2 | `((3, 4), (4, 5), (5, 2))` |

The shape identifies the tensor *space*. Two tensors can only be added if they
share the exact same `tensor_shape`. This is enforced at construction time by
`AlgebraicTensor.__new__` (raises `ShapeMismatchError` on violation).

Shape normalization accepts several input forms:
- `((3,4), (4,5))` — already a tuple of pairs, unchanged
- `(3, 4)` — bare pair, auto-wrapped to `((3, 4),)`
- `[(3,4)]` — list of pairs, converted to tuple of tuples

---

## Tensor Factors

A **tensor factor** is any object that carries a `.shape` attribute returning
a `(rows, cols)` pair. Typical factors are:

- `MatrixSymbol("A", m, n)` — symbolic matrix
- `Matrix` instances from `sympy.matrices`
- 1×1 wrappers around non-commutative symbols, e.g. a `MatrixSymbol` of size
  1×1 whose entry is `Symbol("X", commutative=False)`

Factors are **non-commutative**: `A ⊗ B ≠ B ⊗ A`. Order is always preserved.

---

## AlgebraicZeroTensor

`AlgebraicZeroTensor(shape)` is the additive identity for the tensor space identified
by `shape`. It is a plain Python class (not a SymPy `Basic` subclass) to keep
it lightweight and avoid inheritance conflicts.

Key behaviors:
- `AlgebraicZeroTensor(s) + AlgebraicZeroTensor(s)` → `AlgebraicZeroTensor(s)`
- `AlgebraicZeroTensor(s) + T` where `T.tensor_shape == s` → returns `T` (or an
  `AlgebraicTensor` with the zero anchored, if there are multiple terms)
- `AlgebraicZeroTensor(s1) + AlgebraicZeroTensor(s2)` with `s1 != s2` → `ShapeMismatchError`
- `-AlgebraicZeroTensor(s)` → `AlgebraicZeroTensor(s)`
- `bool(AlgebraicZeroTensor(s))` → `False`

When all terms in an `AlgebraicTensor` cancel to zero, the result is a
`AlgebraicZeroTensor` of the appropriate shape.

---

## AlgebraicPureTensor

`AlgebraicPureTensor(f1, f2, ..., fn)` represents the non-commutative tensor product
`f1 ⊗ f2 ⊗ … ⊗ fn`.

Extends `Mul` with `evaluate=False` so that SymPy's simplification,
differentiation, and pattern-matching machinery is available. The
`_eval_is_commutative` is `False`.

### Coefficient handling

A `AlgebraicPureTensor` may carry a leading coefficient:

```
AlgebraicPureTensor(alpha, A, B)   →  alpha * (A ⊗ B)
AlgebraicPureTensor(2, A)          →  2 * A
AlgebraicPureTensor(A, B)          →  A ⊗ B  (implicit coefficient 1)
```

The coefficient is always the first arg if it is a `Number` or a commutative
symbol. The `.factors` property returns only the tensor factors (strips the
coefficient). Use `._get_coeff()` to extract the coefficient.

### Special cases in `__new__`

- `AlgebraicPureTensor(0, A, B)` → `AlgebraicZeroTensor(((A.shape), (B.shape),))`
- `AlgebraicPureTensor(A)` with no coefficient → returns `A` directly (unwrap)
- `AlgebraicPureTensor(1, A, B)` → `AlgebraicPureTensor(A, B)`

### Arithmetic

| Operation | Result |
|---|---|
| `AlgebraicPureTensor * scalar` | `AlgebraicPureTensor` with absorbed coefficient |
| `AlgebraicPureTensor * AlgebraicPureTensor` | New `AlgebraicPureTensor` with concatenated factors |
| `AlgebraicPureTensor + AlgebraicPureTensor` | `AlgebraicTensor(self, other)` |
| `-AlgebraicPureTensor` | `AlgebraicPureTensor` with negated coefficient |

Scalar multiplication never produces `MatMul` — the scalar is absorbed as a
coefficient into the `AlgebraicPureTensor`.

---

## AlgebraicTensor

`AlgebraicTensor(T1, T2, ...)` represents `T1 + T2 + …` where every `Ti`
shares the same `tensor_shape`.

Built on `Basic` (not `Add`) to avoid the `is_commutative` descriptor conflict
in `AssocOp._from_args`, while still exposing `is_Add = True` so the wider
SymPy ecosystem recognizes it as an additive expression.

### Construction (`__new__`)

1. Sympify all args (except `AlgebraicZeroTensor`, which is left as-is)
2. If single arg is already a `AlgebraicPureTensor` / `AlgebraicZeroTensor` / `AlgebraicTensor`,
   return it directly
3. Flatten nested `AlgebraicTensor` instances (collect all leaf terms)
4. Validate that all terms share the same `tensor_shape` (raise
   `ShapeMismatchError` if not)
5. Remove `S.Zero` identity terms
6. If exactly one non-zero term remains and no `AlgebraicZeroTensor` was given,
   unwrap to that single term
7. If a `AlgebraicZeroTensor` was provided, keep it as an anchor in the args

### Flattening (`_flatten_args`)

The flattening logic walks through all args:
- `AlgebraicZeroTensor` → validates shape, stored as `zero_term`
- `AlgebraicTensor` → pushes its args onto the work stack
- `AlgebraicPureTensor` → validates shape, added to `terms`
- `Mul(coeff, AlgebraicPureTensor)` → validates shape via `_tensor_shape_of`, added to `terms`
- `Number` → added to `terms` (identity coefficients)
- Anything else with `.shape` → validates shape, added to `terms`

### Properties

- `.tensor_shape` — shape shared by every term (derived from first term with
  a recognizable shape)
- `.terms` — non-zero, non-coefficient terms in the sum
- `.has_zero_term()` — whether a `AlgebraicZeroTensor` anchors this sum

### Arithmetic

| Operation | Result |
|---|---|
| `AlgebraicTensor + X` | `AlgebraicTensor(self, other)` |
| `AlgebraicTensor - X` | `AlgebraicTensor(self, -other)` |
| `AlgebraicTensor * scalar` | Distributes scalar to each term |
| `AlgebraicTensor * AlgebraicPureTensor` | Tensor-product: each term gets the factor appended |
| `-AlgebraicTensor` | Negates each term |

### Factorization helpers

- `.as_common_left()` — extracts common leading factors from all AlgebraicPureTensor terms
- `.as_common_right()` — extracts common trailing factors
- `.as_common_factors()` — extracts both left and right simultaneously

Each returns `(left_factors, middle_expr, right_factors)` where `middle_expr`
is the remaining sub-expression after factoring out the common parts.

---

## Add Dispatcher Registration

At module load time, `AlgebraicTensor` is registered with SymPy's `add`
dispatcher:

```python
add.register_handlerclass((AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor)
```

This ensures that expressions like `A + B` where `A` or `B` involves a
`AlgebraicPureTensor` or `AlgebraicTensor` are routed through `AlgebraicTensor.__new__`
instead of `Add.__new__`.

---

## Simplification (`tensorsimplify`)

The `simplify.py` module provides `tensorsimplify(expr)` and `.simplify()`
methods on all three classes. The pipeline for `AlgebraicTensor` is:

1. **Combine like terms** — group terms sharing the same factor sequence
   (by `id(f)` of each factor), sum their coefficients with SymPy's `Add`,
   simplify the sum, drop groups that cancel to zero
2. **Extract common factors** — use `as_common_factors()` to pull out shared
   left/right factors, leaving a middle sub-sum
3. **Recurse** on the middle sub-sum and reassemble with the common factors

For `AlgebraicPureTensor`, simplification simplifies the coefficient and each factor
individually via SymPy's `simplify`, then reconstructs.

`AlgebraicZeroTensor.simplify()` is a no-op (returns self).

---

## Design Principles

1. **Shape is the contract** — `tensor_shape` must be consistent across all
   terms in an `AlgebraicTensor`. Mismatch raises `ShapeMismatchError`.
2. **Non-commutative by default** — tensor product order is always preserved.
3. **Coefficient absorption** — scalars are absorbed into `AlgebraicPureTensor` as
   leading coefficients; they never produce `MatMul`.
4. **Unwrap when possible** — single-term `AlgebraicTensor` returns the bare
   term; single-factor `AlgebraicPureTensor` with coeff=1 returns the bare matrix.
5. **AlgebraicZeroTensor anchors** — when a `AlgebraicZeroTensor` is explicitly provided or
   all terms cancel, the result carries shape information via `AlgebraicZeroTensor`.
6. **Built on SymPy foundations** — `AlgebraicPureTensor` extends `Mul`, `AlgebraicTensor`
   extends `Basic` with `is_Add=True`, so existing SymPy tools work naturally.

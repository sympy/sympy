# Merge Readiness Checklist: `sympy/tensor/algebraic`

Generated: 2026-07-24

This document tracks every action item from the code review of the `sympy/tensor/algebraic` module. Each item includes the precise code location, reasoning, and a detailed implementation plan.

---

## Table of Contents

1. [Critical Issues (Must Fix)](#1-critical-issues-must-fix)
2. [Style and Convention Violations](#2-style-and-convention-violations)
3. [Missing Core Properties and Methods](#3-missing-core-properties-and-methods)
4. [Test Coverage Gaps](#4-test-coverage-gaps)
5. [Documentation](#5-documentation)
6. [Feature Completeness](#6-feature-completeness)
7. [Code Quality Improvements](#7-code-quality-improvements)
8. [Cleanup](#8-cleanup)

---

## 1. Critical Issues (Must Fix)

### 1.1 [ ] (Missing) Add `is_commutative` class attribute to `AlgebraicPureTensor` and `AlgebraicTensor`

**Location:**
- `algebraic_pure_tensor.py:74` — class definition of `AlgebraicPureTensor`
- `algebraic_tensor.py:294` — class definition of `AlgebraicTensor`

**Current state:** Neither class defines `is_commutative`. The `Basic` base class resolves it through the `_ask()` assumption inference system, which is expensive and can trigger cascades of assumption queries on complex expressions. In contrast, `AlgebraicZeroTensor` correctly sets `is_commutative = True` at `algebraic_zero_tensor.py:44`.

The codebase currently works around this absence with 29 scattered `hasattr(other, 'is_commutative')` guard checks (see `algebraic_pure_tensor.py:139,152,234,285,313,512`, `algebraic_tensor.py:827,841,1097,1134`, `algebraic_zero_tensor.py:131,148`, and throughout `simplify.py`). These are defensive checks that exist because the commutativity flag isn't reliable.

**Reasoning:** In SymPy, `is_commutative` is a foundational property that determines how `Add`, `Mul`, and the dispatcher system treat an object. A tensor product is inherently non-commutative (order of factors matters: `A ⊗ B ≠ B ⊗ A`). Not declaring this forces SymPy to fall back to expensive inference, and worse — the `_ask()` system may not return `False` (it may return `None` for "unknown"), which means expressions mixing algebraic tensors with other SymPy objects may not flatten or simplify correctly.

**Implementation plan:**
```python
# In AlgebraicPureTensor (algebraic_pure_tensor.py, after line 127):
is_commutative = False

# In AlgebraicTensor (algebraic_tensor.py, after line 342):
is_commutative = False
```

After adding these, audit all 29 `hasattr(other, 'is_commutative')` checks. Many can be simplified to direct attribute access since the flag will always be defined. For example, `algebraic_pure_tensor.py:285-287`:

```python
# Before:
if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
        other.is_commutative and not (...)):

# After (once our own classes have is_commutative):
if (hasattr(other, 'is_commutative') and other.is_commutative
        and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor))):
```

The `hasattr` guard on `other` should remain (since `other` is an external object), but the logic becomes more reliable once our own classes have the attribute.

---

### 1.2 [ ] (Missing) Add `free_symbols` property to `AlgebraicPureTensor` and `AlgebraicTensor`

**Location:**
- `algebraic_pure_tensor.py` — class `AlgebraicPureTensor`
- `algebraic_tensor.py` — class `AlgebraicTensor`
- `algebraic_zero_tensor.py:77-78` — `AlgebraicZeroTensor` already returns `set()` (correct)

**Current state:** `AlgebraicPureTensor` and `AlgebraicTensor` do not define `free_symbols`. The `Basic` base class has a default `free_symbols` that iterates over `args` and collects `free_symbols` from each arg. For `AlgebraicPureTensor`, `args` includes the coefficient and factors, so the default *may* work. However, this depends on `Basic`'s iteration behavior and may miss symbols if the `args` structure is unusual (e.g., when a single-factor tensor unwraps to a bare matrix).

**Reasoning:** `free_symbols` is used by `diff()`, `has()`, `subs()`, `replace()`, and many other SymPy operations. If it's wrong or missing, these operations silently fail or produce incorrect results. The property must explicitly collect symbols from both the coefficient and all tensor factors.

**Implementation plan:**

```python
# In AlgebraicPureTensor (add as a property after the existing properties):
@property
def free_symbols(self):
    """Union of free symbols from the coefficient and all factors."""
    syms = set()
    if self.coeff.free_symbols:
        syms |= self.coeff.free_symbols
    for f in self.factors:
        if hasattr(f, 'free_symbols'):
            syms |= f.free_symbols
    return syms

# In AlgebraicTensor (add as a property after the existing properties):
@property
def free_symbols(self):
    """Union of free symbols from all terms in this sum."""
    syms = set()
    for a in self.args:
        if hasattr(a, 'free_symbols'):
            syms |= a.free_symbols
    return syms
```

**Verification:** After implementing, test with:
```python
from sympy.abc import x, y
A = MatrixSymbol("A", 3, 4)
T = AlgebraicPureTensor(x**2, A, B)
assert T.free_symbols == {x}  # A has no free_symbols (MatrixSymbol doesn't contribute)
T2 = AlgebraicPureTensor(x, M, N)  # M contains x
assert x in T2.free_symbols
```

---

### 1.3 [ ] (Missing) Verify and test `subs()` behavior

**Location:** Both `AlgebraicPureTensor` and `AlgebraicTensor` inherit `subs()` from `Basic`.

**Current state:** No custom `_eval_subs` or `subs` override. `Basic.subs()` iterates over `args` and calls `subs` on each. This should work for coefficient substitution, but the behavior with pattern matching (using `Wild`) or with complex substitution rules is untested.

**Reasoning:** Users will naturally write `T.subs(x, 2)` to substitute a symbol in a tensor coefficient or factor. This must work correctly. If `Basic.subs()` doesn't recurse into factors properly (e.g., into `ImmutableDenseMatrix` entries), symbols inside matrix entries won't be substituted.

**Implementation plan:**
1. First, test the default `Basic.subs()` behavior:
   - `T.subs(x, 2)` where `x` is in the coefficient
   - `T.subs(x, 2)` where `x` is inside a concrete matrix entry
   - `T.subs(A, C)` where `A` is a `MatrixSymbol` factor
   - `T.subs({x: 2, A: C})` with multiple substitutions

2. If any case fails, implement `_eval_subs`:
```python
def _eval_subs(self, old, new):
    # Substitute in coefficient
    new_coeff = self.coeff.subs(old, new)
    # Substitute in each factor
    new_factors = [f.subs(old, new) for f in self.factors]
    # Rebuild
    if new_coeff is S.One:
        return AlgebraicPureTensor(*new_factors)
    return AlgebraicPureTensor(new_coeff, *new_factors)
```

3. For `AlgebraicTensor`, apply the same pattern to each term in `self.args`.

---

## 2. Style and Convention Violations

### 2.1 [ ] (Missing) Remove `display()` methods from all three classes

**Location:**
- `algebraic_pure_tensor.py:589-607` — `AlgebraicPureTensor.display()`
- `algebraic_tensor.py:946-964` — `AlgebraicTensor.display()`
- `algebraic_zero_tensor.py:188-205` — `AlgebraicZeroTensor.display()`

**Current state:** All three classes have a `display(mode="latex")` method that imports `IPython.display.display` and `IPython.display.Latex`. If IPython is not installed, it falls back to `print()`.

**Reasoning:** SymPy does not include IPython-dependent methods on core classes. SymPy's printing system (`str`, `repr`, `latex`, `pretty`) is the standard interface, and it works through the `printing/` module handlers (which are already correctly implemented for these classes). The `display()` method:
- Introduces a hard dependency check on IPython in a core SymPy module
- Is redundant — users can call `display(Latex(latex(expr)))` themselves in a notebook
- Violates SymPy's separation of concerns (printing is handled by the `printing/` subsystem, not by the classes themselves)
- The docstring uses `Parameters` with `----------` underlines (NumPy style), not the `==========` underlines SymPy requires

**Implementation plan:** Delete the `display()` method from all three classes. No replacement needed — SymPy's existing printers handle all output formats.

---

### 2.2 [ ] (Missing) Remove `tests_old/` directory

**Location:** `tensor/algebraic/tests_old/`

**Current contents:**
- `playground_functions.py` — 607-line interactive script with physics examples (Dirac matrices, representation theory)
- `tempCodeRunnerFile.py` — VS Code temp file (development artifact)
- `Untitled.ipynb` — Jupyter notebook with LaTeX display experiments
- `.ipynb_checkpoints/` — Notebook checkpoint directory
- `__pycache__/` — Python cache
- `.DS_Store` — macOS metadata

**Reasoning:** This directory contains development artifacts that should not be in the repository. A PR with `tests_old/` containing temp files and notebooks would be rejected by reviewers. The `playground_functions.py` contains valuable example code that should be extracted as proper doctests or documentation examples.

**Implementation plan:**
1. Review `playground_functions.py` for useful examples
2. Convert useful examples to doctest examples in the relevant class docstrings or to proper test functions in the `tests/` directory
3. Delete the entire `tests_old/` directory

---

### 2.3 [ ] (Missing) Remove `.DS_Store` files

**Location:**
- `tensor/algebraic/.DS_Store`
- `tensor/algebraic/tests/.DS_Store`
- `tensor/algebraic/tests_old/.DS_Store`

**Reasoning:** `.DS_Store` is a macOS metadata file that should never be committed. These should be in `.gitignore`.

---

### 2.4 [ ] (Missing) Fix docstring parameter section format in `AlgebraicZeroTensor.display()`, `AlgebraicPureTensor.display()`, `AlgebraicTensor.display()`

**Location:** The `display()` method docstrings in all three classes.

**Current state:** Uses `Parameters` with `----------` (dashes) underlines.

**SymPy convention:** Per the docstring style guide, section headings must use `==========` (equals signs) underlines. The `display()` methods will be removed (see 2.1), so this is noted for completeness only.

---

## 3. Missing Core Properties and Methods

### 3.1 [ ] (Done — verified) `__hash__` correctness

**Location:**
- `algebraic_zero_tensor.py:67-68` — `_hashable_content()` defined
- `algebraic_pure_tensor.py` — inherits `Basic.__hash__`
- `algebraic_tensor.py` — inherits `Basic.__hash__`

**Current state:** `AlgebraicZeroTensor` correctly defines `_hashable_content()` returning `(self._shape,)`. The other two classes inherit from `Basic`, which uses `self.args` for hashing. Since `AlgebraicPureTensor.args` contains the coefficient and factors (which are themselves hashable SymPy objects), and `AlgebraicTensor.args` contains the terms, the inherited behavior is correct.

**Verification needed:** Ensure that two equal `AlgebraicPureTensor` objects have the same hash:
```python
T1 = AlgebraicPureTensor(A, B)
T2 = AlgebraicPureTensor(A, B)
assert hash(T1) == hash(T2) and T1 == T2
```

---

### 3.2 [ ] (Missing) Add `__eq__` / equality semantics documentation

**Location:** `algebraic_pure_tensor.py` class definition.

**Current state:** Equality is inherited from `Basic`, which compares `args` tuples. This means two `AlgebraicPureTensor` objects are equal if and only if they have identical args (same coefficient, same factors in same order).

**Reasoning:** This is the correct behavior for a non-commutative tensor product. However, it should be documented in the class docstring that `A ⊗ B == B ⊗ A` is `False` (order matters), to prevent user confusion.

**Implementation plan:** Add to the `AlgebraicPureTensor` class docstring:
```
Equality is structural: two pure tensors are equal if and only if they
have the same coefficient and identical factors in the same order.
Because the tensor product is non-commutative, ``A ⊗ B`` and ``B ⊗ A``
are not equal (assuming ``A`` and ``B`` have different shapes).
```

---

### 3.3 [ ] (Missing) Consider adding `evalf()` support

**Location:** Neither class defines `evalf()` or `_eval_evalf()`.

**Current state:** `Basic.evalf()` will attempt to call `evalf` on each `arg`. For `AlgebraicPureTensor`, this would call `evalf` on the coefficient and each factor. If a factor is a concrete `ImmutableDenseMatrix`, `evalf` will evaluate its entries. This should work by inheritance.

**Reasoning:** Users may want to numerically evaluate tensors with symbolic coefficients and concrete matrix factors. The inherited behavior likely works, but it should be tested and documented.

**Implementation plan:**
1. Test inherited `evalf()`:
```python
T = AlgebraicPureTensor(pi, M, N)  # M is concrete matrix, N is MatrixSymbol
result = T.evalf()
# Expected: coefficient is Float, matrix entries are Float, N stays symbolic
```
2. If it works, add a docstring note and a test. If not, implement `_eval_evalf`.

---

### 3.4 [ ] (Missing) Add `as_coeff_mul` and `as_coeff_Add` support

**Location:** `algebraic_pure_tensor.py`, `algebraic_tensor.py`.

**Current state:** Not implemented.

**Reasoning:** `as_coeff_mul()` is used by SymPy's rewriting system and simplification to separate commutative coefficients from non-commutative parts. For `AlgebraicPureTensor`, the natural implementation would be `(self.coeff, self.factors)`. For `AlgebraicTensor`, `as_coeff_Add()` doesn't apply (it's already the Add form), but `as_coeff_mul()` could return `(S.One, (self,))`.

**Implementation plan:**
```python
# In AlgebraicPureTensor:
def as_coeff_mul(self):
    """Return (coeff, (factor1, factor2, ...))."""
    return (self.coeff, self.factors)

def as_coeff_Mul(self):
    """Return (coeff, rest) where rest is the tensor product of factors."""
    if len(self.factors) == 1:
        return (self.coeff, self.factors[0])
    if self.coeff is S.One:
        return (S.One, self)
    return (self.coeff, AlgebraicPureTensor(*self.factors))
```

---

## 4. Test Coverage Gaps

### 4.1 [ ] (Missing) Fix tautological assertion in `test_algebraic_tensor.py`

**Location:** `tests/test_algebraic_tensor.py:278`

**Current code:**
```python
assert S2_doit == AlgebraicTensor(T3, T4)
```

**Problem:** `S2_doit` IS `AlgebraicTensor(T3, T4)` (it was constructed that way just before). This asserts identity, not that `doit()` performed any transformation.

**Implementation plan:** Replace with a test that verifies `doit()` actually evaluates sub-expressions. For example:
```python
from sympy import UnevaluatedExpr
T3 = AlgebraicPureTensor(UnevaluatedExpr(2 + 3), A, B)
T4 = AlgebraicPureTensor(UnevaluatedExpr(1 * 4), C, D)
S2 = AlgebraicTensor(T3, T4)
S2_doit = S2.doit()
assert S2_doit.coeff == 5  # UnevaluatedExpr(2+3) was evaluated to 5
```

Or, if `doit()` on the coefficient works through the `AlgebraicTensor.__new__` coefficient merging:
```python
T3 = AlgebraicPureTensor(x + y, A, B)
T4 = AlgebraicPureTensor(x, A, B)
S2 = AlgebraicTensor(T3, T4)
assert S2.doit().coeff == 2*x + y  # Coefficients should be merged
```

---

### 4.2 [ ] (Missing) Add tests for `free_symbols`

**Location:** New tests in `tests/test_algebraic_pure_tensor.py` and `tests/test_algebraic_tensor.py`.

**Implementation plan:**
```python
def test_free_symbols_pure_tensor():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    x, y = symbols('x y')
    T = AlgebraicPureTensor(x**2, A, B)
    assert x in T.free_symbols
    assert y not in T.free_symbols

def test_free_symbols_with_matrix_entries():
    x = symbols('x')
    M = ImmutableDenseMatrix([[x, 1], [2, x**2]])
    N = MatrixSymbol("N", 2, 3)
    T = AlgebraicPureTensor(M, N)
    assert x in T.free_symbols

def test_free_symbols_tensor_sum():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    x, y = symbols('x y')
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(y, C, D))
    assert S.free_symbols == {x, y}

def test_free_symbols_zero_tensor():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.free_symbols == set()
```

---

### 4.3 [ ] (Missing) Add tests for `subs()` behavior

**Location:** New tests in `tests/test_algebraic_pure_tensor.py`.

**Implementation plan:**
```python
def test_subs_coefficient():
    x = symbols('x')
    T = AlgebraicPureTensor(x**2, A, B)
    result = T.subs(x, 3)
    assert result.coeff == 9

def test_subs_matrix_factor():
    C = MatrixSymbol("C", 3, 4)
    T = AlgebraicPureTensor(A, B)
    result = T.subs(A, C)
    assert result.factors[0] == C

def test_subs_concrete_matrix_entries():
    x = symbols('x')
    M = ImmutableDenseMatrix([[x, 1], [2, x]])
    N = MatrixSymbol("N", 2, 3)
    T = AlgebraicPureTensor(M, N)
    result = T.subs(x, 5)
    assert result.factors[0][0, 0] == 5
```

---

### 4.4 [ ] (Missing) Add tests for 3+ factor tensors

**Location:** New tests in `tests/test_algebraic_pure_tensor.py`.

**Current state:** All existing tests use 2-factor tensors. The code handles arbitrary numbers of factors, but 3+ factor tensors are untested.

**Implementation plan:**
```python
def test_three_factor_pure_tensor():
    A = MatrixSymbol("A", 2, 3)
    B = MatrixSymbol("B", 3, 4)
    C = MatrixSymbol("C", 4, 5)
    T = AlgebraicPureTensor(A, B, C)
    assert T.num_factors == 3
    assert T.shape == ((2, 3), (3, 4), (4, 5))

def test_compose_three_factor_tensors():
    A = MatrixSymbol("A", 2, 3)
    B = MatrixSymbol("B", 3, 4)
    C = MatrixSymbol("C", 4, 5)
    D = MatrixSymbol("D", 3, 2)
    E = MatrixSymbol("E", 4, 3)
    F = MatrixSymbol("F", 5, 4)
    T1 = AlgebraicPureTensor(A, B, C)
    T2 = AlgebraicPureTensor(D, E, F)
    result = compose_algebraic_pure_tensors(T1, T2)
    assert result.num_factors == 3

def test_expand_three_factors():
    A = MatrixSymbol("A", 2, 3)
    B1 = MatrixSymbol("B1", 3, 4)
    B2 = MatrixSymbol("B2", 3, 4)
    C = MatrixSymbol("C", 4, 5)
    T = AlgebraicPureTensor(A, MatAdd(B1, B2), C)
    expanded = T.expand()
    assert isinstance(expanded, AlgebraicTensor)
    assert len(expanded.args) == 2
```

---

### 4.5 [ ] (Missing) Add tests for `expand(deep=False)`

**Location:** New tests in `tests/test_algebraic_pure_tensor.py`.

**Implementation plan:**
```python
def test_expand_deep_false():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 4, 5)
    T = AlgebraicPureTensor(A, MatAdd(B, C))
    # With deep=False, the MatAdd should NOT be distributed
    result = T.expand(deep=False)
    assert result is T  # Should return self unchanged
```

---

### 4.6 [ ] (Missing) Add test for `add.register_handlerclass` behavior

**Location:** New test in `tests/test_algebraic_tensor.py`.

**Current code:** `algebraic_tensor.py:1193-1195`:
```python
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)
```

**Reasoning:** This registration routes `Add` construction through `AlgebraicTensor.__new__` when any arg is an algebraic tensor type. Without a test, this could break silently if the registration is removed or if `Add`'s dispatcher changes.

**Implementation plan:**
```python
def test_add_dispatcher():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    # Using Add directly should route through AlgebraicTensor
    from sympy.core.add import Add
    result = Add(T1, T2)
    assert isinstance(result, AlgebraicTensor)
    assert not isinstance(result, Add)
```

---

### 4.7 [ ] (Missing) Add tests for top-level `simplify()` interaction

**Location:** New test in `tests/test_simplify.py`.

**Current state:** Only `tensorsimplify()` is tested. Users will naturally try `from sympy import simplify; simplify(tensor_expr)`.

**Reasoning:** SymPy's top-level `simplify()` calls `_simplify` on the expression's args. For `AlgebraicTensor`, this may or may not recurse into the tensor-specific simplification. If it doesn't, users get suboptimal results.

**Implementation plan:**
```python
def test_sympy_simplify_interaction():
    from sympy import simplify
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    T = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(2, A, B))
    # This should at least not crash and should simplify the coefficient
    result = simplify(T)
    assert result.coeff == 3  # If simplify recurses into our types correctly
```

If the above fails (which is likely), document in the module that users should use `tensorsimplify()` for algebraic tensor expressions, and consider registering a handler in `sympy.simplify.simplify`.

---

### 4.8 [ ] (Missing) Strengthen simplify test assertions

**Location:** Multiple tests in `tests/test_simplify.py`.

**Current state:** Many tests use `isinstance(result, AlgebraicPureTensor)` without verifying the result is actually simpler.

**Reasoning:** An `isinstance` check passes even if the function returns the input unchanged. Tests should verify that simplification actually improved the expression.

**Implementation plan:** Replace weak assertions with specific checks. For example:
```python
# Before:
result = tensorsimplify(expr)
assert isinstance(result, AlgebraicPureTensor)

# After:
result = tensorsimplify(expr)
assert isinstance(result, AlgebraicPureTensor)
assert result.coeff == expected_coeff  # Verify the coefficient was simplified
assert len(result.factors) == expected_factors  # Verify factor count
```

---

### 4.9 [ ] (Missing) Add tests for `_compose_with_term` method

**Location:** New test in `tests/test_algebraic_tensor.py`.

**Current code:** `algebraic_tensor.py:749-788` — `AlgebraicTensor._compose_with_term(other)`.

**Reasoning:** This is a private method used by `AlgebraicTensor.__mul__` when composing with a `PureTensor`. It has no direct test.

**Implementation plan:**
```python
def test_compose_with_term():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    E = MatrixSymbol("E", 4, 2)
    F = MatrixSymbol("F", 5, 3)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    T = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T)
    # Should produce the same as compose_algebraic_tensors(S, T)
    expected = compose_algebraic_tensors(S, T)
    assert result == expected
```

---

### 4.10 [ ] (Missing) Add test for `is_Add = True` behavior

**Location:** New test in `tests/test_algebraic_tensor.py`.

**Current code:** `algebraic_tensor.py:343` sets `is_Add = True`.

**Reasoning:** This flag allows SymPy functions to recognize `AlgebraicTensor` as an additive expression. It should work with `Add.flatten()`, `as_add_terms()`, and similar utilities.

**Implementation plan:**
```python
def test_is_add_flag():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    assert S.is_Add is True
    # Verify it works with Add operations
    assert isinstance(S, Basic)
```

---

## 5. Documentation

### 5.1 [ ] (Done — skeleton exists) Sphinx documentation file exists

**Location:** `doc/src/modules/tensor/algebraic.rst` (32 lines).

**Current state:** File exists with `automodule`, `autoclass`, and `autofunction` directives for all public API. Listed in `doc/src/modules/tensor/index.rst:15`.

---

### 5.2 [ ] (Missing) Expand Sphinx documentation with tutorial and examples

**Location:** `doc/src/modules/tensor/algebraic.rst`.

**Current state:** The file is a minimal skeleton (32 lines) with only `autoclass`/`autofunction` directives. It has no introduction, no mathematical background, no usage examples, and no "when to use this" guidance.

**Reasoning:** A new module in SymPy should have documentation that helps users understand:
- What problem algebraic tensors solve
- How they differ from `Indexed`, `tensorproduct`, `Array`, and the `Tensor` class
- A realistic worked example
- The mathematical notation used

**Implementation plan:** Restructure `algebraic.rst` to include:

```rst
.. _tensor-algebraic:

Algebraic Tensor Products
=========================

Overview
--------
Explain what algebraic tensors are, how they model tensor products of
matrix-like objects, and how they differ from other tensor representations
in SymPy (Indexed, Array, Tensor).

Mathematical Background
-----------------------
Brief description of the tensor product notation, composition as
factor-wise matrix multiplication, and the role of commutativity patterns.

Quick Example
-------------
A worked example showing creation, arithmetic, composition, and simplification.

.. jupyter-execute::

   from sympy.matrices.expressions import MatrixSymbol
   from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor
   from sympy.tensor.algebraic.simplify import tensorsimplify

   A = MatrixSymbol("A", 3, 4)
   B = MatrixSymbol("B", 4, 5)
   T = AlgebraicPureTensor(A, B)
   print(T)

Classes
-------
... (existing autoclass directives)

Functions
---------
... (existing autofunction directives)
```

---

### 5.3 [ ] (Missing) Add equality semantics to `AlgebraicPureTensor` docstring

**Location:** `algebraic_pure_tensor.py:74` — class docstring.

See item 3.2 above.

---

### 5.4 [ ] (Missing) Verify all docstring examples pass doctest

**Location:** All files with docstrings.

**Reasoning:** Per SymPy's docstring style guide, every example in a docstring must pass `./bin/doctest`. Some examples in the current code use `print()` (e.g., `algebraic_pure_tensor.py:101`), which produces output that doctest verifies. Other examples use `>>> expr` and rely on the repr, which may change.

**Implementation plan:**
1. Run `./bin/doctest tensor/algebraic` to verify all examples pass
2. Fix any failing examples
3. Ensure examples use `from sympy import ...` style imports (per SymPy convention), not `from sympy.tensor.algebraic import ...` in the outer scope — though this is a gray area for module-local examples

---

## 6. Feature Completeness

### 6.1 [ ] (Missing) Register handler for top-level `simplify()`

**Location:** `sympy/simplify/simplify.py`.

**Current state:** Only `tensorsimplify()` works. Calling `simplify(tensor_expr)` will use SymPy's general simplification, which won't apply tensor-aware rules.

**Reasoning:** Users expect `simplify()` to work on any SymPy expression. If `simplify()` doesn't delegate to `tensorsimplify()` for algebraic tensor types, users will get suboptimal results and may file bug reports.

**Implementation plan:** Option A (preferred) — Add a handler in `sympy/simplify/simplify.py`:
```python
# In the _simplify_functions or similar dispatch mechanism:
def _simplify_algebraic_tensor(expr):
    from sympy.tensor.algebraic.simplify import tensorsimplify
    return tensorsimplify(expr)
```

Option B — Document in the module docstring that `tensorsimplify()` should be used instead of `simplify()` for algebraic tensor expressions.

---

### 6.2 [ ] (Missing) Document scope: tensor contraction is not supported

**Location:** `algebraic/__init__.py` module docstring.

**Reasoning:** A user familiar with tensors will expect contraction (summing over paired indices) as a basic operation. This module does not support contraction. This should be explicitly documented to set correct expectations.

**Implementation plan:** Add a "Limitations" section to the module docstring:
```
Limitations
===========

This module supports tensor products, linear combinations, and
factor-wise composition (matrix multiplication) of tensor factors.
The following operations are not yet supported:

- Tensor contraction (summing over paired indices)
- Tensor transposition that reorders factors (only factor-wise .T)
- Inner products between tensor spaces
```

---

### 6.3 [ ] (Missing) Consider adding `equals()` for numerical comparison

**Location:** `algebraic_pure_tensor.py`, `algebraic_tensor.py`.

**Current state:** Not implemented. Relies on `Basic.equals()`.

**Reasoning:** `Basic.equals()` tries random numerical substitution to check equality. For tensors, this may not work well because tensor expressions can't always be substituted with random numbers (factors are matrices, not scalars). It's worth testing whether `Basic.equals()` works for these types, and if not, consider a custom implementation.

**Implementation plan:** Test `Basic.equals()` on tensor expressions first. If it fails, document the limitation.

---

## 7. Code Quality Improvements

### 7.1 [ ] (Missing) Simplify `hasattr(other, 'is_commutative')` checks after adding `is_commutative`

**Location:** 29 occurrences across all module files (see item 1.1 for the full list).

**Current state:** Every multiplication check uses the pattern:
```python
hasattr(other, 'is_commutative') and other.is_commutative and not (...)
```

**Reasoning:** After adding `is_commutative = False` to our own classes (item 1.1), the `not (hasattr(other, 'is_AlgebraicPureTensor') or ...)` guards become less necessary since our classes will have `is_commutative = False` and won't pass the `other.is_commutative` check. However, the guards also exclude `AlgebraicZeroTensor` (which has `is_commutative = True`), so the `not` clause is still needed to prevent zero tensors from being treated as scalar multipliers.

**Implementation plan:** After fixing item 1.1, simplify the checks where possible. The exclusion of `AlgebraicZeroTensor` must remain since it has `is_commutative = True`.

---

### 7.2 [ ] (Missing) Reduce code duplication in `compose_algebraic_pure_tensors` zero-tensor handling

**Location:** `algebraic_pure_tensor.py:723-793` and `algebraic_tensor.py:169-242`.

**Current state:** The zero-tensor shortcut logic in `compose_algebraic_pure_tensors` (lines 723-793) duplicates similar logic in `compose_algebraic_tensors` (lines 169-242). Both handle the `(Zero, Zero)`, `(Zero, non-zero)`, and `(non-zero, Zero)` cases with nearly identical shape computation.

**Reasoning:** The duplication is understandable (the functions need to be self-contained), but the shape computation logic is identical in both. Extract a helper function to compute the composed shape from two shapes.

**Implementation plan:**
```python
def _compute_composed_shape(left_shape, right_shape):
    """Compute the resulting shape from composing two tensor shapes.

    Raises ShapeMismatchError if shapes are incompatible.
    """
    if len(left_shape) != len(right_shape):
        raise ShapeMismatchError(...)
    for i, (l, r) in enumerate(zip(left_shape, right_shape)):
        if l[1] != r[0]:
            raise ShapeMismatchError(...)
    return tuple((l[0], r[1]) for l, r in zip(left_shape, right_shape))
```

Then replace the duplicated shape computation in both functions with calls to this helper.

---

### 7.3 [ ] (Missing) Single-factor unwrap behavior — document or make configurable

**Location:** `algebraic_pure_tensor.py:254-259`.

**Current code:**
```python
# Single factor with coefficient 1: unwrap to bare factor
if len(processed) == 1 and coeff is S.One:
    return processed[0]

# Single factor with non-trivial coefficient: return coeff * factor
if len(processed) == 1:
    return coeff * processed[0]
```

**Reasoning:** When a user constructs `AlgebraicPureTensor(A)`, they get back `A` (a `MatrixSymbol`), not an `AlgebraicPureTensor`. This is a convenience that breaks type expectations:
- `isinstance(AlgebraicPureTensor(A), AlgebraicPureTensor)` is `False`
- The return type of `AlgebraicPureTensor.__new__` is not `AlgebraicPureTensor`

This is consistent with how SymPy's `Mul` and `Add` unwrap single-argument expressions. But it should be clearly documented. Also, consider whether this is the right default behavior — a user constructing a pure tensor may explicitly want a tensor object, not the bare factor.

**Implementation plan:** At minimum, add to the docstring:
```
Note that a single-factor tensor unwraps to the bare factor:
``AlgebraicPureTensor(A)`` returns ``A`` (a ``MatrixSymbol``), not an
``AlgebraicPureTensor`` instance. Use ``AlgebraicPureTensor(1, A)``
to force tensor wrapping.
```

---

### 7.4 [ ] (Missing) `_is_zero_like` function robustness

**Location:** `algebraic_pure_tensor.py:56-71`.

**Current code:**
```python
def _is_zero_like(expr):
    if expr is S.Zero:
        return True
    if getattr(expr, '_assumptions', None) and expr._assumptions.get('zero') is True:
        return True
    return expr == S.Zero * expr
```

**Reasoning:** The `expr == S.Zero * expr` fallback is clever but has subtle issues:
- For a `MatrixSymbol`, `S.Zero * A` produces `ZeroMatrix` times `A`, which may not equal the original `A` comparison cleanly
- For complex expressions, `S.Zero * expr` triggers `Mul.__new__` which may have side effects
- The comment acknowledges this: "avoiding the descriptor to prevent triggering SymPy's assumption inference"

The `_assumptions` direct access is a SymPy internal API and may break between versions.

**Implementation plan:** This is acceptable as-is but should have a test for edge cases:
```python
def test_is_zero_like():
    from sympy.tensor.algebraic.algebraic_pure_tensor import _is_zero_like
    assert _is_zero_like(S.Zero) is True
    assert _is_zero_like(0) is True
    ZM = ZeroMatrix(3, 4)
    assert _is_zero_like(ZM) is True
    A = MatrixSymbol("A", 3, 4)
    assert _is_zero_like(A) is False
```

---

### 7.5 [ ] (Missing) `conjugate()` fallback to `Conjugate(...)` wrapper

**Location:** `algebraic_pure_tensor.py:428-432` and `algebraic_tensor.py:724-728`.

**Current code:**
```python
def conjugate(self):
    result = self._eval_conjugate()
    if result is not None:
        return result
    from sympy.functions.elementary.complexes import conjugate as c
    return c(self)
```

**Reasoning:** `_eval_conjugate()` is implemented and returns a valid result for all current tensor types. The fallback to `c(self)` would wrap the tensor in `Conjugate(AlgebraicPureTensor(...))`, creating an unevaluated wrapper. This fallback path is dead code if `_eval_conjugate()` always succeeds. If it can fail (return `None`), document when and why.

Looking at `_eval_conjugate()` in both classes, it always returns a result (either a new tensor or `self`). The fallback is unreachable. Either remove it or document it.

**Implementation plan:** Remove the fallback and simplify:
```python
def conjugate(self):
    """Return the complex conjugate of this tensor."""
    return self._eval_conjugate()
```

---

### 7.6 [ ] (Missing) `_coeff_map` in `AlgebraicTensor` — verify it stays in sync

**Location:** `algebraic_tensor.py:340` — `__slots__ = ('_coeff_map',)`, set at `algebraic_tensor.py:388`.

**Current state:** `_coeff_map` is computed in `_collect_coefficients()` and stored on the object. It's used by `_combine_coeff_maps()` for efficient merge/subtract operations. It's excluded from pickle via `__getstate__` (`algebraic_tensor.py:567-572`).

**Reasoning:** The `_coeff_map` is a performance optimization for the `_merge` and `_subtract` fast paths (used when both operands have "simple terms"). If the map ever gets out of sync with `self.args`, the fast path produces incorrect results. The slow path (via `AlgebraicTensor.__new__`) is always correct since it recomputes from scratch.

**Implementation plan:** Add an assertion or invariant check in test mode:
```python
def test_coeff_map_consistency():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    S = AlgebraicTensor(
        AlgebraicPureTensor(2, A, B),
        AlgebraicPureTensor(3, A, B),
        AlgebraicPureTensor(C, D)
    )
    # The coeff_map should reflect the merged coefficients
    assert S._coeff_map is not None
    # After merging, A⊗B should have coefficient 5
    key = AlgebraicPureTensor(A, B)
    assert S._coeff_map.get(key) == 5
```

---

## 8. Cleanup

### 8.1 [ ] (Missing) Delete `tests_old/` directory

See item 2.2 above.

---

### 8.2 [ ] (Missing) Add `tensor/algebraic/` to `.gitignore` for `.DS_Store` and `__pycache__`

**Location:** Root `.gitignore`.

**Implementation plan:** Ensure the global `.gitignore` covers:
```
.DS_Store
__pycache__/
*.pyc
.ipynb_checkpoints/
```

These should already be in the SymPy `.gitignore`, but the files are currently present in the working tree, suggesting they were committed. Run `git rm --cached` on each.

---

### 8.3 [ ] (Missing) Remove `.DS_Store` from all directories

**Location:**
- `tensor/algebraic/.DS_Store`
- `tensor/algebraic/tests/.DS_Store`
- `tensor/algebraic/tests_old/.DS_Store`
- `tensor/.DS_Store`

**Implementation plan:**
```bash
git rm --cached tensor/algebraic/.DS_Store
git rm --cached tensor/algebraic/tests/.DS_Store
git rm --cached tensor/algebraic/tests_old/.DS_Store
git rm --cached tensor/.DS_Store
```

---

## Summary of Priority

| Priority | Items | Count |
|----------|-------|-------|
| **Critical** (block merge) | 1.1 `is_commutative`, 1.2 `free_symbols`, 1.3 `subs()`, 2.1 `display()` removal, 2.2 `tests_old/` cleanup | 5 |
| **High** (fix before PR) | 4.1 tautological test, 4.2-4.10 test gaps, 5.2 Sphinx docs, 8.1-8.3 cleanup | 12 |
| **Medium** (include if time) | 3.3 `evalf()`, 3.4 `as_coeff_mul`, 6.1 `simplify()` handler, 6.2 scope doc, 7.2-7.6 code quality | 7 |
| **Low** (future work) | 6.3 `equals()`, tensor contraction feature | 2 |

**Total actionable items: 33**

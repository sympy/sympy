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
6. [Code Quality Improvements](#6-code-quality-improvements)
7. [Cleanup](#7-cleanup)

---

## 1. Critical Issues (Must Fix)

### 1.1 [x] (Done) Add `is_commutative` class attribute to `AlgebraicPureTensor` and `AlgebraicTensor`

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

### 1.2 [x] (Done) Add `free_symbols` property to `AlgebraicPureTensor` and `AlgebraicTensor`

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

### 1.3 [x] (Done) Verify and test `subs()` behavior

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

### 2.1 [x] (Done) Remove `display()` methods from all three classes

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

### 2.2 [x] (Done) Remove `.DS_Store` files

**Location:**
- `tensor/algebraic/.DS_Store`
- `tensor/algebraic/tests/.DS_Store`
- `tensor/algebraic/tests_old/.DS_Store`

**Reasoning:** `.DS_Store` is a macOS metadata file that should never be committed. These should be in `.gitignore`.

---

### 2.3 [x] (Done — superseded by 2.1) Fix docstring parameter section format in `AlgebraicZeroTensor.display()`, `AlgebraicPureTensor.display()`, `AlgebraicTensor.display()`

**Location:** The `display()` method docstrings in all three classes.

**Current state:** Uses `Parameters` with `----------` (dashes) underlines.

**SymPy convention:** Per the docstring style guide, section headings must use `==========` (equals signs) underlines. The `display()` methods will be removed (see 2.1), so this is noted for completeness only.

---

## 3. Missing Core Properties and Methods

### 3.1 [x] (Done) `__hash__` correctness

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

### 3.2 [x] (Done) Add `__eq__` / equality semantics documentation

**Location:** `algebraic_pure_tensor.py` and `algebraic_tensor.py` class definitions.

**Current state:** Equality is inherited from `Basic`, which compares `args` tuples. This means two `AlgebraicPureTensor` objects are equal if and only if they have identical args (same coefficient, same factors in same order). On the other hand, two `AlgebraicTensor` objects are equal if and only if they are given as additions of identical AlgebraicPureTensors. This also means that the sequence of additions should be the same. To reflect the commutativity of addition, two `AlgebraicTensor` objects should be equal if and only if the sets construced from their args are equal. 

**Reasoning:** This is the correct behavior for a non-commutative tensor product. However, it should be documented in the `AlgebraicPureTensor` class docstring that `A ⊗ B == B ⊗ A` is `False` (order matters), to prevent user confusion. Additionally, it should be documented that different orderings of `AlgebraicPureTensor` object sums are considered equal `AlgebraicTensor` objects

**Implementation plan:** Add to the `AlgebraicPureTensor` class docstring:
```
Equality is structural: two pure tensors are equal if and only if they
have the same coefficient and identical factors in the same order.
Because the tensor product is non-commutative, ``A ⊗ B`` and ``B ⊗ A``
are not equal (assuming ``A`` and ``B`` have different shapes).
```

Add to the `AlgebraicTensor` class docstring:
```
Addition is commutative: two algebraic tensors are equal if and only if the sets
of pure tensors that define them (via addition) are equal. The ordering in which
summands appear is not relevant. E.g., ``A ⊗ B + C ⊗ D`` is equal to ``C ⊗ D + A ⊗ B``
```

**Verification needed** Ensure that two AlgebraicTensor objects are equal in case their args sets are equal:
T1 = AlgebraicPureTensor(A, B)
T2 = AlgebraicPureTensor(A2, B2)
assert (T1 + T2) == (T2 + T1)

---

## 4. Test Coverage Gaps

### 4.1 [x] (Done) Add tests for 3+ factor tensors

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

### 4.2 [x] (Done) Add tests for top-level `simplify()` interaction

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

If the above fails (which is likely), register a handler in `sympy.simplify.simplify` which invokes tensorsimplify and redo the tests.

---

### 4.3 [x] (Done) Strengthen simplify test assertions

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

### 4.4 [x] (Done) Add tests for `_compose_with_term` method

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

### 4.5 [x] (Done) Add test for `is_Add = True` behavior

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

### 5.1 [x] (Done) Sphinx documentation file exists

**Location:** `doc/src/modules/tensor/algebraic.rst` (32 lines).

**Current state:** File exists with `automodule`, `autoclass`, and `autofunction` directives for all public API. Listed in `doc/src/modules/tensor/index.rst:15`.

---

### 5.2 [x] (Missing) Expand Sphinx documentation with overview and worked example

**Location:** `doc/src/modules/tensor/algebraic.rst`.

**Current state:** The file is a minimal skeleton (32 lines) with only `autoclass`/`autofunction` directives. It has no introduction, no usage examples, and no guidance on when to use this module.

**Reasoning:** A new module in SymPy should have documentation that helps users understand:
- What problem algebraic tensors solve
- How they differ from `Indexed`, `tensorproduct`, `Array`, and the `Tensor` class
- A realistic worked example demonstrating the simplification power

This module was used by its author to calculate all results of the
mathematical physics paper `arXiv:2511.08159 <https://arxiv.org/abs/2511.08159>`_.
The documentation should include a representative example from that work:
the commutator of a Dirac operator with an algebra element in a finite
geometry model. The commutator expands to many pure tensor terms but
simplifies to only 8, demonstrating the module's simplification power.

SymPy documentation does **not** use `jupyter-execute`. Examples are written as inline `>>>` doctest-style blocks directly in the RST, following the convention seen in `doc/src/modules/matrices/expressions.rst` and `doc/src/modules/diffgeom.rst`. Mathematical notation uses the `:math:` role.

**Implementation plan:** Restructure `algebraic.rst` to follow the SymPy convention (matching `expressions.rst` and `diffgeom.rst`):

```rst
.. _tensor-algebraic:

Algebraic Tensor Products
=========================

.. module:: sympy.tensor.algebraic

The algebraic tensor module provides support for tensor products of
matrix-like objects and their linear combinations. Unlike
:class:`~sympy.tensor.array.ndim_array.NDimArray` which operates on
concrete indexed components, algebraic tensors work with symbolic
matrix expressions and preserve the tensor product structure.

This module was used to calculate all results of the mathematical
physics paper `arXiv:2511.08159 <https://arxiv.org/abs/2511.08159>`_
on finite noncommutative geometry.

A pure tensor is a tensor product of matrix factors:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> T
    A ⊗ B

Pure tensors support addition, forming algebraic tensors (linear
combinations of pure tensors):

    >>> C = MatrixSymbol("C", 3, 4)
    >>> D = MatrixSymbol("D", 4, 5)
    >>> S = AlgebraicPureTensor(C, D)
    >>> S + T
    C ⊗ D + A ⊗ B

Composition of pure tensors performs factor-wise matrix multiplication:

    >>> from sympy.tensor.algebraic import compose_algebraic_pure_tensors
    >>> X = MatrixSymbol("X", 4, 3)
    >>> Y = MatrixSymbol("Y", 5, 4)
    >>> compose_algebraic_pure_tensors(T, AlgebraicPureTensor(X, Y))
    A*X ⊗ B*Y

The :func:`~sympy.tensor.algebraic.simplify.tensorsimplify` function
can simplify algebraic tensor expressions by combining like terms and
reducing zero factors.

Worked Example: Dirac Commutator in Finite Geometry
----------------------------------------------------

The following example demonstrates the simplification power of this
module. We define a Dirac operator :math:`D` and an algebra element
:math:`a` as algebraic tensors in a finite noncommutative geometry
model (from `arXiv:2511.08159 <https://arxiv.org/abs/2511.08159>`_).
Their commutator :math:`[D, a] = D \cdot a - a \cdot D` expands
to many pure tensor terms but simplifies to only 8.

First, define the fermion mass symbols (noncommutative) and the
complex scalar symbols for the algebra element:

    >>> from sympy import symbols, Matrix, eye, zeros
    >>> from sympy.tensor.algebraic import (
    ...     AlgebraicPureTensor, AlgebraicTensor, tensorsimplify)

    >>> (upsilon_R, upsilonc_nu, cupsilon_nu, upsilon_nu,
    ...  upsilont_nu, upsilonc_R) = symbols(
    ...     r"\Upsilon_R, \Upsilon^*_\nu, \overline{\Upsilon}_\nu, "
    ...     r"\Upsilon_\nu, \Upsilon^t_\nu, \Upsilon^*_R",
    ...     commutative=False)
    >>> (upsilonc_u, cupsilon_u, upsilon_u, upsilont_u) = symbols(
    ...     r"\Upsilon^*_u, \overline{\Upsilon}_u, \Upsilon_u, "
    ...     r"\Upsilon^t_u", commutative=False)
    >>> (upsilonc_e, cupsilon_e, upsilon_e, upsilont_e) = symbols(
    ...     r"\Upsilon^*_e, \overline{\Upsilon}_e, \Upsilon_e, "
    ...     r"\Upsilon^t_e", commutative=False)
    >>> (upsilonc_d, cupsilon_d, upsilon_d, upsilont_d) = symbols(
    ...     r"\Upsilon^*_d, \overline{\Upsilon}_d, \Upsilon_d, "
    ...     r"\Upsilon^t_d", commutative=False)

    >>> (z, w, alpha, beta, gamma, delta) = symbols(
    ...     r"z, w, \alpha, \beta, \gamma, \delta",
    ...     complex=True)
    >>> (m11, m12, m13, m21, m22, m23,
    ...  m31, m32, m33) = symbols(
    ...     r"m_{11}, m_{12}, m_{13}, m_{21}, m_{22}, "
    ...     r"m_{23}, m_{31}, m_{32}, m_{33}", complex=True)

Define the Dirac operator as a sum of four pure tensors, each with
three matrix factors (2x2, 4x4, 4x4):

    >>> proj_L = Matrix([[1, 0], [0, 0]])
    >>> proj_R = Matrix([[0, 0], [0, 1]])
    >>> id4 = eye(4)
    >>> pL_f = Matrix([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    >>> pR_f = Matrix([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    >>> m_nu = Matrix([
    ...     [0, 0, upsilon_R, upsilonc_nu],
    ...     [0, 0, cupsilon_nu, 0],
    ...     [upsilonc_R, upsilont_nu, 0, 0],
    ...     [upsilon_nu, 0, 0, 0]])
    >>> m_u = Matrix([
    ...     [0, 0, 0, upsilonc_u],
    ...     [0, 0, cupsilon_u, 0],
    ...     [0, upsilont_u, 0, 0],
    ...     [upsilon_u, 0, 0, 0]])
    >>> m_e = Matrix([
    ...     [0, 0, 0, upsilonc_e],
    ...     [0, 0, cupsilon_e, 0],
    ...     [0, upsilont_e, 0, 0],
    ...     [upsilon_e, 0, 0, 0]])
    >>> m_d = Matrix([
    ...     [0, 0, 0, upsilonc_d],
    ...     [0, 0, cupsilon_d, 0],
    ...     [0, upsilont_d, 0, 0],
    ...     [upsilon_d, 0, 0, 0]])

    >>> D1 = AlgebraicPureTensor(proj_L, pL_f, m_nu)
    >>> D2 = AlgebraicPureTensor(proj_L, pR_f, m_u)
    >>> D3 = AlgebraicPureTensor(proj_R, pL_f, m_e)
    >>> D4 = AlgebraicPureTensor(proj_R, pR_f, m_d)
    >>> Dirac = D1 + D2 + D3 + D4

Define the algebra element :math:`a` as a sum of three pure tensors
(the representation :math:`\pi(a)`):

    >>> a_f1 = Matrix([[z, 0], [0, w]])
    >>> a_f2 = Matrix([[alpha, beta], [gamma, delta]])
    >>> a_m  = Matrix([
    ...     [z, 0, 0, 0],
    ...     [0, m11, m12, m13],
    ...     [0, m21, m22, m23],
    ...     [0, m31, m32, m33]])
    >>> a_f3 = Matrix([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,0]])

    >>> a_pi1 = AlgebraicPureTensor(a_f1, id4, pL_f)
    >>> a_pi2 = AlgebraicPureTensor(a_f2, id4, pR_f)
    >>> a_pi3 = AlgebraicPureTensor(eye(2), a_m, a_f3)
    >>> a = a_pi1 + a_pi2 + a_pi3

Compute the commutator, expand, and simplify:

    >>> da = (Dirac * a1 - a1 * Dirac).expand()
    >>> len(da.terms)
    24
    >>> da1_simp = tensorsimplify(da)
    >>> len(da_simp.terms)
    8

The simplification reduces the expression from 24 pure tensor terms to
just 8, by combining like terms and eliminating zero factors.

Classes
-------

.. autoclass:: AlgebraicPureTensor
   :members:

.. autoclass:: AlgebraicTensor
   :members:

.. autoclass:: AlgebraicZeroTensor
   :members:

.. autoclass:: ShapeMismatchError
   :members:

Functions
---------

.. autofunction:: algebraic_tensor_product

.. autofunction:: compose_algebraic_pure_tensors

.. autofunction:: compose_algebraic_tensors

.. autofunction:: tensorsimplify
```

Note the following SymPy conventions:
- `.. module::` directive at the top (not `.. automodule::`), which sets the module context for cross-references
- `>>>` doctest-style examples inline in the RST text (NOT `jupyter-execute`)
- `:math:` role for mathematical notation (NOT raw LaTeX or `.. math::` blocks for inline math)
- `:class:~` and `:func:~` cross-reference roles for internal links (the `~` hides the module prefix in the displayed text)
- `autoclass` and `autofunction` directives reference objects without the full module prefix since `.. module::` sets the context
- Brief prose sections with examples, not separate "Overview", "Mathematical Background", or "Quick Example" sections

**Verification required before committing:** The agent must run the full example
(interactive Python session or script) to verify the exact number of terms before
and after simplification. The numbers 24 and 8 above are based on the playground
file (`tests_old/playground_functions.py` lines 567-569) and may differ depending
on how `.expand()` distributes terms. The agent should adjust the numbers in the
RST to match the actual output. Also verify that the `len(da1.terms)` call works
-- the attribute may be `.terms`, `.args`, or something else depending on the
`AlgebraicTensor` API.

---

### 5.3 [x] (Missing) Verify all docstring examples pass doctest

**Location:** All files with docstrings.

**Reasoning:** Per SymPy's docstring style guide, every example in a docstring must pass `./bin/doctest`. Some examples in the current code use `print()` (e.g., `algebraic_pure_tensor.py:101`), which produces output that doctest verifies. Other examples use `>>> expr` and rely on the repr, which may change.

**Implementation plan:**
1. Run `./bin/doctest tensor/algebraic` to verify all examples pass
2. Fix any failing examples
3. Ensure examples use `from sympy import ...` style imports (per SymPy convention), not `from sympy.tensor.algebraic import ...` in the outer scope — though this is a gray area for module-local examples

---

## 6. Code Quality Improvements

### 6.1 [ ] (Missing) Reduce code duplication in `compose_algebraic_pure_tensors` zero-tensor handling

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

## 7. Cleanup

### 7.1 [ ] (Missing) Delete `tests_old/` directory

See item 2.2 above.

---

### 7.2 [ ] (Missing) Add `tensor/algebraic/` to `.gitignore` for `.DS_Store` and `__pycache__`

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

### 7.3 [ ] (Missing) Remove `.DS_Store` from all directories

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

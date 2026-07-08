# SymPy Contribution Checklist — Algebraic Tensor Module

> **How to use this file:** Each section header carries a status flag such as
> `(Missing)` or `(In Progress)`. After every checkbox in a section is marked
> `- [x]`, update that section's flag to `(Done)`. When all sections are
> `(Done)`, the module is ready for PR submission.

## 1. Tests (Done)

5 test modules created in `tensor/algebraic/tests/` — 177 tests, all passing.
SymPy test guide: https://docs.sympy.org/latest/contributing/new-contributors-guide/writing-tests.html

- [x] `tensor/algebraic/tests/test_algebraic_pure_tensor.py` — 51 tests covering AlgebraicPureTensor constructor, properties, arithmetic, expand, and composition
- [x] `tensor/algebraic/tests/test_algebraic_tensor.py` — 35 tests covering AlgebraicTensor constructor, terms, arithmetic, shape mismatch, and composition dispatch
- [x] `tensor/algebraic/tests/test_algebraic_zero_tensor.py` — 33 tests covering AlgebraicZeroTensor singleton, shape, arithmetic identity, and preservation
- [x] `tensor/algebraic/tests/test_scalar_mul.py` — 32 tests covering ScalarMul constructor, scalar/tensor properties, multiplication fusion, and arithmetic
- [x] `tensor/algebraic/tests/test_simplify.py` — 26 tests covering tensorsimplify, proportionality_factoring, and edge cases

Each test file covers (per DESIGN.md §9.10):

- [x] Constructor unwrapping behavior (single term, zero coefficient, zero tensor factor, empty)
- [x] Shape mismatch detection in addition and composition
- [x] Full composition dispatch table (all operand type combinations)
- [x] Scalar multiplication vs composition dispatch (commutative vs non-commutative)
- [x] Negation behavior
- [x] `proportionality_factoring` with all/one/more-than-one non-proportional slots
- [x] `AlgebraicZeroTensor` preservation through operations
- [x] `tensorsimplify` on each type

---

## 2. Documentation (Done)

No Sphinx documentation exists for the algebraic module. Follow the official Sympy guide at https://docs.sympy.org/latest/contributing/documentation-style-guide.html

- [x] Create `doc/src/modules/tensor/algebraic.rst` with `automodule`/`autoclass` directives for all public classes and functions
- [x] Update `doc/src/modules/tensor/index.rst` to add `algebraic.rst` to the toctree

---

## 3. Module-level Docstrings (Missing)

Follow the official Sympy guide at https://docs.sympy.org/latest/contributing/docstring.html

- [ ] Add module-level docstring to `tensor/algebraic/__init__.py` (summary + examples)
- [ ] Add module-level docstring to `tensor/algebraic/algebraic_zero_tensor.py` (summary + examples)
- [ ] Add module-level docstring to `tensor/algebraic/algebraic_pure_tensor.py` (summary + examples)
- [ ] Add module-level docstring to `tensor/algebraic/algebraic_tensor.py` (summary + examples)
- [ ] Add module-level docstring to `tensor/algebraic/simplify.py` (summary + examples)
- [ ] Add module-level docstring to `tensor/algebraic/scalar_mul.py` (summary + examples)

Per SymPy's Docstring Style Guide: every public function/class/method needs a
**Single-Sentence Summary** and an **Examples** section with doctests.

### `algebraic_pure_tensor.py`
- [ ] `AlgebraicPureTensor` class docstring — add `Examples` section with doctests
- [ ] `AlgebraicPureTensor.factors` property — add `Examples`
- [ ] `AlgebraicPureTensor.tensor_shape` property — add `Examples`
- [ ] `AlgebraicPureTensor.commutativity_shape` property — add `Examples`
- [ ] `AlgebraicPureTensor.__new__` — add docstring with `Examples`
- [ ] `AlgebraicPureTensor.__neg__` — add `Examples`
- [ ] `AlgebraicPureTensor.__mul__` — add `Examples`
- [ ] `AlgebraicPureTensor.__rmul__` — add `Examples`
- [ ] `AlgebraicPureTensor.__add__` / `__radd__` — add `Examples`
- [ ] `AlgebraicPureTensor.__sub__` / `__rsub__` — add `Examples`
- [ ] `AlgebraicPureTensor._eval_expand_mul` — add docstring with `Examples`
- [ ] `algebraic_tensor_product` — add `Examples`
- [ ] `compose_algebraic_pure_tensors` — add `Examples`

### `algebraic_tensor.py`
- [ ] `AlgebraicTensor` class docstring — add `Examples` section with doctests
- [ ] `AlgebraicTensor.tensor_shape` property — add `Examples`
- [ ] `AlgebraicTensor.commutativity_shape` property — add `Examples`
- [ ] `AlgebraicTensor.terms` property — add `Examples`
- [ ] `AlgebraicTensor.__neg__` — add `Examples`
- [ ] `AlgebraicTensor.__add__` / `__radd__` — add `Examples`
- [ ] `AlgebraicTensor.__sub__` / `__rsub__` — add `Examples`
- [ ] `AlgebraicTensor.__mul__` — add `Examples`
- [ ] `AlgebraicTensor.__rmul__` — add `Examples`
- [ ] `AlgebraicTensor.expand` — add `Examples`
- [ ] `ShapeMismatchError` — add docstring with `Examples`
- [ ] `compose_algebraic_tensors` — add `Examples`

### `algebraic_zero_tensor.py`
- [ ] `AlgebraicZeroTensor` class docstring — add `Examples` section with doctests
- [ ] `AlgebraicZeroTensor.shape` property — add `Examples`
- [ ] `AlgebraicZeroTensor.tensor_shape` property — add `Examples`
- [ ] `AlgebraicZeroTensor.commutativity_shape` property — add `Examples`
- [ ] `AlgebraicZeroTensor.__neg__` — add `Examples`
- [ ] `AlgebraicZeroTensor.__add__` / `__radd__` — add `Examples`
- [ ] `AlgebraicZeroTensor.__sub__` / `__rsub__` — add `Examples`
- [ ] `AlgebraicZeroTensor.__mul__` / `__rmul__` — add `Examples`
- [ ] `algebraic_zero_tensor` — add `Examples`

### `scalar_mul.py`
- [ ] `ScalarMul` class docstring — add `Examples` section with doctests
- [ ] `ScalarMul.scalar` property — add `Examples`
- [ ] `ScalarMul.tensor` property — add `Examples`
- [ ] `ScalarMul.factors` property — add `Examples`
- [ ] `ScalarMul.tensor_shape` property — add `Examples`
- [ ] `ScalarMul.commutativity_shape` property — add `Examples`
- [ ] `ScalarMul.__neg__` — add `Examples`
- [ ] `ScalarMul.__mul__` — add `Examples`
- [ ] `ScalarMul.__rmul__` — add `Examples`
- [ ] `ScalarMul.__add__` / `__radd__` — add `Examples`
- [ ] `ScalarMul.__sub__` / `__rsub__` — add `Examples`
- [ ] `ScalarMul.expand` — add `Examples`

### `simplify.py`
- [ ] `tensorsimplify` — add `Examples`
- [ ] `proportionality_factoring` — add `Examples`
- [ ] `_extract_pt_and_coeff` — add docstring with `Examples`
- [ ] `_proportionality_ratio` — add `Examples`
- [ ] `_build_pt` — add `Examples`
- [ ] `_matrix_proportionality_ratio` — add docstring with `Examples`
- [ ] `_decompose_commutative_factors` — add `Examples`
- [ ] `_reconstruct_term` — add `Examples`
- [ ] `_extract_commutative_from_factor` — add `Examples`
- [ ] `_extract_commutative_prefactors` — add `Examples`
- [ ] `_normalize_factor_sign` — add docstring with `Examples`
- [ ] `_is_exactly_divisible` — add docstring with `Examples`
- [ ] `_deduplicate_proportional` — add docstring with `Examples`

---

## 4. Fix `_repr_latex_` / `display()` Methods (Missing)

The four `display()` methods call `self._repr_latex_()` which is never defined
on these classes. This will raise `AttributeError`.

- [ ] Implement `_repr_latex_()` on `AlgebraicPureTensor` (use the printer dispatch system: `latex(self)`)
- [ ] Implement `_repr_latex_()` on `AlgebraicTensor`
- [ ] Implement `_repr_latex_()` on `ScalarMul`
- [ ] Implement `_repr_latex_()` on `AlgebraicZeroTensor`

---

## 5. Pretty Printing (Missing)

`printing/pretty.py` has no handlers for the algebraic types.

- [ ] Add `_print_AlgebraicPureTensor` to `sympy/printing/pretty/pretty.py`
- [ ] Add `_print_AlgebraicTensor` to `sympy/printing/pretty/pretty.py`
- [ ] Add `_print_ScalarMul` to `sympy/printing/pretty/pretty.py`
- [ ] Add `_print_AlgebraicZeroTensor` to `sympy/printing/pretty/pretty.py`

---

## 6. Str Printing (Missing)

`printing/str.py` has no handlers for the algebraic types. The `__str__` methods
exist on the classes but should integrate with the printer system.

- [ ] Add `_print_AlgebraicPureTensor` to `sympy/printing/str.py`
- [ ] Add `_print_AlgebraicTensor` to `sympy/printing/str.py`
- [ ] Add `_print_ScalarMul` to `sympy/printing/str.py`
- [ ] Add `_print_AlgebraicZeroTensor` to `sympy/printing/str.py`

---

## 7. Code Quality and Validation (Missing)

- [ ] Run `./bin/doctest` to verify all docstring examples pass
- [ ] Run `./bin/test sympy/tensor/algebraic` to verify all tests pass
- [ ] Check for lines >80 chars in docstrings
- [ ] Check for trailing whitespace in all files
- [ ] Verify `__slots__` usage is consistent across all classes
- [ ] Remove `DESIGN.md` or move it to `doc/src/` (not a standard SymPy file location)
- [ ] Remove `.DS_Store` from `tensor/algebraic/`
- [ ] Remove `__pycache__/` from `tensor/algebraic/`

---

## 8. GitHub PR Preparation (Missing)

- [ ] Create a GitHub issue describing the new module (or link to existing one)
- [ ] Open a PR with a clear description referencing the issue
- [ ] Ensure all CI checks pass (tests, flake8, doctests, docs build)

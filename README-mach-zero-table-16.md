# Round 1

- Project chosen: `sympy`

- Reason: All team members are familiar with the tech stack required, and also we found a lot of open issues that looked challenging, doable, as well as meaningful to the project.

- Issue chosen: [#29714](https://github.com/sympy/sympy/issues/29714): Derivative of function with composite arguments creates non-Dummy symbols with expression names

### Issue

- The bug occurs when differentiating an undefined function with multiple arguments that are compositions of base symbols (e.g. `f(x + y, x * y)`).
- SymPy introduces *temporary symbols* whose names are derived from the string representations of these composite expressions.
- For example, `Symbol('_x + y')` appears during the evaluation of nested `Derivative` and `Subs` expressions to represent `x + y`.
- These symbols are ordinary `Symbol` objects rather than `Dummy` objects, so they are not treated as internal temporary variables. As a result, they are not correctly substituted back with their corresponding expressions (e.g., `x + y`).
- These unintended symbols leak into `free_symbols` and remain in the final expression, breaking subsequent evaluation and `lambdify`.

### Action Plan

- Identify where the unwanted symbols are created.
  - Investigate Derivative.__new__() and its handling of non-symbol differentiation variables.
  - Trace where symbols such as `Symbol('_x + y')` are constructed.
- Trace the nested evaluation process.
  - Investigate how `Subs.doit()` and `Subs._eval_subs()` interact with nested `Derivative` expressions.
  - Identify where the relationship between temporary symbols and their original expressions is lost.
- Identify why the symbols leak.
  - Determine why the generated symbols are not converted back into their original expressions.
  - Check why they remain in `free_symbols` after evaluation.
  - Develop and test potential fixes.
  - Ensure generated intermediate variables are proper `Dummy` objects where appropriate.
  - Preserve the mapping between temporary variables and their original expressions.
  - Add regression tests covering multiple composite function arguments, nested differentiation, and subsequent `doit()` / lambdify() operations.

**Goal:** Ensure that intermediate differentiation variables remain internal and that the final expression contains only the original free symbols (x and y), without unintended symbols such as Symbol('_x + y').

---

# Round 2 & Round 3

### What we tried:

- Started with `sympy/core/function.py`

- We looked `Derivatives` and it did not give any leads.

- `Subs.__new__()`, `Subs._hashable_content()`, and other functions looked relevant to us under `Subs`, especially this snippet because it involved creating a Symbol.

  ```python
  Sols
  ...
  __new__
  ...
  if any(r in expr.free_symbols and
      r in variables and
      Symbol(pre + mystr(point[variables.index(r)])) != r
      ...
  ...
  ```

- We used this test case to represent the issue:

  ```python
  def test_Subs_composite_arguments_derivative():
    x, y = symbols('x y')
    f = Function('f')
  
    result = f(x + y, x*y).diff(y).diff(x)
  
    assert result.free_symbols == {x, y}
  ```

- This is the hotfix we initially tried before understanding the real problem (i.e. replacing `Symbol` with `Dummy`)

  ```python
  Sols
  ...
  __new__
  ...
  if any(r in expr.free_symbols and
      r in variables and
      Dummy(pre + mystr(point[variables.index(r)])) != r
      ...
  ...
  ```

- It failed with `RecursionError`.

- Then we assumed `doit()` was the culprit and tried to find leads in that direction. Then we realized that `doit()` method was not actually causing the problem. 

- The actual problem had to be with the expression evaluation and related to checking whether or not to substitute, and not with the actual substitution part.

### What we came up with:

- Then after that we finally came up with this fix:

  ```python
  if any(r in expr.free_symbols and (
       r not in variables or
       s_pts[point[variable.index(r)]] != r)
  ```

- Changed `hashable_content` as well to match the new implementation:

  ```python
      return (self._expr.xreplace(self.canonical_variables),
          ) + tuple(ordered([(v, p) for v, p in 
              zip(self.variables, self.point) if not self.expr.has(v)]))
  ```

  ```python
      return (
          self._expr.xreplace(self.canonical_variables),
          len(self.variables),
          frozenset(self.free_symbols),
          tuple(sorted(self.point, key=default_sort_key)),
          ) + tuple(ordered([(v, p) for v, p in
              zip(self.variables, self.point) if not self.expr.has(v)]))
  ```

- This passes our test case.

### Why the fix works

- What does `_hashable_content()` do?

  - SymPy uses this method to determine the identity of an expression for hashing and equality comparisons.

  - In particular, two `Subs` objects can have different internal representations but may be considered equivalent after canonicalization.

  - The problem is that comparing only the canonicalized expression and the remaining variable-point pairs may not capture all the information needed to distinguish different `Subs` objects.

  - The new implementation adds:

  | Addition                                          | Purpose                                                                 |
  | ------------------------------------------------- | ----------------------------------------------------------------------- |
  | `len(self.variables)`                             | Records the number of substitution variables.                           |
  | `frozenset(self.free_symbols)`                    | Includes the free-symbol information in the object's identity.          |
  | `tuple(sorted(self.point, key=default_sort_key))` | Includes the substitution points in a canonical, sorted representation. |

- Why is this relevant to our bug?

  - The issue involves generated symbols that represent composite substitution points. During canonicalization, different `Subs` objects can have similar-looking internal expressions even though their substitution information differs.

  - By including the additional information, SymPy can distinguish these objects more reliably rather than accidentally treating different substitutions as equivalent.

- Change in the collision-checking condition

  - This change is most relevant to the original bug.

  ```python
  r in expr.free_symbols and (
  r not in variables or 
  s_pts[point[variables.index(r)]] != r)
  ```

  - Subs.__new__() generates replacement symbols. However, if the generated symbol already exists as a free symbol in the expression, that creates a symbol collision. SymPy needs to ensure that its internally generated replacement symbol does not accidentally refer to an existing user-defined symbol.

  - The old condition only checked whether r was among the substitution variables.

  - The new condition handles two cases:

    1. r is not a substitution variable. This indicates a potential collision that needs attention.

    2. r is a substitution variable. SymPy checks whether the generated replacement symbol for its corresponding substitution point differs from r.

  - It also reuses the already-created replacement symbol from the s_pts mapping. That avoids independently reconstructing a symbol and ensures the collision check uses the same object mapping as the rest of the canonicalization process.

# Round 3 Testing and Validation

## 1. Original Issue Reproduction

We used the following standalone Python script to reproduce the original issue:

```python
import sympy as sp

x, y = sp.symbols('x y')
f = sp.Function('f')

# Reproduce the issue with multiple composite function arguments.
expr = f(x + y, x*y)
result = expr.diff(y).diff(x)

expected_symbols = {x, y}
actual_symbols = result.free_symbols
unexpected_symbols = actual_symbols - expected_symbols

print(f"SymPy version: {sp.__version__}")
print(f"Input:         {expr}")
print(f"Result:        {result}")
print(f"Free symbols:  {actual_symbols}")

print("\n--- Regression check ---")
print(f"Expected:      {expected_symbols}")
print(f"Unexpected:    {unexpected_symbols}")

if unexpected_symbols:
    print("FAIL: Unexpected free symbols detected!")
else:
    print("PASS: No unexpected free symbols.")
```

### Before the fix

The script was executed against the original `master` branch, before applying our changes.

**Output:**

```text
SymPy version: 1.15.0.dev3298+gc508f90ea0.uninstalled
Input:         f(x + y, x*y)

Free symbols:  {x, _x + y, y}

--- Regression check ---
Expected:      {y, x}
Unexpected:    {_x + y}
FAIL: Unexpected free symbols detected!
```

The resulting expression contained an unexpected symbol, `_x + y`, confirming that the original issue was reproducible.

### After the fix

The same script was executed against our fixed branch, after applying our fixes.

**Output:**

```text
SymPy version: 1.15.0.dev3304+gbbff6115cf.uninstalled
Input:         f(x + y, x*y)

Free symbols:  {x, y}

--- Regression check ---
Expected:      {x, y}
Unexpected:    set()
PASS: No unexpected free symbols.
```

The unexpected symbol no longer appears in `free_symbols`. The result contains only the original free symbols, confirming that the original reproducer is fixed.

---

## 2. Expanded Regression Tests

We expanded the original regression test to cover additional combinations of composite arguments and differentiation orders.

```python
def test_Subs_composite_arguments_derivative():
    x, y, z = symbols('x y z')
    f = Function('f')

    # Regression test: differentiating a function with multiple
    # composite arguments must not introduce spurious free symbols.
    result = f(x + y, x*y).diff(y).diff(x)
    assert result.free_symbols == {x, y}

    # Check that repeated differentiation with different composite
    # arguments and differentiation orders does not introduce
    # unexpected free symbols.
    cases = [
        # Different differentiation order.
        (f(x + y, x*y), (x, y), {x, y}),

        # Composite arguments involving powers and constants.
        (f(x**2 + y, x*y + 1), (y, x), {x, y}),

        # Multiple arguments containing the same composite expression.
        (f(x + y, x + y), (y, x), {x, y}),

        # Three composite arguments with overlapping free symbols.
        (f(x + y, x*y, x*z + y), (y, x), {x, y, z}),

        # Nested composite expressions.
        (f((x + y)**2, x*y + x), (y, x), {x, y}),

        # Repeated differentiation of the same expression.
        (f(x + y, x*y), (y, x, y), {x, y}),

        # Different powers and a different differentiation order.
        (f(x**2 + y**2, x*y + y**3), (x, y, x), {x, y}),
    ]

    for expr, variables, allowed_symbols in cases:
        result = expr
        for variable in variables:
            result = result.diff(variable)

        assert result.free_symbols <= allowed_symbols
```

The additional cases cover:

- Different differentiation orders.
- Composite arguments involving powers and constants.
- Repeated composite arguments.
- Functions with three composite arguments.
- Nested composite expressions.
- Multiple successive differentiations.
- Different combinations of powers and differentiation variables.

The original reproducer remains as an explicit assertion checking the exact expected free-symbol set. The additional cases verify that no unexpected free symbols are introduced across the other scenarios.

---

## 3. Test Suite Results

### Targeted test suite

Executed:

```bash
./bin/test test_function
```

**Output:**

```text
==========================================================================================
145 passed, 1 xfailed in 8.56s
==========================================================================================
```

All 145 passing tests succeeded, with one expected failure (`xfail`).

### Full SymPy test suite

Executed:

```bash
./bin/test
```

**Result:** All tests passed.

This provides broader validation that the changes to `Subs` construction and hashing do not introduce regressions elsewhere in SymPy.

---

## 4. Before-and-After Comparison

| Validation | Original version | Fixed version |
|---|---|---|
| Original reproducer | Fails | Passes |
| Free symbols | `{x, y, _x + y}` | `{x, y}` |
| Unexpected symbols | `_x + y` | None |
| Expanded regression cases | Not covered | Pass |
| Targeted test suite | Existing baseline | 145 passed, 1 xfailed |
| Full test suite | — | All tests passed |

## 5. Final Validation

- Successfully reproduced the original issue on the unmodified `master` branch.
- Verified that the same standalone script passes against the fixed branch.
- Added expanded regression coverage for multiple composite argument structures and differentiation sequences.
- Confirmed that the targeted test suite passes.
- Confirmed that the full SymPy test suite passes.

**Conclusion:** The fix resolves the original composite-argument differentiation issue and passes both targeted and codebase-wide testing.

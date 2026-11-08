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
  - Add regression tests covering multiple composite function arguments, nested differentiation, and subsequent `doit()` / `lambdify()` operations.

**Goal:** Ensure that intermediate differentiation variables remain internal and that the final expression contains only the original free symbols (x and y), without unintended symbols such as `Symbol('_x + y')`.

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

- We used this simplified version of the github issue reproducer as the test case to represent the issue:

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

- The tests failed with `RecursionError`.

- Then we assumed `doit()` was the culprit and tried to find leads in that direction. Then we realized that `doit()` method was not actually causing the problem. 

- The actual problem had to be with the expression evaluation and related to checking whether or not to substitute, and not with the actual substitution part.

### What we came up with

After investigating the interaction between `Subs` construction, symbol replacement, and expression evaluation, we identified two changes that resolved the issue.

**1. Improved the symbol collision check in `Subs.__new__()`**

We updated the condition responsible for detecting collisions between generated replacement symbols and existing free symbols.

**Previous implementation:**

```python
if any(
    r in expr.free_symbols
    and r in variables
    and Symbol(pre + mystr(point[variables.index(r)])) != r
    for _, r in reps
):
```

**Updated implementation:**

```python
if any(
    r in expr.free_symbols and (
        r not in variables or
        s_pts[point[variables.index(r)]] != r
    )
    for _, r in reps
):
```

The original condition only checked for collisions when the existing symbol was also one of the substitution variables. The updated condition handles both cases:

- If the symbol is not among the substitution variables, it is treated as a potential collision.
- If it is a substitution variable, SymPy checks whether its corresponding generated replacement symbol differs from it.

We also reused the existing `s_pts` mapping instead of constructing another `Symbol` for the comparison. This keeps the collision check consistent with the replacement symbols already created during substitution.

**2. Expanded `Subs._hashable_content()`**

We also modified `_hashable_content()` to include additional information when determining the structural identity of a `Subs` object.

**Previous implementation:**

```python
def _hashable_content(self):
    return (self._expr.xreplace(self.canonical_variables),) + tuple(
        ordered([
            (v, p)
            for v, p in zip(self.variables, self.point)
            if not self.expr.has(v)
        ])
    )
```

**Updated implementation:**

```python
def _hashable_content(self):
    return (
        self._expr.xreplace(self.canonical_variables),
        len(self.variables),
        frozenset(self.free_symbols),
        tuple(sorted(self.point, key=default_sort_key)),
    ) + tuple(ordered([
        (v, p)
        for v, p in zip(self.variables, self.point)
        if not self.expr.has(v)
    ]))
```

The updated implementation adds three pieces of information:

- `len(self.variables)` records the number of substitution variables.
- `frozenset(self.free_symbols)` includes the free-symbol information in the object's identity.
- `tuple(sorted(self.point, key=default_sort_key))` includes the substitution points in a consistently ordered representation.

These additions provide more information when SymPy compares and hashes `Subs` objects. This helps distinguish substitutions whose canonicalized expressions may look similar but whose substitution details differ.

**Outcome**

Together, these changes resolved the original reproducer. Differentiating a function with multiple composite arguments no longer introduces unintended symbols such as `Symbol('_x + y')` into `free_symbols`.

We subsequently verified the fix using the original reproducer, expanded regression tests, and SymPy's test suite. The results are documented in the [Testing and Validation](./TESTING_Mach-Zero.md) section.

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

  - `Subs.__new__()` generates replacement symbols. However, if the generated symbol already exists as a free symbol in the expression, that creates a symbol collision. SymPy needs to ensure that its internally generated replacement symbol does not accidentally refer to an existing user-defined symbol.

  - The old condition only checked whether r was among the substitution variables.

  - The new condition handles two cases:

    1. r is not a substitution variable. This indicates a potential collision that needs attention.

    2. r is a substitution variable. SymPy checks whether the generated replacement symbol for its corresponding substitution point differs from r.

  - It also reuses the already-created replacement symbol from the s_pts mapping. That avoids independently reconstructing a symbol and ensures the collision check uses the same object mapping as the rest of the canonicalization process.

# Round 3 Testing and Validation

We thoroughly validated that any changes we made were not breaking changes in terms of the general functionality of the SymPy upstream. Refer to [Testing and Validation](./TESTING_Mach-Zero.md) for a detailed record of what we did to ensure this.


---

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

# Round 2

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

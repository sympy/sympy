# Writing Tests

The most important thing for a mathematical library like SymPy is correctness.
Functions should never return mathematically incorrect results. Correctness is
always the top concern, even if it comes at the cost of things like
performance or modularity.

Consequently, all functionality in SymPy is tested extensively. This guide
goes over how tests in SymPy are written.

## Testing Policies

In order to ensure the high standard of correctness, SymPy has the following
rules that apply to all pull requests:

1. All new functionality must be tested. Tests should aim to cover all
   possible cases to best ensure correctness. This means not only maximizing
   code coverage, but also covering all possible corner cases.

2. Every pull request must pass all tests before it can be merged. The tests
   are automatically run by the GitHub Actions CI on every pull request. If
   any tests fail, the CI will fail with a red ❌. These failures must be
   addressed before the pull request can be merged.

3. Bug fixes should be accompanied by a [regression
   test](writing-tests-regression-tests).


## Basics for Writing Tests

Tests are located alongside the code in `tests/` directories, in files named
`test_<thing>.py`. In most cases, if you modified
`sympy/<submodule>/<file>.py` then the test for the functionality will go in
`sympy/<submodule>/tests/test_<file>.py`. For example, the tests for the
functions in `sympy/simplify/sqrtdenest.py` are in
`sympy/simplify/tests/test_sqrtdenest.py`. There are some exceptions to this
rule, so in general try to find where the existing tests are for a function
and add your tests alongside them. If you are adding tests for a new function,
follow the general pattern of tests in the module you are adding to.

Tests follow a simple pattern, which should be apparent from reading the
existing test files. Tests are in functions that start with `test_` and
contain lines like

```py
assert function(arguments) == result
```

For example

```py
# from sympy/functions/elementary/tests/test_trigonometric.py

def test_cos_series():
    assert cos(x).series(x, 0, 9) == \
        1 - x**2/2 + x**4/24 - x**6/720 + x**8/40320 + O(x**9)
```

New test cases can be added to an existing test function if it is relevant, or
you can create a new test function.


(writing-tests-regression-tests)=
## Regression Tests

Regression tests are tests that would fail before the bug fix but now passes.
Often you can use a code example from an issue as a test case, although it is
also OK to simplify such examples or to write your own, so long as it tests
the issue in question.

For example, consider [issue
#21177](https://github.com/sympy/sympy/issues/21177), which identified the
following wrong result:

```py
>>> residue(cot(pi*x)/((x - 1)*(x - 2) + 1), x, S(3)/2 - sqrt(3)*I/2) # doctest: +SKIP
-sqrt(3)*tanh(sqrt(3)*pi/2)/3
>>> residue(cot(pi*x)/(x**2 - 3*x + 3), x, S(3)/2 - sqrt(3)*I/2) # doctest: +SKIP
0
```

Here the first expression was correct but the second is not. In the issue, the
cause of the issue was identified in the `as_leading_term` method, and several
other related issues were also found.

In the corresponding pull request
([#21253](https://github.com/sympy/sympy/pull/21253/files)), several
regression tests were added. For example (from that PR):

```
# In sympy/core/tests/test_expr.py

def test_as_leading_term():
    ...
    # This test was already existing. The following was added to the end:

    # https://github.com/sympy/sympy/issues/21177
    f = -3*x + (x + Rational(3, 2) - sqrt(3)*S.ImaginaryUnit/2)**2\
        - Rational(3, 2) + 3*sqrt(3)*S.ImaginaryUnit/2
    assert f.as_leading_term(x) == \
        (3*sqrt(3)*x - 3*S.ImaginaryUnit*x)/(sqrt(3) + 3*S.ImaginaryUnit)

    # https://github.com/sympy/sympy/issues/21245
    f = 1 - x - x**2
    fi = (1 + sqrt(5))/2
    assert f.subs(x, y + 1/fi).as_leading_term(y) == \
        (-36*sqrt(5)*y - 80*y)/(16*sqrt(5) + 36)
```

```
# In sympy/functions/elementary/tests/test_trigonometric.py

def test_tan():
    ...
    # This test was already existing. The following was added to the end:

    # https://github.com/sympy/sympy/issues/21177
    f = tan(pi*(x + S(3)/2))/(3*x)
    assert f.as_leading_term(x) == -1/(3*pi*x**2)
```

```py
# In sympy/series/tests/test_residues.py

def test_issue_21177():
    r = -sqrt(3)*tanh(sqrt(3)*pi/2)/3
    a = residue(cot(pi*x)/((x - 1)*(x - 2) + 1), x, S(3)/2 - sqrt(3)*I/2)
    b = residue(cot(pi*x)/(x**2 - 3*x + 3), x, S(3)/2 - sqrt(3)*I/2)
    assert a == r
    assert (b - a).cancel() == 0
```

This example shows some important aspects of regression tests:

- Tests should be added for the underlying fix, not just the originally
  reported issue. The originally reported issue in this example was with the `residue()`
  function but the underlying issue was with the `as_leading_term()` method.

- At the same time, it can also be beneficial to add a test for the high-level
  issue as reported. This ensure that `residue` itself won't break in the
  future, even if the implementation details of it change so that it no longer
  uses the same code path that was fixed.

- This example does not show it, but in some cases it may be prudent to
  simplify the originally reported issue for the test case. For example,
  sometimes users will include unnecessary details in the report that don't
  actually matter for the reproduction of the issue (like unnecessary
  assumptions on symbols), or make the input expression to large or have too
  many unnecessary constant symbols. This is especially important to do if the
  originally stated issue is slow to compute. If the same thing can be tested
  with a test that runs more quickly, this should be preferred.

- Regression tests should also be added for additional bugs that are
  identified in the issue.

- It is useful to cross-reference the issue number in a regression test,
  either using a comment or in the test name. A comment is preferred if the
  test is being added to an existing test.

Regression tests aren't just for bug fixes. They should also be used for new
features, to make sure the newly implemented functionality remains implemented
and correct.

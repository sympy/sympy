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

(writing-tests-basics)=
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


## Special Types of Tests

Most tests will be of the form `assert function(input) == output`. However,
there are other types of things that you might want to test that should be
tested in certain ways.

### Testing Exceptions

To test that a function raises a given exception, use
`sympy.testing.pytest.raises`. `raises()` takes an exception class and a
lambda. For example

```py
from sympy.testing.pytest.raises
raises(TypeError, lambda: cos(x, y)
```

Remember to include the `lambda`. Otherwise, the code will be executed
immediately and will raise the exception, causing the test to fail.

```
# BAD
raises(TypeError, cos(x, y)) # This test will fail
```

`raises` can also be used as a context manager, like

```
with raises(TypeError):
    cos(x, y)
```

However, be careful using this form, as it can only check one expression at a
time. If the code under context manager raises multiple exceptions, only the
first one will actually be tested

```
# BAD
with raises(TypeError):
   cos(x, y)
   sin(x, y) # THIS WILL NEVER BE TESTED
```

The `lambda` form is generally better because it avoids this problem, although
if you are testing something that cannot be represented in a `lambda` you will
need to use the context manager form.

(writing-tests-testing-warnings)=
### Testing Warnings

[Warnings](https://docs.python.org/3/library/warnings.html) can be tested with
the [`sympy.testing.pytest.warns()`](warns) context manager. Note that
`SymPyDeprecationWarning` is special and should be tested with
`warns_deprecated_sympy()` instead (see
[below](writing-tests-test-deprecated-functionality)).

The context manager should take a warning class (`warnings.warn()` uses
`UserWarning` by default), and, optionally, a regular expression that the
warning message should match as the `match` keyword argument.

```
from sympy.testing.pytest import warns
with warns(UserWarning):
    function_that_emits_a_warning()

with warns(UserWarning, match=r'warning'):
    function_that_emits_a_warning()

```

**Any test functionality that emits a warning should use `warns()`.** That
way, no warnings are actually emitted during the tests themselves. This
includes warnings coming from external libraries.

Warnings within SymPy itself should be used very sparingly. Aside from
[deprecation warnings](deprecation-policy), warnings are generally not used in
SymPy, as they may be too annoying for users, especially those who use SymPy
as a library, to be warranted.

When you do use them, you must set the `stacklevel` parameter in the warning
so that it shows the user code that called the function that emitted the
warning. If the `stacklevel` parameter is impossible to set correctly, use
`warns(test_stacklevel=False)` to disable the check in `warns` that
`stacklevel` is used properly. `warns(SymPyDeprecationWarning,
test_stacklevel=False)` must be used in place of `warns_deprecated_sympy()` if
this applies to a `SymPyDeprecationWarning`

(writing-tests-test-deprecated-functionality)=
### Test Deprecated Functionality

Deprecated functionality should be tested with the
[`sympy.testing.pytest.warns_deprecated_sympy()`](warns_deprecated_sympy)
context manager.

The only purpose of this context manager is to test that the deprecation
warning itself is functioning correctly. This should be the only place in the
test suite where deprecated functionality is called. All other tests should
use non-deprecated functionality. If it is impossible to avoid deprecated
functionality, this may be a sign that the functionality should not actually
be deprecated.

The [deprecation policy](deprecation-policy) page goes into detail about how
to add a deprecation to a function.

For example,

```
from sympy.testing.pytest import warns_deprecated_sympy
x = symbols('x')

# expr_free_symbols is deprecated
def test_deprecated_expr_free_symbols():
    with warns_deprecated_sympy():
        assert x.expr_free_symbols == {x}
```

If code is using deprecated functionality from another library, this code
should be updated. Until then, the normal
[`warns()`](writing-tests-testing-warnings) context manager should be used in
the corresponding tests to prevent the warning from being emitted.

### Testing that Something is Unchanged

The normal test style of

```py
assert function(input) == output
```

works for most types of tests. However, it doesn't work in the case where a
SymPy object should remain unchanged. Consider the following example:

```py
assert sin(pi) == 0
assert sin(pi/2) == 1
assert sin(1) == sin(1)
```

The first two tests here are fine. The test that `sin` returns the
corresponding special value for the inputs `pi` and `pi/2`. However, the last
test nominally checks that `sin(1)` doesn't return anything. But upon closer
inspection, we see that it doesn't do that at all. `sin(1)` could in fact
return anything. It could return complete nonsense or even a wrong answer like
`0`. The test would still pass, because all it is doing is checking that the
result of `sin(1)` equals the result of `sin(1)`, which it always will so long
as it always returns the same thing.

We really want to check that `sin(1)` remains unevaluated. The
`sympy.core.expr.unchanged` helper will do this.

Use it like

```
from sympy.core.expr import unchanged

def test_sin_1_unevaluated():
    assert unchanged(sin, 1)
```

This test now actually checks the correct thing. If `sin(1)` were made to
return some value, the test would fail.

### Consistency Checks

Checking a set of known inputs and outputs can only get you so far. A test
like

```
assert function(input) == expression
```

will check that `function(input)` returns `expression`, but it doesn't check
that `expression` itself is actually mathematically correct.

However, depending on what `function` is, sometimes a consistency check can be
done to check that `expression` itself is correct. This typically boils down
to "computing `expression` in two different ways". If both ways agree, there
is a pretty high chance it is correct, as it is unlikely that two completely
different methods will produce the same wrong answer.

For example, the inverse of indefinite integration is differentiation. The
tests for integrals can be checked for consistency by seeing if the derivative
of the result produces the original integrand:

```
expr = sin(x)*exp(x)
expected == exp(x)*sin(x)/2 - exp(x)*cos(x)/2

# The test for integrate()
assert integrate(expr, x) == expected
# The consistency check that the test itself is correct
assert diff(expected, x) == expr
```

The implementation for `diff` is very simple compared to `integrate`, and it
is tested separately, so this confirms the answer is correct.

Of course, one could also just confirm the answer by hand, and this is what
most tests in SymPy do. But a consistency check does not hurt, especially when
it is easy to do.

The use of consistency checks in the SymPy test suite is not, itself,
consistent. Some modules make heavy use of them, e.g., every test in the ODE
module checks itself using [`checkodesol()`](sympy.solvers.ode.checkodesol),
for instance.

When making heavy use of consistency checks, it's often a good idea to factor
out the logic into a helper function in the test file to avoid duplication.
Helper functions should start with an underscore so they aren't mistaken for
test functions by the test runner.

### Random Tests

Another way that tests can check themselves for consistency is to check the
expressions on random numerical inputs. The helper functions in
`sympy.core.random` can be used for this. See the tests in
`sympy/functions/special/` which make heavy use of this functionality.

If you add a random test, be sure to run the test multiple times to ensure
that it always passes. Random tests can be reproduced by using the random seed
printed at the top of the tests. For example

```
$./bin/test
========================================================================== test process starts ==========================================================================
executable:         /Users/aaronmeurer/anaconda3/bin/python  (3.9.13-final-0) [CPython]
architecture:       64-bit
cache:              yes
ground types:       gmpy 2.1.2
numpy:              1.22.4
random seed:        7357232
hash randomization: on (PYTHONHASHSEED=3923913114)
```

Here the random seed is `7357232`. It can be reproduced with

```
./bin/test --seed 7357232
```

In general you may need to use the same Python version and architecture as
shown in the test header to reproduce a random test failure. You may also in
some situations, need to run the tests using the exact same input
arguments (i.e., running the full test suite or running only a subset) in
order to reproduce a test that fails randomly.

(writing-tests-skip)=
### Skipping Tests

Tests can be skipped using the `sympy.testing.pytest.SKIP` decorator or using
the `sympy.testing.pytest.skip()` function. Note that tests that are skipped
because they are expected to fail should use the `@XFAIL` decorator instead
(see [below](writing-tests-xfail)). Test that are skipped because they are
too slow should use the [`@slow` decorator instead](writing-tests-slow).

Tests that are skipped unconditionally should be avoided. Such a test is
almost completely useless, as it will never be actually run. The only reason
to skip a test unconditionally is if it would otherwise be `@XFAIL` or `@slow`
but cannot use one of those decorators for some reason.

Both `@SKIP()` and `skip()` should include a message that explains why the
test is being skipped, like `skip('numpy not installed')`.

The typical usage of skipping a test is when a test depends on an [optional
dependency](optional-dependencies).

Such tests are generally written like

```
from sympy.external import import_module

# numpy will be None if NumPy is not installed
numpy = import_module('numpy')

def test_func():
    if not numpy:
       skip('numpy is not installed')

    assert func(...) == ...
```

When the test is written in this way, the test will not fail when NumPy is not
installed, which is important since NumPy is not a hard dependency of SymPy.
See also [](writing-tests-external-dependencies) below.

(writing-tests-xfail)=
### Marking Tests as Expected to Fail

Some tests in SymPy are expected to fail. They are written so that when the
functionality the check is finally implemented, a test is already written for
it.

Tests that are expected to fail are called XFAIL tests. They show up as `f`
in the test runner when they fail as expected and `X` when they pass (or
"XPASS"). A test that XPASSes should have its `@XFAIL` decorator removed so
that it becomes a normal test.

To XFAIL a test, add the `sympy.testing.pytest.XFAIL` decorator to it

```
from sympy.testing.pytest import XFAIL

@XFAIL
def test_failing_integral():
    assert integrate(sqrt(x**2 + 1/x**2), x) == x*sqrt(x**2 + x**(-2))*(sqrt(x**4 + 1) - atanh(sqrt(x**4 + 1)))/(2*sqrt(x**4 + 1))
```

Care should be taken when writing an XFAIL test so that it actually passes
when the functionality starts working. If you mistype the output, for
example, the test may never pass. For example, the integral in the above test
might start working, but return a result in a slightly different form than the
one being checked. A more robust test would be

```
from sympy.testing.pytest import XFAIL

@XFAIL
def test_failing_integral():
    # Should be x*sqrt(x**2 + x**(-2))*(sqrt(x**4 + 1) - atanh(sqrt(x**4 + 1)))/(2*sqrt(x**4 + 1))
    assert not integrate(sqrt(x**2 + 1/x**2), x).has(Integral)
```

This will cause the test to XPASS once the integral starts working, at which
time the test can be updated with the actual output of `integrate()` (which
can be compared against the expected output).

(writing-tests-slow)=
### Marking Tests as Slow

A test that is slow to run should be marked with the `@slow` decorator from
`sympy.testing.pytest.slow`. The `@slow` decorator should be used for tests
that take more than a minute to run. Tests that hang should use `@SKIP`
instead of `@slow`. The slow tests will be run automatically in a separate CI
job, but are skipped by default. You can manually run the slow tests with

```
./bin/test --slow
```

(writing-tests-external-dependencies)=
### Writing Tests with External Dependencies

When writing a test for a function that uses one of SymPy's [optional
dependencies](optional-dependencies), the test should be written in a way that
makes it so that the test does not fail when the module is not installed.

The way to do this is to use `sympy.external.import_module()`.
This will import the module if it is installed and return `None` otherwise.

`sympy.testing.pytest.skip` should be used to skip tests when the module in
question is not installed (see [](writing-tests-skip) above). This can be done
at the module level if the entire test file should be skippped, or in each
individual function.

You should also make sure the test is run in the "Optional Dependencies" CI
run. To do this, edit `bin/test_optional_dependencies.py` and make sure the
test is included (most SymPy submodules that test optional dependencies are
already included automatically).

If the optional dependency is new, add it to the list of packages that are
installed in the optional dependencies build in
`.github/workflows/runtests.yml`, and add it to the optional dependencies
document at `doc/src/contributing/dependencies.md`.

Note that it is not necessary to do any of this when using `mpmath`, as it is
already a [hard dependency](hard-dependencies) of SymPy and will always be
installed.

(writing-tests-doctests)=
## Doctests

Every public function should have a docstring, and every docstring should have
a examples. Code examples are all tested, which is why they are also sometimes
called *doctests*. The [docstring style
guide](style_guide_docstring_examples_section) has more details on how to
format examples in docstrings.

To run the doctests, use the

```
./bin/doctest
```

command. This command can also take arguments to test a specific file or
submodule, similar to `bin/test`.

Doctests should be written in a self-contained manner, with each doctest acting like
a Python session. This means that each doctest must manually import each
function used in the doctest and define the symbols used. For example

```
>>> from sympy import Function, dsolve, cos, sin
>>> from sympy.abc import x
>>> f = Function('f')
>>> dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x),
... f(x), hint='1st_exact')
Eq(x*cos(f(x)) + f(x)**3/3, C1)
```

The doctest output should look exactly as it would in a `python` session, with
`>>>` before the inputs and the outputs after. The doctester tests that the
output string matches, unlike normal tests which typically check that the
Python objects are the same with `==`. Consequently, the output needs to look
*exactly* the same as it does in a Python session.

Like tests, all doctests must pass for a change to be accepted. However, when
writing doctests, it is important to remember that **doctests should not be
thought of as tests. Rather, they are examples that happen to be tested.**

Therefore, you should always think about what will make a good, readable
example when writing doctests. Doctests do not need to extensively cover all
possible inputs, and should not include corner cases unless they are
important for users to be aware of.

Everything that is tested in a doctest should also be tested in a [normal
test](writing-tests-basics). You should always be free to remove or change a
doctest example at any time if it improves the documentation (to contrast, a
normal test should never be changed or removed, except in [certain exceptional
situations](writing-tests-updating-existing-tests)).

This also means that doctests should be written first and foremost in a way
that makes them understandable by someone reading the documentation. It can
sometimes be tempting to write a doctest in some indirect way to please the
doctester, but this should be avoided if it makes the example harder to
understand. For example

```
# BAD
>>> from sympy import sin, cos, trigsimp, symbols
>>> x = symbols('x')
>>> result = trigsimp(sin(x)*cos(x))
>>> result == sin(2*x)/2
True
```

This passes the doctest, and something along these lines would be fine a
normal test. But in a docstring example, it is much clearer to just show the actual output

```
# BETTER
>>> from sympy import sin, cos, trigsimp, symbols
>>> x = symbols('x')
>>> trigsimp(sin(x)*cos(x))
sin(2*x)/2
```

Of course, in some situations, the full output is unwieldy and showing it
would make the example harder to read, so this sort of thing may be
appropriate. Use your best judgment, keeping in mind that the
understandability of the doctest as a *documentation example* is the most
important thing. In some extreme instances, it may be preferable to just skip
a doctest (see below) rather than writing it in a convoluted way that is
difficult to read just to please the doctester.

Here are some additional tips for writing doctests:

- Long input lines can be broken into multiple lines by using `...` as a
  continuation prompt, as in the example above. The doctest runner also allows
  long outputs to be line wrapped (it ignores newlines in the output).

- Common symbol names can be imported from `sympy.abc`. Uncommon symbol names,
  or symbols that use assumptions, should be defined using `symbols`.

  ```py
  >>> from sympy.abc import x, y
  >>> x + y
  x + y
  ```

  ```py
  >>> from sympy import symbols, sqrt
  >>> a, b = symbols('a b', positive=True)
  >>> sqrt((a + b)**2)
  a + b
  ```

- If a test shows a traceback, everything between `Traceback (most recent call
  last):` and the last line with the exception message should be replaced with
  `...`, like

  ```
  >>> from sympy import Integer
  >>> Integer('a')
  Traceback (most recent call last):
  ...
  ValueError: invalid literal for int() with base 10: 'a'
  ```

- `...` is special in that whenever it appears in the output of an example,
  the doctester will allow it to replace any amount of text. It should also be
  used in instances where the exact output differs between runs, like

  ```
  >>> from sympy import simplify
  >>> simplify
  <function simplify at ...>
  ```

  Here the actual output is something like `<function simplify at
  0x10e997790>` but the `0x10e997790` is a memory address which will differ
  with every Python session.

  `...` in outputs should be used sparingly, as it prevents the doctest from
  actually checking that part of the output. It also may not be clear to the
  reader of the documentation what it is meant. Note that it's fine if the
  output of a doctest is updated to something else in the future. `...` should
  not be used in an attempt to "future-proof" doctest output.

- Another option if a doctest cannot pass is to skip it, by adding `#
  doctest:+SKIP` to the end of the input line, like

  <!-- We have to trick it into not removing # doctest:+SKIP in this example -->

  <pre class="highlight">
  >>> import random
  >>> random.random()      # doctest: +SKIP
  0.6868680200532414
  </pre>

  The `# doctest:+SKIP` part will be automatically hidden in the HTML
  documentation. When skipping a doctest, always be sure to test the output
  manually, as the doctester will not check it for you. Doctests should only be
  skipped when it is impossible to run them.

- Doctests that require a dependency to run should not be skipped with `#
  doctest: +SKIP`. Instead, use the
  [`@doctest_depends_on`](doctest_depends_on) decorator on the function to
  indicate which libraries should be installed for the doctest to run.

- Avoid using `pprint()` in doctest examples. If you need to show an
  expression in an easier to read way, you can include it inline as LaTeX math
  using dollar signs. If you absolutely must use `pprint()`, always use
  `pprint(use_unicode=False)` as the Unicode characters used for pretty
  printing do not always render correctly in the HTML documentation.

- If you want to show that something returns `None` use `print`, like

  ```py
  >>> from sympy import Symbol
  >>> x = Symbol('x', positive=True)
  >>> x.is_real
  True
  >>> x = Symbol('x', real=True)
  >>> x.is_positive # Shows nothing, because it is None
  >>> print(x.is_positive)
  None
  ```

- You can add short comments to doctests, either at the end of a line or by
  themselves after `>>>`. However, these should typically be only a few words
  long. Detailed explanations of what is happening in the doctest should go
  in the surrounding text.

- Dictionaries and sets are automatically sorted by the doctester, and any
  expressions are automatically sorted so that the order of terms is always
  printed in the same way. Usually you can just include the output as the
  doctester says it "expects" it and it will always pass subsequently.

  ```py
  >>> {'b': 1, 'a': 2}
  {'a': 2, 'b': 1}
  >>> {'b', 'a'}
  {'a', 'b'}
  >>> y + x
  x + y
  ```

(writing-tests-updating-existing-tests)=
## Updating Existing Tests

Sometimes when you change something or fix a bug, some existing tests will
fail. If this happens, you should check the test to see why it is failing. In
many cases, the test will be checking for something you didn't consider, or
your change has an unexpected side effect that broke something else. When this
happens, you may need to revisit your change. If you are unsure what to do,
you should discuss it on the issue or pull request.

If the test that fails is a [code quality test](code-quality-checks), that
usually means you just need to fix your code so that it satisfies the code
quality check (e.g., remove trailing whitespace).

Occasionally, however, it can happen that the test fails but there is nothing
wrong. In this case, the test should be updated. The most common instance of
this is a test that checks for a specific expression, but the function now
returns a different, but mathematically equivalent expression. This is
especially common with [doctests](writing-tests-doctests), since they check
not just the output expression but the way it is printed.

If a function output is mathematically equivalent, the existing test can be
updated with the new output. However, even when doing this, you should be
careful:

- Carefully check that the new output is indeed the same. Manually check
  something like if the difference of old and new expressions simplifies to 0.
  Sometimes, two expressions are equivalent for some assumptions but not for
  all, so check that the two expressions are really the same for all complex
  numbers. This can particularly happen with expressions involving square
  roots or other radicals. You can check random numbers, or use the `equals()`
  method to do this.

- If the new output is considerably more complicated than the old output, then
  it may not be a good idea to update the test, even if they are
  mathematically equivalent. Instead, you may need to adjust the change so
  that the function still returns the simpler result.

- It's not common, but it can happen that an existing test is itself
  incorrect. If a test is plain wrong, it should just be deleted, and updated.

In any case, when updating an existing test, you should always explain the
rationale for doing so in a commit message or in a pull request comment. Do
not explain the change in a code comment or documentation. Code comments and
documentation should only refer to the code as it is. Discussion of changes
belongs in the commit messages or issue tracker. Code comments that talk about
how the code used to be will only become confusing and won't actually be
relevant anymore once the change is made.

Again, the default should be to not change existing tests. The tests exist for
a reason, and changing them defeats the purpose of having them in the first
place. The exception to this rule is doctests, which are allowed to change or
be removed if they improve the documentation, as the primary purpose of
doctests is to serve as examples for users.

(code-quality-checks)=
## Code Quality Checks

SymPy has several code quality checks that must pass. The first job that is
run on the CI on a pull request is the code quality checks. If this job fails,
none of the other tests are run. Your PR may be ignored by reviewers until
they are fixed.

The code quality checks are all straightforward to fix. You can run the checks
locally using

```
./bin/test quality
```

and

```
flake8 sympy
```

This second command requires you to install `flake8`. Make sure you have the
latest version of flake8 and its dependencies `pycodestyle` and `pyflakes`
installed. Sometimes newer versions of these packages will add new checks and
if you have an older version installed you won't see the checks for them.

The `./bin/test quality` check tests for very basic code quality things. The
most common of these that will cause the test to fail is trailing whitespace.
Trailing whitespace is when a line of code has spaces at the end of it. These
spaces do nothing, and they only cause the code diff to be polluted. The best
way to handle trailing whitespace is to configure your text editor to
automatically strip trailing whitespace when you save. You can also use the
`./bin/strip_whitepace` command in the SymPy repo.

The `flake8` command will check the code for basic code errors like undefined
variables. These are restricted by the configuration in `setup.cfg` to only
check for things that are logical errors. The usual flake8 checks for cosmetic
style errors are disabled. In rare situations, a flake8 warning will be a
false positive. If this happens, add a `# noqa: <CODE>` comment to the
corresponding line, where `<CODE>` is the code for the error from
https://flake8.pycqa.org/en/latest/user/error-codes.html. For example, code
that uses `multipledispatch` will need to use

```py
@dispatch(...)
def funcname(arg1, arg2): # noqa: F811
    ...

@dispatch(...)
def funcname(arg1, arg2): # noqa: F811
    ...
```

to avoid warnings about redefining the same function multiple times.

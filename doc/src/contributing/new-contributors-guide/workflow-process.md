# SymPy Development Workflow Process


```{note}
This guide is intended for people are already familiar with contributing to
open source projects on GitHub. If you are new to GitHub, read the
[](./dev-setup.md) guide first.
```

## Pick an issue to fix

The best way to start with the main code base is to fix some existing bugs. Peruse the [Easy to fix issues] in the issue tracker and see if one interests you. If you'd like to try to fix it, then create a message in the issue saying that you'd like to work on it. If it isn't clear how to fix it, ask for suggestions on how to do it in the issue itself, on the [mailing list], or on [Gitter sympy/sympy].

SymPy's code is organized into Python packages and modules. The core code is
in the `sympy/core` directory and other packages in the sympy directory have
more specific code. For example, `sympy/printing` contains the code that
handles how SymPy objects are printed to the terminal and Jupyter.

If you are making a change that does not yet have an issue, it is not required
to open an issue first. This is only necessary if you feel you want to discuss
the change before you make a pull request, for example, if you're not sure if
something is actually a bug, or if you're not sure if a new feature is in
scope. It's also fine to just open a pull request and discuss it there.
Discussions are easier to have when there is actual code to talk about, so a
pull request is preferable if you have changes, even if they aren't fully
ready to merge yet.

## Create a new branch

The first thing to do before making changes to the code is to make a branch in git.

Remember, **never commit to `master`**. `master` should only be used to pull
in upstream changes from the main sympy/sympy repo. If you commit to `master`,
it will be difficult to pull these changes in, and it will also be difficult
if you wish to make more than one pull request at a time.

First pick a name for your branch. See [](dev-workflow-branch-names) below.
To create and checkout (that is, make it the working branch) a new branch run

```
# Pull any upstream changes from the main SymPy repo first
git checkout master
git pull

git branch <your-branch-name>
git checkout <your-branch-name>
```

The last two commands can also be combined into a single command:

```
git checkout -b <your-branch-name>
```

To view all branches, with your current branch highlighted, type:

```
git branch
```

And remember, **never type the following commands in master**: `git merge`,
`git add`, `git commit`, `git rebase`. If you made some commits to your local
master by accident, you will have to hard reset to drop the commits.

(dev-workflow-branch-names)=
### Branch names

Use a short, easy to type branch name that somehow relates to the changes.
Remember that developers who wish to try out your code will need to type your
branch name in the command line to do so.

Avoid using issue numbers in branch names as these are not easy to type (most
SymPy issue numbers are 5 digits long) and they won't really be indicative of
what the change is about without looking up the issue.

Some examples of good branch names are

```
fix-solve-bug
typo-fix
core-improvements
add-simplify-tests
```

Ultimately the branch name is not that important, so don't spend too much time
thinking about it. It's only function is to distinguish the code for this
contribution from other contributions you may make.

## Modify code

When making your fix, keep in mind there are several requirements that every
contribution should follow:

### Code Quality

SymPy contributions must have sufficient code quality to be accepted. There
are some code quality checks that will run automatically on the CI once you
create a pull request, but you can also run them locally with

```
./bin/test quality
flake8 sympy/
```

Additionally, all tests are required to pass. The CI will automatically run
the tests, but you can also [run them yourself](dev-workflow-run-tests). It is
recommended to run at least the tests that relate to the code you modified
before committing to ensure you did not make any mistakes or accidentally
break something.

Once you submit your pull request, you should check the GitHub Actions checks
once they complete to see if there are any test failures. If there are, you
will need to fix them before the pull request will be accepted.

### Add Tests

All new functionality should be tested. If you are fixing a bug, it should be
accompanied with a regression test. That is, a test that would fail before the
bug fix but now passes. Often you can use a code example from an issue as a
test case, although it is also OK to simplify such examples or to write your
own, so long as it tests the issue in question.

Tests are located alongside the code in `tests/` directories, in files named
`test_<thing>.py`. I most cases, if you modified `sympy/<submodule>/<file>.py`
then the test for the functionality will go in
`sympy/<submodule>/tests/test_<file>.py`. For example, the tests for the
functions in `sympy/simplify/sqrtdenest.py` are in
`sympy/simplify/tests/test_sqrtdenest.py`. There are some exceptions to this
rule so in general try to find where the existing tests are for a function and
add your tests alongside them.

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

### Documentation

All new methods, functions, and classes should have a docstring showing how to
use them. A docstring is a triple quoted string that goes right after the
`def` line that describes the function. Docstrings should follow the format as
outlined in the [](style_guide_docstring_guidelines).

One important thing that should be included in every docstring is examples.
Examples are also called *doctests*, because they are tested with the
`bin/doctest` script to ensure they output is correct.

Doctests need to include imports for every function used define any symbols
that are used. Users should be able to copy and paste the example inputs into
their own Python session and get the exact same outputs. `from sympy import *`
is not allowed in doctests, as this would make it unclear which functions come
from SymPy.

The [docstring style guide](style_guide_docstring_examples_section) has more
details on how to format examples in docstrings.

Keep in mind, doctests are *not* tests. Think of them as examples that happen
to be tested. Some key differences:

- write doctests to be informative; write regular tests to check for regressions and corner cases.
- doctests can be changed at any time; regular tests should not be changed.

In particular, we should be able to change or delete any doctest at any time if it makes the docstring better to understand.

Here is an example docstring with doctests (from `sympy/functions/special/delta_functions.py`).

```py
def fdiff(self, argindex=1):
    """
    Returns the first derivative of a Heaviside Function.

    Examples
    ========

    >>> from sympy import Heaviside, diff
    >>> from sympy.abc import x

    >>> Heaviside(x).fdiff()
    DiracDelta(x)

    >>> Heaviside(x**2 - 1).fdiff()
    DiracDelta(x**2 - 1)

    >>> diff(Heaviside(x)).fdiff()
    DiracDelta(x, 1)

    """
    if argindex == 1:
        return DiracDelta(self.args[0])
    else:
        raise ArgumentIndexError(self, argindex)
```

Additionally, all public function docstrings should be included in the Sphinx
API documentation. Depending on the module, this may mean you need to add an
`.. autofunction::` line to the corresponding `doc/src/modules/<module>.rst`
file. You should [build the documentation](build-the-documentation) and look
at it to ensure there are no markup errors in the rendered HTML.

If you want to write a more long-form guide or tutorial, you may include this
in the Sphinx documentation as a Markdown or RST file instead of putting it in
a docstring. While this is not a requirement for new contributions, we are
always looking to add new well-written long-form guides to our documentation.

Once you have made a pull request on GitHub, the CI will automatically build a
preview of the documentation that you can view. On the pull request page,
scroll to the bottom where the checks are, and find the link that says "Click
here to see a preview of the documentation."

(dev-workflow-run-tests)=
## Run the Tests

There are several ways of running SymPy tests but the easiest is to use the `bin/test` script, consult 'the wiki details on running tests \<<https://github.com/sympy/sympy/wiki/Running-tests>>\`\_.

The script takes a number of options and arguments and then passes them to `sympy.test(*paths, **kwargs)`. Run `bin/test --help` for all supported arguments.

Run all tests by using the command:

```bash
$ bin/test
```

To run tests for a specific file, use:

```bash
$ bin/test test_basic
```

Where `test_basic` is from file `sympy/core/basic.py`.

To run tests for modules, use:

```bash
$  bin/test /core /utilities
```

This will run tests for the `core` and `utilities` modules.

Similary, run quality tests with:

```bash
$ bin/test code_quality
```

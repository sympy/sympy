(workflow-process)=
# Development Workflow Process


```{note}
This guide is intended for people are already familiar with contributing to
open source projects on GitHub. If you are new to GitHub, read the
[](./dev-setup.md) guide first.
```

(workflow-process-checklist)=
## Checklist for Contributions

Here's a checklist of things that need to be done when making a pull request
to SymPy. These things all need to be done before a pull request is merged. It
is not necessary to do all of them before the pull request is opened at all,
but it is often a good idea to check the basic things first before opening a
pull request, or even before committing a change.

- [ ] **Make sure [code quality
  checks](workflow-process-code-quality) pass.**

  ```bash
  ./bin/test quality
  flake8 sympy/
  ```

- [ ] **[Add tests](workflow-process-add-tests).** All new functionality
  should be tested. Bug fixes should add regression tests. Tests are written
  in pytest `assert f(x) == y` style and are included in corresponding `tests`
  directories in the `sympy/` source. See the guide on [writing
  tests](writing-tests).

- [ ] **New public functions and methods should have a docstring.**

- [ ] **Docstrings should include [doctests](writing-tests-doctests).**

- [ ] **Make sure all tests pass.** [You may want to run a relevant subset of
  the test suite locally](workflow-process-run-tests) before committing (e.g.,
  `./bin/test solvers`). When you open a pull request, all tests will be run
  on CI. The CI must be all green before a PR can be merged.

- [ ] **[Write good commit messages](workflow-process-commit-messages).**

- [ ] (First time contributors only) **[Add your name to the `.mailmap`
  file](mailmap-instructions)**. The "test / authors" CI build on GitHub will
  fail if this is not done correctly.

- [ ] **Cross reference relevant issues in the pull request description.** If
  the pull request fixes an issue (i.e., the issue should be closed once the
  PR is merged), use the ["fixes #123" syntax
  ](https://docs.github.com/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue).

- [ ] **Add a comment to the original issue cross-referencing the pull
  request** for visibility. If there is not a corresponding issue, this is OK.
  It is not necessary to open an issue unless there are further improvements
  needed beyond your PR.

- [ ] **[Add a release notes
  entry](https://github.com/sympy/sympy/wiki/Writing-Release-Notes).** This
  should be done when opening the pull request, in the pull request
  description field. It can be edited at any time before the pull request is
  merged.

- [ ] **Respond to review comments.** All SymPy pull requests must be reviewed
  by someone else before they can be merged.


## Pick an issue to fix

The best way to start with the main code base is to fix some existing bugs.
Peruse the ["Easy to fix"
issues](https://github.com/sympy/sympy/issues?q=is%3Aopen+is%3Aissue+label%3A%22Easy+to+Fix%22)
in the issue tracker and see if one interests you. If you'd like to try to fix
it, then create a message in the issue saying that you'd like to work on it.
If it isn't clear how to fix it, ask for suggestions on how to do it in the
issue itself or on the [mailing
list](https://groups.google.com/g/sympy).

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

First pick a name for your branch. See [](workflow-process-branch-names) below.
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

(workflow-process-branch-names)=
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

(workflow-process-code-quality)=
### Code Quality

SymPy contributions must have sufficient code quality to be accepted. There
are some code quality checks that will run automatically on the CI once you
create a pull request, but you can also run them locally with

```
./bin/test quality
flake8 sympy/
```

Additionally, all tests are required to pass. The CI will automatically run
the tests, but you can also [run them yourself](workflow-process-run-tests). It is
recommended to run at least the tests that relate to the code you modified
before committing to ensure you did not make any mistakes or accidentally
break something.

Once you submit your pull request, you should check the GitHub Actions checks
once they complete to see if there are any test failures. If there are, you
will need to fix them before the pull request will be accepted.

(workflow-process-add-tests)=
### Add Tests

All new functionality should be tested. If you are fixing a bug, it should be
accompanied with a regression test. That is, a test that would fail before the
bug fix but now passes. Often you can use a code example from an issue as a
test case, although it is also OK to simplify such examples or to write your
own, so long as it tests the issue in question.

Tests are located alongside the code in `tests/` directories, in files named
`test_<thing>.py`. In most cases, if you modified
`sympy/<submodule>/<file>.py` then the test for the functionality will go in
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

(workflow-process-run-tests)=
## Run the Tests

There are several ways of running SymPy tests but the easiest is to use the
`bin/test` script.

The script takes a number of options and arguments. Run `bin/test --help` for
all supported arguments. Under the hood it uses `pytest`, and you can use that
directly as well if you prefer.

Run all tests by using the command:

```bash
$ ./bin/test
```

To run tests for a specific file, use:

```bash
$ ./bin/test test_basic
```

Where `test_basic` is from file `sympy/core/basic.py`.

To run tests for modules, use:

```bash
$ ./bin/test /core /utilities
```

This will run tests for the `core` and `utilities` modules.

Similary, run quality tests with:

```bash
$ ./bin/test code_quality
```

## Commit the changes

Once the changes are ready, you should commit them. You can check what files are changed:

```
git status
```

Check total changes:

```
git diff
```

If you created any new files, add them with:

```
git add new_file.py
```

You are ready to commit changes locally. A commit also contains a `commit
message` which describes it.  See the next section for guidelines on writing
good commit messages. Type:

```
git commit
```

An editor window will appear automatically in this case. In Linux, this is vim
by default. You can change what editor pops up by changing the `$EDITOR` shell
variable.

Also with the help of option `-a` you can tell the command `commit` to
automatically stage files that have been modified and deleted, but new files
you have not told git about will not be affected, e.g.,:

```
git commit -a
```

If you want to stage only part of your changes, you can use the interactive commit feature.  Just type:

```
git commit --interactive
```

and choose the changes you want in the resulting interface.

### Deleting junk files

A lot of editors can create some configuration files, binary files, or temporary files
in your SymPy directory, which should be removed before merging your commits.

Tracking down individual files can be cumbersome.

You may think of using `.gitignore`, however, editing the `.gitignore` itself
would have the agreement from the community.

Using `.git/info/exclude` would be the best, because it is only applied locally.

<https://stackoverflow.com/questions/22906851/when-would-you-use-git-info-exclude-instead-of-gitignore-to-exclude-files>

<https://docs.github.com/get-started/getting-started-with-git/ignoring-files>

(workflow-process-commit-messages)=
### Writing commit messages

The commit message has two parts: a title (first line) and the body. The two
are separated by a blank line.

Commit messages summarize what the commit does. Just as with the code, your
commit messages will become a permanent part of the project git history. So
you should put some effort into making them high quality. Commit messages are
intended for human readers, both for people who will be reviewing your code
right now, and for people who might come across your commit in the future
while researching some change in the code. Thus, include information that
helps others understand your commit here, if necessary.

Tools like `git shortlog` and the GitHub UI only show the first line of the
commit by default, so it is important to convey the most important aspects of
the commit in the first line.

- Keep the first line 71 characters or less and subsequent lines to 78
  characters or less. This allows the one-line form of the log to display the
  summary without wrapping.

- **Make sure to leave a blank line after the summary**

- Do not end the first line with a period (full stop). Subsequent lines should
  use periods.

- Provide context for the commit if possible,

  e.g. `integrals: Improved speed of heurisch()`
  instead of just `Improved speed of heurisch()`

- Reference any relevant issue numbers. You do not need to reference the pull
  request for the change itself, but issues that are fixed should be
  referenced either by `#12345` or
  `https://github.com/sympy/sympy/issues/12345`. You should also provide a
  brief summary of an issue rather than just referring to the issue number so
  that people don't have to look around for context.

- A commit won't always be seen in the context of your branch, so it is often
  helpful to give each commit some context. This is not required, though, as
  it is not hard to look at the commit metadata to see what files were
  modified or at the commit history to see the nearby related commits.

- Use plain English. Write in complete sentences.

- Describe what actually changed. Don't just write something like `Modified
  solvers.py`. People can already see what files are modified from the commit
  diff. What the message is there for is to tell them what the diff actually
  does, so they don't have to try to figure it out. Similarly, although
  relevant issues should be cross-referenced as noted above, the message
  should contain enough of a basic summary that people can understand what is
  going on without having to look up the issue. The issue can provide more
  detailed context for people who are interested.

- Try to avoid short commit messages, like "Fix", and commit messages that
  give no context, like "Found the bug". When in doubt, a longer commit
  message is probably better than a short one. Avoid using the `-m` switch to
  `git commit` to write a commit message on the command line. Rather, let it
  open your editor so you can write a longer commit message.

- Give an overview of what the commit does if it is difficult to figure out
  just from looking at the diff.

- Include other relevant information, e.g.

  - Known issues
  - A concrete example (for commits that add new features/improve performance etc.)

- Use bullet lists when suitable

- Feel free to use Unicode characters, such as output from the SymPy Unicode
  pretty printer.


#### Example of a good commit message

Here is an example commit message from the commit
[bf0e81e12a2f75711c30f0788daf4e58f72b2a41](https://github.com/sympy/sympy/commit/bf0e81e12a2f75711c30f0788daf4e58f72b2a41),
which is part of the SymPy history:

```
integrals: Improved speed of heurisch() and revised tests

Improved speed of anti-derivative candidate expansion and solution
phases using explicit domains and solve_lin_sys(). The upside of
this change is that large integrals (those that generate lots of
monomials) are now computed *much* faster. The downside is that
integrals involving Derivative() don't work anymore. I'm not sure
if they really used to work properly or it was just a coincidence
and/or bad implementation. This needs further investigation.

Example:

In [1]: from sympy.integrals.heurisch import heurisch

In [2]: f = (1 + x + x*exp(x))*(x + log(x) + exp(x) - 1)/(x + log(x) + exp(x))**2/x

In [3]: %time ratsimp(heurisch(f, x))
CPU times: user 7.27 s, sys: 0.04 s, total: 7.31 s
Wall time: 7.32 s
Out[3]:
   ⎛ 2        x                 2⋅x      x             2   ⎞
log⎝x  + 2⋅x⋅ℯ  + 2⋅x⋅log(x) + ℯ    + 2⋅ℯ ⋅log(x) + log (x)⎠          1
──────────────────────────────────────────────────────────── + ───────────────
                             2                                      x
                                                               x + ℯ  + log(x)

Previously it took 450 seconds and 4 GB of RAM to compute.
```

#### Co-Author

Occasionally, there can be multiple people working as a team for one PR,
or you have applied some suggestions from the community.

For these cases, you may use co-author feature of GitHub by adding

```
Co-authored-by: NAME NAME@EXAMPLE.COM
Co-authored-by: AUTHOR-NAME ANOTHER-NAME@EXAMPLE.COM
```

to the bottom of the commit message. See
https://docs.github.com/pull-requests/committing-changes-to-your-project/creating-and-editing-commits/creating-a-commit-with-multiple-authors.

## Make a Pull Request

Once your changes are ready for review, push them up to GitHub, and make a
pull request.

It is also OK to make a pull request before the changes are completely ready,
to get some early feedback. It is better to get feedback early before you
spend too much time on it. If your pull request is not completely ready for
merging, open it in the "DRAFT" state on GitHub. You can also add "\[WIP\]"
(which stands for "work in progress") to the beginning of the pull request
title to indicate this. Just be sure to remove the DRAFT state or \[WIP\] when
your PR is ready for final review.

### Writing pull request title and description

When you make a pull request, be sure to fill out the pull request description
template. This includes adding cross-references to any relevant issues (with
"fixes" if appropriate), and adding a release notes entry.

- **Descriptive titles are very important.** The pull request title should
  indicate what is fixed. Pull requests with undescriptive titles will tend to
  be ignored by reviewers.

  Examples of bad pull request titles are

  - "Modified solvers\.py"
  - "Fix issue #12345"

  These do indicate to reviewers what is actually changed in the pull request,
  and hence, they are likely to just pass it over instead of reviewing it.
  An example of a better pull request title is

  - "Fix a bug with solve() on transcendental functions"

- Do not put issue numbers or file names in the pull request title. Issue
  numbers should go in the description.

- Use the DRAFT status or include the prefix "\[WIP\]" in the title if you
  aren't ready to have the pull request merged and remove the status/prefix
  when you are ready.

The description is a good place to:

- Show what you have done, perhaps comparing output from master with the output after your changes

- Refer to the issue that was addressed like "#1234"; that format will
  automatically create a link to the corresponding issue or pull request, e.g.
  "This is similar to the problem in issue #1234...". This format also works
  in the discussion section of the pull request.

- Use phrases like "closes #1234" or "fixed #1234" (or similar that [follow
  the auto-close
  syntax](https://docs.github.com/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue)
  and are also [discussed
  here](https://github.blog/2013-05-14-closing-issues-via-pull-requests/)). Then
  those other issues or pull requests will be closed automatically when your
  pull request is merged. Note: this syntax does not work in the discussion of
  the pull request. See this [quick
  guide](https://github.com/sympy/sympy/wiki/Issue-PR-Autoclosing-syntax) for
  the valid and invalid syntax for automatically closing issues from pull
  requests.

- the pull request needs a release notes entry. See
  <https://github.com/sympy/sympy/wiki/Writing-Release-Notes> on how to write
  release notes in the pull request description. The SymPy Bot will check that
  your PR has release notes automatically.

It is best to just fill out the pull request template (the text that is there
when you open a pull request). If you fill out all the sections in the
template, you will have a good pull request description.

(mailmap-instructions)=
### Add your name and email address to the .mailmap file.

Every author's name and email address is stored in the
[AUTHORS](https://github.com/sympy/sympy/blob/master/AUTHORS) file, but this
file should not be edited directly. The AUTHORS file is updated automatically
when a new version of SymPy is released based on the name and email addresses
that are recorded in the commits. Every commit made with git stores the name
and email address that git is configured with (see [Configure git settings](#configure-git-settings)). The [.mailmap](https://github.com/sympy/sympy/blob/master/.mailmap)
file is used to associate the name/email recorded in the commits with an
author name and email address that will be listed in the AUTHORS file.

The first time you make a pull request you will need to add your name and email address to the .mailmap file by adding a line like

```
Joe Bloggs <joe@bloggs.com>  joeb <joe@bloggs.com>
```

This line in the .mailmap file associates the author name with the corresponding commits. The first name and email address is what will eventually go in the AUTHORS file. The second entry is what is recorded in the commit metadata. (see [Mapping user names to AUTHORS file entry](#mailmap-mapping-names))

The commit metadata name and email should exactly match the name and email that you have configured with git before making the commits (see [Configure git settings](#configure-git-settings)). The `bin/mailmap_check.py` script can check that this has been done correctly. If you have made a commit but not yet added yourself to the .mailmap file then you will see this:

```bash
$ python bin/mailmap_check.py
This author is not included in the .mailmap file:
Joe Bloggs <joe@bloggs.com>

The .mailmap file needs to be updated because there are commits with
unrecognised author/email metadata.


For instructions on updating the .mailmap file see:
https://docs.sympy.org/dev/contributing/new-contributors-guide/workflow-process.html#mailmap-instructions

The following authors will be added to the AUTHORS file at the
time of the next SymPy release.
```

This means that you should add your name and email address to the .mailmap file. If you add this at the end of the file then `git diff` will show:

```bash
$ git diff
diff --git a/.mailmap b/.mailmap
index 3af6dc1..7fa63b1 100644
--- a/.mailmap
+++ b/.mailmap
@@ -1307,3 +1307,4 @@ zsc347 <zsc347@gmail.com>
 Øyvind Jensen <jensen.oyvind@gmail.com>
 Łukasz Pankowski <lukpank@o2.pl>
+Joe Bloggs <joe@bloggs.com>
```

Now you can rerun the `bin/mailmap_check.py` script and you should see:

```bash
$ python bin/mailmap_check.py
The mailmap file was reordered

For instructions on updating the .mailmap file see:
https://docs.sympy.org/dev/contributing/new-contributors-guide/workflow-process.html#mailmap-instructions

The following authors will be added to the AUTHORS file at the
time of the next SymPy release.

Joe Bloggs <joe@bloggs.com>
```

The first line their says that the .mailmap file was "reordered". This is because the file should be in alphabetical order. The script will have moved your name into the correct position so now you can see the change as:

```bash
$ git diff
diff --git a/.mailmap b/.mailmap
index 3af6dc1..7598d94 100644
--- a/.mailmap
+++ b/.mailmap
@@ -562,6 +562,7 @@ Joannah Nanjekye <joannah.nanjekye@ibm.com> Joannah Nanjekye <jnanjekye@python.o
 Joannah Nanjekye <joannah.nanjekye@ibm.com> nanjekyejoannah <joannah.nanjekye@ibm.com>
 Joaquim Monserrat <qmonserrat@mailoo.org>
 Jochen Voss <voss@seehuhn.de>
+Joe Bloggs <joe@bloggs.com>
 Jogi Miglani <jmig5776@gmail.com> jmig5776 <jmig5776@gmail.com>
 Johan Blåbäck <johan_bluecreek@riseup.net> <johan.blaback@cea.fr>
 Johan Guzman <jguzm022@ucr.edu>
```

Now if you rerun the script you will see:

```bash
$ python bin/mailmap_check.py
No changes needed in .mailmap

The following authors will be added to the AUTHORS file at the
time of the next SymPy release.

Joe Bloggs <joe@bloggs.com>
```

The key information here is "No changes needed in .mailmap" which means that you have correctly updated the .mailmap file. You should now add and commit these changes as well:

```bash
git add .mailmap
git commit -m 'author: add Joe Bloggs to .mailmap'
```

(mailmap-mapping-names)=
### Mapping user names to AUTHORS file entry
Sometimes a commit will be made with an incorrect name or email address or an author will make multiple commits with different names and email addresses or an author wishes to use a proper name that differs from their github name. In this case a line should be added to the .mailmap file where the first name and email address is what should be recorded in the AUTHORS file and the others are the name and email address that was incorrectly used in the other commits. For example if the commit was recorded with the name `joeb` and the email address `wrong@email.com` but the AUTHORS file should show `Joe Bloggs` as above then there should be a line in the .mailmap file like:

```
Joe Bloggs <joe@bloggs.com> joeb <wrong@email.com>
```

A common reason that this can happen is if making commits with the GitHub web UI which always recorded the name as github username and the email as something like `1785690389+joeb@users.noreply.github.com`. In this case a line will need to be added to .mailmap like:

```
Joe Bloggs <joe@bloggs.com> joeb <1785690389+joeb@users.noreply.github.com>
```

Multiple lines like this can be added to the .mailmap file. They should record all of the different name and email address combinations that have been used by an author and map all of them to a single author name that will show in the AUTHORS file.

If your pull request is merged and you have not previously been added to the AUTHORS file then your name will be added at the time of the next release of SymPy.

(deprecation-policy)=
# Deprecation Policy

This page outlines SymPy's policy on doing deprecations, and describes the
steps developers should take to properly deprecate code.

## What is a deprecation?

A deprecation is a way to make backwards incompatible changes in a way that
allows users to update their code. Deprecated code continues to work as it
used to, but whenever someone uses it, it prints a warning to the screen
indicating that it will be removed in a future version of SymPy, and
indicating what the user should be using instead.

This allows users a chance to update their code without it completely
breaking. It also gives SymPy an opportunity to give users an informative
message on how to update their code, rather than making their code simply
error or start giving wrong answers.

## Try to avoid backwards incompatible changes in the first place

Backwards incompatible API changes should not be made lightly. Any backwards
compatibility break means that users will need to fix their code. Whenever you
want to make a change, you should consider whether this is worth the pain for
users. Users that have to update their code to match new APIs with every SymPy
release will become frustrated with the library, and may go seek a more stable
alternative. Consider whether the behavior you want can be done in a way that
is compatible with existing APIs. New APIs do not necessarily need to
completely supplant old ones. It is possible for old APIs to exist alongside
newer, better designed ones without removing them. For example, the newer
{ref}`solveset <solveset>` API was designed to be a superior replacement for
the older {ref}`solve <solvers>` API. But the older `solve()` function remains
intact and is still supported.

It is important to try to be cognizant of API design whenever adding new
functionality. Try to consider what a function may do in the future, and
design the API in a way that it can do so without having to make a breaking
change. For example, if you add a property to an object `A.attr`, it is
impossible to later convert that property into a method `A.attr()` so that it
can take arguments except by doing so in a backwards incompatible way. If you
are unsure about your API design for a new functionality, one option is to
mark the new functionality as explicitly private or as experimental.

With that being said, it may be decided that the API of SymPy must change in
some incompatible way. Some reasons APIs are changed can include:

- The existing API is confusing.
- There is unnecessary redundancy in the API.
- The existing API limits what is possible.

Because one of the core use-cases of SymPy is to be usable as a library, we
take API breakage very seriously. Whenever an API breakage is necessary, the
following steps should be taken:

- Discuss the API change with the community. Be sure that the improved API is
  indeed better, and worth a breakage. It is important to get API right so
  that we will not need to break the API again to "fix" it a second time.
- If it is possible, deprecate the old API. The technical steps for doing this
  are described [below](deprecation-how-to).
- Document the change so that users know how to update their code. The
  documentation that should added is described
  [below](deprecation-documentation).

## When does a change require deprecation?

When considering whether a change requires a deprecation, two things must be
considered:

- Is the change backwards incompatible?
- Is the behavior being changed public API?

A change is backwards incompatible if user code making use of it would stop
working after the change.

What counts as "public API" needs to be considered on a case-by-case basis.
The exact rules for what does and doesn't constitute public API for SymPy are
still not yet fully codified. Here are some thing that constitute public API.
*Note: these are just general guidelines. This list is not exhaustive, and
there are always exceptions to the rules.*

```{admonition} Public API
- Function names.
- Keyword argument names.
- Keyword argument default values.
- Submodule names.
- The mathematical conventions used to define a function.
```

And here are some things that generally aren't public API, and therefore don't
require deprecations to change (again, this list is only a general set of
guidelines).

```{admonition} Not Public API
- The precise form of an expression. In general, functions may be changed to
  return a different but mathematically equivalent form of the same
  expression. This includes a function returning a value which it was not able
  to compute previously.
- Positional argument names.
- Functions and methods that are private, i.e., for internal use only. Such
  things should generally be prefixed with an underscore `_`, although this
  convention is not currently universally adhered to in the SymPy codebase.
- Anything explicitly marked as "experimental".
- Changes to behavior that were mathematically incorrect previously (in
  general, bug fixes are not considered breaking changes, because despite the
  saying, bugs in SymPy are not features).
- Anything that was added before the most recent release. Code that has not
  yet made it into a release does not need to be deprecated. If you are going
  to change the API of new code, it is best to do it before a release is made
  so that no deprecations are necessary for future releases.
```

Note: both public and private API functions are included in the [reference
documentation](reference), and many functions are not included there which
should be, so this should not be used to determine whether something public or
not. Cleaning up the distinction between public and private APIs, as well as
the categorization in the reference documentation is currently an open issue
for SymPy.

If you're unsure, there is no harm in deprecating something even if it might
not actually be "public API".

## The purpose of deprecation

Deprecation has several purposes:

- To allow existing code to continue to work for a while, giving people a
  chance to upgrade SymPy without fixing all deprecation issues immediately.
- To warn users that their code will break in a future version.
- To inform users how to fix their code so that it will continue work in
  future versions.

All deprecation warnings should be something that users can remove by updating
their code. Deprecation warnings that fire unconditionally, even when using
the "correct" newer APIs, should be avoided.

This also means all deprecated code must have a completely functional
replacement. If there is no way for users to update their code, then this
means API in question is not yet read to be deprecated. The deprecation
warning should inform users of a way to change their code so that it works in
the same version of SymPy, as well as all future versions (see
[below](deprecation-documentation)).

Deprecations should always

1. Allow users to continue to use the existing APIs unchanged during the
   deprecation period (with a warning, which can be silenced with
   `warnings.filterwarnings`).
2. Allow users to always fix their code so that it stops giving the warning.
3. After users fix their code, it should continue to work after the deprecated
   code is removed.

The third point is important. We do not want the "new" method to itself cause
another API break when the deprecation period is over. Doing this would
completely defeat the purpose of doing a deprecation.

### When it is not technically possible to deprecate

In some cases, this is not technically possible to make a deprecation that
follows the above three rules. API changes of this nature should be considered
the most heavily, as they will break people's code immediately without
warning. Consideration into how easy it will be for users to support multiple
versions of SymPy, one with the change and one without, should also be taken
into account.

If you decide that the change is nonetheless worth making, there are two
options:

- Make the non-deprecatable change immediately, with no warnings. This will
  break user code.
- Warn that the code will change in the future. There won't be a way for users
  to fix their code until a version is released with the breaking change, but
  they will at least be aware that changes are coming.

Which of the two to make should decided on a case-by-case basis.

## How long should deprecations last?

Deprecations should remain intact for **at least 1 year** after the first
major release is made with the deprecation. This is only a minimum period:
deprecations are allowed to remain intact for longer than this. If a change is
especially hard for users to migrate, the deprecation period should be
lengthened. The period may also be lengthened for deprecated features that do
not impose a significant maintenance burden to keep around.

The deprecation period policy is time-based rather than release-based for a
few reasons. Firstly, SymPy does not have a regular release schedule.
Sometimes multiple releases will be made in a year, and sometimes only a
single release will be made. Being time-based assures that users have
sufficient opportunity to update their code regardless of how often releases
happen to be made.

Secondly, SymPy does not make use of a rigorous versioning scheme like
semantic versioning. The API surface of SymPy and number of contributions are
both large enough that virtually every major release has some deprecations and
backwards incompatible changes made in some submodule. Encoding this into the
version number would be virtually impossible. The development team also does
not backport changes to prior major releases, except in extreme cases. Thus a
time-based deprecation scheme is more accurate to SymPy's releasing model than
a version-based one would be.

Finally, a time-based scheme removes any temptation to "fudge" a deprecation
period down by releasing early. The best way for the developers to accelerate
the removal of deprecated functionality is to make a release containing the
deprecation as early as possible.

(deprecation-how-to)=
## How to deprecate code

All deprecations should use
{class}`sympy.utilities.exceptions.SymPyDeprecationWarning`. If a function or
method is deprecated, you can use the
{func}`sympy.utilities.decorator.deprecated` decorator. The
`deprecated_since_version` and `active_deprecations_target` flags are required. Please see the docstring of `SymPyDeprecationWarning` for more
information.

```{note}
`warns_deprecated_sympy` is only intended to be used internally by the SymPy
test suite. Users of SymPy should use the
[warnings](https://docs.python.org/3/library/warnings.html) module directly to
filter SymPy deprecation warnings. See {ref}`silencing-sympy-deprecation-warnings`.
```

Add a test for the deprecated behavior. You can use the
{func}`sympy.testing.pytest.warns_deprecated_sympy` context manager.

```py
from sympy.testing.pytest import warns_deprecated_sympy

with warns_deprecated_sympy():
    <deprecated behavior>
```

Note that this is used to assert that a warning is raised by the block as well
as ensuring that any tests in the block pass. If you want to test multiple
things and assert that each emits a warning then use separate with blocks for
each:

```py
with warns_deprecated_sympy():
    <deprecated behavior1>
with warns_deprecated_sympy():
    <deprecated behavior2>
```

This should be the only part of the codebase and test suite that uses the
deprecated behavior. Everything else should be changed to use the new,
non-deprecated behavior. The SymPy test suite is configured to fail if a
`SymPyDeprecationWarning` is issued anywhere except in a
`warns_deprecated_sympy()` block (see below). You should not use this function
or a `warnings.filterwarnings(SymPyDeprecationWarning)` anywhere except in the
test for the deprecation. This includes the documentation examples. The
documentation for a deprecated function should just have a note pointing to
the non-deprecated alternative.

If it is not possible to remove the deprecated behavior somewhere, that is a
sign that it is not ready to be deprecated yet. Consider that users may not be
able to replace the deprecated behavior for exact same reason.

(deprecation-documentation)=
## Documenting a deprecation

All deprecations should be documented. Every deprecation needs to be
documented in at three primary places:

- The `SymPyDeprecationWarning` warning text. This text is allowed to be long
  enough to describe the deprecation, but it should not be more than one
  paragraph. The primary purpose of the warning text should be *to inform
  users how to update their code*. The warning text should *not* discuss why a
  feature was deprecated or unnecessary internal technical details. This
  discussion can go in the other sections mentioned below. Do not include
  information in the message that is already part of the metadata provided to
  the keyword arguments to `SymPyDeprecationWarning`, like the version number.
  Remember that the warning text will be shown in plain-text, so do not use
  RST or Markdown markup in the text. Code blocks should be clearly delineated
  with newlines so that they are easy to read. All text in the warning message
  should be wrapped to 80 characters, except for code examples that cannot be
  wrapped.

  Always include full context of what is deprecated in the message. For
  example, write "the abc keyword to func() is deprecated" instead of just
  "the abc keyword is deprecated". That way if a user has a larger line of
  code that is using the deprecated functionality, it will be easier for them
  to see exactly which part is causing it.

- A deprecation note in the relevant docstring(s). This should use the
  [`deprecated`](https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-deprecated)
  Sphinx directive. This uses the syntax `.. deprecated:: <version>`. If the
  entire function is deprecated, this should be placed at the top of the
  docstring, right below the first line. Otherwise, if only part of a function
  is deprecated (e.g., a single keyword argument), it should be placed near
  the part of the docstring that discusses that feature, e.g., in the
  parameters list.

  The text in the deprecation should be short (no more than a paragraph),
  explaining what is deprecated and what users should use instead. If you
  want, you may use the same text here as in the `SymPyDeprecationWarning`. Be
  sure to use RST formatting, including cross-references to the new function
  if relevant, and a cross-reference to the longer description (see
  [below](deprecations-longer-description)).

  If the documentation for the feature is otherwise the same as the replaced
  feature (i.e., the deprecation is just a renaming of a function or
  argument), you may replace the rest of the documentation with a note like
  "see the documentation for \<new feature\>". Otherwise, the documentation
  for the deprecated functionality should be left intact.

  Consider these (imaginary) examples:

  ```py
  @deprecated("""\
  The simplify_this(expr) function is deprecated. Use simplify(expr)
  instead.""", deprecated_since_version="1.1",
  active_deprecations_target='simplify-this-deprecation')
  def simplify_this(expr):
      """
      Simplify ``expr``.

      .. deprecated:: 1.1

         The ``simplify_this`` function is deprecated. Use :func:`simplify`
         instead. See its documentation for more information. See
         :ref:`simplify-this-deprecation` for details.

     """
     return simplify(expr)
  ```

  ```py
  def is_this_zero(x, y=0):
      """
      Determine if x = 0.

      Parameters
      ==========

      x : Expr
        The expression to check.

      y : Expr, optional
        If provided, check if x = y.

        .. deprecated:: 1.1

           The ``y`` argument to ``is_this_zero`` is deprecated. Use
           ``is_this_zero(x - y)`` instead. See
           :ref:`is-this-zero-y-deprecated` for more details.

      """
      if y != 0:
          sympy_deprecation_warning("""\
  The y argument to is_zero() is deprecated. Use is_zero(x - y) instead.""",
              deprecated_since_version="1.1",
              active_deprecations_target='is-this-zero-y-deprecation')
      return simplify(x - y) == 0
  ```

(deprecations-longer-description)=
- A longer description of the deprecation should be added to [the page listing
  all currently active deprecations](active-deprecations) in the
  documentation (in `doc/src/explaination/active-deprecations.md`).

  This page is where you can go into more detail about the technical details
  of a deprecation. Here you should also list *why* a feature was deprecated.
  You should link to relevant issues, pull requests, and mailing list
  discussions about the deprecation. However, these discussion should be
  summarized so that users can get the basic idea of why the deprecation
  without having to read through pages of old discussions. You may also here
  give longer examples that would not fit in the `SymPyDeprecationWarning`
  message or `.. deprecated::` text.

  Every deprecation should have a cross-reference target (using
  `(target-name)=` above the section header) so that the `.. deprecated::`
  note in the relevant docstring can refer to it. This target should also be
  passed to the `active_deprecations_target` option of
  {class}`SymPyDeprecationWarning
  <sympy.utilities.exceptions.SymPyDeprecationWarning>` or {func}`@deprecated
  <sympy.utilities.decorator.deprecated>`. This will automatically put a link
  to the page in the documentation in the warning message. The target name
  should include the word "deprecation" or "deprecated" (target names are
  global in Sphinx, so the target name needs to be unique across the entire
  documentation).

  If multiple deprecations are related to one another, they can all share a
  single section on this page.

  If the deprecated function is not included in the top-level
  `sympy/__init__.py` be sure to clearly indicate which submodule the object
  is referring to.

  Note that examples here are helpful, but you generally should not use
  doctests to show the deprecated functionality, as this will itself raise the
  deprecation warning and fail the doctest. Instead you may use `# doctest:
  +SKIP`, or just show the example as a code block instead of a doctest.

  ```md
  (simplify-this-deprecation)=
  ## This is an example deprecation description

  **Note: this section is just an example. I will remove it once the *real*
  deprecations are added to this page (it is also duplicated in the deprecations
  guide.**

  The `simplify_this` function is deprecated. It has been replaced with the
  `simplify` function. Code using `simplify_this` can be fixed by replacing
  `simplfiy_this` with `simplify`. The behavior of the two functions is
  otherwise identical.

  The change was made because `simplify` is a much more Pythonic name than
  `simplify_this`.
  ```


## Release notes entry

The deprecation should be added to the "Backwards compatibility breaks and
deprecations" section of the release notes for the upcoming release. This
needs to be done manually, in addition to the change from the bot. See
https://github.com/sympy/sympy/wiki/Writing-Release-Notes#backwards-compatibility-breaks-and-deprecations

Whenever a deprecated functionality is removed entirely after its deprecation
period, this also needs to be added to the "Backwards compatibility breaks and
deprecations" section of the release notes.

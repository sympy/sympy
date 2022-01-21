# Deprecation Policy

## What is a deprecation?

Sometimes the API of SymPy must change in an incompatible way. We try to avoid this, because changing API breaks people's code, but it is often deemed worthy to do so. Some reasons APIs are changed can include:

- The existing API is confusing.
- There is unnecessary redundancy in the API.
- The existing API limits what is possible.

Because one of the core use-cases of SymPy is to be usable as a library, we take API breakage very seriously. Whenever an API breakage is necessary, the following steps should be taken:

- Discuss the API change with the community. Be sure that the improved API is indeed better, and worth a breakage. It is important to get API right so that things will not need to break again to "fix" it a second time.
- Document any API breakages in the relevant section of the release notes for the relevant release. This section is typically near the top of the release notes.
- If it is possible, deprecate the old API.

Deprecations are not always possible. Non-deprecatable API changes should be considered the most heavily, as they will break people's code immediately without warning.  Consideration into how easy it will be for users to support multiple versions of SymPy, one with the change and one without, should also be taken into account.

Note: Code that has not yet made it into a release does not need to be
deprecated. If you are going to change the API of new code, it is best to do
it before a release is made so that no deprecations are necessary for future
releases

## The purpose of deprecation

Deprecation has several purposes:

- To allow existing code to continue to work for a while, giving people a chance to upgrade SymPy without fixing all deprecation issues immediately.
- To warn users that their code will break in a future version.
- To inform users how to fix their code so that it will continue work in future versions.

To this end, a deprecation warning should be something that users can remove by fixing their code. Deprecation warnings that fire unconditionally, or which can otherwise not effectively be silenced by users other than by using a warnings filter should be avoided.

Furthermore, deprecation warnings should inform users of a way to change their code so that it works in the same version of SymPy, as well as all future versions. Features that do not have a suitable replacement should not be deprecated, as there will be no way for users to fix their code to continue to work.

## How to deprecate code

All deprecations should use {class}`sympy.utilities.exceptions.SymPyDeprecationWarning`. If a function or method is deprecated, you can use the {decorator}`sympy.utilities.decorator.deprecated` decorator. There are useful flags to this function that will generate warning messages automatically. The `issue`, `feature`, and `deprecated_since_version` flags are required (the `@deprecated` decorator sets `feature` automatically). The other flags, such as `use_instead` are useful, but you can also write your own custom message instead.  Please see the docstring of `SymPyDeprecationWarning` for more information.

Add a test for the deprecated behavior. You can use the
{func}`sympy.testing.pytest.warns_deprecated_sympy` context manager.

```py
from sympy.testing.pytest import warns_deprecated_sympy

with warns_deprecated_sympy():
    <deprecated behavior>
```

Note that this is used to assert that a warning is raised by the block as well as ensuring that any tests in the block pass. If you want to test multiple things and assert that each emits a warning then use separate with blocks for each:

```py
with warns_deprecated_sympy():
    <deprecated behavior1>
with warns_deprecated_sympy():
    <deprecated behavior2>
```

This should be the only part of the test suite that uses the deprecated behavior. Everything else should be changed to use the new, undeprecated behavior.

## Deprecation removal issues

Every deprecation should have a deprecation removal issue. This issue's number goes in the `issue` field of `SymPyDeprecationWarning`. This will generate a link automatically to the issue on GitHub when the warning is displayed.

The purpose of these issues is twofold:

1. So that we remember to eventually remove the deprecated behavior.
2. To give people more information why the feature is deprecated.

As such, the issue should give a description of why the behavior was removed
and what to use instead. This can be used to give a longer explanation than
what is suitable in the warning text itself (also, reasoning behind the
removal does not belong in the warning text, but it should be in the
issue). Please include a summary of these issues on the issue itself, rather
than relying on cross-references to long discussions, as users will want to
know why the API was changed and how to fix their code without reading through
pages of discussion, much of which may be irrelevant or contradictory to the
final decision.

Deprecation removal issues should have the "Deprecation removal" tag. These issues should not be closed until the deprecated feature is removed entirely.

Finally, make sure to add the deprecation removal issue to the "backwards
compatibility breaks and deprecations" section of the release notes (the bot
does not do this automatically).

## How long should deprecations last?

We have not decided on a formal policy for this yet. For now, it should be considered case-by-case.  Both the number of releases since the deprecation and the length of time since the first release with the deprecation should be considered.

## Release note entry

The deprecation should be added to the "Backwards compatibility breaks and deprecations" section of the release notes for the upcoming release. This needs to be done manually, in addition to the change from the bot. See https://github.com/sympy/sympy/wiki/Writing-Release-Notes#backwards-compatibility-breaks-and-deprecations

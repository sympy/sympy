(active-deprecations)=
# List of active deprecations

This pages lists all active deprecations in the SymPy codebase. See the
{ref}`deprecation-policy` page for a description of SymPy's deprecation
policy, as well as instructions on how for contributors to deprecate things.

Of particular, the deprecation policy for SymPy is for deprecations to last at
least **1 year** after the first major release that includes the deprecation.
After that period, deprecated code may be removed from SymPy, and code will
need to be updated to use the replacement feature to continue working. Prior
to this, any code using deprecated functionality will have a
`SymPyDeprecationWarning` message printed.

(simplify-this-deprecation)=
## This is an example deprecation description

```{note}
This section is just an example. I will remove it once the *real*
deprecations are added to this page (it is also duplicated as an example in
the deprecations guide).
```

The `simplify_this` function is deprecated. It has been replaced with the
`simplify` function. Code using `simplify_this` can be fixed by replacing
`simplfiy_this` with `simplify`. The behavior of the two functions is
otherwise identical.

The change was made because `simplify` is a much more Pythonic name than
`simplify_this`.

This feature has been deprecated since SymPy version 1.1.

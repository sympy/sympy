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

(silencing-sympy-deprecation-warnings)=
## Silencing SymPy Deprecation Warnings

To silence SymPy deprecation warnings, add a filter using the
[warnings](https://docs.python.org/3/library/warnings.html) module. For example:

```py
import warnings
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings(
    # replace "ignore" with "error" to make the warning raise an exception.
    # This useful if you want to test you aren't using deprecated code.
    "ignore",
    message=r"(?s).*<regex matching the warning message>",
    category=SymPyDeprecationWarning,
    module=r"<regex matching your module>"
)
```

Here `(?s).*<regex matching the warning message>` is a regular expression
matching the warning message. For example, to filter a warning about
`sympy.printing`, you might use `message=r"(?s).*sympy\.printing"`. The
leading `(?s).*` is there because the warnings module matches `message`
against the start of the warning message.

`<regex matching your module>` should be a regular expression matching your
module that uses the deprecated code. It is recommended to include this so
that you don't also silence the same warning for unrelated modules.

This same pattern may be used to instead turn `SymPyDeprecationWarning` into
an error so that you can test that you aren't using deprecated code. To do
this, just replace `"ignore"` with `"error"` in the above example. You may
also omit `message` to make this apply to all `SymPyDeprecationWarning`
warnings.

If you are using pytest, you can use the [pytest warnings filtering
capabilities](https://docs.pytest.org/en/latest/how-to/capture-warnings.html)
to either ignore `SymPyDeprecationWarning` or turn them into errors.

```{note}
The Python [`-W`
flag](https://docs.python.org/3/using/cmdline.html#cmdoption-W) and
[`PYTHONWARNINGS` environment
variable](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONWARNINGS)
will NOT work to filter SymPy deprecation warnings (see [this blog
post](https://nedbatchelder.com/blog/201810/why_warnings_is_mysterious.html)
by Ned Batchelder and [this SymPy
issue](https://github.com/sympy/sympy/issues/15130) for details on why). You
will need to either add a `warnings` filter as above or use pytest to filter
SymPy deprecation warnings.
```

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

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
[`warnings`](https://docs.python.org/3/library/warnings.html) module. For
example:

```py
import warnings
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings(
    # replace "ignore" with "error" to make the warning raise an exception.
    # This useful if you want to test you aren't using deprecated code.
    "ignore",

    # message may be omitted to filter all SymPyDeprecationWarnings
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
this, replace `"ignore"` with `"error"` in the above example. You may
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

## Version 1.10

(deprecated-traversal-functions-moved)=
### Some traversal functions have been moved

Some traversal functions have moved. Specifically, the functions

- `bottom_up`
- `interactive_traversal`
- `postorder_traversal`
- `preorder_traversal`
- `use`

have moved to different SymPy submodules.

These functions should be used from the top-level `sympy` namespace, like

```py
sympy.preorder_traversal
```

or

```py
from sympy import preorder_traversal
```

In general, end-users should use the top-level `sympy` namespace for any
functions present there. If a name is in the top-level namespace, its specific
SymPy submodule should not be relied on, as functions may move around due to
internal refactorings.

### `sympy.core.trace.Tr`

(deprecated-sympy-core-compatibility)=
### The `sympy.core.compatibility` submodule

The `sympy.core.compatibility` submodule is deprecated.

This submodule was only ever intended for internal use. Now that SymPy no
longer supports Python 2, this module is no longer necessary, and the
remaining helper functions have been moved to more convenient places in the
SymPy codebase.

Some of the functions that were in this module are available from the
top-level SymPy namespace, i.e.,

```py
from sympy import ordered, default_sort_key
```

The remaining were only intended for internal SymPy use and should not be used
by user code.

## Version 1.9

### `expr_free_symbols`

(deprecated-sympy-stats-numsamples)=
### `sympy.stats.sample(numsamples=n)`

The `numsamples` parameter to {func}`sympy.stats.sample` is deprecated.

`numsamples` makes `sample()` return a list of size `numsamples`, like

```py
>>> from sympy.stats import Die, sample
>>> X = Die('X', 6)
>>> sample(X, numsamples=3) # doctest: +SKIP
[3, 2, 3]
```

However, this functionality can be easily implemented by the user with a list
comprehension

```py
>>> [sample(X) for i in range(3)] # doctest: +SKIP
[5, 4, 3]
```

Additionally, it is redundant with the `size` parameter, which makes `sample`
return a NumPy array with the given shape.

```py
>>> sample(X, size=(3,)) # doctest: +SKIP
array([6, 6, 1])
```

Historically, `sample` was changed in SymPy 1.7 so it returned an iterator
instead of sample value. Since an iterator was returned, a numsamples
parameter was added to specify the length of the iterator.

However, this new behavior was considered confusing, as discussed in issue
[#21563](https://github.com/sympy/sympy/issues/21563), so it was reverted.
Now, `sample_iter` should be used if a iterator is needed. Consequently, the
`numsamples` parameter is no longer needed for `sample()`.

### `sympy.polys.solvers.RawMatrix`

### The `get_segments` attribute of plotting objects

### The `mdft` function in `sympy.physics.matrices`

### The private `SparseMatrix._smat` and `DenseMatrix._mat` attributes

### Non-`Expr` objects in a Matrix

### laplace_transform of a Matrix with noconds=False

## Version 1.8

(theanocode-deprecated)=
### `sympy.printing.theanocode`

[Theano](https://github.com/Theano/Theano) has been discontinued, and forked
into a new project called [Aesara](https://github.com/aesara-devs/aesara). The
`sympy.printing.theanocode` module has been renamed to
`sympy.printing.aesaracode`, and all the corresponding functions have been
renamed (e.g., `theano_code` has been renamed to `aesara_code`,
`TheanoPrinter` has been renamed to `AesaraPrinter`, and so on).

(deprecated-askhandler)=
### `sympy.assumptions.handlers.AskHandler` and related methods

`Predicate` has experienced a big design change. Previously, its handler was a
list of `AskHandler` classes and registration was done by `add_handler()` and
`remove_handler()` functions. Now, its handler is a multipledispatch instance
and registration is done by `register()` or `register_many()` methods. User
must define predicate class to introduce a new one.

Previously, handlers were defined and registered this way:

```python
class AskPrimeHandler(AskHandler):
    @staticmethod
    def Integer(expr, assumptions):
        return expr.is_prime

register_handler('prime', AskPrimeHandler)
```

It should be changed to this:

```python
# Predicate definition.
# Not needed if you are registering the handler to existing predicate.
class PrimePredicate(Predicate):
    name = 'prime'
Q.prime = PrimePredicate()

# Handler registration
@Q.prime.register(Integer)
def _(expr, assumptions):
    return expr.is_prime
```

See GitHub issue [#20209](https://github.com/sympy/sympy/issues/20209).

## Version 1.7.1

(deprecated-distribution-randomindexedsymbol)=
### Calling `sympy.stats.StochasticProcess.distribution` with
`RandomIndexedSymbol`

The `distribution` method of `sympy.stats` [stochastic
processes](sympy-stats-stochastic-processes) used to accept a
`RandomIndexedSymbol` (that is, a stochastic process indexed with a
timestamp), but should now only be called with the time stamp.

For example, if we have

```py
>>> from sympy import symbols
>>> from sympy.stats import WienerProcess
>>> W = WienerProcess('W')
>>> t = symbols('t', positive=True)
```

Previously this would work

```py
W.distribution(W(t)) # DEPRECATED
```

It should now be called like

```py
>>> W.distribution(t)
NormalDistribution(0, sqrt(t))
```

This was change was made as part of a change to store only `Basic` objects in
`sympy.stats` `.args`. See issue
[#20078](https://github.com/sympy/sympy/issues/20078) for details.

## Version 1.7

(deprecated-absorbing_probabilites)=
### `sympy.stats.DiscreteMarkovChain.absorbing_probabilites()`

The `absorbing_probabilites` method name was misspelled. The correct spelling
`absorbing_probabilities` ("absorbing probabilit*i*es") should be used
instead.

### `sympy.utilities.misc.find_executable()`

### Mutable attributes in `sympy.diffgeom`

### The `unicode` argument and attribute to `sympy.printing.pretty.stringpict.prettyForm`

(deprecated-printing-code-submodules)=
### The `sympy.printing.fcode`, `sympy.printing.ccode`, and
`sympy.printing.cxxcode` modules

The submodules `sympy.printing.ccode`, `sympy.printing.fcode`, and
`sympy.printing.cxxcode` were renamed to `sympy.printing.c`,
`sympy.printing.fortran`, and `sympy.printing.cxx`, respectively. These
modules were renamed because they conflict with the corresponding function
names. This causes issues because `from sympy.printing import ccode` can give
the function or the module, depending on whether the `ccode` submodule has
been imported yet or not. See [this
comment](https://github.com/sympy/sympy/issues/20234#issuecomment-707574283)
for a technical discussion on why this happens.

(deprecated-lambdify-arguments-set)
### Passing the arguments to `lambdify` as a `set`

Passing the function arguments to lambdify as a set is deprecated. Instead
pass them as a list or tuple. For
example, instead of

```py
lambdify({x, y}, x + 2*y) # WRONG
```

use

```py
lambdify((x, y), x + 2*y) # RIGHT
```

This is because sets are unordered. For instance, in the above example it
would be impossible for `lambidfy` to know if it was called with `{x, y}` or
`{y, x}`. Thus, when passed the arguments as a set `lambdify` would have to
guess their order, which would lead to an incorrect function if it guessed
incorrectly.

(non-expr-args-deprecated)=
### Core operators no longer accept non-Expr args

The core operator classes {class}`~.Add`, {class}`~.Mul`, and {class}`~.Pow`
can no longer be constructed directly with objects that are not subclasses of
{class}`~.Expr`.

{class}`~.Expr` is the superclass of all SymPy classes that represent scalar
numeric quantities. For example, {class}`~.sin`, {class}`~.Symbol`, and
{class}`~.Add` are all subclasses of {class}`~.Expr`. However, may objects in
SymPy are not {class}`~.Expr` because they represent some other type of
mathematical object. For example, {class}`~.Set`, {class}`~.Poly`, and
{class}`~.Boolean` are all non-`Expr`. These do not make direct mathematical
sense inside of Add, Mul, and Pow, which are designed specifically to
represent the addition, multiplication, and exponentiation of scalar complex
numbers.

Manually constructing one of these classes with such an object is possible,
but it will generally create something that will then break. For example

```py
Mul(1, Tuple(2)) # This is deprecated
```

works and creates `Tuple(2)`, but only because `Mul` is "tricked" by always
treating $1 \cdot x = x$. If instead you try

```py
Mul(2, Tuple(2)) # This is deprecated
```

it fails with an exception

```pytb
AttributeError: 'Tuple' object has no attribute 'as_coeff_Mul'
```

because it tries to call a method of `Expr` on the `Tuple` object, which does
not have all the `Expr` methods (because it is not a subclass of `Expr`).

If you want to use the `+`, `*` or `**` operation on an object, use it
directly rather than using `Mul`, `Add` or `Pow`. If functional versions of
these are desired, you can use a `lambda` or the
[`operator`](https://docs.python.org/3/library/operator.html) module.

## Version 1.6

### `sympy.utilities.tmpfiles`

### `sympy.utilities.runtests`

### `sympy.utilities.randtest`

### `sympy.utilities.pytest`

### `sympy.utilities.benchmarking`

### `sympy.testing.randtest` functions

(deprecated-poly-nonpoly-binary-operations)=
### Mixing `Poly` and non-polynomial expressions in binary operations

In previous versions of SymPy, {class}`~.Poly` was a subclass of
{class}`~.Expr`, but it has been changed to only be a subclass of
{class}`~.Basic`. This means that some things that used to work with `Poly`
are now deprecated because they are only designed to work with {class}`~.Expr`
objects.

This includes combining `Poly` with `Expr` objects using binary operations,
for example

```py
Poly(x)*sin(x) # DEPRECATED
```

To do this, either explicitly convert the non-`Poly` operand to a `Poly` using
{meth}`.Expr.as_poly` or convert the `Poly` operand to an :class:`~.Expr`
using {meth}`.Poly.as_expr`, depending on which type you want the result to
be.

### The `print_cyclic` flag of `sympy.combinatorics.Permutation`

(deprecated-integrate-poly)=
### Using `integrate` with `Poly`

In previous versions of SymPy, {class}`~.Poly` was a subclass of
{class}`~.Expr`, but it has been changed to only be a subclass of
{class}`~.Basic`. This means that some things that used to work with `Poly`
are now deprecated because they are only designed to work with {class}`~.Expr`
objects.

This includes calling {func}`~.integrate` or {class}`~.Integral` with `Poly`.

To integrate a `Poly`, use the {meth}`.Poly.integrate` method. To compute the
integral as an {class}`~.Expr` object, call the {meth}`.Poly.as_expr` method
first.

(deprecated-sympify-string-fallback)=
### The string fallback in `sympify()`

The current behavior of {func}`~.sympify` is that `sympify(expr)` tries
various methods to try to convert `expr` into a SymPy objects. If all these
methods fail, it takes `str(expr)` and tries to parse it. This string fallback
feature is deprecated. It is problematic for a few reasons:

- It can affect performance in major ways. See for instance
  https://github.com/sympy/sympy/issues/18056 and
  https://github.com/sympy/sympy/issues/15416 where it's caused up to 100x
  slowdowns. The issue is that whenever a function, and hence `sympify` is
  passed something it doesn't know how to `sympify`, for instance, a Python
  function type, it passes the string to {func}`~.parse_expr`. This is
  significantly slower than the direct conversions that happen by default.
  This occurs specifically whenever `sympify()` is used in library code
  instead of `_sympify()` (or equivalently `sympify(strict=True)`), but
  presently this is done a lot. Using `strict=True` will at some point be the
  default for all library code, but this is a [harder change to
  make](https://github.com/sympy/sympy/issues/11003).

- It can cause security issues, since strings are evaled, and objects can
  return whatever string they want in their `__repr__`. See also
  https://github.com/sympy/sympy/pull/12524.

- It really isn't very useful to begin with. Just because an object's string
  form can be parsed into a SymPy expression doesn't mean it should be parsed
  that way. This is usually correct for custom numeric types, but an object's
  repr could be anything. For instance, if the string form of an object looks
  like a valid Python identifier, it will parse as a `Symbol`.

There are plenty of ways to make custom objects work inside of
{func}`~.sympify`.

- Firstly, if an object is intended to work alongside other SymPy expressions,
  it should subclass from {class}`~.Basic` (or {class}`~.Expr`). If it does,
  {func}`~.sympify` will just return it unchanged because it will already be a
  valid SymPy object.

- For objects that you control, you can add the `_sympy_` method. The [sympify
  docstring](sympy.core.sympify.sympify) has an example of this.

- For objects that you don't control, you can add a custom converter to the
  `sympy.core.sympify.converter` dictionary. The {func}`~.sympify` docstring
  also has an example of this.

To silence this deprecation warning in all cases, you can pass `strict=True`
to `sympify()`. However, note that this will also disable some other
conversions such as conversion of strings (for converting strings to SymPy
types, you can explicitly use {func}`~.parse_expr`).

(deprecated-indefinite-integral-eq)=
### Creating an indefinite `Integral` with an `Eq` argument

Passing an [`Eq()`](sympy.core.relational.Equality) object to
{func}`~.integrate` is deprecated in the case where the integral is
indefinite. This is because if $a = b$, $\int a\,dx = \int b\,dx$ is not true
in general, due to the arbitrary constants (which `integrate` does not
include).

If you want to make an equality of indefinite integrals, use `Eq(integrate(a, x),
integrate(b, x))` explicitly.

## Version 1.5

### `TensExpr.fun_eval` and `Tensor.__call__`

### `TensorType`

### The `set_dimension`, `set_scale_factor`, `get_dimensional_expr`, and `_collect_factor_and_dimension` methods to `sympy.physics.units.Quantity`

### The `is_EmptySet` attribute of sets

### `ProductSet(iterable)`

### The `set_potential_energy` method in `sympy.physics.mechanics`

### Using a set for the condition in `ConditionSet`

### The `sympy.polys.multivariate_resultants.DixonResultant.max_degree` property

### The `sympy.polys.multivariate_resultants.DixonResultant.get_upper_degree` property

### `Eq(expr)` with the rhs defaulting to 0

### Non-tuple iterable for the `symbols` argument to `Lambda`

(deprecated-differentiate_finite-evaluate)=
### The `evaluate` flag to `differentiate_finite`

The evaluate flag to {func}`~.differentiate_finite` is deprecated.

`differentiate_finite(expr, x, evaluate=True)` expands the intermediate
derivatives before computing differences, but this usually not what you want,
as it does not satisfy the product rule.

If you really do want this behavior, you can emulate it with

```py
diff(expr, x).replace(
    lambda arg: arg.is_Derivative,
    lambda arg: arg.as_finite_difference())
```

See the discussion on issue [#17881](https://github.com/sympy/sympy/pull/17881).

## Version 1.4

### `TensorIndexType.data`, `TensorIndexType.get_matrix()`,
`TensorIndexType.__getitem__`, and `TensExpr.__pow__`

### The `clear_cache` and `clear_subproducts` keywords to
`Matrix.is_diagonalizable`

### The `rows` and `cols` keyword arguments to `Matrix.jordan_block`

## Version 1.3

### The `source()` function

### The `dimension` and `scale_factor` arguments to `sympy.physics.units.Quanitity`

### Importing `classof` and `a2idx` from `sympy.matrices.matrices`

## Version 1.2

### Dot product of non row/column vectors

### `sympy.geometry.Line3D.equation` no longer needs the `k` argument

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

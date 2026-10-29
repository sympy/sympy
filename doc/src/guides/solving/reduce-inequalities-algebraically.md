(solving-guide-inequalities)=
# Reduce Inequalities Algebraically

Use SymPy's {func}`~.reduce_inequalities` to reduce one inequality or a
system of inequalities with respect to a selected symbol. For example,
reducing $x^2 < \pi$ together with $x > 0$ gives
$0 < x < \sqrt{\pi}$.

{func}`~.reduce_inequalities` performs transformations that preserve the
meaning of the original relations under the assumptions on their symbols.
For that reason, the requested symbol is not always isolated in the result.

```{note}
{func}`~.solve` currently uses {func}`~.reduce_inequalities` when solving
inequalities. When working explicitly with inequalities, call
{func}`~.reduce_inequalities` directly.
```

## Basic Usage

### Reduce One Inequality

Pass a single inequality and the symbol of interest:

```py
>>> from sympy import pi, reduce_inequalities, symbols
>>> x = symbols("x")
>>> reduce_inequalities(x**2 <= pi, x)
(x <= sqrt(pi)) & (-sqrt(pi) <= x)
```

### Reduce a System of Inequalities

Pass several simultaneous conditions in a list or tuple:

```py
>>> from sympy import pi, reduce_inequalities, symbols
>>> x = symbols("x")
>>> reduce_inequalities([x >= 0, x**2 <= pi], x)
(0 <= x) & (x <= sqrt(pi))
```

### Reduction Can Lead to Resolution

The result of {func}`~.reduce_inequalities` is a Boolean expression. If the
inequality is always true under the assumptions on its symbols, the result can
be the Boolean atom `S.true`:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x = symbols("x", positive=True)
>>> 2*x + 3 > x*(1 + 1/x)
2*x + 3 > x*(1 + 1/x)
>>> reduce_inequalities(_)
True
```

Similarly, an incompatible system of inequalities reduces to `S.false`.

```py
>>> from sympy import reduce_inequalities, symbols
>>> x = symbols("x", real=True)
>>> reduce_inequalities([x < 0, x > pi], x)
False
```

A system should be passed as a collection of relations rather than as a
Boolean combination such as an {class}`~.And` or {class}`~.Or`.

### Specify the Symbol of Interest

It is generally best to supply the symbol to reduce for explicitly:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y", real=True)
>>> reduce_inequalities(x + y < 1, x)
x < 1 - y
```

If the symbol is omitted, the free symbols are used as symbols of interest.
This can be useful when the inequalities are independent and univariate:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y, z = symbols("x y z", real=True)
>>> reduce_inequalities([x - 3 > 0, y > 0, z - 2 < 0])
(0 < y) & (3 < x) & (z < 2)
```

When an inequality contains more than one symbol, explicitly specifying the
symbol of interest is particularly important.

## Work with the Result

The result of {func}`~.reduce_inequalities` is a symbolic Boolean expression
built from relational objects. It can therefore be inspected or transformed
like other SymPy expressions.

### Recognize Boolean Atoms

When a result resolves completely to true or false, it is returned as the
SymPy Boolean atom S.true or S.false, which prints as True or False.
Do not test its identity against the Python objects True or False:

```py
>>> from sympy import reduce_inequalities, pi, symbols, S
>>> x = symbols("x", real=True)
>>> ans = reduce_inequalities([x < 0, x > pi], x)
>>> ans == False == S.false
True
>>> ans is S.false
True
>>> ans is False  # watch out for doing this
False
```

### Convert a Univariate Result to a Set

For a result with less than two symbols, {meth}`~.Boolean.as_set` can
often provide a convenient set representation:

```py
>>> from sympy import pi, reduce_inequalities, symbols, S
>>> x = symbols("x")
>>> result = reduce_inequalities([3*x >= 1, x**2 <= pi], x); result
(1/3 <= x) & (x <= sqrt(pi))
>>> result.as_set()
Interval(1/3, sqrt(pi))
>>> S.false.as_set()
EmptySet
>>> S.true.as_set()
UniversalSet
```

Likewise, a univariate Boolean expression can be converted directly to a set:

```py
>>> from sympy.abc import x
>>> condition = (x > 0) & (x < 2)
>>> condition.as_set()
Interval.open(0, 2)
```

### Inspect the Individual Relations

Use relational atoms to obtain the individual relations contained in a result.
The {any}`canonical <sympy.core.relational.Relational.canonical>` form places
the symbol on the left when possible:

```py
>>> from sympy.core.relational import Relational
>>> from sympy.abc import x
>>> result = (3 <= x) & (pi >= x)
>>> relations = [r.canonical for r in result.atoms(Relational)]
>>> sorted([(r.lhs, r.rel_op, r.rhs) for r in relations],
...        key=lambda item: float(item[2]))
[(x, '>=', 3), (x, '<=', pi)]
```

When a result is an {class}`~.And`, its
{any}`args <sympy.core.basic.Basic.args>` are the conjuncts:

```py
>>> from sympy.abc import x
>>> result = (3 <= x) & (x <= pi)
>>> result.args
(x >= 3, x <= pi)
```

This is a property of the particular result, not a general rule for all
outputs. For example, the arguments of a single relational are its two sides.

## Behavior and Limitations

### More Than One Symbol of Interest in an Inequality

{func}`~.reduce_inequalities` can currently reduce an individual inequality
with respect to only one symbol of interest.

For example, asking it to reduce the same inequality simultaneously for both
$x$ and $y$ is not supported:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> reduce_inequalities([x + y > 1, y > 0], [x, y])
Traceback (most recent call last):
...
NotImplementedError: inequality has more than one symbol of interest.
```

A system may nevertheless contain other symbols. If only one symbol of
interest occurs in a multivariate inequality, that inequality can be reduced
with respect to that symbol. Inequalities not containing the selected symbol
will also be reduced independently:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y, z = symbols("x y z", real=True)
>>> reduce_inequalities([x + y < 1, x - z < 3, 2*y > 4], x)
(2 < y) & (x < z + 3) & (x < 1 - y)
```

### Reduction Does Not Always Isolate the Symbol

Reducing an inequality does not necessarily mean isolating its symbol of
interest. Algebraic operations are performed only when they preserve the
original relation under the known assumptions.

If an operation would require an assumption that has not been made,
{func}`~.reduce_inequalities` can return a valid partially reduced relation
rather than raise an exception.

#### Dividing by a Factor Requires Enough Information

Dividing a relation by a factor is valid only when the factor is known to be
finite and nonzero.

For the ordered relations `<`, `<=`, `>`, and `>=`, the sign of the factor
must additionally be known to be positive or negative, because division by a
negative quantity reverses the direction of the inequality.

For `Eq` and `Ne`, the sign is irrelevant, but the factor must still be known
to be finite and nonzero.

Thus, an unknown factor is not divided out:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> reduce_inequalities(x*y < 1, x)
x*y < 1
```

If its sign is known, isolation can proceed:

```py
>>> from sympy import Symbol
>>> p = Symbol("p", positive=True)  # finite, not extended_positive
>>> reduce_inequalities(x*p < 1, x)
x < 1/p
```

#### Additive Rearrangement Must Preserve Indeterminate Expressions

Moving terms from one side of a relation to the other can also be unsafe when
nonfinite values are possible.

For example, rearranging terms that can contain opposing infinities can create
or remove an indeterminate expression such as $\infty-\infty$. In that case,
SymPy keeps the relevant additive terms together rather than returning a
relation that is only conditionally equivalent to the original one:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> reduce_inequalities(x + y < 1, x)
x + y < 1
```

Declaring the symbols real makes them finite, so this particular rearrangement
is safe:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y", real=True)
>>> reduce_inequalities(x + y < 1, x)
x < 1 - y
```

Finiteness is sufficient but is not always necessary. The important question
is whether the rearrangement can create or conceal an indeterminate
combination of nonfinite terms.

### Nonlinear Polynomial Inequalities Must Be Univariate

SymPy can solve many nonlinear polynomial inequalities when they are
univariate in the symbol of interest:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x = symbols("x")
>>> reduce_inequalities([x**2 - 12 < 4, x > 0], x)
(0 < x) & (x < 4)
```

With additional symbolic parameters, the result depends both on which
transformations are safe and on what the current nonlinear inequality solver
can handle.

For symbols without a finiteness assumption, SymPy may be able to perform only
part of the reduction:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> reduce_inequalities([x**2 - 4 < y, x > 0], x)
(0 < x) & (x < oo) & (x**2 < y + 4)
```

Here `oo` denotes positive infinity.

If the symbols are known to be finite, further algebraic rearrangement is safe,
but the resulting nonlinear polynomial still contains a symbolic parameter.
The current univariate solver does not handle that case:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y", real=True)
>>> reduce_inequalities([x**2 - 4 < y, x > 0], x)
Traceback (most recent call last):
...
NotImplementedError:
The inequality, x**2 - y - 4 < 0, cannot be solved using
solve_univariate_inequality.
```

This illustrates an important distinction:

- missing assumptions can prevent further reduction without being an error;
- an unsupported mathematical form can still cause
  `NotImplementedError`.

### Periodic Inequalities Return a Representative Interval

For periodic functions, {func}`~.reduce_inequalities` returns solutions over a
representative period rather than explicitly listing all infinitely many
periodic copies.

For example:

```py
>>> from sympy import cos, reduce_inequalities
>>> from sympy.abc import x
>>> from sympy.calculus.util import periodicity
>>> reduce_inequalities(2*cos(x) < 1, x)
(pi/3 < x) & (x < 5*pi/3)
>>> periodicity(2*cos(x), x)
2*pi
```

The remaining solutions are obtained by translating the returned interval by
integer multiples of the period, here $2\pi$.

### Some Inequalities Cannot Be Reduced Symbolically

Some inequalities are mathematically well defined but are not supported by
the current symbolic inequality solvers.

For example:

```py
>>> from sympy import cos, reduce_inequalities, symbols
>>> x = symbols("x")
>>> reduce_inequalities([cos(x) - x > 0, x > 0], x)
Traceback (most recent call last):
...
NotImplementedError:
The inequality, -x + cos(x) > 0, cannot be solved using
solve_univariate_inequality.
```

This is different from a reduction that stops because an algebraic operation
would require an unknown assumption. In that case a valid, possibly
unisolated relation can be returned. `NotImplementedError` indicates that the
current solver does not have a method for the mathematical form it has
reached.

For problems that cannot be solved symbolically, numerical root-finding or
other problem-specific numerical methods may be appropriate.

## Lower-Level Inequality Functions

{func}`~.reduce_inequalities` is the top-level inequality-reduction function.
It calls lower-level routines such as
{func}`~.reduce_rational_inequalities`,
{func}`~.solve_univariate_inequality`, and the absolute-value inequality
routines as needed.

Use those lower-level functions directly when their more specialized interface
is useful.

## Report a Bug

If you find a bug in {func}`~.reduce_inequalities` or another SymPy inequality
solver, please report it on the
[SymPy issue tracker](https://github.com/sympy/sympy/issues).
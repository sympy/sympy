(solving-guide-inequalities)=
# Reduce One or a System of Inequalities for a Single Variable Algebraically

Use SymPy to reduce one or a system of inequalities for a single variable
algebraically. For example, reducing $x^2 < \pi$, $x > 0$ yields $0 < x <
\sqrt{\pi}$.

```{note}
SymPy can currently reduce for only one symbol (variable) in an inequality.
```

SymPy can reduce a system containing more than one symbol, if there is only one
symbol per inequality.

## Alternatives to Consider
- To reduce for more than one symbol in an inequality, try SciPy's
  {external:func}`~scipy.optimize.linprog`
- To reduce Boolean expressions, use {func}`as_set
  <sympy.logic.boolalg.Boolean.as_set>`

## Examples

### Reducing a System of Inequalities for a Single Variable Algebraically
{func}`~.reduce_inequalities` accepts a list or tuple of inequalities to be
reduced as a system:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> reduce_inequalities([x >= 0, x**2 <= pi], x)
(0 <= x) & (x <= sqrt(pi))
```

```{note}
While {func}`~.solve` currently accomplishes the same thing (by calling
{func}`~.reduce_inequalities` internally), that functionality may be
deprecated or removed from {func}`~.solve`. We thus recommend using
{func}`~.reduce_inequalities`.
```

{func}`~.reduce_inequalities` is the top-level inequality-reducing function
which will internally call any other lower-level [inequality-reducing
functions](../../modules/solvers/inequalities.rst) as needed.

### Reducing One Inequality for a Single Variable Algebraically
If you have only one inequality, you can optionally exclude the list construct
and simply pass {func}`~.reduce_inequalities` the inequality as an expression:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> reduce_inequalities(x**2 <= pi, x)
(x <= sqrt(pi)) & (-sqrt(pi) <= x)
```

## Guidance

### Include the Variable to Be Reduced for in the Function Call

We recommend you include the variable to be reduced for as the second argument
for {func}`~.reduce_inequalities` to ensure that it reduces for the desired
variable.

## Reduce a System of Inequalities Algebraically

You can create your inequalities, then reduce the system as a list:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> reduce_inequalities([3*x >= 1, x**2 <= pi], x)
(1/3 <= x) & (x <= sqrt(pi))
```

## Use the Result

A common way to use the result is to extract the bounds for the symbol
(variable). For example, for a solution of $0 < x < \sqrt{\pi}$, you might want
to extract $0$ and $\sqrt{\pi}$.

### Extract a List of Decomposed Relations

You can decompose a set of relations which is joined by `^` ({class}`~.Or`) or
`&` ({class}`~.And`) into individual relations using relational atoms. Using
{any}`canonical <sympy.core.relational.Relational.canonical>` will put order
each relation so the symbol is on the left, so you can take the right-hand side
{any}`rhs <sympy.core.relational.Relational.rhs>` to extract the constants:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> from sympy.core.relational import Relational
>>> x = symbols('x')
>>> eq = reduce_inequalities([3*x >= 1, x**2 <= pi], x); eq
(1/3 <= x) & (x <= sqrt(pi))
>>> relations = [(i.lhs, i.rel_op, i.rhs) for i in [i.canonical for i in eq.atoms(Relational)]]
>>> relations_sorted = sorted(relations, key=lambda x: float(x[2])) # Sorting relations just to ensure consistent list order for docstring testing
>>> relations_sorted
[(x, '>=', 1/3), (x, '<=', sqrt(pi))]
```

### Extract a Tuple of Relations

The {any}`args <sympy.core.basic.Basic.args>` (arguments) of reduced relations
are the individual relations, so you can extract the constants from the left- or
right-hand side of the `args`:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> eq = reduce_inequalities([3*x >= 1, x**2 <= pi], x); eq
(1/3 <= x) & (x <= sqrt(pi))
>>> eq.args
(1/3 <= x, x <= sqrt(pi))
>>> constants = []
>>> for arg in eq.args:
...     if arg.lhs == x:
...         constants.append(arg.rhs)
...     else:
...         constants.append(arg.lhs)
>>> constants
[1/3, sqrt(pi)]
```

## Limitations of Inequality Reduction Using SymPy

### SymPy Can Reduce for Only One Symbol of Interest Per Inequality

SymPy can currently reduce for only one symbol (variable) of interest in a given
inequality.

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> reduce_inequalities([x + y > 1, y > 0], [x, y])
Traceback (most recent call last):
...
NotImplementedError: inequality has more than one symbol of interest.
```

You can use SciPy's {external:func}`~scipy.optimize.linprog` to reduce this
system of inequalities.

SymPy can reduce for more than one symbol in a system, if there is only one
symbol of interest per inequality. For example, the following system of
inequalities has two variables, $x$ and $y$. SymPy can reduce for $x$, and gives
the constraints on $y$.

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> reduce_inequalities([x + y > 1, y > 0], x)
(0 < y) & (y < oo) & (x > 1 - y)
```

(`oo` is {class}`~.Infinity`.)

If each inequality contains only one symbol to be reduced for, SymPy can reduce
the set of inequalities for multiple symbols:

```py
>>> from sympy import reduce_inequalities, symbols
>>> x, y = symbols("x y")
>>> x_y_reduced = reduce_inequalities([x > 1, y > 0], [x, y]); x_y_reduced
(0 < y) & (1 < x) & (x < oo) & (y < oo)
```

Note that this provides no mathematical insight beyond reducing the inequalities
separately:

```py
>>> from sympy import And
>>> x_reduced = reduce_inequalities(x > 1, x); x_reduced
(1 < x) & (x < oo)
>>> y_reduced = reduce_inequalities(y > 0, y); y_reduced
(0 < y) & (y < oo)
>>> And(x_reduced, y_reduced) == x_y_reduced
True
```

so the benefit of solving such inequalities as a set maybe only convenience.

### Limitations on Types of Inequalities That SymPy Can Solve

{func}`~.reduce_inequalities` can solve a system of inequalities involving a
power of the symbol to be reduced for, or involving another symbol, but not
both:

```py
>>> from sympy import reduce_inequalities
>>> from sympy.abc import x, y
>>> reduce_inequalities([x ** 2 < 4, x > 0], x)
(0 < x) & (x < 2)
>>> reduce_inequalities([x < y, x > 0], x)
(0 < x) & (x < oo) & (x < y)
>>> reduce_inequalities([x ** 2 - y < 4, x > 0], x)
Traceback (most recent call last):
...
NotImplementedError: The inequality, -_y + x**2 - 4 < 0, cannot be solved using
solve_univariate_inequality.
```

### Not All Results Are Returned for Periodic Functions

The results returned for trigonometric inequalities are restricted in its
periodic interval. {func}`~.reduce_inequalities` tries to return just enough
solutions so that all (infinitely many) solutions can generated from the
returned solutions by adding integer multiples of the {func}`~.periodicity` of
the equation, here $2\pi$.

```py
>>> from sympy import reduce_inequalities, cos
>>> from sympy.abc import x, y
>>> from sympy.calculus.util import periodicity
>>> reduce_inequalities([2*cos(x) < 1, x > 0], x)
(0 < x) & (x < oo) & (pi/3 < x) & (x < 5*pi/3)
>>> periodicity(2*cos(x), x)
2*pi
```

## Not All Systems of Inequalities Can Be Reduced

### Systems of Inequalities Which Cannot Be Satisfied

If the system of inequalities has incompatible conditions, for example $x < 0$
and $x > \pi$, SymPy will return `False`:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> reduce_inequalities([x < 0, x > pi], x)
False
```

### Systems of Inequalities That Cannot Be Reduced Analytically

SymPy may reflect that your system of inequalities has no solutions that can be
expressed algebraically (symbolically) by returning an error such as
`NotImplementedError`:

```py
>>> from sympy import symbols, reduce_inequalities, cos
>>> x = symbols('x')
>>> reduce_inequalities([cos(x) - x > 0, x > 0], x)
Traceback (most recent call last):
...
NotImplementedError: The inequality, -x + cos(x) > 0, cannot be solved using solve_univariate_inequality.
```

so you may have to reduce your inequalities numerically instead using SciPy's
{external:func}`~scipy.optimize.linprog`.

### Inequalities Which Can Be Reduced Analytically, and SymPy Cannot Reduce
Refer to [](#limitations-of-inequality-reduction-using-sympy) above.

## Report a Bug

If you find a bug with {func}`~.diophantine`, please post the problem on the
[SymPy mailing list](https://groups.google.com/g/sympy). Until the issue is
resolved, you can use a different method listed in
[](#alternatives-to-consider).

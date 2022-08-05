# Reduce a system of inequalities of a single variable algebraically

Use SymPy to reduce a system of univariate inequalities algebraically. For
example, solving $x^2 < \pi$, $x > 0$ yields $0 < x < \sqrt{\pi}$.

```{note}
SymPy can currently reduce inequalities involving only one variable (symbol).
```

Alternatives to consider:
- For multivariate systems (more than one symbol), try SciPy's
  {external:func}`~scipy.optimize.linprog`
- To reduce Boolean expressions, use {func}`sympy.logic.boolalg.Boolean.as_set`

Here is a simple example of reducing a system of inequalities of a single
variable algebraically. {func}`~.reduce_inequalities` accepts a list or tuple of
inequalities to be reduced as a system:

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

## Guidance

### Include the variable to be reduced for in the function call

We recommend you include the variable to be reduced for as the second argument
for {func}`~.reduce_inequalities` to ensure that it reduced for the desired
variable.

## Reduce a system of inequalities algebraically

You can create your inequalities, then reduce the system as a list:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> reduce_inequalities([3*x >= 1, x**2 <= pi], x)
(1/3 <= x) & (x <= sqrt(pi))
```

## Use the solution result

A common way to use the solution result is to extract the bounds for the symbol
(variable). For example, for a solution of $0 < x < \sqrt{\pi}$, you might want
to extract $0$ and $\sqrt{\pi}$.

### Extract a list of decomposed relations

You can decompose a set of relations which is joined by `^` (or) or `&` (and)
into individual relations using relational atoms. Using {any}`canonical
<sympy.core.relational.Relational.canonical>` will put order each relation so
the symbol is on the left, so you can take the right-hand side {any}`rhs
<sympy.core.relational.Relational.lhs>` to extract the constants:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> from sympy.core.relational import Relational
>>> x = symbols('x')
>>> eq = reduce_inequalities([3*x >= 1, x**2 <= pi], x); eq
(1/3 <= x) & (x <= sqrt(pi))
>>> relations = [(i.lhs, i.rel_op, i.rhs) for i in [i.canonical for i in eq.atoms(Relational)]]
>>> # Sorting relations just to ensure consistent list order for docstring testing
>>> relations_sorted = sorted(relations, key=lambda x: float(x[2]))
>>> print(relations_sorted)
[(x, '>=', 1/3), (x, '<=', sqrt(pi))]
```

### Extract a tuple of relations

The {any}`args <sympy.core.basic.Basic.args>` of a solution set are the
individual relations, so you can extract the constants from the left- or
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

## Not all systems of inequalities can be reduced

### Systems of inequalities with no solution

If the system of inequalities has incompatible conditions, for example $x < 0$
and $x > \pi$, SymPy will return `False`:

```py
>>> from sympy import symbols, reduce_inequalities, pi
>>> x = symbols('x')
>>> reduce_inequalities([x < 0, x > pi], x)
False
```

### Systems of inequalities that cannot be reduced analytically

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

### Inequalities which can be reduced analytically, and SymPy cannot reduce

SymPy has implemented algorithms to reduce inequalities involving only one
symbol (variable), so it cannot reduce a set of inequalities for more than one
symbol:

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

### Report a problem with SymPy

If you encounter a problem with SymPy, please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
may be able to use SciPy's {external:func}`~scipy.optimize.linprog` to reduce
the system of inequalities.

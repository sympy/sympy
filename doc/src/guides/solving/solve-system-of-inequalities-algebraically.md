# Solve a system of inequalities algebraically

Use SymPy to solve a system of inequalities algebraically. For example, solving
$x^2 < \pi$, $x > 0$ yields $0 < x < \sqrt{pi}$.

Alternatives to consider:
- For systems with more than one symbol (variable), try [Wolfram
Alpha](https://www.wolframalpha.com/)

Here is a simple example of solving a system of inequalities algebraically.
solve {func}`~.solve` accepts a list or tuple of inequalities to be solved as a system:

```py
>>> from sympy import symbols, solve, pi
>>> x = symbols('x')
>>> solve([x >= 0, x**2 <= pi], x)
(0 <= x) & (x <= sqrt(pi))
```

## Guidance

### Include the variable to be solved for in the function call

We recommend you include the variable to be solved for as the second argument 
for {func}`~.solve`. While {func}`~.solve` can currently solve systems of equations with only one symbol, it is a good practice in case that capability is expanded, and because {func}`~.solve` can solve an equation with more than one symbol.

### Use exact values

If you want to preserve the exact mathematical values of symbols such as
[fractions](tutorial-gotchas-final-notes) and {any}`square roots
<sympy.core.basic.Basic.args>`, define them so that SymPy can
interpret them symbolically, for example define one-third as
`Integer(1)/Integer(3)`:

```py
>>> from sympy import symbols, solve, pi, Integer
>>> x = symbols('x')
>>> one_third = Integer(1)/Integer(3)
>>> solve((x >= one_third, x**2 <= pi), x)
(1/3 <= x) & (x <= sqrt(pi))
```

If you use the division operator between standard Python numbers, for example
`1/3`, [Python will perform the division numerically](python-vs-sympy-numbers) and pass that inexact value to SymPy,
leading to an inexact, numerical expression for all relations:

```py
>>> from sympy import symbols, solve, pi, Integer
>>> x = symbols('x')
>>> one_third = 1/3
>>> solve([x >= one_third, x**2 <= pi], x)
(0.333333333333333 <= x) & (x <= 1.77245385090552)
```

## Solve a system of inequalities algebraically

You can solve a system of inequalities algebraically in several ways.

### Enter your inequalities directly

You can create your inequalities directly, then solve the system as a list:

```py
>>> from sympy import symbols, solve, pi, Integer
>>> x = symbols('x')
>>> one_third = Integer(1)/Integer(3)
>>> solve([x >= one_third, x**2 <= pi], x)
(1/3 <= x) & (x <= sqrt(pi))
```

### Put your inequalities into the Relational class

You can create each inequality using the {class}`Relational class
<sympy.core.relational.Relational>` by specifying
the left-hand side, the right-hand side, and then a relational operator such
strict greater than (`gt` or `>`), or less than or equal to (`le` or `<=`):

```py
>>> from sympy.core import Rel
>>> from sympy import symbols, Integer, pi, solve
>>> x = symbols('x')
>>> inequality1 = Rel(x, Integer(1)/Integer(3), 'gt')
>>> inequality1
x > 1/3
>>> inequality2 = Rel(x**2, pi, '<=')
>>> inequality2
x**2 <= pi
>>> solve([inequality1, inequality2], x)
(1/3 < x) & (x <= sqrt(pi))
```

### Parse a string representing each inequality

Parse a string representing each inequality into a form that SymPy can
understand, then apply solve() to the list of parsed inequalities. This approach
is convenient if you are programmatically reading in a string for each
inequality. If you are creating the expression yourself, we [recommend against
parsing a
string](https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#strings-as-input).

```py
>>> from sympy import parse_expr, pi, solve
>>> from sympy.abc import x
>>> inequality1 = 'x >= 0'
>>> inequality2 = 'x**2 <= pi'
>>> inequalities = [inequality1, inequality2]
>>> inequalities_parsed = [parse_expr(inequality) for inequality in inequalities]
>>> inequalities_parsed
[x >= 0, x**2 <= pi]
>>> solve(inequalities_parsed, x)
(0 <= x) & (x <= sqrt(pi))
```

## Use the solution result

A common way to use the solution result is to extract the bounds for the symbol
(variable). For example, for a solution of $0 < x < \sqrt{pi}$, you might want
to extract $0$ and $\sqrt{pi}$.

### Extract relational atoms

You can decompose a set of relations which is joined by `^` (or) or `&` (and)
into individual relations using relational atoms. Using {any}`canonical
<sympy.core.relational.Relational.canonical>` will put order each relation so
the symbol is on the left, so you can take the right-hand side {any}`rhs
<sympy.core.relational.Relational.lhs>` to extract the constants:

```py
>>> from sympy import symbols, solve, Integer, pi
>>> from sympy.core.relational import Relational
>>> x = symbols('x')
>>> one_third = Integer(1)/Integer(3)
>>> eq = solve([x >= one_third, x**2 <= pi], x)
>>> relations = [(i.lhs, i.rel_op, i.rhs) for i in [i.canonical for i in eq.atoms(Relational)]]
>>> # Sorting relations just to ensure consistent list order for docstring testing
>>> relations_sorted = sorted(relations, key=lambda x: float(x[2]))
>>> print(relations_sorted)
[(x, '>=', 1/3), (x, '<=', sqrt(pi))]
```

### Extract arguments

The {any}`args <sympy.core.basic.Basic.args>` of a solution set are the
individual relations, so you can extract the constants from the left- or
right-hand side of the `args`:

```py
>>> from sympy import symbols, solve, Integer, pi
>>> x = symbols('x')
>>> eq = solve([x >= Integer(1)/Integer(3), x**2 <= pi], x)
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

## *Tradeoffs (speed vs. accuracy, etc.) for function*

### *Tradeoff 1*

*Speed-up option 1 content*

### *Speed-up option 2*

*Speed-up option 2 content*

## Not all systems of inequalities can be solved

### Systems of inequalities with no solution

If the system of inequalities has incompatible conditions, for example $x < 0$
and $x > \pi$, SymPy will return `False`:

```py
>>> from sympy import symbols, solve, pi
>>> x = symbols('x')
>>> solve([x < 0, x > pi])
False
```

### Equations with no analytical solution

SymPy may reflect that your equation has no solutions that can be expressed
algebraically (symbolically) by returning an error such as
`NotImplementedError`:

```py
>>> from sympy import symbols, solve, cos
>>> x = symbols('x')
>>> solve([cos(x) - x > 0, x > 0], x)
Traceback (most recent call last):
    ...
NotImplementedError: The inequality, -x + cos(x) > 0, cannot be solved using solve_univariate_inequality.
```

so you may have to solve your equation numerically instead using [Wolfram
Alpha](https://www.wolframalpha.com/input?i2d=true&i=solve%5C%2840%29cos%5C%2840%29x%5C%2841%29+-+x+%3E+0+and+++x%3E0%5C%2844%29x%5C%2841%29).

### Equations which have an analytical solution, and SymPy cannot solve

SymPy has implemented algorithms to solve inequalities involving only one symbol
(variable), so it cannot solve a set of inequalities involving more than one
symbol:

```py
>>> from sympy import solve, symbols
>>> x, y = symbols('x y')
>>> from sympy.abc import x, y
>>> solve([x**2 < y, x > 0], x)
Traceback (most recent call last):
    ...
NotImplementedError: inequality has more than one symbol of interest.
```

You can use [Wolfram
Alpha](https://www.wolframalpha.com/input?i2d=true&i=solve%5C%2840%29Power%5Bx%2C2%5D+%3C+y+and+++x%3E0%5C%2844%29x%5C%2841%29)
to solve this problem.

If you encounter a problem with SymPy, please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
may be able to use [Wolfram
Alpha](https://www.wolframalpha.com/input?i2d=true&i=solve%5C%2840%29Power%5Bx%2C2%5D+%3C+y+and+++x%3E0%5C%2844%29x%5C%2841%29)
to solve the problem.

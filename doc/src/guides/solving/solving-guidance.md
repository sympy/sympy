# Solving Guidance

These guidelines apply to many types of solving.

## Numeric solutions

### Numeric solutions by choice
Solving functions such as {func}`~.solve` and {func}`~.solveset` will not try to
find a numeric solution, only a mathematically-exact symbolic solution. So if
you want a numeric solution, consider {func}`~.nsolve`.

In some situations, even though a closed-form solution is available, it may be
too cumbersome to be desirable. In that case, you can use {func}`evalf n
<sympy.core.evalf>` instead if a numerical solution is acceptable:

```py
>>> from sympy import symbols, solve
>>> x = symbols('x')
>>> solutions = solve(x**4 + 10*x**2 + x + 1, x, dict=True); solutions
[{x: -sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2 - sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) + 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2}, {x: sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2 - sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) - 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2}, {x: sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) - 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2 + sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2}, {x: sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) + 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2 - sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2}]
>>> for solution in solutions:
...     print(solution[x].n())
-0.0509758447494279 + 0.313552108895239*I
0.0509758447494279 + 3.14751999969868*I
0.0509758447494279 - 3.14751999969868*I
-0.0509758447494279 - 0.313552108895239*I
```

### Equations with no analytical solution

The vast majority of arbitrary nonlinear equations are not analytically
solvable. The classes of equations that are solvable are basically:
1. Linear equations
2. Polynomials, except where limited by the Abel-Ruffini theorem
3. Equations that can be solved by inverting some transcendental functions
4. Problems that can be transformed into the cases above (e.g., by turning
trigonometric functions into polynomials)
5. A few other special cases that can be solved with something like the
{class}`Lambert W function <sympy.functions.elementary.exponential.LambertW>` 

SymPy may reflect that your equation has no solutions that can be expressed
algebraically (symbolically) by returning an error such as
`NotImplementedError`:

```py
>>> from sympy import solve, cos
>>> from sympy.abc import x
>>> solve(cos(x) - x, x)
Traceback (most recent call last):
  ...
NotImplementedError: multiple generators [x, cos(x)]
No algorithms are implemented to solve equation -x + cos(x)
```

so you may have to {func}`solve your equation numerically
<sympy.solvers.solvers.nsolve>` instead, for example

```py
>>> from sympy import nsolve, cos
>>> from sympy.abc import x
>>> nsolve(cos(x) - x, x, 2)
0.739085133215161
```

If you receive non-closed-form solutions such as {class}`CRootOf
<sympy.polys.rootoftools.ComplexRootOf>`, you can evaluate them numerically
using {func}`evalf n <sympy.core.evalf.EvalfMixin.n>`:

```py
>>> from sympy import solve
>>> from sympy.abc import x
>>> solutions = solve(x**5 - x - 1, x, dict=True)
>>> print(solutions)
[{x: CRootOf(x**5 - x - 1, 0)}, {x: CRootOf(x**5 - x - 1, 1)}, {x: CRootOf(x**5 - x - 1, 2)}, {x: CRootOf(x**5 - x - 1, 3)}, {x: CRootOf(x**5 - x - 1, 4)}]
>>> [solution[x].n(3) for solution in solutions]
[1.17, -0.765 - 0.352*I, -0.765 + 0.352*I, 0.181 - 1.08*I, 0.181 + 1.08*I]
```

({class}`CRootOf <sympy.polys.rootoftools.ComplexRootOf>` represents an indexed
complex root of a polynomial.)

## Use exact values

If you want to preserve the exact mathematical values of symbols such as
[fractions](tutorial-gotchas-final-notes) and [square
roots](symbolic-computation), define them so that SymPy can interpret them
symbolically, for example use SymPy's {class}`.Pi`:

```py
>>> from sympy import symbols, solve, pi
>>> x = symbols('x')
>>> solve(x**2 - pi, x, dict=True)
[{x: -sqrt(pi)}, {x: sqrt(pi)}]
```

If you use the standard Python math version of $\pi$, Python will pass that
inexact value to SymPy, leading to an inexact, numerical solution:

```py
>>> from sympy import symbols, solve
>>> from math import pi
>>> x = symbols('x')
>>> solve(x**2 - pi, x, dict=True)
[{x: -1.77245385090552}, {x: 1.77245385090552}]
```

In certain cases, using an inexact value will prevent SymPy from finding a
result. For example, this exact equation can be solved:

```py
>>> from sympy import symbols, solve, sqrt
>>> x = symbols('x')
>>> eq = x**sqrt(2) - 2
>>> solve(eq, x, dict=True)
[{x: 2**(sqrt(2)/2)}]
```

but if you use the inexact equation `eq = x**1.4142135623730951 - 2`, SymPy will
not return a result despite attempting for a long time. 

(speed_up_solve)=
## Options that can speed up {func}`~.solve`

### Include solutions making any denominator zero by using `check=False`

Normally, {func}`~.solve` checks whether any solutions make any denominator
zero, and automatically excludes them. If you want to include those solutions,
and speed up {func}`~.solve` (at the risk of obtaining invalid solutions), set
`check=False`:

```py
>>> from sympy import Symbol, sin, solve
>>> x = Symbol("x")
>>> solve(sin(x)/x)  # 0 is excluded
[pi]
>>> solve(sin(x)/x, check=False)  # 0 is not excluded
[0, pi]
```

### Do not simplify solutions by using `simplify=False`

Normally, {func}`~.solve` simplifies all but polynomials of order 3 or greater
before returning them and (if `check` is not False) uses the general
{func}`simplify <sympy.simplify.simplify.simplify>` function on the solutions
and the expression obtained when they are substituted into the function which
should be zero. If you do not want the solutions simplified, and want to speed
up {func}`~.solve`, use `simplify=False`.

```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> expr = x**2 - (y**5 - 3*y**3 + y**2 - 3)
>>> solve(expr, x, dict=True)
[{x: -sqrt(y**5 - 3*y**3 + y**2 - 3)}, {x: sqrt(y**5 - 3*y**3 + y**2 - 3)}]
>>> solve(expr, x, dict=True, simplify=False)
[{x: -sqrt((y + 1)*(y**2 - 3)*(y**2 - y + 1))}, {x: sqrt((y + 1)*(y**2 - 3)*(y**2 - y + 1))}]
```

## Parse a string representing the equation

If you are creating the expression yourself, we [recommend against using parsing
a string](
https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#strings-as-input).
But if you are programmatically reading in a string, this approach is
convenient.

You can parse a string representing the equation into a form that SymPy can
understand (for example, {any}`Eq <sympy.core.relational.Eq>` form), then apply
solve the parsed expression. Parsing an equation from a string requires you to
use {func}`transformations <sympy.parsing.sympy_parser.parse_expr>` for SymPy to
- interpret equals signs
- create symbols from your variables
-  use more mathematical (rather than standard Python) notation, for example the
exponent operator can be parsed from `^` rather than having to use Python's
`**`.

To extract the solutions, you can iterate through the list of dictionaries:  
    
```py
>>> from sympy import parse_expr, solve, solveset
>>> from sympy.abc import x
>>> expr = "x^2 = y"
>>> parsed = parse_expr(expr, transformations="all")
>>> print(parsed)
Eq(x**2, y)
>>> solutions = solve(parsed, x, dict=True)
>>> print(solutions)
[{x: -sqrt(y)}, {x: sqrt(y)}]
>>> for solution in solutions:
...     for key, val in solution.items():
...         print(val)
-sqrt(y)
sqrt(y)
>>> solveset(parsed, x)
{-sqrt(y), sqrt(y)}
```

If you already have the equation in {any}`Eq <sympy.core.relational.Eq>` form,
you can parse that string:

```py
>>> from sympy import parse_expr, solve, solveset
>>> from sympy.abc import x
>>> expr = "Eq(x^2, y)"
>>> parsed = parse_expr(expr, transformations="all")
>>> print(parsed)
Eq(x**2, y)
>>> solutions = solve(parsed, x, dict=True)
>>> print(solutions)
[{x: -sqrt(y)}, {x: sqrt(y)}]
>>> for solution in solutions:
...     for key, val in solution.items():
...         print(val)
-sqrt(y)
sqrt(y)
>>> solutions_set = solveset(parsed, x)
>>> print(solutions_set)
{-sqrt(y), sqrt(y)}
```

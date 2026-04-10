# Solving Guidance

These guidelines apply to many types of solving.

## Numeric Solutions

### Equations With no Closed-Form Solution

The vast majority of arbitrary nonlinear equations have no closed-form solution.
The classes of equations that are solvable are basically:
1. Linear equations
2. Polynomials, except where limited by the [Abel-Ruffini
   theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem) (learn
   more about solving polynomials using a {class}`~.GroebnerBasis`)
3. Equations that can be solved by inverting some transcendental functions
4. Problems that can be transformed into the cases above (e.g., by turning
trigonometric functions into polynomials)
5. A few other special cases that can be solved with something like the
{class}`Lambert W function <sympy.functions.elementary.exponential.LambertW>`
6. Equations that you can {func}`~sympy.polys.polytools.decompose` via any of
   the above

SymPy may reflect that your equation has no solutions that can be expressed
algebraically (symbolically), or that SymPy lacks an algorithm to find a
closed-form solution that does exist, by returning an error such as
`NotImplementedError`:

```py
>>> from sympy import solve, cos
>>> from sympy.abc import x
>>> solve(cos(x) - x, x, dict=True)
Traceback (most recent call last):
  ...
NotImplementedError: multiple generators [x, cos(x)]
No algorithms are implemented to solve equation -x + cos(x)
```

so you may have to solve your equation numerically instead, for example using
{func}`~.nsolve`

```py
>>> from sympy import nsolve, cos
>>> from sympy.abc import x
>>> nsolve(cos(x) - x, x, 2)
0.739085133215161
```

If you receive non-closed-form solutions such as {class}`CRootOf()
<sympy.polys.rootoftools.ComplexRootOf>` (which represents an indexed complex
root of a polynomial), you can evaluate them numerically using
{func}`~sympy.core.evalf`:

```py
>>> from sympy import solve
>>> from sympy.abc import x
>>> solutions = solve(x**5 - x - 1, x, dict=True)
>>> solutions
[{x: CRootOf(x**5 - x - 1, 0)}, {x: CRootOf(x**5 - x - 1, 1)}, {x: CRootOf(x**5 - x - 1, 2)}, {x: CRootOf(x**5 - x - 1, 3)}, {x: CRootOf(x**5 - x - 1, 4)}]
>>> [solution[x].evalf(3) for solution in solutions]
[1.17, -0.765 - 0.352*I, -0.765 + 0.352*I, 0.181 - 1.08*I, 0.181 + 1.08*I]
```

### When You Might Prefer a Numeric Solution
Even if your problem has a closed-form solution, you might prefer a numeric
solution.

Solving functions such as {func}`~.solve` and {func}`~.solveset` will not try to
find a numeric solution, only a mathematically-exact symbolic solution. So if
you want a numeric solution, consider {func}`~.nsolve`.

In some situations, even though a closed-form solution is available, it may be
too cumbersome to be desirable. In that case, you can use
{func}`~sympy.core.evalf` instead if a numerical solution is acceptable. For
example, the following solution set contains more than 40 terms total when
expressed exactly (scroll horizontally in the code block below if you want to
view them all), compared to eight when expressed numerically:

```py
>>> from sympy import symbols, solve
>>> x = symbols('x')
>>> solutions = solve(x**4 + 10*x**2 + x + 1, x, dict=True)
>>> solutions
[{x: -sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2 - sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) + 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2}, {x: sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2 - sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) - 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2}, {x: sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) - 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2 + sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2}, {x: sqrt(-40/3 - 2*(1307/432 + sqrt(434607)*I/144)**(1/3) + 2/sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3)) - 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)))/2 - sqrt(-20/3 + 56/(9*(1307/432 + sqrt(434607)*I/144)**(1/3)) + 2*(1307/432 + sqrt(434607)*I/144)**(1/3))/2}]
>>> for solution in solutions:
...     solution[x].evalf()
-0.0509758447494279 + 0.313552108895239*I
0.0509758447494279 + 3.14751999969868*I
0.0509758447494279 - 3.14751999969868*I
-0.0509758447494279 - 0.313552108895239*I
```

In other situations, even if the exact solution has few terms, you may want a
numeric solution so you know its approximate numerical value. For example, it
may be difficult to estimate that $\sqrt{2} e^{\pi}/2$ is approximately $16$:

```py
>>> from sympy import pi, sqrt, exp, solve, evalf
>>> shorter = solve(sqrt(2)*x - exp(pi), x, dict=True)
>>> shorter
[{x: sqrt(2)*exp(pi)/2}]
>>> [solution[x].evalf(3) for solution in shorter]
[16.4]
```

## Use Exact Values

If you want to preserve the exact mathematical values of symbols such as
transcendental numbers and [square roots](symbolic-computation), define them so
that SymPy can interpret them symbolically, for example use SymPy's
{class}`.Pi`:

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

To use exact values for numbers such as $6.2$ or $1/2$, refer to
[](python-vs-sympy-numbers).

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

## Include the Variable to be Solved for in the Function Call

We recommend you include the variable to be solved for as the second argument
for solving functions including {func}`~.solve` and {func}`~.solveset`. While
this is optional for univariate equations, it is a good practice because it
ensures SymPy will solve for the desired symbol. For example, you might be
interested in a solution for $x$, but SymPy solves for $y$:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x**2 - y, dict=True)
[{y: x**2}]
```

Specifying the variable to solve for ensures that SymPy solves for it:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x**2 - y, x, dict=True)
[{x: -sqrt(y)}, {x: sqrt(y)}]
```

(ensure-consistent-formatting-from-solve)=
## Ensure Consistent Formatting From {func}`~.solve`

{func}`~.solve` produces a variety of output as explained in
{ref}`solve_output`. Using `dict=True` will give a consistent output format
which is especially important when extracting information about the solution
programmatically.

To extract the solutions, you can iterate through the list of dictionaries:

```py
>>> from sympy import parse_expr, solve, solveset
>>> from sympy.abc import x
>>> expr = "x^2 = y"
>>> parsed = parse_expr(expr, transformations="all")
>>> parsed
Eq(x**2, y)
>>> solutions = solve(parsed, x, dict=True)
>>> [solution[x] for solution in solutions]
[-sqrt(y), sqrt(y)]
>>> solveset(parsed, x)
{-sqrt(y), sqrt(y)}
```

(options-that-can-speed-up-solve)=
## Options That Can Speed up {func}`~.solve`

### Include Solutions Making Any Denominator Zero

Normally, {func}`~.solve` checks whether any solutions make any denominator
zero, and automatically excludes them. If you want to include those solutions,
and speed up {func}`~.solve` (at the risk of obtaining invalid solutions), set
`check=False`:

```py
>>> from sympy import Symbol, sin, solve
>>> x = Symbol("x")
>>> solve(sin(x)/x, x, dict=True) # 0 is excluded
[{x: pi}]
>>> solve(sin(x)/x, x, dict=True, check=False) # 0 is not excluded
[{x: 0}, {x: pi}]
```

### Do Not Simplify Solutions

Normally, {func}`~.solve` simplifies many results before returning them and (if
`check` is not False) uses the general {func}`~sympy.simplify.simplify.simplify`
function on the solutions and the expression obtained when they are substituted
into the function which should be zero. If you do not want the solutions
simplified, and want to speed up {func}`~.solve`, use `simplify=False`.

```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> expr = x**2 - (y**5 - 3*y**3 + y**2 - 3)
>>> solve(expr, x, dict=True)
[{x: -sqrt(y**5 - 3*y**3 + y**2 - 3)}, {x: sqrt(y**5 - 3*y**3 + y**2 - 3)}]
>>> solve(expr, x, dict=True, simplify=False)
[{x: -sqrt((y + 1)*(y**2 - 3)*(y**2 - y + 1))}, {x: sqrt((y + 1)*(y**2 - 3)*(y**2 - y + 1))}]
```

## Parse a String Representing the Equation

If you are creating the expression yourself, we advise [against using string
parsing to create expressions](
https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#user-content-strings-as-input).
But if you are programmatically reading in a string, this approach is
convenient.

You can parse a string representing the equation into a form that SymPy can
understand (for example, {class}`~sympy.core.relational.Eq` form), then solve
the parsed expression. Parsing an equation from a string requires you to use
{func}`transformations <sympy.parsing.sympy_parser.parse_expr>` for SymPy to
- interpret equals signs
- create symbols from your variables
-  use more mathematical (rather than standard Python) notation, for example the
exponent operator can be parsed from `^` rather than having to use Python's
`**`.

If you already have the equation in {class}`~sympy.core.relational.Eq`
(equation) form, you can parse that string:

```py
>>> from sympy import parse_expr, solve, solveset
>>> from sympy.abc import x
>>> expr = "Eq(x^2, y)"
>>> parsed = parse_expr(expr, transformations="all")
>>> parsed
Eq(x**2, y)
```

SymPy can also parse [LaTeX](https://www.latex-project.org/) into expressions
using {func}`~.parse_latex`.

## Report a Bug

If you find a bug with these commands, please post the problem on the [SymPy mailing
list](https://groups.google.com/g/sympy).

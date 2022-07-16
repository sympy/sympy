# Solving Guidance

Here are guidelines that apply to many types of solving.

## Guidance

### Equations with no analytical solution

The vast majority of arbitrary nonlinear equations are not analytically 
solvable. The classes of equations that are solvable are basically:
1. Linear equations
2. Polynomials, except where limited by the Abel-Ruffini theorem
3. Equations that can be solved by inverting some transcendental functions
4. Problems that can be transformed into the cases above 
(e.g., by turning trigonometric functions into polynomials)
5. A few other special cases that can be solved with something like the
Lambert W function

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

If you receive non-closed-form solutions, you can evaluate them numerically 
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

({class}`CRootOf <sympy.polys.rootoftools.ComplexRootOf>` represents an 
indexed complex root of a polynomial.)

### Use exact values

If you want to preserve the exact mathematical values of symbols such as
[fractions](tutorial-gotchas-final-notes) and [square
roots](symbolic-computation), define them so that SymPy can interpret them
symbolically, for example define one-third as a {class}`rational object <sympy.core.numbers.Rational>` using `Rational(1, 3)`:

```py
>>> from sympy import symbols, solve, pi, Rational
>>> x = symbols('x')
>>> one_third = Rational(1, 3)
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

### How to parse a string representing the equation

If you are creating the expression yourself, we 
[recommend against using parsing a string](
https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#strings-as-input). But if you are programmatically reading in a string, this approach is convenient.

You can parse a string representing the equation into a form that SymPy can
understand (for example, `Eq` form), then apply solve the parsed expression.
Parsing an equation from a string requires you to use {func}`transformations
<sympy.parsing.sympy_parser.parse_expr>` for SymPy to
- interpret equals signs
- create symbols from your variables
-  use more mathematical (rather than standard Python) notation, 
for example the exponent operator can be parsed from `^` rather than having 
to use Python's `**`.

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

If you already have the equation in `Eq` form, you can parse that string:

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

## *Tradeoffs (speed vs. accuracy, etc.) for function*

### *Tradeoff 1*

*Speed-up option 1 content*

### *Speed-up option 2*

*Speed-up option 2 content*

## Not all equations can be solved

### Equations with no solution

*Equations with no solution content*

### Equations with no analytical solution

*Equations with no analytical solution content*

### Equations which have an analytical solution, and SymPy cannot solve

*Equations which have an analytical solution, and SymPy cannot solve content*

Please post the problem on the 
[mailing list](https://groups.google.com/g/sympy), or open an issue on 
[SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue 
is resolved, you can *workaround*.

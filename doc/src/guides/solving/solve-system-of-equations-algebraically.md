# Solve a System of Equations Algebraically

Use SymPy to solve a system of equations algebraically. For example, $x^2 + y =
2z, y = -4z$ yields $\{(x = -\sqrt{6z}, y = -4z),$ ${(x = \sqrt{6z}, y =
-4z)\}}$.

Alternatives to consider:
- Some systems of equations cannot be solved algebraically (either at all or by
SymPy), so you may have to {func}`solve your system of equations numerically
using nsolve() <sympy.solvers.solvers.nsolve>` instead.

Here is an example of solving a system of equations algebraically:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x**2 + y - 2*z, y + 4*z], x, y, dict=True)
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
```

## Guidance

*Just link to solving guidance page?*

### *Guidance 1*

*Guidance 1 content*

### *Guidance 2*

*Guidance 2 content*


## Examples

### Return a Dictionary of Solutions

You can solve a system of equations leaving another symbol as a variable. You
can specify the variables to solve for as multiple arguments, or as a list (or
tuple):

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, y + 4*z]
>>> solve(equations, x, y, dict=True)
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
>>> solve(equations, [x, y], dict=True)
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
```

### Return a Set of Solution(s)

To get a list of symbols and set of solution(s) use `set=True` instead of
`dict=True`:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x**2 + y - 2*z, y + 4*z], x, y, set=True)
([x, y], {(-sqrt(6)*sqrt(z), -4*z), (sqrt(6)*sqrt(z), -4*z)})
```

## Use the Solution Result

### Use a Dictionary of Solutions

You can extract solutions by specifying in brackets ("slicing") the solution
number, and then the symbol, for example `solution[0][x]` gives the result for
`x` in the first solution:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, y + 4*z]
>>> solution = solve([x**2 + y - 2*z, y + 4*z], x, y, dict=True)
>>> solution
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
>>> solution[0][x]
-sqrt(6)*sqrt(z)
>>> solution[0][y]
-4*z
```

## Options That Can Speed up {func}`~.solve`

*After solving-guidance merged, link to target in it speed_up_solve as:* `Refer
to [solving guidance](speed_up_solve).`

## Not All Systems of Equations Can be Solved

### Systems of Equations With no Solution

Some systems of equations have no solution. For example, the following system
reduces to $z = 2z$, which has no general solution:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x + y - z, x + y - 2*z], [x, y], dict=True)
[]
```

The following system is overconstrained, meaning there are more equations (3)
than unknowns to be solved for (2, namely $x$ and $y$):

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x + y - z, x - 2*z, 2*x - y], [x, y], dict=True)
[]
```

### Systems of Equations With no Analytical Solution

Some systems of equations cannot be solved algebraically, for example:

```py
>>> from sympy import cos, solve
>>> from sympy.abc import x, y, z
>>> solve([x - y, cos(x) - y], [x, y], dict=True)
Traceback (most recent call last):
    ...
NotImplementedError: could not solve -y + cos(y)
```

So you can use {func}`~.nsolve` to find a numerical solution.

### Equations Which Have an Analytical Solution, and SymPy Cannot Solve

It is also possible that there is an algebraic solution to your equation, and
SymPy has not implemented an appropriate algorithm. If that happens, or SymPy
returns an empty set or list when there is a mathematical solution (indicating a
bug in SymPy), please post it on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>`
instead.

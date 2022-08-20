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

Refer to
[](solving-guidance.md#include-the-variable-to-be-solved-for-in-the-function-call)
and [](solving-guidance.md#ensure-consistent-formatting-from).

There are two methods below for containing solution results: dictionary or set.
A dictionary is easier to interrogate programmatically, so if you need to
extract solutions using code, we recommend the dictionary approach.

## Solve and Use Results in a Dictionary

### Solve Into a Dictionary of Solutions

You can solve a system of equations for some variables (for example, $x$ and
$y$) leaving another symbol as a variable (for example, $z$). You can specify
the variables to solve for as multiple separate arguments, or as a list (or
tuple):

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, y + 4*z]
>>> solutions = solve(equations, x, y, dict=True)
>>> solutions
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
```

### Use a Dictionary of Solutions

You can then extract solutions by specifying in brackets ("slicing") the
solution number, and then the symbol. For example `solutions[0][x]` gives the
result for `x` in the first solution:

```py
>>> solutions[0][x]
-sqrt(6)*sqrt(z)
>>> solutions[0][y]
-4*z
```

## Solve Results in a Set

To get a list of symbols and set of solutions, use `set=True` instead of
`dict=True`:

```py
from sympy import solve
from sympy.abc import x, y, z
solve([x**2 + y - 2*z, y + 4*z], [x, y], set=True)
([x, y], {(-sqrt(6)*sqrt(z), -4*z), (sqrt(6)*sqrt(z), -4*z)})
```

## Options That Can Speed up {func}`~.solve`

Refer to [](solving-guidance.md#options-that-can-speed-up).

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

The following system is overconstrained, meaning there are more equations
(three) than unknowns to be solved for (two, namely $x$ and $y$):

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x + y - z, x - 2*z, 2*x - y], [x, y], dict=True)
[]
```

### Systems of Equations With no Closed-Form Solution

Some systems of equations cannot be solved algebraically, for example those
containing [transcendental
equations](https://en.wikipedia.org/wiki/Transcendental_equation):

```py
>>> from sympy import cos, solve
>>> from sympy.abc import x, y, z
>>> solve([x - y, cos(x) - y], [x, y], dict=True)
Traceback (most recent call last):
    ...
NotImplementedError: could not solve -y + cos(y)
```

So you can use {func}`~.nsolve` to find a numerical solution.

### Equations Which Have a Closed-Form Solution, and SymPy Cannot Solve

It is also possible that there is an algebraic solution to your equation, and
SymPy has not implemented an appropriate algorithm. If SymPy returns an empty
set or list when you know there is a closed-form solution (indicating a bug in
SymPy), please post it on the [mailing list](https://groups.google.com/g/sympy),
or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>`
instead.

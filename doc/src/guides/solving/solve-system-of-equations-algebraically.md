# Solve a System of Equations Algebraically

Use SymPy to algebraically solve a system of equations, whether linear or
nonlinear. For example, solving $x^2 + y = 2z, y = -4z$ for x and y (assuming z
is a constant or parameter) yields $\{(x = -\sqrt{6z}, y = -4z),$ ${(x =
\sqrt{6z}, y = -4z)\}}$.

## Alternatives to Consider
- Some systems of equations cannot be solved algebraically (either at all or by
SymPy), so you may have to [solve your system of equations
numerically](solve-numerically.md) using {func}`~.nsolve` instead.

## Examples of Solving a System of Equations Algebraically

Whether your equations are linear or nonlinear, you can use {func}`~.solve`:

### Solve a System of Linear Equations Algebraically

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x + y - 2*z, y + 4*z], [x, y], dict=True)
[{x: 6*z, y: -4*z}]
```

### Solve a System of Nonlinear Equations Algebraically

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x**2 + y - 2*z, y + 4*z], x, y, dict=True)
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
```

## Guidance

Refer to
[](solving-guidance.md#include-the-variable-to-be-solved-for-in-the-function-call)
and [](ensure-consistent-formatting-from-solve).

There are two methods below for containing solution results:
[dictionary](#solve-and-use-results-in-a-dictionary) or
[set](#solve-results-in-a-set). A dictionary is easier to interrogate
programmatically, so if you need to extract solutions using code, we recommend
the dictionary approach.

## Solve and Use Results in a Dictionary

### Solve Into a Solution Given as a Dictionary

You can solve a system of equations for some variables (for example, $x$ and
$y$) leaving another symbol as a constant or parameter (for example, $z$). You
can specify the variables to solve for as multiple separate arguments, or as a
list (or tuple):

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, y + 4*z]
>>> solutions = solve(equations, x, y, dict=True)
>>> solutions
[{x: -sqrt(6)*sqrt(z), y: -4*z}, {x: sqrt(6)*sqrt(z), y: -4*z}]
```

### Use a Solution Given as a Dictionary

You can then extract solutions by indexing (specifying in brackets) the solution
number, and then the symbol. For example `solutions[0][x]` gives the result for
`x` in the first solution:

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

Refer to [](options-that-can-speed-up-solve).

## Not All Systems of Equations Can be Solved

### Systems of Equations With no Solution

Some systems of equations have no solution. For example, the following two
systems have no solution because they reduce to `1 == 0`, so SymPy returns an
empty list:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> solve([x + y - 1, x + y], [x, y], dict=True)
[]
```

```py
from sympy import solve
from sympy.abc import x, y, z
solve([x + y - (z + 1), x + y - z)], [x, y], dict=True)
[]
```

The following system reduces to $z = 2z$, so it has no general solution, but it
could be satisfied if $z=0$. Note that {func}`~.solve` will not assume that
$z=0$, even though that is the only value of $z$ that makes the system of
equations consistent, because $z$ is a parameter rather than an unknown. That
is, {func}`~.solve` does not treat $z$ as an unknown because it is not in the
list of symbols specified as unknowns (`[x, y]`) and all such symbols are
treated like parameters with arbitrary value. Whether a symbol is treated as a
variable or a parameter is determined only by whether it is specified as a
symbol to solve for in {func}`~.solve`. There is no such distinction made when
creating the symbol using {func}`~.symbols` (or importing from {mod}`~.abc`).

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x + y - z, x + y - 2*z], [x, y], dict=True)
[]
```

The following system is
[overconstrained](https://en.wikipedia.org/wiki/Overdetermined_system), meaning
there are more equations (three) than unknowns to be solved for (two, namely $x$
and $y$). It has no solution:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> solve([x + y - z, x - (z + 1), 2*x - y], [x, y], dict=True)
[]
```

Note that some overconstrained systems do have solutions (for example, if an
equation is a linear combination of the others), in which case SymPy can solve
the overconstrained system.

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

So you can use {func}`~.nsolve` to [find a numerical
solution](solve-numerically.md):

```py
>>> from sympy import cos, nsolve
>>> from sympy.abc import x, y, z
>>> nsolve([x - y, cos(x) - y], [x, y], [1,1])
    Matrix([
    [0.739085133215161],
    [0.739085133215161]])
```

### Equations Which Have a Closed-Form Solution, and SymPy Cannot Solve

It is also possible that there is an algebraic solution to your equation, and
SymPy has not implemented an appropriate algorithm. If SymPy returns an empty
set or list when you know there is a closed-form solution (indicating a bug in
SymPy), please post it on the [mailing list](https://groups.google.com/g/sympy),
or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can use a different method listed in [](#alternatives-to-consider).

## Report a Bug

If you find a bug with {func}`~.solve`, please post the problem on the [SymPy
mailing list](https://groups.google.com/g/sympy). Until the issue is resolved,
you can use a different method listed in [](#alternatives-to-consider).

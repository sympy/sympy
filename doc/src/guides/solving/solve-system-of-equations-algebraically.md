# Solve a system of equations algebraically

Use SymPy to solve a system of equations algebraically. For example, $x^2 + y =
2, x - y = 4$ yields $\{(x = -3, y = -7), (x = 2, y = 2)\}$.

Alternatives to consider:
- *alternative 1*
- *alternative 2*

Here is a simple example of solving a system of equations algebraically:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> solve([x**2 + y - 2, x - y - 4], x, y, dict=True)
[{x: -3, y: -7}, {x: 2, y: -2}]
```

## Guidance

*Just link to solving guidance page?*

### *Guidance 1*

*Guidance 1 content*

### *Guidance 2*

*Guidance 2 content*


## Examples

### Return a dictionary of solutions

You can solve a system of equations leaving another symbol as a variable. You
can specify the variables to solve for as multiple arguments, or as a list (or
tuple):

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, x - y - 4*z]
>>> solve(equations, x, y, dict=True)
[{x: -sqrt(24*z + 1)/2 - 1/2, y: -4*z - sqrt(24*z + 1)/2 - 1/2}, {x: sqrt(24*z + 1)/2 - 1/2, y: -4*z + sqrt(24*z + 1)/2 - 1/2}]
>>> solve(equations, [x, y], dict=True)
[{x: -sqrt(24*z + 1)/2 - 1/2, y: -4*z - sqrt(24*z + 1)/2 - 1/2}, {x: sqrt(24*z + 1)/2 - 1/2, y: -4*z + sqrt(24*z + 1)/2 - 1/2}]
```

### Return a set of solution(s)

To get a list of symbols and set of solution(s) use `set=True` instead of
`dict=True`:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, x - y - 4*z]
>>> solve(equations, x, y, set=True)
([x, y], {(-sqrt(24*z + 1)/2 - 1/2, -4*z - sqrt(24*z + 1)/2 - 1/2), (sqrt(24*z + 1)/2 - 1/2, -4*z + sqrt(24*z + 1)/2 - 1/2)})
```

## Use the solution result

### Use a dictionary of solutions

You can extract solutions by specifying in brackets ("slicing") the solution number, and then the symbol, for example `solution[0][x]` gives the result for `x` in the first solution:

```py
>>> from sympy import solve
>>> from sympy.abc import x, y, z
>>> equations = [x**2 + y - 2*z, x - y - 4*z]
>>> solution = solve(equations, x, y, dict=True)
>>> solution
[{x: -sqrt(24*z + 1)/2 - 1/2, y: -4*z - sqrt(24*z + 1)/2 - 1/2}, {x: sqrt(24*z + 1)/2 - 1/2, y: -4*z + sqrt(24*z + 1)/2 - 1/2}]
>>> solution[0][x]
-sqrt(24*z + 1)/2 - 1/2
>>> solution[0][y]
-4*z - sqrt(24*z + 1)/2 - 1/2
```

### Use a set of solutions

*Usage method 2 content*

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

Please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can *workaround*.

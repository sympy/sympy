# Solve an equation algebraically

Use SymPy to solve an equation algebraically. For example, solving $x^2 = 4$ yields $x \in \{-2,2\}$.

There are two high-level functions to solve equations, [`solve`](#) and [`solveset`](#). Here are their advantages and disadvantages: ... *table?*

## Using {func}`~.solve`

{func}`~.solve`
- is best for solving equations that...
- produces outputs in the format of a Python list.

You can solve an equation using {func}`~.solve` in several ways.

1. Use the fact that any expression not in an `Eq` (equation) is automatically assumed to equal 0 by the solving functions. You can rearrange the equation $x^2 = 4$ to $x^2 - 4 = 0$, and {func}`~.solve` that expression. This approach is convenient if you are interactively solving an equation which can easily be rearranged to $expression = 0$.

```
>>> from sympy import solve
>>> from sympy.abc import x
>>> solve(x**2 - 4)
[-2, 2]
```

2. Put your equation into `Eq` form, then apply {func}`~.solve` to the `Eq`. This approach is convenient if you are interactively solving an equation which cannot easily be rearranged to $expression = 0$.

```
>>> from sympy import solve, Eq
>>> from sympy.abc import x
>>> eqn = Eq(x**2, 4)
>>> solve(eqn)
[-2, 2]
```

3. Parse a string representing the equation into a form that SymPy can understand, then apply {func}`~.solve` to the parsed expression. 

*But this doesn't work for equality; does for inequality. Is there some way to do this for equality?*

```
>>> from sympy import solve, parse_expr
>>> expr = "x**2 == 4"
>>> parsed = parse_expr(expr)
>>> solve(parsed)
[-2, 2]
```

## Using {func}`~.solveset`

{func}`~.solveset`
- is best for solving equations that...
- produces outputs in the format of [SymPy mathematical Sets](https://docs.sympy.org/dev/modules/sets.html?highlight=sets#module-sympy.sets.sets) rather than [Python sets](https://docs.python.org/3/library/stdtypes.html#set)
- can return infinitely many solutions
- clearly separates the complex and real domains

```
>>> from sympy import solveset
>>> from sympy.abc import x
>>> solveset(x**2 - 4)
FiniteSet(-2, 2)
```

## Not all equations can be solved algebraically

Some equations have no algebraic solution, in which case SymPy may return an empty set:

```
>>> from sympy import solve, Eq
>>> from sympy.abc import x
>>> eqn = Eq(3 * x - 7, 3 * x + 2)
>>> solve(eqn)
[]
```

You may want to check whether there is a mistake in the equation.

If SymPy returns an error such as `NotImplementedError`, there may be no way to solve the equation algebraically:

```
>>> from sympy import solve, cos
>>> from sympy.abc import x
>>> solve(cos(x) - x)
...
NotImplementedError: multiple generators [x, cos(x)]
No algorithms are implemented to solve equation -x + cos(x)
```

so you may have to {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` instead.

It is also possible that there is a way to solve your equation algebraically, and SymPy has not have implemented an appropriate algorithm. You can ask about this on the mailing list, or open an issue on GitHub. 

*Is it appropriate to include such instructions on a page like this?*

# Solve an equation algebraically

Use SymPy to solve an equation algebraically (symbolically). For example, solving $x^2 = y$ yields $x \in \{-\sqrt{y},\sqrt{y}\}$.

There are two high-level functions to solve equations, {func}`~.solve` and {func}`~.solveset`. Here are their advantages and disadvantages: ... *table?*

We recommend you include the variable to be solved for, as the second argument for either function. While this is optional for equations with a single symbol, it is a good practice because it ensures SymPy will solve for the desired symbol.

## Using {func}`~.solve`

{func}`~.solve`
- is best for solving equations that...
- produces explicit solution expressions that you can use to substitute into something else
- produces various output formats depending on the answer
    - unless you use `dict=True` to ensure the result will be formatted as a dictionary, which we recommend if you want to extract information from the result programmatically

You can solve an equation using {func}`~.solve` in several ways.

*Mention assumptions: can set to real for a cubic equation to remove complex roots, if only care about real roots*

1. Use the fact that any expression not in an `Eq` (equation) is automatically assumed to equal 0 by the solving functions. You can rearrange the equation $x^2 = y$ to $x^2 - y = 0$, and {func}`~.solve` that expression. This approach is convenient if you are interactively solving an equation which can easily be rearranged to $expression = 0$.

```
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> solve(x**2 - y, x)
[-sqrt(y), sqrt(y)]
```

2. Put your equation into `Eq` form, then apply {func}`~.solve` to the `Eq`. This approach is convenient if you are interactively solving an equation which cannot easily be rearranged to $expression = 0$.

```
>>> from sympy import solve, Eq
>>> from sympy.abc import x, y
>>> eqn = Eq(x**2, y)
>>> solve(eqn, x)
[-sqrt(y), sqrt(y)]
```

3. Parse a string representing the equation into a form that SymPy can understand, then apply {func}`~.solve` to the parsed expression.  This approach is convenient if you are programmatically reading in a string; we recommend against using it if you are creating the expression yourself.

Should always include the variable to solve for if want to extract results programmatically

*But this doesn't work for equality; does for inequality. Is there some way to do this for equality? -- transformations https://docs.sympy.org/dev/modules/parsing.html?highlight=parse_expr#sympy.parsing.sympy_parser.parse_expr convert equal sign transformations. Should be in a different guide?*

https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#strings-as-input

```
>>> from sympy import solve, parse_expr
>>> from sympy.abc import x
>>> expr = "Eq(x**2, y)"
>>> parsed = parse_expr(expr)
>>> solve(parsed, x)
[-sqrt(y), sqrt(y)]
```

## Using {func}`~.solveset`

{func}`~.solveset`
- is best for solving equations that...
- produces outputs in the format of [SymPy mathematical Sets](https://docs.sympy.org/dev/modules/sets.html?highlight=sets#module-sympy.sets.sets) rather than [Python sets](https://docs.python.org/3/library/stdtypes.html#set)
- can return infinitely many solutions
- clearly separates the complex and real domains
- the solution set can be more difficult to parse programmatically (trig function as example)

```
>>> from sympy import solveset
>>> from sympy.abc import x, y
>>> solveset(x**2 - y, x)
FiniteSet(sqrt(y), -sqrt(y))
```

*Discuss domain argument--e.g. complex or real numbers: domain=reals*

## Not all equations can be solved

### Equations with no solution

Some equations have no solution, in which case SymPy may return an empty set. For example, the following equation reduces to -7 = 2, which has no solution because no value of $x$ will make it true:

```
>>> from sympy import solve, Eq
>>> from sympy.abc import x
>>> eqn = Eq(3 * x - 7, 3 * x + 2)
>>> solve(eqn)
[]
```

You may want to check whether there is a mistake in the equation.

### Equations with no analytical solution

The vast majority of arbitrary nonlinear equations are not analytically solvable. The classes of equations that are solvable are basically:
1. Linear equations.
2. Polynomials (except where limited by the Abel-Ruffini theorem).
3. Equations that can be solved by inverting some transcendental functions.
4. Problems that can be transformed into the cases above (e.g. by turning trig into polynomials).
5. A few other special cases that can be solved with something like the Lambert W function.

SymPy may reflect that your equation has no solutions that can be expressed algebraically (symbolically) by returning an error such as `NotImplementedError`:

```
>>> from sympy import solve, cos
>>> from sympy.abc import x
>>> solve(cos(x) - x, x)
...
NotImplementedError: multiple generators [x, cos(x)]
No algorithms are implemented to solve equation -x + cos(x)
```

so you may have to {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` instead.

### Equations which have an analytical solution, and SymPy cannot solve

It is also possible that there is a way to solve your equation algebraically, and SymPy has not have implemented an appropriate algorithm. You can ask about this on the [mailing list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub page](https://github.com/sympy/sympy/issues).

*Ask if someone has an example where there is a mathematical solution but SymPy returns an empty list. Related: sin(x) = 0 will return two solutions but there are infinitely many solutions*
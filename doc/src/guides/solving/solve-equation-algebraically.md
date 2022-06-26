# Solve an equation algebraically

Use SymPy to solve an equation algebraically (symbolically). For example, solving $x^2 = y$ yields $x \in \{-\sqrt{y},\sqrt{y}\}$.

There are two high-level functions to solve equations, {func}`~.solve` and {func}`~.solveset`. Here are recommendations on when to use:

- {func}`~.solve`
    - You need to programmatically extract components (expressions, or individual symbols or constants) from the output by traversing it.
    - You want an explicit solution expression that you can use to substitute into something else.

- {func}`~.solveset`
    - You want a consistent input and output interface.
    - You want to get all the solutions, including if there are infinitely many.
    - You want a clear separation between equations in the complex domain and the real domain.

We recommend you include the variable to be solved for as the second argument for either function. While this is optional for equations with a single symbol, it is a good practice because it ensures SymPy will solve for the desired symbol. For example, you may expect the following to solve for $x$, and SymPy will solve for $y$:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x ** 2 - y)
[{y: x**2}]
```

Specifying the variable to solve for ensures that SymPy solves for it:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x ** 2 - y, x)
[-sqrt(y), sqrt(y)]
```

## Using {func}`~.solve`

- Produces various output formats depending on the answer
    - unless you use `dict=True` to ensure the result will be formatted as a dictionary, which we recommend if you want to extract information from the result programmatically

You can solve an equation using {func}`~.solve` in several ways.

### Expression equals zero

Use the fact that any expression not in an `Eq` (equation) is automatically assumed to equal zero (0) by the solving functions. You can rearrange the equation $x^2 = y$ to $x^2 - y = 0$, and {func}`~.solve` that expression. This approach is convenient if you are interactively solving an equation which can easily be rearranged to $expression = 0$.
    
```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> solve(x**2 - y, x)
[-sqrt(y), sqrt(y)]
```
    
### Put your equation into `Eq` form

Put your equation into `Eq` form, then apply {func}`~.solve` to the `Eq`. This approach is convenient if you are interactively solving an equation which cannot easily be rearranged to $expression = 0$.
    
    ```py
    >>> from sympy import solve, Eq
    >>> from sympy.abc import x, y
    >>> eqn = Eq(x**2, y)
    >>> solve(eqn, x)
    [-sqrt(y), sqrt(y)]
    ```
    
### Parse a string representing the equation

Parse a string representing the equation into a form that SymPy can understand (`Eq` form), then apply {func}`~.solve` to the parsed expression.  This approach is convenient if you are programmatically reading in a string. We [recommend against using parsing a string if you are creating the expression yourself](https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#strings-as-input).  
You should always include the variable to solve for if you want to extract results programmatically, to ensure that SymPy solves for the desired variable. To ensure SymPy will produce results in a consistent format, use `dict=True`. To extract the solutions, you can iterate through the list of dictionaries:  
    
```py
>>> from sympy import parse_expr, solve
>>> from sympy.abc import x
>>> from sympy.parsing.sympy_parser import convert_equals_signs, standard_transformations
>>> expr = "x ** 2 = y"
>>> parsed = parse_expr(expr, transformations=standard_transformations + (convert_equals_signs,))
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
```

If you already have the equation in `Eq` form, you can parse that string:

```py
>>> from sympy import parse_expr, solve
>>> from sympy.abc import x
>>> expr = "Eq(x**2, y)"
>>> parsed = parse_expr(expr)
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
```

### Restricting the domain of solutions using {func}`~.solve`

By default, SymPy will return solutions in the complex domain, which also includes purely real and imaginary values. Here, the first two solutions are real, and the last two are imaginary:

```py
>>> from sympy import Symbol, solve
>>> x = Symbol('x')
>>> solution = solve(x ** 4 - 256, x)
>>> print(solution)
[-4, 4, -4*I, 4*I]
```

If you want to restrict returned solutions to real numbers, you can place an assumption on the symbol to be solved for, $x$:

```py
>>> from sympy import Symbol, solve
>>> x = Symbol('x', real=True)
>>> solution = solve(x ** 4 - 256, x)
>>> print(solution)
[-4, 4]
```

## Using {func}`~.solveset`

- Produces outputs in the format of [SymPy mathematical Sets](https://docs.sympy.org/dev/modules/sets.html?highlight=sets#module-sympy.sets.sets) rather than [Python sets](https://docs.python.org/3/library/stdtypes.html#set)
- can return infinitely many solutions
- the solution set can be more difficult to parse programmatically *trig function as example*

```py
>>> from sympy import solveset
>>> from sympy.abc import x, y
>>> solution = solveset(x**2 - y, x)
>>> print(solution)
{-sqrt(y), sqrt(y)}
```

### Restricting the domain of solutions using {func}`~.solveset`

By default, SymPy will return solutions in the complex domain, which also includes purely real and imaginary values. Here, the first two solutions are real, and the last two are imaginary:

```py
>>> from sympy import S, solveset
>>> from sympy.abc import x
>>> solution = solveset(x**4 - 256, x)
>>> print(solution)
{-4, 4, -4*I, 4*I}
```

If you want to restrict returned solutions to real numbers, you can specify the [domain to solve in](https://docs.sympy.org/dev/modules/solvers/solveset.html?highlight=solveset#what-is-this-domain-argument-about) as `S.Reals`:

```py
>>> from sympy import S, solveset
>>> from sympy.abc import x
>>> solution = solveset(x**4 - 256, x, domain=S.Reals)
>>> print(solution)
{-4, 4}
```

## Not all equations can be solved

### Equations with no solution

Some equations have no solution, in which case SymPy may return an empty set. For example, the following equation reduces to -7 = 2, which has no solution because no value of $x$ will make it true:

```py
>>> from sympy import solve, Eq
>>> from sympy.abc import x
>>> eqn = Eq(3 * x - 7, 3 * x + 2)
>>> solve(eqn)
[]
```

So if SymPy returns an empty set, you may want to check whether there is a mistake in the equation.

### Equations with no analytical solution

The vast majority of arbitrary nonlinear equations are not analytically solvable. The classes of equations that are solvable are basically:
1. Linear equations
2. Polynomials, except where limited by the Abel-Ruffini theorem
3. Equations that can be solved by inverting some transcendental functions
4. Problems that can be transformed into the cases above (e.g., by turning trigonometric functions into polynomials)
5. A few other special cases that can be solved with something like the Lambert W function

SymPy may reflect that your equation has no solutions that can be expressed algebraically (symbolically) by returning an error such as `NotImplementedError`:

```py
>>> from sympy import solve, cos
>>> from sympy.abc import x
>>> solve(cos(x) - x, x)
Traceback (most recent call last):
  ...
NotImplementedError: multiple generators [x, cos(x)]
No algorithms are implemented to solve equation -x + cos(x)
```

so you may have to {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` instead.

### Equations which have an analytical solution, and SymPy cannot solve

It is also possible that there is a way to solve your equation algebraically, and SymPy has not have implemented an appropriate algorithm. If that happens, or SymPy returns an empty set when there is a mathematical solution (indicating a bug in SymPy), please post it on the [mailing list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you can {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` instead.

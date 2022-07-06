# Solve an equation algebraically

Use SymPy to solve an equation algebraically (symbolically). For example, solving $x^2 = y$ yields $x \in \{-\sqrt{y},\sqrt{y}\}$.

Alternatives to consider:
- SymPy can also [solve many other types of problems including sets of equations](index.md).
- Some equations cannot be solved algebraically (either at all or by SymPy), 
so you may have to {func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` instead.

There are two high-level functions to solve equations, {func}`~.solve` and {func}`~.solveset`.
Here is a simple example of each:

{func}`~.solve`

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solution = solve(x ** 2 - y, x, dict=True)
>>> print(solution)
[{x: -sqrt(y)}, {x: sqrt(y)}]
```

{func}`~.solveset`

```py
>>> from sympy import solveset
>>> from sympy.abc import x, y
>>> solution = solveset(x**2 - y, x)
>>> print(solution)
{-sqrt(y), sqrt(y)}
```

Here are recommendations on when to use:

- {func}`~.solve`
    - You want to get explicit symbolic representations of the different values
    a variable could take that would satisfy the equation.
    - You want to substitute those explicit solution values into other equations
    or expressions involving the same variable using {meth}`~sympy.core.basic.Basic.subs`

- {func}`~.solveset`
    - You want to represent the solutions in a mathematically precise way, using [mathematical sets](../../modules/sets.rst).
    - You want a representation of all the solutions, including if there are infinitely many.
    - You want a consistent input interface.
    - You want to limit the domain of the solutions to any arbitrary set.
    - You do not need to programmatically extract solutions from the solution set:
    solution sets cannot necessarily be interrogated programmatically.

## Guidance

### Include the variable to be solved for in the function call

We recommend you include the variable to be solved for as the second argument for either function. 
While this is optional for equations with a single symbol, it is a good practice because it ensures 
SymPy will solve for the desired symbol. For example, you may expect the following to solve for $x$, 
and SymPy will solve for $y$:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solution = solve(x ** 2 - y, dict=True)
>>> print(solution)
[{y: x**2}]
```

Specifying the variable to solve for ensures that SymPy solves for it:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solution = solve(x ** 2 - y, x, dict=True)
>>> print(solution)
[{x: -sqrt(y)}, {x: sqrt(y)}]
```

### Ensure consistent formatting from {func}`~.solve` by using `dict=True`

{func}`~.solve` produces various output formats depending on the answer, 
unless you use `dict=True` to ensure the result will be formatted as a dictionary. 
We recommend using `dict=True`, especially if you want to 
extract information from the result programmatically.

## Solve an equation using {func}`~.solve` or {func}`~.solveset`

You can solve an equation using in several ways. 
The examples below demonstrate using both {func}`~.solve` and {func}`~.solveset` where applicable. 
You can choose the function best suited to your equation.

### Make your equation into an expression that equals zero

Use the fact that any expression not in an `Eq` (equation) is automatically assumed to equal zero (0) 
by the solving functions. You can rearrange the equation $x^2 = y$ to $x^2 - y = 0$, and 
solve that expression. This approach is convenient if you are interactively solving 
an expression which already equals zero, or an equation that you do not mind rearranging to 
$expression = 0$.

```py
>>> from sympy import solve, solveset
>>> from sympy.abc import x, y
>>> solution = solve(x ** 2 - y, x, dict=True)
>>> print(solution)
[{x: -sqrt(y)}, {x: sqrt(y)}]
>>> solution_set = solveset(x**2 - y, x)
>>> print(solution_set)
{-sqrt(y), sqrt(y)}
```
    
### Put your equation into `Eq` form

Put your equation into `Eq` form, then solve the `Eq`. 
This approach is convenient if you are interactively solving an equation 
which you already have in the form of an equation,
or which you think of as an equality.
    
```py
>>> from sympy import Eq, solve, solveset
>>> from sympy.abc import x, y

>>> eqn = Eq(x**2, y)
>>> solutions = solve(eqn, x, dict=True)
>>> print(solutions)
[{x: -sqrt(y)}, {x: sqrt(y)}]
>>> solutions_set = solveset(eqn, x)
>>> print(solutions_set)
{-sqrt(y), sqrt(y)}
>>> for solution_set in solutions_set:
...     print(solution_set)
sqrt(y)
-sqrt(y)
```
    
### Parse a string representing the equation

Parse a string representing the equation into a form that SymPy can understand (`Eq` form), 
then apply {func}`~.solve` to the parsed expression.  
This approach is convenient if you are programmatically reading in a string. 
We [recommend against using parsing a string if you are creating the expression yourself](https://github.com/sympy/sympy/wiki/Idioms-and-Antipatterns#strings-as-input). 
Parsing an equation from a string requires you to use 
{func}`transformations <sympy.parsing.sympy_parser.parse_expr>`
for SymPy to interpret equals signs and create symbols from your variables.

You should always include the variable to solve for if you want to extract results programmatically, 
to ensure that SymPy solves for the desired variable. 
To ensure SymPy will produce results in a consistent format, use `dict=True`. 
To extract the solutions, you can iterate through the list of dictionaries:  
    
```py
>>> from sympy import parse_expr, solve, solveset
>>> from sympy.abc import x
>>> expr = "x ** 2 = y"
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

If you already have the equation in `Eq` form, you can parse that string:

```py
>>> from sympy import parse_expr, solve, solveset
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
>>> solutions_set = solveset(parsed, x)
>>> print(solutions_set)
{-sqrt(y), sqrt(y)}
```

### Restrict the domain of solutions

By default, SymPy will return solutions in the complex domain, which also includes 
purely real and imaginary values. Here, the first two solutions are real, 
and the last two are imaginary:

```py
>>> from sympy import Symbol, solve, solveset
>>> x = Symbol('x')
>>> solution = solve(x ** 4 - 256, x, dict=True)
>>> print(solution)
[{x: -4}, {x: 4}, {x: -4*I}, {x: 4*I}]
>>> solution_set = solveset(x ** 4 - 256, x)
>>> print(solution_set)
{-4, 4, -4*I, 4*I}
```

If you want to restrict returned solutions to real numbers, you can
- for {func}`~.solve`, place an assumption on the symbol to be solved for, $x$

```py
from sympy import Symbol, solve, solveset
x = Symbol('x', real=True)
solution = solve(x ** 4 - 256, x, dict=True)
print(solution)
[{x: -4}, {x: 4}]
```

- for {func}`~.solveset`, limit the output domain in the function call
```py
>>> from sympy import S, solveset
>>> from sympy.abc import x
>>> solution = solveset(x**4 - 256, x, domain=S.Reals)
>>> print(solution)
{-4, 4}
```

If you restrict the solutions to a domain in which there are no solutions, {func}`~.solveset` will return the empty set, [EmptySet](../../modules/sets.rst):

```py
>>> from sympy import solveset, S
>>> from sympy.abc import x
>>> solution = solveset(x**2 + 1, x, domain=S.Reals)
>>> print(solution)
EmptySet
```

Using {func}`~.solveset`, you can restrict returned solutions to any arbitrary set, including an interval:

```py
>>> from sympy import Interval, pi, sin, solveset
>>> from sympy.abc import x
>>> solution = solveset(sin(x), x, Interval(-pi, pi))
>>> print(solution)
{0, -pi, pi}
```

### Explicitly represent infinite sets of possible solutions using {func}`~.solveset`

{func}`~.solveset` 
[can represent infinite sets of possible solutions](why-solveset)
and express them in standard mathematical notation, 
for example $\sin(x) = 0$ for $x = n * \pi$ for every integer value of $n$:

```py
>>> from sympy import pprint, sin, solveset
>>> from sympy.abc import x
>>> solution = solveset(sin(x), x)
>>> pprint(solution)
{2*n*pi | n in Integers} U {2*n*pi + pi | n in Integers}
```

However, {func}`~.solve` will return only a finite number of solutions:

```py
>>> from sympy import sin, solve
>>> from sympy.calculus.util import periodicity
>>> from sympy.abc import x
>>> f = sin(x)
>>> solution = solve(f, x)
>>> print(solution)
[0, pi]
>>> print(periodicity(f, x))
2*pi
```

{func}`~.solve` tries to return just enough solutions so that all (infinitely many) solutions can generated 
from the returned solutions by adding integer multiples of the {func}`~.periodicity` of the equation, here $2\pi$.

## Use the solution result

### Substitute solutions from {func}`~.solve` into an expression

You can substitute solutions from {func}`~.solve` into an expression.

A common use case is finding the critical points and values for a function $f$. 
At the critical points, the {func}`~.Derivative` equals zero (or is undefined). 
You can then obtain the function values at those critical points
by substituting the critical points back into the function using {meth}`~sympy.core.basic.Basic.subs`. You can also tell if the critical point is a maxima or minima by substituting the values
into the expression for the second derivative: a negative value indicates a maximum, and a positive value indicates a minimum.

```py
>>> from sympy.abc import x
>>> from sympy import solve, diff
>>> f = x**3 + x**2 - x
>>> derivative = diff(f, x)
>>> critical_points = solve(derivative, x, dict=True)
>>> print(critical_points)
[{x: -1}, {x: 1/3}]
>>> point1, point2 = critical_points
>>> print(f.subs(point1))
1
>>> print(f.subs(point2))
-5/27
curvature = diff(f, x, 2)
print(curvature.subs(point1))
print(curvature.subs(point2))
-4
4
```

### {func}`~.solveset` solution sets cannot necessarily be interrogated programmatically

If {func}`~.solveset` returns a finite set (class {class}`~.FiniteSet`), you can iterate through the solutions:

```py
>>> from sympy import solveset
>>> from sympy.abc import x, y
>>> solution_set = solveset(x**2 - y, x)
>>> print(solution_set)
{-sqrt(y), sqrt(y)}
>>> solution_list = list(solution_set)
>>> print(solution_list)
[sqrt(y), -sqrt(y)]
```

However, for more complex results, it may not be possible to list the solutions:

```py
>>> from sympy import S, solveset, symbols
>>> x, y = symbols('x, y')
>>> solution_set = solveset(x**2 - y, x, domain=S.Reals)
>>> print(solution_set)
Intersection({-sqrt(y), sqrt(y)}, Reals)
>>> solution_list = list(solution_set)
Traceback (most recent call last):
    ...
TypeError: The computation had not completed because of the undecidable set membership is found in every candidates.
```

## Not all equations can be solved

### Equations with no solution

Some equations have no solution, in which case SymPy may return an empty set. For example, the equation 
$x - 7 = x + 2$ reduces to $-7 = 2$, which has no solution because no value of $x$ will make it true:

```py
>>> from sympy import solve, Eq
>>> from sympy.abc import x
>>> eqn = Eq(x - 7, x + 2)
>>> solve(eqn, x)
[]
```

So if SymPy returns an empty list, you may want to check whether there is a mistake in the equation.

### Equations with no analytical solution

The vast majority of arbitrary nonlinear equations are not analytically solvable. 
The classes of equations that are solvable are basically:
1. Linear equations
2. Polynomials, except where limited by the Abel-Ruffini theorem
3. Equations that can be solved by inverting some transcendental functions
4. Problems that can be transformed into the cases above 
(e.g., by turning trigonometric functions into polynomials)
5. A few other special cases that can be solved with something like the Lambert W function

SymPy may reflect that your equation has no solutions that can be expressed 
algebraically (symbolically) by returning an error such as `NotImplementedError`:

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

It is also possible that there is a way to solve your equation algebraically, and SymPy has not have 
implemented an appropriate algorithm. If that happens, or SymPy returns an empty set when there is a 
mathematical solution (indicating a bug in SymPy), please post it on the 
[mailing list](https://groups.google.com/g/sympy), or open an issue on 
[SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you can 
{func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` instead.

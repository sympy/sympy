# Solve an equation algebraically

Use SymPy to solve an equation algebraically (symbolically). For example, solving $x^2 = y$ for $x$ yields $x \in \{-\sqrt{y},\sqrt{y}\}$.

Alternatives to consider:
- SymPy can also 
[solve many other types of problems including sets of equations](index.md).
- Some equations cannot be solved algebraically (either at all or by SymPy), 
so you may have to 
{func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` 
instead.

There are two high-level functions to solve equations, {func}`~.solve` and 
{func}`~.solveset`. Here is a simple example of each:

{func}`~.solve`

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x ** 2 - y, x, dict=True)
[{x: -sqrt(y)}, {x: sqrt(y)}]
```

{func}`~.solveset`

```py
>>> from sympy import solveset
>>> from sympy.abc import x, y
>>> solveset(x**2 - y, x)
{-sqrt(y), sqrt(y)}
```

Here are recommendations on when to use:

- {func}`~.solve`
    - You want to get explicit symbolic representations of the different 
    values a variable could take that would satisfy the equation.
    - You want to substitute those explicit solution values into other 
    equations or expressions involving the same variable using 
    {meth}`~sympy.core.basic.Basic.subs`

- {func}`~.solveset`
    - You want to represent the solutions in a mathematically precise way, 
    using [mathematical sets](../../modules/sets.rst).
    - You want a representation of all the solutions, including if there are 
    infinitely many.
    - You want a consistent input interface.
    - You want to limit the domain of the solutions to any arbitrary set.
    - You do not need to programmatically extract solutions from the solution 
    set: solution sets cannot necessarily be interrogated programmatically.

## Guidance

### Include the variable to be solved for in the function call

We recommend you include the variable to be solved for as the second argument 
for either function. While this is optional for equations with a single 
symbol, it is a good practice because it ensures 
SymPy will solve for the desired symbol. For example, you may expect the 
following to solve for $x$, 
and SymPy will solve for $y$:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x ** 2 - y, dict=True)
[{y: x**2}]
```

Specifying the variable to solve for ensures that SymPy solves for it:

```py
>>> from sympy.abc import x, y
>>> from sympy import solve
>>> solve(x ** 2 - y, x, dict=True)
[{x: -sqrt(y)}, {x: sqrt(y)}]
```

### Ensure consistent formatting from {func}`~.solve` by using `dict=True`

{func}`~.solve` produces various output formats depending on the answer, 
unless you use `dict=True` to ensure the result will be formatted as a 
dictionary. We recommend using `dict=True`, especially if you want to 
extract information from the result programmatically.

## Solve an equation using {func}`~.solve` or {func}`~.solveset`

You can solve an equation in several ways. 
The examples below demonstrate using both {func}`~.solve` and 
{func}`~.solveset` where applicable. 
You can choose the function best suited to your equation.

### Make your equation into an expression that equals zero

Use the fact that any expression not in an `Eq` (equation) is automatically 
assumed to equal zero (0) by the solving functions. You can rearrange the 
equation $x^2 = y$ to $x^2 - y = 0$, and solve that expression. This approach 
is convenient if you are interactively solving an expression which already 
equals zero, or an equation that you do not mind rearranging to 
$expression = 0$.

```py
>>> from sympy import solve, solveset
>>> from sympy.abc import x, y
>>> solve(x ** 2 - y, x, dict=True)
[{x: -sqrt(y)}, {x: sqrt(y)}]
>>> solveset(x**2 - y, x)
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
>>> eqn
Eq(x**2, y)
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

### Restrict the domain of solutions

By default, SymPy will return solutions in the complex domain, which also 
includes purely real and imaginary values. Here, the first two solutions are 
real, and the last two are imaginary:

```py
>>> from sympy import Symbol, solve, solveset
>>> x = Symbol('x')
>>> solve(x ** 4 - 256, x, dict=True)
[{x: -4}, {x: 4}, {x: -4*I}, {x: 4*I}]
>>> solveset(x ** 4 - 256, x)
{-4, 4, -4*I, 4*I}
```

To restrict returned solutions to real numbers, or another domain or range, 
the different solving functions use different methods.

For {func}`~.solve`, place an assumption on the symbol to be solved for, $x$

```py
>>> from sympy import Symbol, solve
>>> x = Symbol('x', real=True)
>>> solve(x ** 4 - 256, x, dict=True)
[{x: -4}, {x: 4}]
```

or restrict the solutions with standard Python techniques for filtering a 
list such as a list comprehension, or 
by adding inequalities to {func}`~.solve` (but the range must be continuous):

```py
>>> from sympy import Or, Symbol, solve
>>> x = Symbol('x', real=True)
>>> expr = (x-4)*(x-3)*(x-2)*(x-1)
>>> solution = solve(expr, x)
>>> print(solution)
[1, 2, 3, 4]
>>> solution_outside_2_3 = [v for v in solution if (v.is_real and Or(v<2,v>3))]
>>> print(solution_outside_2_3)
[1, 4]
>>> solution_2_3 = solve((expr,x>=2,x<=3), x)
>>> print(solution_2_3)
Eq(x, 2) | Eq(x, 3)
```

For {func}`~.solveset`, limit the output domain in the function call by 
setting a domain

```py
>>> from sympy import S, solveset
>>> from sympy.abc import x
>>> solveset(x**4 - 256, x, domain=S.Reals)
{-4, 4}
```

or by restricting returned solutions to any arbitrary set, including an 
interval:

```py
>>> from sympy import Interval, pi, sin, solveset
>>> from sympy.abc import x
>>> solveset(sin(x), x, Interval(-pi, pi))
{0, -pi, pi}
```

and if you restrict the solutions to a domain in which there are no solutions,
{func}`~.solveset` will return the empty set,
[EmptySet](../../modules/sets.rst):

```py
>>> from sympy import solveset, S
>>> from sympy.abc import x
>>> solveset(x**2 + 1, x, domain=S.Reals)
EmptySet
```

### Explicitly represent infinite sets of possible solutions using {func}`~.solveset`

{func}`~.solveset` [can represent infinite sets of possible
solutions](why-solveset) and express them in standard mathematical notation, for
example $\sin(x) = 0$ for $x = n * \pi$ for every integer value of $n$:

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
>>> solve(f, x)
[0, pi]
>>> periodicity(f, x)
2*pi
```

{func}`~.solve` tries to return just enough solutions so that all 
(infinitely many) solutions can generated from the returned solutions 
by adding integer multiples of the {func}`~.periodicity` of the 
equation, here $2\pi$.

## Use the solution result

### Substitute solutions from {func}`~.solve` into an expression

You can substitute solutions from {func}`~.solve` into an expression.

A common use case is finding the critical points and values for a function 
$f$. At the critical points, the {class}`~.Derivative` equals zero (or is 
undefined). You can then obtain the function values at those critical points
by substituting the critical points back into the function using 
{meth}`~sympy.core.basic.Basic.subs`. You can also tell if the critical point 
is a maxima or minima by substituting the values into the expression for the 
second derivative: a negative value indicates a maximum, and a positive value 
indicates a minimum.

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
>>> curvature = diff(f, x, 2)
>>> print(curvature.subs(point1))
-4
>>> print(curvature.subs(point2))
4
```

### {func}`~.solveset` solution sets cannot necessarily be interrogated programmatically

If {func}`~.solveset` returns a finite set (class {class}`~.FiniteSet`), you can
iterate through the solutions:

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

However, for more complex results, it may not be possible to list the 
solutions:

```py
>>> from sympy import S, solveset, symbols
>>> x, y = symbols('x, y')
>>> solution_set = solveset(x**2 - y, x, domain=S.Reals)
>>> print(solution_set)
Intersection({-sqrt(y), sqrt(y)}, Reals)
>>> list(solution_set)
Traceback (most recent call last):
    ...
TypeError: The computation had not completed because of the undecidable set 
membership is found in every candidates.
```

In this case, it is because, if $y$ is negative, its square root would be
imaginary rather than real and therefore outside the declared domain of the
solution set. By declaring $y$ to be real and positive, SymPy can determine that
its square root is real, and thus resolve the intersection between the solutions
and the set of real numbers:

```py
>>> from sympy import S, Symbol, solveset
>>> x = Symbol('x')
>>> y = Symbol('y', real=True, positive=True)
>>> solution_set = solveset(x**2 - y, x, domain=S.Reals)
>>> print(solution_set)
{-sqrt(y), sqrt(y)}
>>> list(solution_set)
[sqrt(y), -sqrt(y)]
```

Alternatively, you can extract the sets from the solution set using {any}`args
<sympy.core.basic.Basic.args>`, then create a list from the set containing the
symbolic solutions:

```py
>>> from sympy import S, solveset, symbols
>>> x, y = symbols('x, y')
>>> solution_set = solveset(x**2 - y, x, domain=S.Reals)
>>> print(solution_set)
Intersection({-sqrt(y), sqrt(y)}, Reals)
>>> solution_set_args = solution_set.args
>>> print(solution_set.args)
(Reals, {-sqrt(y), sqrt(y)})
>>> list(solution_set_args[1])
[sqrt(y), -sqrt(y)]
```

## Options that can speed up {func}`~.solve`

### Include solutions making any denominator zero by using `check=False`

Normally, {func}`~.solve` checks whether any solutions make any denominator 
zero, and automatically excludes them. If you want to include those 
solutions, and speed up {func}`~.solve` (at the risk of obtaining invalid 
solutions), set `check=False`:

```py
>>> from sympy import Symbol, sin, solve
>>> x = Symbol("x")
>>> solve(sin(x)/x)  # 0 is excluded
[pi]
>>> solve(sin(x)/x, check=False)  # 0 is not excluded
[0, pi]
```

### Do not simplify solutions by using `simplify=False`

Normally, {func}`~.solve` simplifies all but polynomials of order 3 or 
greater before returning them and (if `check` is not False) uses the general 
{func}`simplify <sympy.simplify.simplify.simplify>` function on the solutions 
and the expression obtained when they are substituted into the function which 
should be zero. If you do not want the solutions simplified, and want to 
speed up {func}`~.solve`, use `simplify=False`.

```py
>>> from sympy import solve
>>> from sympy.abc import x, y
>>> expr = x**2 - (y**5 - 3*y**3 + y**2 - 3)
>>> solve(expr, x, dict=True)
[{x: -sqrt(y**5 - 3*y**3 + y**2 - 3)}, {x: sqrt(y**5 - 3*y**3 + y**2 - 3)}]
>>> solve(expr, x, dict=True, simplify=False)
[{x: -sqrt((y + 1)*(y**2 - 3)*(y**2 - y + 1))}, {x: sqrt((y + 1)*(y**2 - 3)*(y**2 - y + 1))}]
```

## Not all equations can be solved

### Equations with no solution

Some equations have no solution, in which case SymPy may return an empty set. 
For example, the equation $x - 7 = x + 2$ reduces to $-7 = 2$, which has no 
solution because no value of $x$ will make it true:

```py
>>> from sympy import solve, Eq
>>> from sympy.abc import x
>>> eqn = Eq(x - 7, x + 2)
>>> solve(eqn, x)
[]
```

So if SymPy returns an empty list, you may want to check whether there is a 
mistake in the equation.

### Equations which have an analytical solution, and SymPy cannot solve

It is also possible that there is an algebraic solution to your equation, 
and SymPy has not implemented an appropriate algorithm. 
If that happens, or SymPy returns an empty set or list when there is a 
mathematical solution (indicating a bug in SymPy), please post it on the 
[mailing list](https://groups.google.com/g/sympy), or open an issue on 
[SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the issue 
is resolved, you can 
{func}`solve your equation numerically <sympy.solvers.solvers.nsolve>` 
instead.

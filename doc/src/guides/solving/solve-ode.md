# Solve an Ordinary Differential Equation Algebraically

Use SymPy to solve an ordinary differential equation algebraically. For example,
solving $y''(x) + 9y(x)=0 $ yields $ y(x)=C_{1} \sin(3x)+ C_{2} \cos(3x)$.

Alternatives to consider:
- *which SciPy functions? from
  https://docs.scipy.org/doc/scipy/reference/integrate.html?highlight=ode?*

Here is an example of solving an ordinary differential equation algebraically
using {func}`~.dsolve`:

```py
>>> from sympy import Function, dsolve, Derivative
>>> from sympy.abc import x
>>> y = Function('y')
>>> result = dsolve(Derivative(y(x), x, x) + 9*y(x), y(x))
>>> result
Eq(y(x), C1*sin(3*x) + C2*cos(3*x))
```

You can then use SymPy to verify that the solution is correct:

```py
>>> from sympy import checkodesol
>>> solution = result.rhs
>>> solution
C1*sin(3*x) + C2*cos(3*x)
>>> checkodesol(Derivative(y(x), x, x) + 9*y(x), solution)
(True, 0)
```

The output of {func}`~.checkodesol` is a tuple where the first item, a boolean,
tells whether the substitution results in `0`.

## Guidance

### Input Format

There are many ways to express derivatives of functions. For an undefined
function, both {class}`~.Derivative` and {func}`~.diff` represent the undefined
derivative. Thus, all of the following `ypp` represent $y''$ ("y prime prime"),
the second derivative with respect to $x$ of a function $y(x)$:

```py
ypp = y(x).diff(x, x)
ypp = y(x).diff(x, 2)
ypp = y(x).diff((x, 2))
ypp = diff(y(x), x, x)
ypp = diff(y(x), x, 2)
ypp = Derivative(y(x), x, x)
ypp = Derivative(y(x), x, 2)
ypp = Derivative(Derivative(y(y(x), x), x)
ypp = diff(diff(y(x), x), x)
yp = y(x).diff(x)
ypp = yp.diff(x)
```

We recommend specifying the function to be solved for, as the second argument to
{func}`~.dsolve`. Note that it must be a function rather than a variable
(symbol). SymPy will give an error if you specify a variable ($x$) rather than a
function ($f(x)$):

```py
>>> dsolve(Derivative(y(x), x, x) + 9*y(x), x)
Traceback (most recent call last):
    ...
ValueError: dsolve() and classify_ode() only work with functions of one variable, not x
```

You can define the function to be solved for in two ways. The subsequent syntax
for specifying initial conditions depends on your choice.

## Define a Function Without Including Its Independent Variable

As in the example above, you can define a function without including its
independent variable:

```py
>>> from sympy import symbols, Eq, Function, dsolve
>>> f, g = symbols("f g", cls=Function)
>>> x = symbols("x")
>>> eqs = [Eq(f(x).diff(x), g(x)), Eq(g(x).diff(x), f(x))]
>>> dsolve(eqs, [f(x), g(x)])
[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]
```

Note that you supply the functions to be solved for as a list as the second
argument of {func}`~.dsolve`, here `[f(x), g(x)]`.

### Specify Initial (Boundary) Conditions

If your differential equation(s) have initial or boundary conditions, specify
them with the {func}`~.dsolve` optional argument `ics`. It should be given in
the form of `{f(x0): x1, f(x).diff(x).subs(x, x2): x3}` and so on. For power
series solutions, if no initial conditions are specified `f(0)` is assumed to be
`C0` and the power series solution is calculated about 0.

Here is an example of setting the initial values for functions:

```py
>>> from sympy import symbols, Eq, Function, dsolve
>>> f, g = symbols("f g", cls=Function)
>>> x = symbols("x")
>>> eqs = [Eq(f(x).diff(x), g(x)), Eq(g(x).diff(x), f(x))]
>>> dsolve(eqs, [f(x), g(x)])
[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]
>>> dsolve(eqs, [f(x), g(x)], ics={f(0): 1, g(0): 0})
[Eq(f(x), exp(x)/2 + exp(-x)/2), Eq(g(x), exp(x)/2 - exp(-x)/2)]
```

## Define a Function of an Independent Variable

You may prefer to specify a function (for example $x$) of its independent
variable (for example $t$):

```py
>>> from sympy import symbols, Function, dsolve
>>> t = symbols('t')
>>> x = Function('x')(t)
>>> x
x(t)
>>> xp = x.diff()
>>> xpp = xp.diff()
>>> eq = xpp + 2*xp + x
>>> eq
x(t) + 2*Derivative(x(t), t) + Derivative(x(t), (t, 2))
>>> dsolve(eq, x)
Eq(x(t), (C1 + C2*t)*exp(-t))
```

Using this convention, the second argument of {func}`~.dsolve`, `x`, represents
`x(t)`, so SymPy recognizes it as a valid function to solve for.

### Specify Initial (Boundary) Conditions Using {func}`~sympy.core.basic.Basic.subs`

Using that syntax, you specify initial conditions by substituting in values of
the independent variable using {func}`~sympy.core.basic.Basic.subs` because the
function $x$ already has its independent variable as an argument $t$:

```py
>>> dsolve(eq, x, ics={x.subs(t, 0): 0})
Eq(x(t), C2*t*exp(-t))
```

### Beware Copying and Pasting Results

If you choose to define a function of an independent variable, note that copying
a result and pasting it into subsequent code may cause an error because `x` is
already defined as `x(t)`, so if you paste in `x(t)` it is interpreted as
`x(t)(t)`:

```py
>>> dsolve(x(t).diff(), x)
Traceback (most recent call last):
    ...
TypeError: 'x' object is not callable
```

So remember to exclude the independent variable call `(t)`:

```py
>>> dsolve(x.diff(), x)
Eq(x(t), C1)
```

## Use the Solution Result

Unlike other solving functions, {func}`~.dsolve` returns an {class}`~.Equality`
(equation) formatted as, for example, `Eq(y(x), C1*sin(3*x) + C2*cos(3*x))`
which is equivalent to the mathematical notation $y(x) = C1 \sin(3x) + C2
\cos(3x)$.

### Extract the Result for One Solution and Function

You can extract the result from an {class}`~.Equality` using the right-hand side
property {any}`rhs <sympy.core.relational.Relational.rhs>`:

```py
>>> from sympy import Function, dsolve, Derivative
>>> from sympy.abc import x
>>> y = Function('y')
>>> result = dsolve(Derivative(y(x), x, x) + 9*y(x), y(x))
>>> result
Eq(y(x), C1*sin(3*x) + C2*cos(3*x))
>>> result.rhs
C1*sin(3*x) + C2*cos(3*x)
```

### Extract the Result for Multiple Function-Solution Pairs

If you are solving a system of equations with multiple unknown functions,
{func}`~.dsolve` will return a nested list of equalities, the outer list
representing each solution and the inner list representing each function. While
you can extract results by specifying the index of each function ("slicing" each
solution), we recommend an approach which is robust with respect to function
ordering. The following converts each solution into a dictionary so you can
easily extract the result for the desired function. It uses standard Python
techniques such as a loops or comprehensions, in a nested fashion.

```py
>>> from sympy import symbols, Eq, Function, dsolve
>>> y, z = symbols("y z", cls=Function)
>>> x = symbols("x")
>>> eqs = [Eq(y(x).diff(x)**2, z(x)**2), Eq(z(x).diff(x), z(x))]
>>> solutions = dsolve(eqs, [y(x), z(x)])
>>> solutions
[[Eq(y(x), C1 - C2*exp(x)), Eq(z(x), C2*exp(x))], [Eq(y(x), C1 + C2*exp(x)), Eq(z(x), C2*exp(x))]]
>>> solutions_list = [] # nested list approach
>>> for solution in solutions:
...     solution_dict = {}
...     for fn in solution:
...             solution_dict.update({fn.lhs: fn.rhs})
...     solutions_list.append(solution_dict)
>>> solutions_list
[{y(x): C1 - C2*exp(x), z(x): C2*exp(x)}, {y(x): C1 + C2*exp(x), z(x): C2*exp(x)}]
>>> solutions_list = [{fn.lhs:fn.rhs for fn in solution} for solution in solutions]
>>> solutions_list # nested comprehension approach
[{y(x): C1 - C2*exp(x), z(x): C2*exp(x)}, {y(x): C1 + C2*exp(x), z(x): C2*exp(x)}]
>>> solutions_list[0][y(x)]
C1 - C2*exp(x)
```

### Work With Arbitrary Constants

You can manipulate arbitrary constants such as `C1`, `C2`, and `C3`, which are
generated automatically by {func}`~.dsolve`, by creating them as symbols. For
example, if you want to assign values to arbitrary constants, you can create
them as symbols and then substitute in their values using
{meth}`~sympy.core.basic.Basic.subs`:

```py
>>> from sympy import Function, dsolve, Derivative, symbols, pi
>>> x, C1, C2 = symbols("x, C1, C2")
>>> y = Function('y')
>>> result = dsolve(Derivative(y(x), x, x) + 9*y(x), y(x)).rhs
>>> result
C1*sin(3*x) + C2*cos(3*x)
>>> result.subs({C1: 7, C2: pi})
7*sin(3*x) + pi*cos(3*x)
```

## Ordinary Differential Equation Solving Hints

### Return Unevaluated Integrals

By default, {func}`~.dsolve` attempts to evaluate the integrals it produces to
solve your ordinary differential equation. You can disable evaluation of the
integrals by using {ref}`hints` ending with `_Integral`, for example
`separable_Integral`. This is useful because
{func}`~sympy.core.expr.Expr.integrate` is an expensive routine. SymPy may hang
(appear to never complete the operation) because of a difficult or impossible
integral, so using an `_Integral` hint will at least return an (unintegrated)
result, which you can then consider. The simplest way to disable integration is
with the `all_Integral` hint because you do not need to know which hint to
supply: for any hint with a corresponding `_Integral` hint, `all_Integral` only
returns the `_Integral` hint.

### Select a Specific Solver

You may wish to select a specific solver using a hint for a couple of reasons:
- educational purposes: for example if you are learning about a specific method
  to solve ODEs and want to get a result that exactly matches that method
- form of the result: sometimes an ODE can be solved by many different solvers,
  and they can return different results. They will be mathematically equivalent,
  though the arbitrary constants may not be. {func}`~.dsolve` by default tries
  to use the "best" solvers first, which are most likely to return the most
  usable output, but it is not a perfect heuristic. For example, the "best"
  solver may produce a result with an integral that SymPy cannot solve, but
  another solver may produce a different integral that SymPy can solve. So if
  the solution isn't in a form you like, you can try other hints to check
  whether they give a preferable result.

## Not All Equations Can Be Solved

### Equations With No Solution

Not all differential equations can be solved, for example:

```py
>>> from sympy import Function, dsolve, Derivative, symbols
>>> x, C1, C2 = symbols("x, C1, C2")
>>> y = Function('y')
>>> dsolve(Derivative(y(x), x, 3) - (y(x)**2), y(x)).rhs
Traceback (most recent call last):
    ...
NotImplementedError: solve: Cannot solve -y(x)**2 + Derivative(y(x), (x, 3))
```

### Equations With No Analytical Solution

*Equations with no analytical solution content*

### Equations Which Have An Analytical Solution, and SymPy Cannot Solve

*Equations which have an analytical solution, and SymPy cannot solve content*

Please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can *workaround*.

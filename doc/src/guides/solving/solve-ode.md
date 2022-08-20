# Solve an Ordinary Differential Equation Algebraically

Use SymPy to solve an ordinary differential equation algebraically. For example,
solving $y''(x) + 9y(x)=0 $ yields $ y(x)=C_{1} \sin(3x)+ C_{2} \cos(3x)$.

Alternatives to consider:
- To solve a system of ordinary differential equations, use
  {func}`~.dsolve_system` *any reason to recommend dsolve_system instead of
  dsolve?*
- *which SciPy functions? from
  https://docs.scipy.org/doc/scipy/reference/integrate.html?highlight=ode?*

Here is an example of solving an ordinary differential equation algebraically:

```py
>>> from sympy import Function, dsolve, Derivative
>>> from sympy.abc import x
>>> y = Function('y')
>>> result = dsolve(Derivative(y(x), x, x) + 9*y(x), y(x)); result
Eq(y(x), C1*sin(3*x) + C2*cos(3*x))
```

You can then use SymPy to verify that the solution is correct:

```py
>>> from sympy import checkodesol
>>> solution = result.rhs; solution
C1*sin(3*x) + C2*cos(3*x)
>>> checkodesol(Derivative(y(x), x, x) + 9*y(x), solution)
(True, 0)
```

The output of {func}`~.checkodesol` is a tuple where the first item, a boolean,
tells whether the substitution results in `0`.

## Guidance

### Input Format

*Make a recommendation on Derivative() vs. diff()?* The one required input is
the differential equation(s). The derivatives, such as $y''(x)$, should be
expressed using {class}`~.Derivative` rather than {func}`~.diff`. Using
{func}`~.diff` is a common mistake 

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
>>> dsolve(eqs)
[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]
```

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
>>> dsolve(eqs)
[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]
>>> dsolve(eqs, ics={f(0): 1, g(0): 0})
[Eq(f(x), exp(x)/2 + exp(-x)/2), Eq(g(x), exp(x)/2 - exp(-x)/2)]
```

## Define a Function of an Independent Variable

You may prefer to specify a function (for example $x$) of its independent
variable (for example $t$):

```py
>>> from sympy import symbols, Function, dsolve
>>> t = symbols('t')
>>> x = Function('x')(t); x
x(t)
>>> xp = x.diff()
>>> xpp = xp.diff()
>>> eq = xpp + 2*xp + x; eq
x(t) + 2*Derivative(x(t), t) + Derivative(x(t), (t, 2))
>>> dsolve(eq, x)
Eq(x(t), (C1 + C2*t)*exp(-t))
```

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
>>> result = dsolve(Derivative(y(x), x, x) + 9*y(x), y(x)); result
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
>>> solutions = dsolve(eqs); solutions
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

### Work With Arbitrary Constants *include such a section?*

## *Tradeoffs (speed vs. accuracy, etc.) for function*

### *Tradeoff 1*

*Speed-up option 1 content*

### *Speed-up option 2*

*Speed-up option 2 content*

## Ordinary Differential Equation Type and Solving Strategy

*`hint` and `classifiy_ode`*

## Not All Equations Can Be Solved

### Equations With No Solution

*Equations with no solution content*

### Equations With No Analytical Solution

*Equations with no analytical solution content*

### Equations Which Have An Analytical Solution, and SymPy Cannot Solve

*Equations which have an analytical solution, and SymPy cannot solve content*

Please post the problem on the [mailing
list](https://groups.google.com/g/sympy), or open an issue on [SymPy's GitHub
page](https://github.com/sympy/sympy/issues). Until the issue is resolved, you
can *workaround*.

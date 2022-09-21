# Solve an Ordinary Differential Equation (ODE) Algebraically

Use SymPy to solve an ordinary differential equation (ODE) algebraically. For
example, solving $y''(x) + 9y(x)=0 $ yields $ y(x)=C_{1} \sin(3x)+ C_{2}
\cos(3x)$.

## Alternatives to Consider
- To numerically solve a system of ODEs, use a [SciPy ODE
  solver](https://docs.scipy.org/doc/scipy/reference/integrate.html#solving-initial-value-problems-for-ode-systems)
  such as `solve_ivp`. You can also use SymPy to create and then
  {func}`~.lambdify` an ODE to be solved numerically using SciPy's as
  `solve_ivp` as described below in [](#numerically-solve-an-ode-in-scipy).

## Solve an Ordinary Differential Equation (ODE)

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
tells whether substituting the solution into the ODE results in `0`, indicating
the solution is correct.

## Guidance

### Input Format

There are many ways to express derivatives of functions. For an undefined
function, both {class}`~.Derivative` and {func}`~.diff` represent the undefined
derivative. Thus, all of the following `ypp` ("y prime prime") represent $y''$,
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

Similarly, you must specify the argument of the function: $f(x)$, not just $f$.

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
the form of `{f(x0): y0, f(x).diff(x).subs(x, x1): y1}` and so on where, for
example, the value of `f(x)` at $x = $ `x0` is `y0`. For power series solutions,
if no initial conditions are specified `f(0)` is assumed to be `C0` and the
power series solution is calculated about $0$.

Here is an example of setting the initial values for functions, namely namely
$f(0) = 1$ and $g(2) = 3$:

```py
>>> from sympy import symbols, Eq, Function, dsolve
>>> f, g = symbols("f g", cls=Function)
>>> x = symbols("x")
>>> eqs = [Eq(f(x).diff(x), g(x)), Eq(g(x).diff(x), f(x))]
>>> dsolve(eqs, [f(x), g(x)])
[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]
>>> dsolve(eqs, [f(x), g(x)], ics={f(0): 1, g(2): 3})
[Eq(f(x), (1 + 3*exp(2))*exp(x)/(1 + exp(4)) - (-exp(4) + 3*exp(2))*exp(-x)/(1 + exp(4))), Eq(g(x), (1 + 3*exp(2))*exp(x)/(1 + exp(4)) + (-exp(4) + 3*exp(2))*exp(-x)/(1 + exp(4)))]
```

Here is an example of setting the initial value for the derivative of a
function, namely $f'(1) = 2$:

```py
>>> dsolve(eqs, [f(x), g(x)], ics={f(x).diff(x).subs(x, 1): 2})
[Eq(f(x), C2*exp(x) + (C2*exp(2) - 2*E)*exp(-x)), Eq(g(x), C2*exp(x) - (C2*exp(2) - 2*E)*exp(-x))]
```

## Define a Function of an Independent Variable

You may prefer to specify a function (for example $x$) of its independent
variable (for example $t$), so that `x`, represents `x(t)`:

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

### Specify Initial (Boundary) Conditions

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

### Numerically Solve an ODE in SciPy

A common workflow which leverages
[SciPy's](https://docs.scipy.org/doc/scipy/index.html) fast numerical ODE
solving is
1. set up an ODE in SymPy
2. convert it to a numerical function using {func}`~.lambdify`
3. solve the initial value problem by numerically integrating the ODE using
   SciPy's `solve_ivp`.

```{warning}
{func}`~.lambdify` uses {external:func}`~.exec` to dynamically execute Python code, and thus should not be used on unsanitized input.
```

Here is an [example from the field of chemical
kinetics](https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/25-chemical-kinetics-intro.html)
where the nonlinear ordinary differential equations take this form:

$$ r_f = & k_f y_0(t)^2  y_1(t) \\
r_b = & k_b y_2(t)^2 \\
\frac{d y_0(t)}{dt} = & 2(r_b - r_f) \\
\frac{d y_1(t)}{dt} = & r_b - r_f \\
\frac{d y_2(t)}{dt} = & 2(r_f - r_b) $$

and

$$

\vec{y}(t) = \begin{bmatrix} y_0(t) \\
y_1(t) \\
y_2(t) \end{bmatrix} $$

```py
>>> from sympy import symbols, lambdify
>>> import numpy as np
>>> import scipy.integrate
>>> import matplotlib.pyplot as plt
>>> def rhs(t, y, kf, kb):
...     rf = kf * y[0]**2 * y[1]
...     rb = kb * y[2]**2
...     return [2*(rb - rf), rb - rf, 2*(rf - rb)]
>>> y, (kf, kb) = symbols('y:3'), symbols('kf kb')
>>> ydot = rhs(None, y, kf, kb) # derivative of the function y(t); value given by function rhs (right-hand side) for input values y, kf, and kb
>>> ydot
[2*kb*y2**2 - 2*kf*y0**2*y1, kb*y2**2 - kf*y0**2*y1, -2*kb*y2**2 + 2*kf*y0**2*y1]
>>> t = symbols('t') # not used in this case
>>> f = lambdify((t, y, kf, kb), ydot) # convert the SymPy symbolic expression for ydot into a form that SciPy can evaluate numerically, f
>>> k_vals = np.array([0.42, 0.17]) # arbitrary in this case
>>> y0 = [1, 1, 0] # initial condition (initial values)
>>> t_eval = np.linspace(0, 10, 50) # evaluate integral from t = 0-10 for 50 points
>>> solution = scipy.integrate.solve_ivp(f, (0, 10), y0, t_eval=t_eval, args=k_vals) # call SciPy's solve_ivp by passing it the function f, the interval of integration, the initial state, and the arguments to pass to the function f
>>> solution
    {'message': 'The solver successfully reached the end of the integration interval.', 'nfev': 68, 'njev': 0, 'nlu': 0, 'sol': None, 'status': 0, 'success': True, 't': [ 0.          0.20408163  0.40816327  0.6122449   0.81632653  1.02040816
      1.2244898   1.42857143  1.63265306  1.83673469  2.04081633  2.24489796
      2.44897959  2.65306122  2.85714286  3.06122449  3.26530612  3.46938776
      3.67346939  3.87755102  4.08163265  4.28571429  4.48979592  4.69387755
      4.89795918  5.10204082  5.30612245  5.51020408  5.71428571  5.91836735
      6.12244898  6.32653061  6.53061224  6.73469388  6.93877551  7.14285714
      7.34693878  7.55102041  7.75510204  7.95918367  8.16326531  8.36734694
      8.57142857  8.7755102   8.97959184  9.18367347  9.3877551   9.59183673
      9.79591837 10.        ], 't_events': None, 'y': [[1.         0.8588458  0.76037109 0.68865638 0.63448858 0.59298188
      0.56070716 0.53512394 0.51478622 0.49863364 0.48571594 0.47520885
      0.4666383  0.45971742 0.45416064 0.44968802 0.44602514 0.4429604
      0.44043698 0.4383935  0.43676122 0.43546445 0.43442049 0.43353967
      0.43277066 0.43212494 0.43159916 0.43118497 0.43086892 0.43063257
      0.43045243 0.43029995 0.43014456 0.42999498 0.42986356 0.42975504
      0.42967192 0.42961451 0.42958085 0.42956676 0.42956583 0.42956944
      0.42956671 0.42954545 0.42951959 0.4294984  0.42948121 0.42946735
      0.42945611 0.4294468 ]
     [1.         0.9294229  0.88018554 0.84432819 0.81724429 0.79649094
      0.78035358 0.76756197 0.75739311 0.74931682 0.74285797 0.73760443
      0.73331915 0.72985871 0.72708032 0.72484401 0.72301257 0.7214802
      0.72021849 0.71919675 0.71838061 0.71773223 0.71721025 0.71676983
      0.71638533 0.71606247 0.71579958 0.71559248 0.71543446 0.71531629
      0.71522621 0.71514997 0.71507228 0.71499749 0.71493178 0.71487752
      0.71483596 0.71480726 0.71479042 0.71478338 0.71478292 0.71478472
      0.71478335 0.71477272 0.71475979 0.7147492  0.71474061 0.71473367
      0.71472805 0.7147234 ]
     [0.         0.1411542  0.23962891 0.31134362 0.36551142 0.40701812
      0.43929284 0.46487606 0.48521378 0.50136636 0.51428406 0.52479115
      0.5333617  0.54028258 0.54583936 0.55031198 0.55397486 0.5570396
      0.55956302 0.5616065  0.56323878 0.56453555 0.56557951 0.56646033
      0.56722934 0.56787506 0.56840084 0.56881503 0.56913108 0.56936743
      0.56954757 0.56970005 0.56985544 0.57000502 0.57013644 0.57024496
      0.57032808 0.57038549 0.57041915 0.57043324 0.57043417 0.57043056
      0.57043329 0.57045455 0.57048041 0.5705016  0.57051879 0.57053265
      0.57054389 0.5705532 ]], 'y_events': None}
>>> y = solution.y # Extract the y (concentration) values from SciPy solution result
>>> plt.plot(t_eval, y.T) # doctest: +SKIP
>>> plt.xlabel('t') # doctest: +SKIP
>>> plt.ylabel('Concentration') # doctest: +SKIP
>>> plt.legend(['NO', 'Br$_2$', 'NOBr'], shadow=True) # doctest: +SKIP
>>> plt.title('Chemical Kinetics') # doctest: +SKIP
>>> plt.show()
```

SciPy's `solve_ivp` returns a result containing `t` (time) points and
corresponding `y` (numerical function result) values for each initial point in
`y0`.

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

### Equations With No Closed-Form Solution

Some systems of differential equations have no closed-form solution because they
are chaotic, for example the [Lorenz
system](https://en.wikipedia.org/wiki/Lorenz_system#Overview) or a double
pendulum described by these two differential equations (simplified from
[ScienceWorld](https://scienceworld.wolfram.com/physics/DoublePendulum.html)):

$$ 2 \theta_1''(t) + \theta_2''(t) \cos(\theta_1-\theta_2) + \theta_2'^2(t)
\sin(\theta_1 - \theta_2) + 2g \sin(\theta_1) = 0 $$

$$ \theta_2''(t) + \theta_1''(t) \cos(\theta_1-\theta_2) - \theta_1'^2(t)
\sin(\theta_1 - \theta_2) + g \sin(\theta_2) = 0 $$

```py
>>> from sympy import symbols, Function, cos, sin, dsolve
>>> theta1, theta2 = symbols('theta1 theta2', cls=Function)
>>> g, t = symbols('g t')
>>> eq1 = 2*theta1(t).diff(t, t) + theta2(t).diff(t, t)*cos(theta1(t) - theta2(t)) + theta2(t).diff(t)**2*sin(theta1(t) - theta2(t)) + 2*g*sin(theta1(t))
>>> eq2 = theta2(t).diff(t, t) + theta1(t).diff(t, t)*cos(theta1(t) - theta2(t)) - theta1(t).diff(t)**2*sin(theta1(t) - theta2(t)) + g*sin(theta2(t))
>>> dsolve([eq1, eq2], [theta1(t), theta2(t)])
Traceback (most recent call last):
...
NotImplementedError
```

For such cases, you can solve the equations numerically as mentioned in
[](#alternatives-to-consider).

## Report a Problem

If you know your ODE has a solution, and SymPy cannot find it, please post the
problem on the [mailing list](https://groups.google.com/g/sympy), or open an
issue on [SymPy's GitHub page](https://github.com/sympy/sympy/issues). Until the
issue is resolved, you can try one of the [](#alternatives-to-consider).

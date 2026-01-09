(solving-guide-de-laplace)=
# Solve Ordinary Differential Equations (ODEs) and Partial Differential Equations (PDEs) Algebraically with the Laplace Transform

This guide shows ho to solve ODEs and PDEs using the Laplace Transform.
Using the Laplace transform implies that there is a boundary at $0$
with given boundary conditions.

For example, solving $y''(t) + 9y(t)=0$ for initial conditions $y(0)=y_0$ and
$y'(0)=v_0$ yields $y(t)=\theta(t)\cdot\left(\frac{v_0}{3} \sin(3t) + y_0 \cos(3t)\right)$.
All solutions will have a factor like $\theta(t)$ (the Heaviside Theta function, or unit step) in it, signifying that the use of the Laplace
transform implies that the solution is $0$ for $t<0$.

## Solve an Ordinary Differential Equation (ODE)

Here is an example of solving the above ordinary differential equation
algebraically using the Laplace transform. First, the ODE is written
using {func}`~.diff`. Then {func}`~.laplace_transform` is used
to transform `e1` from the $t$ domain into the $s$ domain. This can
be simplified to `e3` by replacing abstract transform expressions by
functions using {func}`~.laplace_correspondence` and by inserting
initial conditions with {func}`~.laplace_initial_conds`. The resulting
equation can then be transformed back to the $t$ domain with
{func}`~.inverse_laplace_transform`.

```py
>>> from sympy import (
...     diff, Function, inverse_laplace_transform, laplace_correspondence,
...     laplace_initial_conds, laplace_transform, solve, symbols)
>>> y, Y = symbols('y, Y', cls=Function)
>>> s = symbols('s')
>>> t, y0, v0 = symbols('t, y0, v0', real=True)
>>> e1 = diff(y(t), t, 2) + 9*y(t)
>>> e1
9*y(t) + Derivative(y(t), (t, 2))
>>> e2 = laplace_transform(e1, t, s, noconds=True)
>>> e2
s**2*LaplaceTransform(y(t), t, s) - s*y(0) + 9*LaplaceTransform(y(t), t, s) - Subs(Derivative(y(t), t), t, 0)
>>> e3 = laplace_initial_conds(laplace_correspondence(e2, {y: Y}), t, {y: [y0, v0]})
>>> e3
s**2*Y(s) - s*y0 - v0 + 9*Y(s)
>>> e4 = solve(e3, Y(s))
>>> e4
[(s*y0 + v0)/(s**2 + 9)]
>>> e5 = inverse_laplace_transform(e4[0], s, t)
>>> e5
(v0*sin(3*t)/3 + y0*cos(3*t))*Heaviside(t)

```

## Solve a Partial Differential Equation (PDE)

Next we want to solve the PDE $y{\left(x,t \right)} + \frac{\partial}{\partial t}
y{\left(x,t \right)} + \frac{\partial}{\partial x} y{\left(x,t \right)} = 0$
for $t\geq 0$, $x\geq 0$, and the initial condition $y(0, t) = f(t)$ with an
unknown function $f$ for which $f(0)=0$. The solution of this problem is
$f{\left(t - x \right)} e^{- x} \theta\left(t - x\right)$.

This is not as straightforward as solving ODEs because
{func}`~.laplace_correspondence` does not work with multivariate functions,
and because simplifying all results is not possible automatically due to
valid constraints which SymPy does not see.

The plan to solve this PDE is to use the Laplace transform twice, once from
the $t$ to the $s$ domain, and a second time from the $x$ to the $p$ domain.
The equation can then be solved directly with {func}`~.solve`, and the
solution in $x$ and $t$ can be obtained by two inverse Laplace transforms.

The example solution below works as follows: `e1` is the PDE above. `e2`
is the PDE with the $t$ axis transformed to the $s$ domain. In this step
we call the Laplace transform of $y(t, x)$ from $t$ to $s$ simply $Y(x)$,
omitting the $s$ variable in the function argument. Then `e3` is the same
with the $x$ axis transformed to the $p$ domain. After this step,
{func}`~.laplace_correspondence` can be used to rewrite the Laplace
transform of $Y(x)$ to $Z(p)$ (both with the $s$ variable not explicitely
written out), resulting in `e4`, which can be solved for $Z$, giving `e5`.

The first {func}`~.inverse_laplace_transform` from $p$ to $x$ gives the
result `e6`. It contains the expression `exp(-x*(s**2 + 2*s + 1)/(s + 1))`
which can be simplified if $s \neq 1$. The {func}`~.inverse_laplace_transform`
from $s$ to $t$ to be done later can easily be calculted with a convergence
plane that excludes the point $s=1$, but SymPy does not know this, so we
rewrite that expression `p1` to `p2` and substitute it, giving `e7`.

Then we apply the initial conditions: $Y(0)$ is the Laplace transform of $y(0, t) = f(t)$, and $f(0)=0$, which we substitute to obtain `e8`. This can then be transformed back to the $t$ domain, resulting in `e9`.

It is posible to simplify this further, but only by giving SymPy the
information that $x$ is strictly positive. So we define a variable `xp` that
is strictly positive, and reach the goal `e11` using {func}`~.expand` and
two substitutions.

```py
>>> from sympy import (
...     Derivative, diff, expand, Function, inverse_laplace_transform,
...     LaplaceTransform, laplace_correspondence, laplace_transform,
...     solve, symbols
... )
>>> f, F, y, Y, Z = symbols('f, F, y, Y, Z', cls=Function)
>>> s, p = symbols('s, p')
>>> t, x = symbols('t, x', real=True)
>>> xp = symbols('x_{positive}', positive=True)
>>> pde = diff(y(x, t), t) + diff(y(x, t), x) + y(x, t)
>>> pde
y(x, t) + Derivative(y(x, t), t) + Derivative(y(x, t), x)
>>> e1 = laplace_transform(pde, t, s, noconds=True)
>>> e1
s*LaplaceTransform(y(x, t), t, s) + LaplaceTransform(y(x, t), t, s) + LaplaceTransform(Derivative(y(x, t), x), t, s) - y(x, 0)
>>> e2 = (
...     e1.subs(LaplaceTransform(y(x, t), t, s), Y(x))
...         .subs(LaplaceTransform(Derivative(y(x, t), x), t, s), diff(Y(x), x))
...         .subs(y(x, 0), f(0)))
>>> e2
s*Y(x) + Y(x) - f(0) + Derivative(Y(x), x)
>>> e3 = laplace_transform(e2, x, p, noconds=True)
>>> e3
p*LaplaceTransform(Y(x), x, p) + s*LaplaceTransform(Y(x), x, p) + LaplaceTransform(Y(x), x, p) - Y(0) - f(0)/p
>>> e4 = laplace_correspondence(e3, {Y: Z})
>>> e4
p*Z(p) + s*Z(p) - Y(0) + Z(p) - f(0)/p
>>> e5 = solve(e4, Z(p))
>>> e5
[(p*Y(0) + f(0))/(p*(p + s + 1))]
>>> e6 = inverse_laplace_transform(e5[0], p, x)
>>> e6
(s*Y(0) + Y(0) - f(0))*exp(-x*(s**2 + 2*s + 1)/(s + 1))*Heaviside(x)/(s + 1) + f(0)*Heaviside(x)/(s + 1)
>>> p1 = e6.args[1].args[3]
>>> p1
exp(-x*(s**2 + 2*s + 1)/(s + 1))
>>> newargs = [a.factor() for a in p1.args]
>>> p2 = p1.func(*newargs)
>>> p2
exp(-x*(s + 1))
>>> e7 = e6.subs(p1, p2)
>>> e7
(s*Y(0) + Y(0) - f(0))*exp(-x*(s + 1))*Heaviside(x)/(s + 1) + f(0)*Heaviside(x)/(s + 1)
>>> e8 = e7.subs(Y(0), laplace_transform(f(t), t, s, noconds=True)).subs(f(0), 0)
>>> e8
(s*LaplaceTransform(f(t), t, s) + LaplaceTransform(f(t), t, s))*exp(-x*(s + 1))*Heaviside(x)/(s + 1)
>>> e9 = inverse_laplace_transform(e8, s, t)
>>> e9
InverseLaplaceTransform(LaplaceTransform(f(t), t, s)*exp(-x*(s + 1)), s, t, _None)*Heaviside(x)
>>> e10 = expand(e9.subs(x, xp))
>>> e10
InverseLaplaceTransform(LaplaceTransform(f(t), t, s)*exp(-x_{positive})*exp(-s*x_{positive}), s, t, _None)
>>> e11 = e10.doit().subs(xp, x)
>>> e11
f(t - x)*exp(-x)*Heaviside(t - x)
```

(solving-guide)=
# Solve Equations

The Python package SymPy can symbolically solve equations, differential
equations, linear equations, nonlinear equations, matrix problems, inequalities,
Diophantine equations, and evaluate integrals. SymPy can also solve numerically.

The [solving guide](solving-guidance.md) provides suggestions for many types of
solving tasks.

Learn how to use SymPy computer algebra system to:

| Description                                                  | Example                                                                                                                     | Solution |
|--------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------|
| [ Solve an equation algebraically ](solve-equation-algebraically.md)                        | $x^2 = y$ | $x \in \{-\sqrt{y},\sqrt{y}\}$                                                                                                |
| [ Solve a system of equations (linear or nonlinear) algebraically ](solve-system-of-equations-algebraically.md)              | $x^2 + y = 2z, y = -4z$ | $\{(x = -\sqrt{6z}, y = -4z),$ ${(x = \sqrt{6z}, y = -4z)\}}$                                                                                        |
|  [Solve one or a system of equations numerically](solve-numerically.md)                           | $\cos(x) = x $ | $ x \approx 0.739085133215161$                                                                                           |
|  {func}`Solve an ordinary differential equation algebraically <sympy.solvers.ode.dsolve>`   | $y''(x) + 9y(x)=0 $ | $ y(x)=C_{1} \sin(3x)+ C_{2} \cos(3x)$                                                    |
|  {func}`Solve a matrix problem algebraically <sympy.matrices.matrices.MatrixBase.solve>`                    | $ \left[\begin{array}{cc} 1 & 1\\1 & -1\end{array}\right] \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} 2\\0\end{array}\right] $ | $ \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} 1\\1\end{array}\right]$  |
|  [ Reduce one or a system of inequalities for a single variable algebraically ](reduce-inequalities-algebraically.md)           | $ x^2 < \pi, x > 0 $ | $ 0 < x < \sqrt{\pi} $                                                                                                  |
| [ Solve (find the roots of) a polynomial algebraically ](../../modules/polys/basics.rst)                       | $ x^2 - x = 0 $ | $ x \in \{0, 1\} $                                                                                                |
|  [ Solve a Diophantine equation (find integer solutions to a polynomial equation) algebraically ](../../modules/solvers/diophantine.rst)             | $x^2 - 4xy + 8y^2 - 3x + 7y - 5 = 0 $ | $ \{(x = 2, y = 1), (x = 5, y = 1)\}$                                                                                  |

Note: SymPy has a function called {func}`~.solve` which is designed to find the
roots of an equation or system of equations. SymPy {func}`~.solve` may or may
not be what you need for a particular problem, so we recommend you use the links
on this page to learn how to "solve" your problem. And while a common,
colloquial expression is, for example, ["solve an
integral"](../../modules/integrals/integrals.rst), in SymPy's terminology it
would be ["evaluate an integral."](../../modules/integrals/integrals.rst)

```{toctree}
:hidden: true

solving-guidance.md
solve-equation-algebraically.md
solve-system-of-equations-algebraically.md
solve-numerically.md
reduce-inequalities-algebraically.md
```

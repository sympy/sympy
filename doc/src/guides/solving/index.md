(solving-guide)=
# Solve Equations

The Python package SymPy can symbolically solve equations, differential
equations, linear equations, nonlinear equations, matrix problems, inequalities,
Diophantine equations, and evaluate integrals. SymPy can also solve numerically.

The [](solving-guidance.md) page provides recommendations applicable to many
types of solving tasks.

Learn how to use SymPy computer algebra system to:

| Description                                                  | Example                                                                                                                     | Solution |
|--------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------|
| [ Solve an equation algebraically ](solve-equation-algebraically.md)                        | $x^2 = y$ | $x \in \{-\sqrt{y},\sqrt{y}\}$                                                                                                |
| [ Solve a system of equations algebraically ](solve-system-of-equations-algebraically.md)              | $x^2 + y = 2z, y = -4z$ | $\{(x = -\sqrt{6z}, y = -4z),$ ${(x = \sqrt{6z}, y = -4z)\}}$                                                                                        |
|  [Solve one or a system of equations numerically](solve-numerically.md)                           | $\cos(x) = x $ | $ x \approx 0.739085133215161$                                                                                           |
|  [Solve an ordinary differential equation algebraically](solve-ode.md)                           | $y''(x) + 9y(x)=0 $ | $ y(x)=C_{1} \sin(3x)+ C_{2} \cos(3x)$                                                                                           |
| [ Find the roots of a polynomial algebraically or numerically ](find-roots-polynomial.md)                       | $ ax^2 + bx + c = 0 $ | $ x = \frac{-b\pm\sqrt{b^2 - 4ac}}{2a} $                                                                                                |
|  [ Solve a matrix equation algebraically ](solve-matrix-equation.md)                    | $ \left[\begin{array}{cc} c & d\\1 & -e\end{array}\right] \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} 2\\0\end{array}\right] $ | $ \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} \frac{2e}{ce+d}\\\frac{2}{ce+d}\end{array}\right]$  |
|  [ Reduce one or a system of inequalities for a single variable algebraically ](reduce-inequalities-algebraically.md)           | $ x^2 < \pi, x > 0 $ | $ 0 < x < \sqrt{\pi} $                                                                                                  |
|  [ Solve a Diophantine equation algebraically ](solve-diophantine-equation.md)             | $a^2 + b^2 = c^2$ | $(a=2pq, b=p^2-q^2, c=p^2+q^2)$                                                                                  |

Notes:
- SymPy has a function called {func}`~.solve` which is designed to find the
solutions of an equation or system of equations, or the roots of a function.
SymPy {func}`~.solve` may or may not be what you need for a particular problem,
so we recommend you use the links on this page to learn how to "solve" your
problem.
- While a common, colloquial expression is, for example, "[solve an
integral](../../modules/integrals/integrals.rst)," in SymPy's terminology it
would be "[evaluate an integral](../../modules/integrals/integrals.rst)." This
page does not provide guidance for such tasks. Please search the documentation
for the type of expression you want to evaluate.

```{toctree}
:hidden: true

solving-guidance.md
solve-equation-algebraically.md
solve-system-of-equations-algebraically.md
solve-numerically.md
solve-ode.md
find-roots-polynomial.md
solve-matrix-equation.md
reduce-inequalities-algebraically.md
solve-diophantine-equation.md
```

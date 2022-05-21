# Solve Using SymPy

The Python package SymPy can symbolically solve equations, differential equations, 
linear equations, nonlinear equations, matrix problems, inequalities, 
Diophantine equations, and evaluate integrals. SymPy can also solve numerically.

Learn how to use SymPy computer algebra system to:


| Description                                                  | Example                                                                                                                     |
|--------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| [ Solve an equation algebraically ](http://docs.sympy.org)                        | $x^2 = 4 \to x =\{-2, 2\} $                                                                                                |
| [ Solve a system of equations algebraically ](http://docs.sympy.org)              | $x^2 + y = 2, x - y = 4 \to \{(x = -3, y = -7), (x = 2, y = 2)\}$                                                                                        |
| [ Solve an equation numerically ](http://docs.sympy.org)                          | $\cos(x) = x \to x \approx 0.739085133215161$                                                                                           |
| [ Solve a partial differential equation algebraically ](http://docs.sympy.org)    | ${\frac {\partial ^{2}(x^3 y^3)}{\partial x\partial y}}=9 \to x \ne 0, y = \pm\frac{1}{x} $                                                                     |
| [ Solve an ordinary differential equation algebraically  ](http://docs.sympy.org) | $y''(x) + 9f(x)=0 \to f(x)=C_{1} \sin(3x)+ C_{2} \cos(3x)$                                                    |
| [ Solve a system of linear equations algebraically ](http://docs.sympy.org)       | $x + y = 2, x - y = 0 \to x = 1, y = 1$                                                                                           |
| [ Solve a system of nonlinear equations algebraically ](http://docs.sympy.org)    | $x^2 + y^3 = 1, x^3 - y^2 = 0 \to x = 1, y = 0$                                                                                       |
| [ Solve a matrix problem algebraically ](http://docs.sympy.org)                   | $ \left[\begin{array}{cc} 1 & 1\\1 & -1\end{array}\right] \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} 2\\0\end{array}\right] \to \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} 1\\1\end{array}\right]$  |
| [ Solve an inequality algebraically ](http://docs.sympy.org)                      | $ x^2 < 4 \to -2 < x < 2 $                                                                                                        |
| [ Solve a system of inequalities algebraically ](http://docs.sympy.org)           | $ x^2 < 4, x > 0 \to 0 < x < 2 $                                                                                                  |
| [ Solve (find the roots of) a polynomial algebraically ](http://docs.sympy.org)                       | $ x^2 - x = 0 \to x = \{0, 1\} $                                                                                                |
| [ Solve a Diophantine equation (find integer solutions to a polynomial equation) algebraically ](http://docs.sympy.org)             | $x^2 - 4xy + 8y^2 - 3x + 7y - 5 = 0 \to \{(x = 2, y = 1), (x = 5, y = 1)\}$                                                                                  |

| Description                                                  | Example                                                                                                                     |
|--------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| [ Solve an equation algebraically ](http://docs.sympy.org)                        | $x^2 = 4 $                                                                                                        |
| [ Solve a system of equations algebraically ](http://docs.sympy.org)              | $x^2 + y = 2, x - y = 4$                                                                                               |
| [ Solve an equation numerically ](http://docs.sympy.org)                          | $\cos(x) = x $                                                                                                      |
| [ Solve a partial differential equation algebraically ](http://docs.sympy.org)    | ${\frac {\partial ^{2}(x^3 y^3)}{\partial x\partial y}}=9 $                                                                     |
| [ Solve an ordinary differential equation algebraically  ](http://docs.sympy.org) | $y''(x) + 9f(x)=0 $                                                    |
| [ Solve a system of linear equations algebraically ](http://docs.sympy.org)       | $x + y = 2, x - y = 0$                                                                                            |
| [ Solve a system of nonlinear equations algebraically ](http://docs.sympy.org)    | $x^2 + y^3 = 1, x^3 - y^2 = 0 $                                                                                       |
| [ Solve a matrix problem algebraically ](http://docs.sympy.org)                   | $ \left[\begin{array}{cc} 1 & 1\\1 & -1\end{array}\right] \left[\begin{array}{cc} x\\y\end{array}\right] = \left[\begin{array}{cc} 2\\0\end{array}\right] $ |
| [ Solve an inequality algebraically ](http://docs.sympy.org)                      | $ x^2 < 4 $                                                                                                         |
| [ Solve a system of inequalities algebraically ](http://docs.sympy.org)           | $ x^2 < 4, x > 0 $                                                                                                  |
| [ Solve (find the roots of) a polynomial algebraically ](http://docs.sympy.org)                       | $ x^2 - x = 0 $                                                                                                 |
| [ Solve a Diophantine equation (find integer solutions to a polynomial equation) algebraically ](http://docs.sympy.org)             | $x^2 - 4xy + 8y^2 - 3x + 7y - 5 = 0 $                                                                                  |

Note: SymPy has a function called 
[`solve`](https://docs.sympy.org/dev/modules/solvers/solvers.html?highlight=solve#sympy.solvers.solvers.solve) 
which is designed to find the roots of an equation or system of equations.
SymPy `solve` may or may not be what you need for a particular problem, 
so we recommend you use the links on this page to learn how to "solve" your problem.
And while a common, colloquial expression is, for example, ["solve an integral"](http://docs.sympy.org), 
in SymPy's terminology it would be ["evaluate an integral."](http://docs.sympy.org)
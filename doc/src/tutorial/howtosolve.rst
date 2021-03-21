===========================================
 How to solve algebraic equations in Sympy
===========================================


Sympy solves algebraic equations with the help of the following functions

* `solve()`

* `solveset()`

* `linsolve()`

* `nonlinsolve()`

* `solve\_decomposition()`

* `solvify()`


Solve()
=============

This is the oldest and as of now, the most versatile solver in Sympy. It can solve polynomial, transcendental, piecewise combination of polynomial and transcendental, systems of linear and polynomial equations, systems containing relational expressions and expressions containing derivatives. Unlike `solveset()`, `solve()` can handle multivariate equations.

The return type is a Sympy expression and is, therefore, easy to handle. However, this comes at a cost. `Solve()` does not guarantee if it has returned all the solutions. It cannot handle cases with infinite solutions.

code::

    >>> from sympy.solvers import solve
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> solve(sin(x),x)
    [0, pi]
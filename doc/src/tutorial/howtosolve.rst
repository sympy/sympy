===================================
 How to solve equations in Sympy
===================================


Sympy solves equations with the help of the following functions

* `solve()`

* `solveset()`

* `nsolve()`

* `linsolve()`

* `nonlinsolve()`

* `solve\_decomposition()`

* `solvify()`

* `rsolve()`

* `diophantine`

* `solve\_univariate\_inequality`

* `roots`

* `Poly.nroots`

* `RootOf`


Numeric and symbolic solutions
+++++++++++++++++++++++++++++++

Sympy has solver functions to return both numeric and symbolic solutions for equations. ``nsolve`` is the main solver to return numeric solutions, whereas ``nroots`` computes numerical approximations of the roots. The rest of the solver functions are for returning symbolic equations.

Since ``nsolve`` cannot solve multivariate equations, ``scipy.optimize.fsolve``, which can handle such cases, is a possible alternative to it.

The solver functions are discussed in detail below.


solve()
=============

This is the oldest and as of now, the most versatile solver in Sympy. It can solve polynomial, transcendental, piecewise combination of polynomial and transcendental, systems of linear and polynomial equations, systems containing relational expressions and expressions containing derivatives. Unlike ``solveset``, ``solve`` can handle multivariate equations.

The return type is a Sympy expression and is, therefore, easy to handle. However, this comes at a cost. ``solve`` does not guarantee if it has returned all the solutions. It cannot handle cases with infinite solutions.

    >>> from sympy.solvers import solve
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> solve(sin(x),x)
    [0, pi]


solveset()
===========

``solveset`` was added in `this project <https://github.com/sympy/sympy/wiki/GSoC-2014-Application-Harsh-Gupta:-Solvers>`_ for GSoC 2014. It has a cleaner input and output interface as compared to ``solve``. The return type is a set object. Unlike ``solve``, which does not guarantee all solutions, ``solveset`` returns a ``ConditionSet`` object with partial solution for cases it doesn't "know". It can also return infinite solutions.

The disadvantage with ``solveset`` is that it can handle only univariate equations.

    >>> from sympy import solveset, Symbol, cos
    >>> x=Symbol('x')
    >>> solveset(cos(x),x)
    Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers)


nsolve()
=========

``nsolve`` uses ``mpmath.findroot`` to figure out the numeric solutions of univariate equations. If the user wants to numerical solutions, then ``nsolve`` should be used. Refer to its `documentation <https://docs.sympy.org/latest/modules/solvers/solvers.html#sympy.solvers.solvers.nsolve>`_ where it has been discussed in detail.


linsolve()
===========

``linsolve`` is able to solve both overdetermined and underdetermined systems of linear equation. However, the return types are incompatible for set operations unlike ``solveset``. For example, in the following code, the set object returned, a tuple wrapped inside a ``FiniteSet``, cannot be worked with.

    >>> from sympy.solvers import linsolve
    >>> from sympy import Symbol
    >>> from sympy.abc import x,y
    >>> linsolve([x + y], [x, y])
    FiniteSet((-y, y))


nonlinsolve()
==============

``nonlinsolve`` is able to solve both overdetermined and underdetermined systems of nonlinear equation. However, the return types are more awkward than ``linsolve``. In the following code, the returned solution is a ``FiniteSet`` containing two tuples, of the two elements of which one is an ``Expr`` and the other is an ``ImageSet``.

    >>> from sympy.solvers import linsolve
    >>> from sympy import Symbol
    >>> from sympy.abc import x,y
    >>> nonlinsolve([sin(y)], [x, y]) 
    FiniteSet((x, ImageSet(Lambda(_n, 2*_n*pi), Integers)), (x, ImageSet(Lambda(_n, 2*_n*pi + pi), Integers)))

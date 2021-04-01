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


solve_decomposition()
======================

``solve_decomposition`` uses the principle of rewriting and decomposition to solve equations.

For instance, if we have to solve an equation ``f(x) = 0`` where ``f(x) = 4**x - 3(2**x) + 2``, by the method of decomposition ``f(x)`` is broken converted to ``g(x)**2 - 3*g(x) + 2`` where ``g(x) = 2 **x``.

Substituting ``g(x)`` with a dummy variable ``t``, we get a polynomial in ``f(t)`` given as

``f(t) = (t - 1)*(t - 2)``

Resubstituting the value of ``t`` as ``g(x)`` and further composing ``f(x)`` from ``g(x)``, we get

``f(x) = (2**x - 1)*(2**x - 2)``

We can now represent ``f(x)`` as the product of two simpler functions by the idea of rewriting

``f(x) = k(x) * l(x)``

Calculating the solutions of these functions is comparatively easier with respect to the original equation. We solve for these functions separately and get individual solutions as

``2**x - 1 = 0 or 2**x - 2 = 0``

that is, ``x = 0 or 1``

The disadvantage with ``solve_decomposition`` is that it cannot solve multivariate equations.

    >>> from sympy.solvers.solveset import solve_decomposition
    >>> from sympy import Symbol
    >>> from sympy.abc import x,y
    >>> solve_decomposition(4**x-3*(2**x)+2, x,S.Reals)
    FiniteSet(0, 1)


solvify()
==========

``solvify`` is essentially a wrapper around ``solveset`` such that the return value is the same as the ``solve`` API.

    >>> from sympy.solvers.solveset import solvify
    >>> from sympy import Symbol
    >>> from sympy.abc import x,y
    >>> solvify(x**2 - 9, x, S.Reals)
    [-3, 3]


rsolve()
=========

``rsolve()`` solves univariate recurrence relations.

    >>> from sympy import Function, rsolve
    >>> from sympy.abc import n
    >>> y = Function('y')
    >>> f = (n - 1)*y(n + 2) - (n**2 + 3*n - 2)*y(n + 1) + 2*n*(n + 1)*y(n)
    >>> rsolve(f, y(n))
    2**n*C0 + C1*factorial(n)


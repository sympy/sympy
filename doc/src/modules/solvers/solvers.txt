Solvers
==========

.. module:: sympy.solvers

The *solvers* module in SymPy implements methods for solving equations.

Algebraic equations
--------------------

Use :func:`solve` to solve algebraic equations. We suppose all equations are equaled to 0,
so solving x**2 == 1 translates into the following code::

    >>> from sympy.solvers import solve
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> solve(x**2 - 1, x)
    [-1, 1]

The first argument for :func:`solve` is an equation (equaled to zero) and the second argument
is the symbol that we want to solve the equation for.

.. autofunction:: sympy.solvers.solvers.solve

.. autofunction:: sympy.solvers.solvers.solve_linear

.. autofunction:: sympy.solvers.solvers.solve_linear_system

.. autofunction:: sympy.solvers.solvers.solve_linear_system_LU

.. autofunction:: sympy.solvers.solvers.solve_undetermined_coeffs

.. autofunction:: sympy.solvers.solvers.tsolve

.. autofunction:: sympy.solvers.solvers.nsolve

.. autofunction:: sympy.solvers.solvers.check_assumptions

.. autofunction:: sympy.solvers.solvers.checksol

Ordinary Differential equations (ODEs)
--------------------------------------

See :ref:`ode-docs`.

Partial Differential Equations (PDEs)
-------------------------------------

.. autofunction:: sympy.solvers.pde.pde_separate

.. autofunction:: sympy.solvers.pde.pde_separate_add

.. autofunction:: sympy.solvers.pde.pde_separate_mul

Recurrence Equtions
-------------------

.. autofunction:: sympy.solvers.recurr.rsolve

.. autofunction:: sympy.solvers.recurr.rsolve_poly

.. autofunction:: sympy.solvers.recurr.rsolve_ratio

.. autofunction:: sympy.solvers.recurr.rsolve_hyper

Systems of Polynomial Equations
-------------------------------

.. autofunction:: sympy.solvers.polysys.solve_poly_system

.. autofunction:: sympy.solvers.polysys.solve_triangulated


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

.. autofunction:: sympy.solvers.solvers.nsolve

.. autofunction:: sympy.solvers.solvers.check_assumptions

.. autofunction:: sympy.solvers.solvers.checksol

Ordinary Differential equations (ODEs)
--------------------------------------

See :ref:`ode-docs`.

Partial Differential Equations (PDEs)
-------------------------------------

See :ref:`pde-docs`.

Deutils (Utilities for solving ODE's and PDE's)
-----------------------------------------------

.. autofunction:: sympy.solvers.deutils.ode_order

Recurrence Equations
--------------------

.. module:: sympy.solvers.recurr

.. autofunction:: rsolve

.. autofunction:: rsolve_poly

.. autofunction:: rsolve_ratio

.. autofunction:: rsolve_hyper

Systems of Polynomial Equations
-------------------------------

.. autofunction:: sympy.solvers.polysys.solve_poly_system

.. autofunction:: sympy.solvers.polysys.solve_triangulated

Diophantine Equations (DEs)
---------------------------

See :ref:`diophantine-docs`

Inequalities
------------

See :ref:`inequality-docs`

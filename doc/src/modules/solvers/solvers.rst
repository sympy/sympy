.. _solvers:

Solvers
=======

.. module:: sympy.solvers

The *solvers* module in SymPy implements methods for solving equations.

.. note::

   It is recommended to use :func:`solveset` to solve univariate equations,
   :func:`~.linsolve` to solve system of linear equations
   instead of :func:`~sympy.solvers.solvers.solve` and :func:`~.nonlinsolve` to
   solve system of non linear equations since sooner or later the :func:`~.solveset`
   will take over :func:`~sympy.solvers.solvers.solve` either internally or externally.


Algebraic equations
--------------------

Use :func:`~sympy.solvers.solvers.solve` to solve algebraic equations. We suppose all equations are equaled to 0,
so solving x**2 == 1 translates into the following code::

    >>> from sympy.solvers import solve
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> solve(x**2 - 1, x)
    [-1, 1]

The first argument for :func:`~sympy.solvers.solvers.solve` is an equation (equaled to zero) and the second argument
is the symbol that we want to solve the equation for.

.. autofunction:: sympy.solvers.solvers::solve

.. autofunction:: sympy.solvers.solvers::solve_linear

.. autofunction:: sympy.solvers.solvers::solve_linear_system

.. autofunction:: sympy.solvers.solvers::solve_linear_system_LU

.. autofunction:: sympy.solvers.solvers::solve_undetermined_coeffs

.. autofunction:: sympy.solvers.solvers::nsolve

.. autofunction:: sympy.solvers.solvers::checksol

.. autofunction:: sympy.solvers.solvers::unrad

Ordinary Differential equations (ODEs)
--------------------------------------

See :ref:`ode-docs`.

Partial Differential Equations (PDEs)
-------------------------------------

See :ref:`pde-docs`.

Deutils (Utilities for solving ODE's and PDE's)
-----------------------------------------------

.. autofunction:: sympy.solvers.deutils::ode_order

Recurrence Equations
--------------------

.. module:: sympy.solvers.recurr

.. autofunction:: rsolve

.. autofunction:: rsolve_poly

.. autofunction:: rsolve_ratio

.. autofunction:: rsolve_hyper

Systems of Polynomial Equations
-------------------------------

.. autofunction:: sympy.solvers.polysys::solve_poly_system

.. autofunction:: sympy.solvers.polysys::solve_triangulated

Diophantine Equations (DEs)
---------------------------

See :ref:`diophantine-docs`

Inequalities
------------

See :ref:`inequality-docs`

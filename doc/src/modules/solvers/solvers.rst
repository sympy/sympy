.. _solvers-docs:

Solvers
=======

.. module:: sympy.solvers

The *solvers* module in SymPy implements methods for solving equations.

.. note::

   For a beginner-friendly guide focused on solving common types of equations,
   refer to :ref:`solving-guide`.

.. note::

   :func:`~sympy.solvers.solvers.solve` is an older more mature general
   function for solving many types of equations.
   :func:`~sympy.solvers.solvers.solve` has many options and uses different
   methods internally to determine what type of equations you pass it, so if
   you know what type of equation you are dealing with you may want to use the
   newer :func:`solveset` which solves univariate equations, :func:`~.linsolve`
   which solves system of linear equations, and :func:`~.nonlinsolve` which
   solves systems of non linear equations.

.. _solvers-algebraic-equations:

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

Linear Programming (Optimization)
---------------------------------

.. module:: sympy.solvers.simplex

.. autofunction:: lpmax

.. autofunction:: lpmin

.. autofunction:: linprog

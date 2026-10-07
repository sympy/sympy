======================
Solver Functions
======================

.. currentmodule:: sympy.solvers.solvers

This page provides a curated list of functions available in SymPy's ``solvers.solvers`` module.

.. note::
   For a beginner-friendly guide focused on solving common types of equations,
   refer to :ref:`solving-guide`.

Functions
---------

.. autofunction:: sympy.solvers.solvers::solve

.. autofunction:: sympy.solvers.solvers::solve_linear

.. autofunction:: sympy.solvers.solvers::solve_linear_system

.. autofunction:: sympy.solvers.solvers::solve_linear_system_LU

.. autofunction:: sympy.solvers.solvers::minsolve_linear_system

.. autofunction:: sympy.solvers.solvers::checksol

.. autofunction:: sympy.solvers.solvers::check_assumptions

.. autofunction:: sympy.solvers.solvers::nsolve

.. autofunction:: sympy.solvers.solvers::solve_undetermined_coeffs

Related Functions
----------------

For other types of solving functionality, see:

* :func:`~sympy.solvers.solveset.solveset` - An improved version of ``solve`` that returns a set of solutions
* :func:`~sympy.solvers.inequalities.solve_univariate_inequality` - For solving inequalities
* :func:`~sympy.solvers.ode.dsolve` - For solving differential equations
* :func:`~sympy.solvers.pde.pdsolve` - For solving partial differential equations


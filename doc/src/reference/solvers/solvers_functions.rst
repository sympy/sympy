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

.. autosummary::
   :toctree: ./generated/
   :template: function.rst

   solve
   solve_linear
   solve_linear_system
   solve_linear_system_LU
   minsolve_linear_system
   checksol
   check_assumptions
   check_solutions
   nsolve
   solve_undetermined_coeffs
   
Related Functions
----------------

For other types of solving functionality, see:

* :func:`~sympy.solvers.solveset.solveset` - An improved version of ``solve`` that returns a set of solutions
* :func:`~sympy.solvers.inequalities.solve_univariate_inequality` - For solving inequalities
* :func:`~sympy.solvers.ode.dsolve` - For solving differential equations
* :func:`~sympy.solvers.pde.pdsolve` - For solving partial differential equations


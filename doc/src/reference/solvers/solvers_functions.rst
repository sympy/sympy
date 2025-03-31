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

.. toctree::
   :hidden:

   generated/sympy.solvers.solvers.solve.rst
   generated/sympy.solvers.solvers.solve_linear.rst
   generated/sympy.solvers.solvers.solve_linear_system.rst
   generated/sympy.solvers.solvers.solve_linear_system_LU.rst
   generated/sympy.solvers.solvers.minsolve_linear_system.rst
   generated/sympy.solvers.solvers.checksol.rst
   generated/sympy.solvers.solvers.check_assumptions.rst
   generated/sympy.solvers.solvers.check_solutions.rst
   generated/sympy.solvers.solvers.nsolve.rst
   generated/sympy.solvers.solvers.solve_undetermined_coeffs.rst
   generated/sympy.solvers.solvers.nsolve_examples.rst
   generated/sympy.solvers.solvers.solve_linear_system_examples.rst
   generated/sympy.solvers.solvers.solve_examples.rst


========================
Solveset Functions
========================

.. currentmodule:: sympy.solvers.solveset

This page provides a curated list of functions available in SymPy's ``solvers.solveset`` module, which introduces a newer solving framework.

.. note::
   For a beginner-friendly guide focused on solving common types of equations,
   refer to :ref:`solving-guide`.

Functions
---------

.. autosummary::
   :toctree: ./generated/
   :template: function.rst

   solveset
   linsolve
   nonlinsolve
   substitution
   linear_eq_to_matrix
   
Related Functions
----------------

For other types of solving functionality, see:

* :func:`~sympy.solvers.solvers.solve` - The original solve function (use solveset when possible)
* :func:`~sympy.solvers.inequalities.solve_univariate_inequality` - For solving inequalities
* :func:`~sympy.solvers.ode.dsolve` - For solving differential equations
* :func:`~sympy.solvers.pde.pdsolve` - For solving partial differential equations 
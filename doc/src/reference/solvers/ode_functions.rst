====================
ODE Functions
====================

.. currentmodule:: sympy.solvers.ode

This page provides a curated list of functions available in SymPy's ``solvers.ode`` module.

.. note::
   For a beginner-friendly guide focused on solving ordinary differential equations,
   refer to :ref:`ode-docs`.

Functions
---------

.. autosummary::
   :toctree: ./generated/
   :template: function.rst

   dsolve
   classify_ode
   checkodesol
   homogeneous_order
   infinitesimals
   checkinfsol

Related Functions
----------------

For other types of solving functionality, see:

* :func:`~sympy.solvers.solvers.solve` - For solving algebraic equations
* :func:`~sympy.solvers.solveset.solveset` - An improved version of ``solve`` that returns a set of solutions
* :func:`~sympy.solvers.inequalities.solve_univariate_inequality` - For solving inequalities
* :func:`~sympy.solvers.pde.pdsolve` - For solving partial differential equations

.. toctree::
   :hidden:

   generated/sympy.solvers.ode.dsolve.rst
   generated/sympy.solvers.ode.classify_ode.rst
   generated/sympy.solvers.ode.checkodesol.rst
   generated/sympy.solvers.ode.homogeneous_order.rst
   generated/sympy.solvers.ode.infinitesimals.rst
   generated/sympy.solvers.ode.checkinfsol.rst
   generated/sympy.solvers.ode.dsolve_examples.rst


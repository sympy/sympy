.. _pde-docs:

PDE
===

.. module::sympy.solvers.pde

User Functions
--------------
These are functions that are imported into the global namespace with ``from sympy import *``.  They are intended for user use.

.. autofunction:: sympy.solvers.pde::pde_separate

.. autofunction:: sympy.solvers.pde::pde_separate_add

.. autofunction:: sympy.solvers.pde::pde_separate_mul

.. autofunction:: sympy.solvers.pde::pdsolve

.. autofunction:: sympy.solvers.pde::classify_pde

.. autofunction:: sympy.solvers.pde::checkpdesol

Hint Methods
------------
These functions are meant for internal use. However they contain useful information on the various solving methods.

.. autofunction:: sympy.solvers.pde::pde_1st_linear_constant_coeff_homogeneous

.. autofunction:: sympy.solvers.pde::pde_1st_linear_constant_coeff

.. autofunction:: sympy.solvers.pde::pde_1st_linear_variable_coeff

Information on the pde module
-----------------------------

.. automodule:: sympy.solvers.pde

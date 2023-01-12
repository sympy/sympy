.. _inequality-docs:

Inequality Solvers
==================

For general cases :func:`~.reduce_inequalities` should be used. Other functions
are the subcategories useful for special dedicated operations, and will be
called internally as needed by ``reduce_inequalities``.

.. note::

   For a beginner-friendly guide focused on solving inequalities, refer to
   :ref:`solving-guide-inequalities`.

.. note::

   Some of the examples below use :func:`~.poly`, which simply transforms an
   expression into a polynomial; it does not change the mathematical meaning of
   the expression.

.. module:: sympy.solvers.inequalities

.. autofunction:: solve_rational_inequalities

.. autofunction:: solve_poly_inequality

.. autofunction:: solve_poly_inequalities

.. autofunction:: reduce_rational_inequalities

.. autofunction:: reduce_abs_inequality

.. autofunction:: reduce_abs_inequalities

.. autofunction:: reduce_inequalities

.. autofunction:: solve_univariate_inequality

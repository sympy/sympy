.. _ode-docs:

ODE
===

.. module::sympy.solvers.ode

User Functions
--------------
These are functions that are imported into the global namespace with ``from sympy import *``.  They are intended for user use.

:func:`dsolve`
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.dsolve

:func:`classify_ode`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.classify_ode

:func:`ode_order`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_order

:func:`checkodesol`
^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.checkodesol

:func:`homogeneous_order`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.homogeneous_order

Hint Methods
------------
These functions are intended for internal use by :func:`dsolve` and others.  Nonetheless, they contain useful information in their docstrings on the various ODE solving methods.

:obj:`preprocess`
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.preprocess

:obj:`odesimp`
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.odesimp

:obj:`constant_renumber`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.constant_renumber

:obj:`constantsimp`
^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.constantsimp

:obj:`sol_simplicity`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_sol_simplicity

:obj:`1st_exact`
^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_exact

:obj:`1st_homogeneous_coeff_best`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_homogeneous_coeff_best

:obj:`1st_homogeneous_coeff_subs_dep_div_indep`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_homogeneous_coeff_subs_dep_div_indep

:obj:`1st_homogeneous_coeff_subs_indep_div_dep`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_homogeneous_coeff_subs_indep_div_dep

:obj:`1st_linear`
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_linear

:obj:`Bernoulli`
^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_Bernoulli

:obj:`Liouville`
^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_Liouville

:obj:`Riccati_special_minus2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_Riccati_special_minus2

:obj:`nth_linear_constant_coeff_homogeneous`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_nth_linear_constant_coeff_homogeneous

:obj:`nth_linear_constant_coeff_undetermined_coefficients`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_nth_linear_constant_coeff_undetermined_coefficients

:obj:`nth_linear_constant_coeff_variation_of_parameters`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_nth_linear_constant_coeff_variation_of_parameters

:obj:`separable`
^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_separable

Information on the ode module
-----------------------------

.. automodule:: sympy.solvers.ode

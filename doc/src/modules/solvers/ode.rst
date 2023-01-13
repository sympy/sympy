.. _ode-docs:

ODE
===

.. note::

   For a beginner-friendly guide focused on solving ODEs, refer to
   :ref:`solving-guide-ode`.

.. module::sympy.solvers.ode

.. automodule:: sympy.solvers.ode


User Functions
--------------
These are functions that are imported into the global namespace with ``from
sympy import *``.  These functions (unlike `Hint Functions`_, below) are
intended for use by ordinary users of SymPy.

.. autofunction:: sympy.solvers.ode::dsolve

.. autofunction:: sympy.solvers.ode.systems::dsolve_system

.. autofunction:: sympy.solvers.ode::classify_ode

.. autofunction:: sympy.solvers.ode::checkodesol

.. autofunction:: sympy.solvers.ode::homogeneous_order

.. autofunction:: sympy.solvers.ode::infinitesimals

.. autofunction:: sympy.solvers.ode::checkinfsol

.. autofunction:: sympy.solvers.ode::constantsimp

.. _hints:

Hint Functions
--------------
These functions are intended for internal use by
:py:meth:`~sympy.solvers.ode.dsolve` and others.  Unlike `User Functions`_,
above, these are not intended for every-day use by ordinary SymPy users.
Instead, functions such as :py:meth:`~sympy.solvers.ode.dsolve` should be used.
Nonetheless, these functions contain useful information in their docstrings on
the various ODE solving methods. For this reason, they are documented here.

.. autodata:: sympy.solvers.ode::allhints

.. autofunction:: sympy.solvers.ode.ode::odesimp

.. autofunction:: sympy.solvers.ode.ode::constant_renumber

.. autofunction:: sympy.solvers.ode.ode::ode_sol_simplicity

.. autoclass:: sympy.solvers.ode.single::Factorable
   :members:

.. autoclass:: sympy.solvers.ode.single.FirstExact
   :members:

.. autoclass:: sympy.solvers.ode.single::HomogeneousCoeffBest
   :members:

.. autoclass:: sympy.solvers.ode.single::HomogeneousCoeffSubsDepDivIndep
   :members:

.. autoclass:: sympy.solvers.ode.single::HomogeneousCoeffSubsIndepDivDep
   :members:

.. autoclass:: sympy.solvers.ode.single::FirstLinear
   :members:

.. autoclass:: sympy.solvers.ode.single::RationalRiccati
   :members:

.. autoclass:: sympy.solvers.ode.single::SecondLinearAiry
   :members:

.. autoclass:: sympy.solvers.ode.single::SecondLinearBessel
   :members:

.. autoclass:: sympy.solvers.ode.single::Bernoulli
   :members:

.. autoclass:: sympy.solvers.ode.single::Liouville
   :members:

.. autoclass:: sympy.solvers.ode.single::RiccatiSpecial
   :members:

.. autoclass:: sympy.solvers.ode.single::NthLinearConstantCoeffHomogeneous
   :members:

.. autoclass:: sympy.solvers.ode.single::NthLinearConstantCoeffUndeterminedCoefficients
   :members:

.. autoclass:: sympy.solvers.ode.single::NthLinearConstantCoeffVariationOfParameters
   :members:

.. autoclass:: sympy.solvers.ode.single::NthLinearEulerEqHomogeneous
   :members:

.. autoclass:: sympy.solvers.ode.single::NthLinearEulerEqNonhomogeneousVariationOfParameters
   :members:

.. autoclass:: sympy.solvers.ode.single::NthLinearEulerEqNonhomogeneousUndeterminedCoefficients
   :members:

.. autoclass:: sympy.solvers.ode.single::NthAlgebraic
   :members:

.. autoclass:: sympy.solvers.ode.single::NthOrderReducible
   :members:

.. autoclass:: sympy.solvers.ode.single::Separable
   :members:

.. autoclass:: sympy.solvers.ode.single::AlmostLinear
   :members:

.. autoclass:: sympy.solvers.ode.single::LinearCoefficients
   :members:

.. autoclass:: sympy.solvers.ode.single::SeparableReduced
   :members:

.. autoclass:: sympy.solvers.ode.single::LieGroup
   :members:

.. autoclass:: sympy.solvers.ode.single::SecondHypergeometric
   :members:

.. autofunction:: sympy.solvers.ode.ode::ode_1st_power_series

.. autofunction:: sympy.solvers.ode.ode::ode_2nd_power_series_ordinary

.. autofunction:: sympy.solvers.ode.ode::ode_2nd_power_series_regular

Lie heuristics
--------------
These functions are intended for internal use of the Lie Group Solver.
Nonetheless, they contain useful information in their docstrings on the
algorithms implemented for the various heuristics.

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco1_simple

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco1_product

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_bivariate

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_chi

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco2_similar

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_function_sum

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco2_unique_unknown

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco2_unique_general

.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_linear

Rational Riccati Solver
-----------------------
These functions are intended for internal use to solve a first order Riccati
differential equation with atleast one rational particular solution.

.. autofunction:: sympy.solvers.ode.riccati::riccati_normal

.. autofunction:: sympy.solvers.ode.riccati::riccati_inverse_normal

.. autofunction:: sympy.solvers.ode.riccati::riccati_reduced

.. autofunction:: sympy.solvers.ode.riccati::construct_c

.. autofunction:: sympy.solvers.ode.riccati::construct_d

.. autofunction:: sympy.solvers.ode.riccati::rational_laurent_series

.. autofunction:: sympy.solvers.ode.riccati::compute_m_ybar

.. autofunction:: sympy.solvers.ode.riccati::solve_aux_eq

.. autofunction:: sympy.solvers.ode.riccati::remove_redundant_sols

.. autofunction:: sympy.solvers.ode.riccati::get_gen_sol_from_part_sol

.. autofunction:: sympy.solvers.ode.riccati::solve_riccati

System of ODEs
--------------
These functions are intended for internal use by
:py:meth:`~sympy.solvers.ode.dsolve` for system of differential equations.

.. autofunction:: sympy.solvers.ode.ode::_linear_2eq_order1_type6

.. autofunction:: sympy.solvers.ode.ode::_linear_2eq_order1_type7

.. autofunction:: sympy.solvers.ode.systems::linear_ode_to_matrix

.. autofunction:: sympy.solvers.ode.systems::canonical_odes

.. autofunction:: sympy.solvers.ode.systems::linodesolve_type

.. autofunction:: sympy.solvers.ode.systems::matrix_exp_jordan_form

.. autofunction:: sympy.solvers.ode.systems::matrix_exp

.. autofunction:: sympy.solvers.ode.systems::linodesolve

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type1

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type2

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type3

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type4

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type5

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type1

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type2

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type3

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type4

.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type5

Information on the ode module
-----------------------------

.. automodule:: sympy.solvers.ode.ode

These functions are not intended for end-user use.

.. autofunction:: sympy.solvers.ode.ode::_handle_Integral

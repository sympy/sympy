.. _ode-docs:

ODE
===

.. module::sympy.solvers.ode

.. automodule:: sympy.solvers.ode

User Functions
--------------
These are functions that are imported into the global namespace with ``from
sympy import *``.  These functions (unlike `Hint Functions`_, below) are
intended for use by ordinary users of SymPy.

dsolve
^^^^^^
.. autofunction:: sympy.solvers.ode::dsolve

dsolve_system
^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::dsolve_system

classify_ode
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode::classify_ode

checkodesol
^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode::checkodesol

homogeneous_order
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode::homogeneous_order

infinitesimals
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode::infinitesimals

checkinfsol
^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode::checkinfsol

constantsimp
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode::constantsimp

Hint Functions
--------------
These functions are intended for internal use by
:py:meth:`~sympy.solvers.ode.dsolve` and others.  Unlike `User Functions`_,
above, these are not intended for every-day use by ordinary SymPy users.
Instead, functions such as :py:meth:`~sympy.solvers.ode.dsolve` should be used.
Nonetheless, these functions contain useful information in their docstrings on
the various ODE solving methods. For this reason, they are documented here.

allhints
^^^^^^^^
.. autodata:: sympy.solvers.ode::allhints

odesimp
^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::odesimp

constant_renumber
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::constant_renumber

sol_simplicity
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::ode_sol_simplicity

factorable
^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::Factorable
   :members:

1st_exact
^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single.FirstExact
   :members:

1st_homogeneous_coeff_best
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::HomogeneousCoeffBest
   :members:

1st_homogeneous_coeff_subs_dep_div_indep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::HomogeneousCoeffSubsDepDivIndep
   :members:

1st_homogeneous_coeff_subs_indep_div_dep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::HomogeneousCoeffSubsIndepDivDep
   :members:

1st_linear
^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::FirstLinear
   :members:

1st_rational_riccati
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::RationalRiccati
   :members:

2nd_linear_airy
^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::SecondLinearAiry
   :members:

2nd_linear_bessel
^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::SecondLinearBessel
   :members:

Bernoulli
^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::Bernoulli
   :members:

Liouville
^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::Liouville
   :members:

Riccati_special_minus2
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::RiccatiSpecial
   :members:

nth_linear_constant_coeff_homogeneous
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthLinearConstantCoeffHomogeneous
   :members:

nth_linear_constant_coeff_undetermined_coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthLinearConstantCoeffUndeterminedCoefficients
   :members:

nth_linear_constant_coeff_variation_of_parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthLinearConstantCoeffVariationOfParameters
   :members:

nth_linear_euler_eq_homogeneous
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthLinearEulerEqHomogeneous
   :members:

nth_linear_euler_eq_nonhomogeneous_variation_of_parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthLinearEulerEqNonhomogeneousVariationOfParameters
   :members:

nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthLinearEulerEqNonhomogeneousUndeterminedCoefficients
   :members:

nth_algebraic
^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthAlgebraic
   :members:

nth_order_reducible
^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::NthOrderReducible
   :members:

separable
^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::Separable
   :members:

almost_linear
^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::AlmostLinear
   :members:

linear_coefficients
^^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::LinearCoefficients
   :members:

separable_reduced
^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::SeparableReduced
   :members:

lie_group
^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::LieGroup
   :members:

2nd_hypergeometric
^^^^^^^^^^^^^^^^^^
.. autoclass:: sympy.solvers.ode.single::SecondHypergeometric
   :members:

1st_power_series
^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::ode_1st_power_series

2nd_power_series_ordinary
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::ode_2nd_power_series_ordinary

2nd_power_series_regular
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::ode_2nd_power_series_regular

Lie heuristics
--------------
These functions are intended for internal use of the Lie Group Solver.
Nonetheless, they contain useful information in their docstrings on the algorithms
implemented for the various heuristics.

abaco1_simple
^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco1_simple

abaco1_product
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco1_product

bivariate
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_bivariate

chi
^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_chi

abaco2_similar
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco2_similar

function_sum
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_function_sum

abaco2_unique_unknown
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco2_unique_unknown

abaco2_unique_general
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_abaco2_unique_general

linear
^^^^^^
.. autofunction:: sympy.solvers.ode.lie_group::lie_heuristic_linear

Rational Riccati Solver
-----------------------
These functions are intended for internal use to solve
a first order Riccati differential equation with atleast
one rational particular solution.

riccati_normal
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::riccati_normal

riccati_inverse_normal
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::riccati_inverse_normal

riccati_reduced
^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::riccati_reduced

construct_c
^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::construct_c

construct_d
^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::construct_d

rational_laurent_series
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::rational_laurent_series

compute_m_ybar
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::compute_m_ybar

solve_aux_eq
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::solve_aux_eq

remove_redundant_sols
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::remove_redundant_sols

get_gen_sol_from_part_sol
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::get_gen_sol_from_part_sol

solve_riccati
^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.riccati::solve_riccati

System of ODEs
--------------
These functions are intended for internal use by
:py:meth:`~sympy.solvers.ode.dsolve` for system of differential equations.

Linear, 2 equations, Order 1, Type 6
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_linear_2eq_order1_type6

Linear, 2 equations, Order 1, Type 7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_linear_2eq_order1_type7

Linear ODE to matrix
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::linear_ode_to_matrix

Canonical Equations Converter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::canonical_odes

LinODESolve Systems Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::linodesolve_type

Matrix Exponential Jordan Form
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::matrix_exp_jordan_form

Matrix Exponential
^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::matrix_exp

Linear, n equations, Order 1 Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.systems::linodesolve

Nonlinear, 2 equations, Order 1, Type 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type1

Nonlinear, 2 equations, Order 1, Type 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type2

Nonlinear, 2 equations, Order 1, Type 3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type3

Nonlinear, 2 equations, Order 1, Type 4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type4

Nonlinear, 2 equations, Order 1, Type 5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_2eq_order1_type5

Nonlinear, 3 equations, Order 1, Type 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type1

Nonlinear, 3 equations, Order 1, Type 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type2

Nonlinear, 3 equations, Order 1, Type 3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type3

Nonlinear, 3 equations, Order 1, Type 4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type4

Nonlinear, 3 equations, Order 1, Type 5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode::_nonlinear_3eq_order1_type5

Information on the ode module
-----------------------------

.. automodule:: sympy.solvers.ode.ode

Internal functions
^^^^^^^^^^^^^^^^^^

These functions are not intended for end-user use.

.. autofunction:: sympy.solvers.ode.ode::_handle_Integral

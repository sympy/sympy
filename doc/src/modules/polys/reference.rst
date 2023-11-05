.. _polys-reference:

=========================================
Polynomials Manipulation Module Reference
=========================================

.. automodule:: sympy.polys

See :ref:`polys-docs` for an index of documentation for the polys module and
:ref:`polys-basics` for an introductory explanation.

Basic polynomial manipulation functions
=======================================

.. currentmodule:: sympy.polys.polytools

.. autofunction:: poly
.. autofunction:: poly_from_expr
.. autofunction:: parallel_poly_from_expr
.. autofunction:: degree
.. autofunction:: degree_list
.. autofunction:: LC
.. autofunction:: LM
.. autofunction:: LT
.. autofunction:: pdiv
.. autofunction:: prem
.. autofunction:: pquo
.. autofunction:: pexquo
.. autofunction:: div
.. autofunction:: rem
.. autofunction:: quo
.. autofunction:: exquo
.. autofunction:: half_gcdex
.. autofunction:: gcdex
.. autofunction:: invert
.. autofunction:: subresultants
.. autofunction:: resultant
.. autofunction:: discriminant
.. autofunction:: terms_gcd
.. autofunction:: cofactors
.. autofunction:: gcd
.. autofunction:: gcd_list
.. autofunction:: lcm
.. autofunction:: lcm_list
.. autofunction:: trunc
.. autofunction:: monic
.. autofunction:: content
.. autofunction:: primitive
.. autofunction:: compose
.. autofunction:: decompose
.. autofunction:: sturm
.. autofunction:: gff_list
.. autofunction:: gff
.. autofunction:: sqf_norm
.. autofunction:: sqf_part
.. autofunction:: sqf_list
.. autofunction:: sqf
.. autofunction:: factor_list
.. autofunction:: factor
.. autofunction:: intervals
.. autofunction:: refine_root
.. autofunction:: count_roots
.. autofunction:: all_roots
.. autofunction:: real_roots
.. autofunction:: nroots
.. autofunction:: ground_roots
.. autofunction:: nth_power_roots_poly
.. autofunction:: cancel
.. autofunction:: reduced
.. autofunction:: groebner
.. autofunction:: is_zero_dimensional

.. autoclass:: Poly
   :members:

.. autoclass:: PurePoly
   :members:

.. autoclass:: GroebnerBasis
   :members:

Extra polynomial manipulation functions
=======================================

.. currentmodule:: sympy.polys.polyfuncs

.. autofunction:: symmetrize
.. autofunction:: horner
.. autofunction:: interpolate
.. autofunction:: viete

Domain constructors
===================

.. currentmodule:: sympy.polys.constructor

.. autofunction:: construct_domain

Monomials encoded as tuples
===========================

.. currentmodule:: sympy.polys.monomials

.. autoclass:: Monomial
   :members:
.. autofunction:: itermonomials
.. autofunction:: monomial_count

Orderings of monomials
======================

.. currentmodule:: sympy.polys.orderings

.. autoclass:: MonomialOrder
   :members:
.. autoclass:: LexOrder
   :members:
.. autoclass:: GradedLexOrder
   :members:
.. autoclass:: ReversedGradedLexOrder
   :members:

Formal manipulation of roots of polynomials
===========================================

.. currentmodule:: sympy.polys.rootoftools

.. autofunction:: rootof
.. autoclass:: RootOf
   :members:
.. autoclass:: ComplexRootOf
   :members:
   :private-members:
.. autoclass:: RootSum
   :members:

Symbolic root-finding algorithms
================================

.. currentmodule:: sympy.polys.polyroots

.. autofunction:: roots

Special polynomials
===================

.. currentmodule:: sympy.polys.specialpolys

.. autofunction:: swinnerton_dyer_poly
.. autofunction:: interpolating_poly
.. autofunction:: cyclotomic_poly
.. autofunction:: symmetric_poly
.. autofunction:: random_poly

Orthogonal polynomials
======================

.. currentmodule:: sympy.polys.orthopolys

.. autofunction:: chebyshevt_poly
.. autofunction:: chebyshevu_poly
.. autofunction:: gegenbauer_poly
.. autofunction:: hermite_poly
.. autofunction:: hermite_prob_poly
.. autofunction:: jacobi_poly
.. autofunction:: legendre_poly
.. autofunction:: laguerre_poly
.. autofunction:: spherical_bessel_fn

Appell sequences
================

.. currentmodule:: sympy.polys.appellseqs

.. autofunction:: bernoulli_poly
.. autofunction:: bernoulli_c_poly
.. autofunction:: genocchi_poly
.. autofunction:: euler_poly
.. autofunction:: andre_poly

Manipulation of rational functions
==================================

.. currentmodule:: sympy.polys.rationaltools

.. autofunction:: together

Partial fraction decomposition
==============================

.. currentmodule:: sympy.polys.partfrac

.. autofunction:: apart
.. autofunction:: apart_list
.. autofunction:: assemble_partfrac_list

Dispersion of Polynomials
=========================

.. currentmodule:: sympy.polys.dispersion

.. autofunction:: dispersionset
.. autofunction:: dispersion

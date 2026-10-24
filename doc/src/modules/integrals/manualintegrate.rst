.. _manualintegrate-rules:

=========================
Rules of manualintegrate
=========================

:func:`~sympy.integrals.manualintegrate.integral_steps` returns a tree of
rules. Each rule either integrates its integrand directly or rewrites the
integral into simpler integrals, which are solved by the rules of its
substeps. Every rule below is described by the integral it applies to, the
result it produces and the conditions under which it applies.

.. autoclass:: sympy.integrals.manualintegrate.Rule

Basic rules
===========

.. autoclass:: sympy.integrals.manualintegrate.ConstantRule
.. autoclass:: sympy.integrals.manualintegrate.ConstantTimesRule
.. autoclass:: sympy.integrals.manualintegrate.AddRule
.. autoclass:: sympy.integrals.manualintegrate.PowerRule
.. autoclass:: sympy.integrals.manualintegrate.NestedPowRule
.. autoclass:: sympy.integrals.manualintegrate.ReciprocalRule
.. autoclass:: sympy.integrals.manualintegrate.ExpRule

Substitution and integration by parts
=====================================

.. autoclass:: sympy.integrals.manualintegrate.URule
.. autoclass:: sympy.integrals.manualintegrate.ReparameterizationRule
.. autoclass:: sympy.integrals.manualintegrate.PartsRule
.. autoclass:: sympy.integrals.manualintegrate.CyclicPartsRule

Trigonometric and hyperbolic functions
======================================

.. autoclass:: sympy.integrals.manualintegrate.SinRule
.. autoclass:: sympy.integrals.manualintegrate.CosRule
.. autoclass:: sympy.integrals.manualintegrate.SinhRule
.. autoclass:: sympy.integrals.manualintegrate.CoshRule

Inverse trigonometric functions and square roots of quadratics
==============================================================

.. autoclass:: sympy.integrals.manualintegrate.ArcsinRule
.. autoclass:: sympy.integrals.manualintegrate.ArcsinhRule
.. autoclass:: sympy.integrals.manualintegrate.ArctanRule
.. autoclass:: sympy.integrals.manualintegrate.CompleteSquareRule
.. autoclass:: sympy.integrals.manualintegrate.ReciprocalSqrtQuadraticRule
.. autoclass:: sympy.integrals.manualintegrate.SqrtQuadraticDenomRule
.. autoclass:: sympy.integrals.manualintegrate.SqrtQuadraticRule

Rational functions
==================

.. autoclass:: sympy.integrals.manualintegrate.RatintRule

Rewriting and combination of results
====================================

.. autoclass:: sympy.integrals.manualintegrate.RewriteRule
.. autoclass:: sympy.integrals.manualintegrate.PiecewiseRule
.. autoclass:: sympy.integrals.manualintegrate.AlternativeRule
.. autoclass:: sympy.integrals.manualintegrate.DerivativeRule
.. autoclass:: sympy.integrals.manualintegrate.DontKnowRule
.. autoclass:: sympy.integrals.manualintegrate.PendingRule

Distributions
=============

.. autoclass:: sympy.integrals.manualintegrate.HeavisideRule
.. autoclass:: sympy.integrals.manualintegrate.DiracDeltaRule

Orthogonal polynomials
======================

.. autoclass:: sympy.integrals.manualintegrate.JacobiRule
.. autoclass:: sympy.integrals.manualintegrate.GegenbauerRule
.. autoclass:: sympy.integrals.manualintegrate.ChebyshevTRule
.. autoclass:: sympy.integrals.manualintegrate.ChebyshevURule
.. autoclass:: sympy.integrals.manualintegrate.LegendreRule
.. autoclass:: sympy.integrals.manualintegrate.HermiteRule
.. autoclass:: sympy.integrals.manualintegrate.LaguerreRule
.. autoclass:: sympy.integrals.manualintegrate.AssocLaguerreRule

Special functions
=================

.. autoclass:: sympy.integrals.manualintegrate.EiRule
.. autoclass:: sympy.integrals.manualintegrate.CiRule
.. autoclass:: sympy.integrals.manualintegrate.SiRule
.. autoclass:: sympy.integrals.manualintegrate.ChiRule
.. autoclass:: sympy.integrals.manualintegrate.ShiRule
.. autoclass:: sympy.integrals.manualintegrate.LiRule
.. autoclass:: sympy.integrals.manualintegrate.ErfRule
.. autoclass:: sympy.integrals.manualintegrate.OwensTRule
.. autoclass:: sympy.integrals.manualintegrate.FresnelCRule
.. autoclass:: sympy.integrals.manualintegrate.FresnelSRule
.. autoclass:: sympy.integrals.manualintegrate.PolylogRule
.. autoclass:: sympy.integrals.manualintegrate.UpperGammaRule
.. autoclass:: sympy.integrals.manualintegrate.EllipticFRule
.. autoclass:: sympy.integrals.manualintegrate.EllipticERule

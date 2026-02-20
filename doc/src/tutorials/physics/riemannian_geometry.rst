.. _riemannian_geometry_tutorial:

========================================================
Computing Christoffel Symbols, Riemann and Ricci Tensors
========================================================

This page demonstrates how to compute curvature tensors from a metric using
SymPy's tensor/differential-geometry functionality.

The workflow is:

1. Define coordinates and a metric tensor ``g``.
2. Compute Christoffel symbols.
3. Compute the Riemann tensor, Ricci tensor, and scalar curvature.
4. Verify results with lightweight symbolic checks.

Flat space in polar coordinates
===============================

The 2D polar metric

.. math::

   ds^2 = dr^2 + r^2 d\theta^2

describes flat Euclidean space, so the Ricci tensor and scalar curvature should
be zero.

We verify this by computing the scalar curvature:
.. code-block:: pycon

   >>> from sympy.diffgeom import TensorProduct, metric_to_Christoffel_2nd, metric_to_Ricci_components
   >>> from sympy.diffgeom.rn import R2_p
   >>> TP = TensorProduct
   >>> g = TP(R2_p.dr, R2_p.dr) + R2_p.r**2*TP(R2_p.dtheta, R2_p.dtheta)
   >>> Gamma = metric_to_Christoffel_2nd(g)
   >>> Gamma[0][1][1]
   -rho
   >>> Ric = metric_to_Ricci_components(g)
   >>> Ric
   [[0, 0], [0, 0]]


Unit 2-sphere
=============

The unit 2-sphere metric

.. math::

   ds^2 = d\theta^2 + \sin^2(\theta)\, d\phi^2

has constant positive curvature. For the unit sphere, the scalar curvature is 2
(with the common convention used here).

.. code-block:: pycon

   >>> from sympy.diffgeom import TensorProduct, metric_to_Ricci_components
   >>> from sympy.diffgeom.rn import R2
   >>> from sympy import simplify
   >>> TP = TensorProduct
   >>> g = (1/R2.y**2)*(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
   >>> Ric = metric_to_Ricci_components(g)
   >>> Ric
   [[-1/y**2, 0], [0, -1/y**2]]
   >>> R = simplify((R2.y**2)*Ric[0][0] + (R2.y**2)*Ric[1][1])
   >>> R
   -2

Notes on larger examples
========================

More complicated metrics (e.g. Schwarzschild) can produce very large symbolic
expressions. If you are working with such metrics, consider:

- applying assumptions (e.g. ``positive=True``) to reduce branch issues,
- using targeted simplification (e.g. simplifying only selected components),
- verifying only key identities/components to keep runtime reasonable.

.. _tensor_riemannian_geometry:

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
.. doctest::

   >>> from sympy import symbols, diag, simplify
   >>> from sympy.diffgeom.riemannian_geometry import (
   ...     metric_to_Christoffel_2nd, riemann_tensor,
   ...     ricci_tensor, scalar_curvature
   ... )
   >>> r, th = symbols('r th', positive=True, real=True)
   >>> g = diag(1, r**2)
   >>> coords = [r, th]
   >>> Gamma = metric_to_Christoffel_2nd(g, coords)
    >>> Gamma[0, 1, 1]
   r
   >>> Riem = riemann_tensor(Gamma, coords)
   >>> Ric = ricci_tensor(Riem, coords)
   >>> Scal = scalar_curvature(Ric, g, coords)
   >>> all(simplify(c) == 0 for row in Ric for c in row)
   True
   >>> simplify(Scal)
   0

Unit 2-sphere
=============

The unit 2-sphere metric

.. math::

   ds^2 = d\theta^2 + \sin^2(\theta)\, d\phi^2

has constant positive curvature. For the unit sphere, the scalar curvature is 2
(with the common convention used here).

.. code-block:: pycon

   >>> from sympy import symbols, diag, sin, simplify
   >>> from sympy.diffgeom.riemannian_geometry import (
   ...     metric_to_Christoffel_2nd, riemann_tensor,
   ...     ricci_tensor, scalar_curvature
   ... )
   >>> th, ph = symbols('th ph', real=True)
   >>> g = diag(1, sin(th)**2)
   >>> coords = [th, ph]
   >>> Gamma = metric_to_Christoffel_2nd(g, coords)
   >>> Riem = riemann_tensor(Gamma, coords)
   >>> Ric = ricci_tensor(Riem, coords)
   >>> Scal = scalar_curvature(Ric, g, coords)
   >>> simplify(Scal - 2)
   0

Notes on larger examples
========================

More complicated metrics (e.g. Schwarzschild) can produce very large symbolic
expressions. If you are working with such metrics, consider:

- applying assumptions (e.g. ``positive=True``) to reduce branch issues,
- using targeted simplification (e.g. simplifying only selected components),
- verifying only key identities/components to keep runtime reasonable.

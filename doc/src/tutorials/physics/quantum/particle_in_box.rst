.. -*- coding: utf-8 -*-
.. _particle_in_box_tutorial:

=====================================
Particle in a Box (Infinite Well)
=====================================

Problem
=======

We consider a quantum particle confined to a one-dimensional infinite potential
well of length ``L`` (often called the "particle in a box" model). The potential
is zero inside the box (``0 < x < L``) and infinite outside, which forces the
wavefunction to vanish at the boundaries.

Input variables
===============

- ``x`` – position along the box (0 to L).
- ``L`` – the length of the box (positive real constant).
- ``n`` – quantum number (positive integer, here we will show the first three).

Solution approach
=================

The time-independent Schrödinger equation in this potential reduces to a simple
boundary value problem with sinusoidal solutions. The normalized eigenfunctions
are:

.. math::

   \psi_n(x) = \sqrt{\frac{2}{L}} \sin\!\left(\frac{n \pi x}{L}\right), \quad n = 1, 2, 3, \dots

We will define these eigenfunctions symbolically using SymPy, convert them to
numerical functions with ``lambdify``, and finally plot them with Matplotlib.

Energy eigenvalues
==================

The corresponding energy spectrum is discrete and scales as :math:`n^2` and
:math:`L^{-2}`:

.. math::

   E_n \;=\; \frac{n^2 \pi^2 \hbar^2}{2 m L^2}\,, \qquad n=1,2,3,\dots

Larger boxes (larger :math:`L`) lower all energies; higher quantum numbers grow
quadratically.

Symbolic eigenfunctions
=======================

Below we define the first three eigenfunctions symbolically with SymPy and then
create a callable function ``f`` using ``lambdify`` so they can be evaluated
numerically for plotting.

.. plot::
   :include-source: True
   :context: reset

   import numpy as np
   import matplotlib.pyplot as plt
   from sympy import symbols, sin, pi, sqrt, lambdify
   from sympy.abc import x

   L = symbols('L', positive=True)
   psi1 = sqrt(2/L) * sin(pi*x/L)
   psi2 = sqrt(2/L) * sin(2*pi*x/L)
   psi3 = sqrt(2/L) * sin(3*pi*x/L)
   print(psi1)

   f = lambdify((x, L), (psi1, psi2, psi3))

Wavefunction plots (SymPy variables → NumPy)
============================================

.. plot::
   :include-source: True
   :context:

   L_val = 1.0
   xs = np.linspace(0.0, L_val, 400)

   plt.figure(figsize=(6,3.5))
   plt.plot(xs, f(xs, L_val)[0], label="psi_1(x)")
   plt.plot(xs, f(xs, L_val)[1], label="psi_2(x)")
   plt.plot(xs, f(xs, L_val)[2], label="psi_3(x)")
   plt.xlabel("x")
   plt.ylabel("psi_n(x)")
   plt.title("Particle in a box: eigenfunctions (driven by SymPy variable L)")
   plt.legend()
   plt.tight_layout()
   plt.show()

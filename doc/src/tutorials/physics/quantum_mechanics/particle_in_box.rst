.. -*- coding: utf-8 -*-
.. _particle_in_box_tutorial:

=====================================
Particle in a Box (Infinite Well)
=====================================

We define the normalized eigenfunctions using SymPy symbols first and only then
evaluate them for plotting.

Symbolic eigenfunctions
=======================

.. code-block:: python

   from sympy import symbols, sin, pi, sqrt
   from sympy.abc import x

   L = symbols('L', positive=True)
   psi1 = sqrt(2/L) * sin(pi*x/L)
   psi2 = sqrt(2/L) * sin(2*pi*x/L)
   psi3 = sqrt(2/L) * sin(3*pi*x/L)
   print(psi1)

Wavefunction plots (SymPy variables → NumPy)
============================================

.. plot::
   :include-source: True

   import numpy as np
   import matplotlib.pyplot as plt
   from sympy import symbols, sin, pi, sqrt, lambdify
   from sympy.abc import x

   L = symbols('L', positive=True)
   psi1 = sqrt(2/L) * sin(pi*x/L)
   psi2 = sqrt(2/L) * sin(2*pi*x/L)
   psi3 = sqrt(2/L) * sin(3*pi*x/L)

   f1 = lambdify((x, L), psi1, modules="numpy")
   f2 = lambdify((x, L), psi2, modules="numpy")
   f3 = lambdify((x, L), psi3, modules="numpy")

   L_val = 1.0
   xs = np.linspace(0.0, L_val, 400)

   plt.figure(figsize=(6,3.5))
   plt.plot(xs, f1(xs, L_val), label="psi_1(x)")
   plt.plot(xs, f2(xs, L_val), label="psi_2(x)")
   plt.plot(xs, f3(xs, L_val), label="psi_3(x)")
   plt.xlabel("x")
   plt.ylabel("psi_n(x)")
   plt.title("Particle in a box: eigenfunctions (driven by SymPy variable L)")
   plt.legend()
   plt.tight_layout()
   plt.show()

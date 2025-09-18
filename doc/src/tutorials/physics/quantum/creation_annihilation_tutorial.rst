.. -*- coding: utf-8 -*-
.. _creation_annihilation_tutorial:

======================================
Creation and Annihilation Operators
======================================

We use :mod:`~sympy` symbols and quantum operators; the visualization uses a SymPy
expression for the energy ladder.

Afterwards we evaluate it via :func:`sympy.utilities.lambdify.lambdify`.

Symbolic setup
==============

.. code-block:: python

   from sympy.physics.quantum import Dagger, Operator, Commutator
   from sympy import symbols

   a = Operator('a')
   adag = Dagger(a)
   print("Commutator [a, a^†] =", Commutator(a, adag))

   n = symbols('n', integer=True, nonnegative=True)
   hw = symbols('hbar_omega', positive=True)  # parameter
   E_n = hw*(n + 1/2)
   print("E_n =", E_n)

Energy ladder plot (SymPy variables → NumPy)
============================================

.. plot::
   :include-source: True

   import numpy as np
   import matplotlib.pyplot as plt
   from sympy import symbols, lambdify

   # SymPy variables & expression
   n = symbols('n', integer=True, nonnegative=True)
   hw = symbols('hbar_omega', positive=True)
   E_n = hw*(n + 1/2)

   # Turn into a numeric callable E(n, hw)
   E = lambdify((n, hw), E_n, modules="numpy")

   ns = np.arange(0, 7)
   hw_val = 1.0
   En = E(ns, hw_val)

   plt.figure(figsize=(6,3.5))
   markerline, stemlines, baseline = plt.stem(ns, En)
   plt.setp(markerline, markersize=6)
   plt.xlabel("n")
   plt.ylabel("E_n")
   plt.title("Harmonic-oscillator energy ladder (driven by SymPy variable hbar_omega)")
   plt.tight_layout()
   plt.show()

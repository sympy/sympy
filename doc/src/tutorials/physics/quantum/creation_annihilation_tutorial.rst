.. -*- coding: utf-8 -*-
.. _creation_annihilation_tutorial:

======================================
Creation and Annihilation Operators
======================================

Creation (``a†``) and annihilation (``a``) operators are fundamental tools in
quantum mechanics, especially in the study of the harmonic oscillator and
quantum field theory. They provide a convenient algebraic way to move between
quantum states:

- The annihilation operator ``a`` lowers the state (removes one quantum of
  energy).
- The creation operator ``a†`` raises the state (adds one quantum of energy).

These operators obey the commutation relation :math:`[a, a^\dagger] = 1`, which
encodes the quantization of the system. With them we can systematically build
the ladder of energy eigenstates and compute physical quantities such as the
spectrum of the harmonic oscillator.

We will demonstrate this algebra in SymPy, symbolically construct the energy
levels :math:`E_n`, and then plot the discrete spectrum.

Symbolic setup
==============

We introduce the following variables and operators:

- ``a`` – the annihilation operator.
- ``a†`` – the creation operator, defined as the Hermitian conjugate of ``a``.
- ``n`` – the quantum number, a nonnegative integer representing the energy level.
- ``hbar_omega`` – shorthand for :math:`\hbar \omega`, the characteristic energy scale.
- ``E_n`` – the symbolic expression for the nth energy level,
  :math:`E_n = \hbar \omega \left(n + \tfrac{1}{2}\right)`.

.. plot::
   :include-source: True
   :context: reset

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

We now convert the symbolic expression for :math:`E_n` into a numerical function
using ``lambdify`` and evaluate it for several values of ``n`` to plot the
discrete energy ladder of the harmonic oscillator.

.. plot::
   :include-source: True
   :context:

   import numpy as np
   import matplotlib.pyplot as plt
   from sympy import lambdify

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

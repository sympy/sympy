.. -*- coding: utf-8 -*-
.. _pauli_commutators_tutorial:

===============================
Pauli Matrices and Commutators
===============================

This tutorial defines Pauli matrices with :mod:`sympy` and computes commutators.
The visualization is **driven by SymPy variables** and only converted for plotting.

Symbolic setup
==============

.. code-block:: python

   from sympy import Matrix, I

   sigma_x = Matrix([[0, 1], [1, 0]])
   sigma_y = Matrix([[0, -I], [I, 0]])
   sigma_z = Matrix([[1, 0], [0, -1]])

   comm_xy = sigma_x * sigma_y - sigma_y * sigma_x
   comm_yz = sigma_y * sigma_z - sigma_z * sigma_y
   comm_zx = sigma_z * sigma_x - sigma_x * sigma_z

   print("[sx, sy] =", comm_xy)
   print("[sy, sz] =", comm_yz)
   print("[sz, sx] =", comm_zx)

Visualization (SymPy → NumPy via ``lambdify``)
===============================================

.. plot::
   :include-source: True

   import numpy as np
   import matplotlib.pyplot as plt
   from sympy import Matrix, I, lambdify, symbols

   # SymPy variables to scale/shift (just to demonstrate true SymPy control)
   s = symbols('s', real=True)
   offset = symbols('o', real=True)

   # SymPy matrices
   sigma_x = Matrix([[0, 1], [1, 0]])
   sigma_y = Matrix([[0, -I], [I, 0]])
   sigma_z = Matrix([[1, 0], [0, -1]])

   # Build SymPy expressions that depend on variables
   Mx = s * sigma_x + offset*Matrix([[1,0],[0,1]])
   My = s * sigma_y.as_real_imag()[0]  # take real part symbolically
   Mz = s * sigma_z

   # Lambdify to get numeric arrays once we pick s, offset
   Mx_fn = lambdify((s, offset), Mx, modules="numpy")
   My_fn = lambdify((s,), My, modules="numpy")
   Mz_fn = lambdify((s,), Mz, modules="numpy")

   # Choose numeric values for the SymPy variables
   sval, oval = 1.0, 0.0
   Mx_np = np.array(Mx_fn(sval, oval), dtype=float)
   My_np = np.array(My_fn(sval), dtype=float)
   Mz_np = np.array(Mz_fn(sval), dtype=float)

   mats = [Mx_np, My_np, Mz_np]
   titles = ["sigma_x (scaled by s)", "Re(sigma_y) (scaled by s)", "sigma_z (scaled by s)"]

   fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.5))
   for ax, M, title in zip(axs, mats, titles):
       im = ax.imshow(M, vmin=-1, vmax=1)
       ax.set_title(title)
       ax.set_xticks([]); ax.set_yticks([])
   #fig.colorbar(im, ax=axs.ravel().tolist(), orientation="vertical", shrink=0.7, label="value",pad=0.05 )
   cax = fig.add_axes([0.95, 0.05, 0.01, 0.7])
   cbar=fig.colorbar(im, cax=cax, orientation="vertical", label="value")
   cbar.ax.tick_params(labelsize=8)
   fig.suptitle("Pauli matrices (driven by SymPy variables s, o)")
   fig.tight_layout()
   plt.show()

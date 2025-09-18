.. -*- coding: utf-8 -*-
.. _pauli_commutators_tutorial:

===============================
Pauli Matrices and Commutators
===============================

Problem
-------

Explore the algebra of the **Pauli matrices**, compute their commutators symbolically using SymPy, and visualize them as scaled numeric matrices controlled by symbolic parameters.

Input variables
---------------

- :math:`\sigma_x, \sigma_y, \sigma_z` — the standard 2×2 Pauli matrices.  
- :math:`s` — a real-valued scaling parameter applied to the matrices in the visualization.  
- :math:`o` — an offset parameter applied to :math:`\sigma_x` as a scalar multiple of the identity.  

Example problem
---------------

**Question.** Verify that the Pauli matrices satisfy the SU(2) commutation relations  

.. math::

   [\sigma_i, \sigma_j] = 2i \,\epsilon_{ijk}\,\sigma_k,

where :math:`\epsilon_{ijk}` is the Levi-Civita symbol. Then, visualize the effect of applying a scaling parameter :math:`s` and an offset :math:`o` to these matrices.

**Approach.**  
1. Define the Pauli matrices symbolically in SymPy.  
2. Compute their commutators and check if they match the SU(2) algebra.  
3. Introduce symbolic variables :math:`s,o`, lambdify the modified matrices, and plot the results as numeric arrays.

Symbolic setup
==============

We begin by defining the three Pauli matrices:

.. math::

   \sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
   \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
   \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.

We then compute their commutators:

.. math::

   [\sigma_x,\sigma_y]=\sigma_x\sigma_y-\sigma_y\sigma_x,\quad
   [\sigma_y,\sigma_z]=\sigma_y\sigma_z-\sigma_z\sigma_y,\quad
   [\sigma_z,\sigma_x]=\sigma_z\sigma_x-\sigma_x\sigma_z.

These reproduce the well-known relations
:math:`[\sigma_i,\sigma_j]=2i\epsilon_{ijk}\sigma_k`.

.. plot::
   :context: reset
   :include-source: True
   :format: python

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

For visualization, we introduce symbolic variables:

- ``s`` is a real-valued scaling factor applied to each matrix.  
- ``offset`` is a real-valued parameter that shifts :math:`\sigma_x` by adding ``offset`` times the identity matrix.  
- ``Mx``, ``My``, ``Mz`` are symbolic matrix expressions that depend on ``s`` and ``offset``.  
- Each is lambdified into a NumPy function for fast numeric evaluation.

After choosing numeric values (``sval=1.0``, ``oval=0.0``), we evaluate these matrices and display them. Each subplot corresponds to one Pauli matrix (or its real part, in the case of :math:`\sigma_y`).

.. plot::
   :context:
   :include-source: True
   :format: python

   import numpy as np
   import matplotlib.pyplot as plt
   from sympy import Matrix, I, lambdify, symbols

   s = symbols('s', real=True)
   offset = symbols('o', real=True)

   sigma_x = Matrix([[0, 1], [1, 0]])
   sigma_y = Matrix([[0, -I], [I, 0]])
   sigma_z = Matrix([[1, 0], [0, -1]])

   Mx = s * sigma_x + offset*Matrix([[1,0],[0,1]])
   My = s * sigma_y.as_real_imag()[0]
   Mz = s * sigma_z

   fns = lambdify((s, offset), [Mx, My, Mz])

   sval, oval = 1.0, 0.0
   Mx_np, My_np, Mz_np = [np.array(M, dtype=float) for M in fns(sval, oval)]

   mats = [Mx_np, My_np, Mz_np]
   titles = ["sigma_x (scaled by s, shifted by o)",
             "Re(sigma_y) (scaled by s)",
             "sigma_z (scaled by s)"]

   plt.close('all')
   fig, axs = plt.subplots(1, 3, figsize=(8.0, 2.5))

   for ax, M, title in zip(axs, mats, titles):
       im = ax.imshow(M, vmin=-1, vmax=1)
       ax.set_title(title, fontsize=8)
       ax.set_xticks([]); ax.set_yticks([])


   cax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  
   cbar = fig.colorbar(im, cax=cax, orientation="vertical", label="value")
   cbar.ax.tick_params(labelsize=8)


   fig.subplots_adjust(left=0.05, right=0.85)

   fig.suptitle("Pauli matrices (driven by SymPy variables s, o)")
   fig.tight_layout(rect=[0, 0, 0.85, 1])  # reserve space for colorbar
   plt.show()

Summary
--------------

- **Algebra.** The commutators confirm the SU(2) structure: :math:`[\sigma_x,\sigma_y]=2i\sigma_z`, and cyclic permutations.  
- **Parameters.** The scaling parameter ``s`` magnifies the matrices, while ``offset`` shifts :math:`\sigma_x` by the identity, demonstrating symbolic control.  
- **Visualization.** Displaying the matrices as color maps emphasizes their entries’ structure and how symbolic parameters affect them.  
- **Symbolic-to-numeric pipeline.** Using :func:`sympy.utilities.lambdify` to convert symbolic matrices into NumPy arrays cleanly separates the algebraic derivation from the plotting stage.

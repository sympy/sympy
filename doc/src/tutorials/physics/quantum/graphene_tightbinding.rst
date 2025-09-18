.. -*- coding: utf-8 -*-
.. _graphene_tightbinding:

==========================================
Graphene Band Structure with Tight Binding
==========================================

Problem
-------

Compute the graphene :math:`\pi`-band dispersion from the nearest-neighbor tight-binding model and visualize (i) the band structure along a high-symmetry path and (ii) the density of states (DOS). Everything is first defined **symbolically** in SymPy, then evaluated numerically via :func:`sympy.utilities.lambdify.lambdify`.

Input variables
---------------

- :math:`k_x, k_y` — components of the crystal momentum vector :math:`\mathbf{k}`.
- :math:`a` — lattice constant setting real-space neighbor distances.
- :math:`t` — nearest-neighbor hopping amplitude (energy scale in eV).

Derived symbolic objects
------------------------

We introduce the **nearest-neighbor displacement vectors** :math:`\mathbf{d}_1,\mathbf{d}_2,\mathbf{d}_3`, which connect A–B sublattices in the honeycomb lattice. Using :math:`k_x,k_y,a`, we define the **structure factor**
:math:`f(\mathbf{k})=\sum_{j=1}^3 e^{\,i\mathbf{k}\cdot\mathbf{d}_j}` and its **magnitude**
:math:`|f(\mathbf{k})|`. The tight-binding eigenvalues are
:math:`E_\pm(\mathbf{k})=\pm\,t\,|f(\mathbf{k})|`. These symbolic quantities will be lambdified into fast NumPy callables for plotting.

.. math::
   \mathbf{d}_1=(0,\,a),\quad
   \mathbf{d}_2=\Bigl(\tfrac{\sqrt{3}}{2}a,\,-\tfrac{1}{2}a\Bigr),\quad
   \mathbf{d}_3=\Bigl(-\tfrac{\sqrt{3}}{2}a,\,-\tfrac{1}{2}a\Bigr),\qquad
   f(\mathbf{k})=\sum_{j=1}^3 e^{i\mathbf{k}\cdot\mathbf{d}_j},\quad
   E_\pm(\mathbf{k})=\pm\,t\,|f(\mathbf{k})|.

Symbolic setup
------------------------------------

This block defines and explains each symbolic variable used later:

- ``kx, ky`` are the symbolic components of :math:`\mathbf{k}`.
- ``a`` is the symbolic lattice constant; it parametrizes the neighbor vectors and rescales the Brillouin-zone coordinates.
- ``t`` is the symbolic hopping that multiplies :math:`|f|` to produce the band energies.
- ``d1, d2, d3`` are the symbolic nearest-neighbor vectors.
- ``f`` is the structure factor; ``fabs`` is its magnitude :math:`|f|`.
- ``fabs_np`` is a NumPy-callable function of ``(kx, ky, a)`` produced by :func:`sympy.utilities.lambdify`.
- ``a_val, t_val`` are **numeric** values assigned once and reused throughout.

.. plot::
   :context: reset
   :include-source: True
   :nofigs:
   :format: python

   from sympy import symbols, Matrix, I, exp, Abs, sqrt, pi, lambdify

   kx, ky, a, t = symbols('kx ky a t', positive=True)
   d1 = Matrix([0, a])
   d2 = Matrix([sqrt(3)*a/2, -a/2])
   d3 = Matrix([-sqrt(3)*a/2, -a/2])

   f = exp(I*(kx*d1[0] + ky*d1[1])) \
     + exp(I*(kx*d2[0] + ky*d2[1])) \
     + exp(I*(kx*d3[0] + ky*d3[1]))
   fabs = Abs(f)

   fabs_np = lambdify((kx, ky, a), fabs, modules=["numpy"])

   a_val = 1.0   # numeric lattice constant
   t_val = 2.7   # numeric hopping (eV)

Band path, sampling variables, and evaluation
-------------------------------------------------------------------

We next define all numeric variables that govern the band-structure path and evaluation:

- ``G``, ``K``, ``M`` are the **high-symmetry points** in the hexagonal Brillouin zone, expressed as numeric 2-vectors in units of :math:`1/a`.
- ``path`` is the ordered list :math:`\Gamma\!\to\!K\!\to\!M\!\to\!\Gamma`.
- ``labels`` are the tick labels corresponding to segment boundaries.
- ``Nk`` is the **number of k-points per segment**; larger values give smoother curves.
- ``kpts`` is the **stacked array of sampled momenta** obtained by linear interpolation on each segment.
- ``magf`` is the **numeric array** of :math:`|f(\mathbf{k}_i)|` values evaluated via ``fabs_np``.
- ``Eplus`` and ``Eminus`` are the **upper and lower band energies** :math:`\pm t\,|f|` along the path.

We close any existing figures and open a **fresh figure** before plotting to avoid state bleed.

.. plot::
   :context:
   :include-source: True
   :format: python

   import numpy as np
   import matplotlib.pyplot as plt

   G = np.array([0.0, 0.0])
   K = np.array([2*np.pi/(3*np.sqrt(3)*a_val), 2*np.pi/(3*a_val)])
   M = np.array([np.pi/(np.sqrt(3)*a_val), np.pi/(3*a_val)])
   path = [G, K, M, G]
   labels = [r"$\Gamma$", "K", "M", r"$\Gamma$"]

   Nk = 200
   kpts = np.vstack([np.linspace(path[i], path[i+1], Nk) for i in range(len(path)-1)])

   magf = fabs_np(kpts[:,0], kpts[:,1], a_val)
   Eplus  = +t_val * magf
   Eminus = -t_val * magf

   plt.close('all')
   plt.figure(figsize=(6,4))
   plt.plot(Eplus,  lw=1.6)
   plt.plot(Eminus, lw=1.6)
   plt.axhline(0.0, linestyle="--", linewidth=1)
   plt.xticks([0, Nk, 2*Nk, 3*Nk], labels)
   plt.ylabel("Energy (eV)")
   plt.title("Graphene Band Structure (driven by SymPy variables a, t)")
   plt.grid(True)
   plt.tight_layout()
   plt.show()

DOS strategy, random-k variables, and histogram
-------------------------------------------------------------------

For the DOS we introduce and discuss the sampling variables:

- ``N`` is the number of **random k-samples** used for Monte Carlo estimation.
- ``kx_s, ky_s`` are arrays of random :math:`k_x,k_y` values uniformly drawn from :math:`[-\pi/a,\pi/a]`. This rectangle covers the first Brillouin zone and is convenient for sampling; importance sampling over the true hexagonal BZ could refine results.
- ``magf`` is again :math:`|f(\mathbf{k})|` evaluated at the sampled points.
- ``E`` is the **combined energy array** containing both bands, i.e., concatenating :math:`+t|f|` and :math:`-t|f|`.
- The **histogram** with ``density=True`` returns a normalized DOS estimate; the visible features include the linear vanishing near the Dirac point and van Hove peaks near :math:`|E|\approx t`.


.. plot::
   :context:
   :include-source: True
   :format: python

   import numpy as np
   import matplotlib.pyplot as plt

   N = 30000
   kx_s = np.random.uniform(-np.pi/a_val, np.pi/a_val, N)
   ky_s = np.random.uniform(-np.pi/a_val, np.pi/a_val, N)

   magf = fabs_np(kx_s, ky_s, a_val)
   E = np.concatenate([+t_val*magf, -t_val*magf])

   plt.close('all')
   plt.figure(figsize=(6,4))
   plt.hist(E, bins=200, density=True, alpha=0.75)
   plt.xlabel("Energy (eV)")
   plt.ylabel("DOS (a.u.)")
   plt.title("Graphene DOS (driven by SymPy variable a)")
   plt.tight_layout()
   plt.show()

Interpretation of variables and results
---------------------------------------

- Role of :math:`a` and :math:`t`.
  The lattice constant :math:`a` sets the real-space scale, rescaling momenta and fixing the positions of the high-symmetry points :math:`K` and :math:`M` in the Brillouin zone.
  The hopping amplitude :math:`t` fixes the overall **energy scale**, stretching or compressing the bands linearly as :math:`E_\pm=\pm\,t\,|f(\mathbf{k})|`.

- Meaning of :math:`f(\mathbf{k})`.
  The complex structure factor :math:`f(\mathbf{k})` represents the **interference** of contributions from the three nearest neighbors. Its zeros identify the **Dirac points**, where conduction and valence bands meet at :math:`E_\pm=0`.

- Sampling sizes (``Nk``, ``N``).
  The variable ``Nk`` controls the number of interpolation points along the high-symmetry path; larger values refine the band-structure curves.
  The variable ``N`` sets the number of random momentum samples used for the DOS; larger values smooth the histogram and make the van Hove singularities more distinct.

- Numerical stability.
  The use of :func:`sympy.utilities.lambdify.lambdify` to create the callable ``fabs_np`` ensures a clean separation between **symbolic definitions** and **numeric evaluation**, preserving exactness in the derivation while enabling efficient computation and plotting.
  
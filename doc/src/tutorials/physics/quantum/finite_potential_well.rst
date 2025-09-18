.. -*- coding: utf-8 -*-
.. _finite_square_well_sympy:

Finite Square Well with SymPy (Analytic → Numeric)
==================================================

Problem
-------

Find bound-state energies and normalized wavefunctions for the symmetric finite
square well with depth :math:`V_0>0` and half-width :math:`a`, where
:math:`V(x)=-V_0` for :math:`|x|\le a` and :math:`V(x)=0` otherwise.

Input variables
---------------

- :math:`m` – particle mass (positive).
- :math:`\hbar` – reduced Planck constant (positive).
- :math:`a` – half-width of the well (positive).
- :math:`V_0` – well depth (positive).
- :math:`\varepsilon` – binding energy with physical energy :math:`E=-\varepsilon` and :math:`0<\varepsilon<V_0`.
- :math:`x` – position along the 1D line (real).

Solution approach
-----------------

1. Derive even/odd transcendental equations and analytic normalization constants
   symbolically in SymPy.
2. Use :func:`sympy.solvers.solvers.nsolve` to find bound-state energies.
3. Plot zero-crossing functions, potential, and analytically
   normalized wavefunctions.

Analytic derivation (all SymPy)
-------------------------------

**Step 1 — Parameterization of bound states.**  
Bound states satisfy :math:`E<0`.
Introduce the *binding energy*
:math:`\varepsilon>0` via :math:`E=-\varepsilon` with :math:`0<\varepsilon<V_0`.
Define inside/outside wave numbers
:math:`k=\sqrt{2m(V_0-\varepsilon)}/\hbar` and
:math:`\alpha=\sqrt{2m\varepsilon}/\hbar`.

**Step 2 — Parity decomposition and boundary conditions.**  
For a symmetric well, eigenfunctions have definite parity.
Imposing continuity at :math:`x=\pm a` and exponential decay outside yields the
bounded transcendental conditions

.. math::
   \begin{aligned}
   k\sin(ka)-\alpha\cos(ka) &= 0 \quad (\text{even}),\\
   k\cos(ka)+\alpha\sin(ka) &= 0 \quad (\text{odd}).
   \end{aligned}

**Step 3 — Piecewise eigenfunctions.**  
Using the parity-appropriate trigonometric form in :math:`|x|\le a` and
exponential tails in :math:`|x|>a`, construct

.. math::
   \psi_{\text{even}}(x)=
   \begin{cases}
   A\cos(kx), & |x|\le a,\\[2pt]
   A\cos(ka)\,e^{-\alpha(|x|-a)}, & |x|>a,
   \end{cases}
   \qquad
   \psi_{\text{odd}}(x)=
   \begin{cases}
   A\sin(kx), & |x|\le a,\\[2pt]
   A\,\operatorname{sgn}(x)\sin(ka)\,e^{-\alpha(|x|-a)}, & |x|>a.
   \end{cases}

**Step 4 — Analytic normalization.**  
By symmetry, integrate over :math:`x\ge 0` and double:

.. math::
   \begin{aligned}
   \|\psi_{\text{even}}\|^2
   &= 2\!\left[\int_0^a \cos^2(kx)\,dx
   + \cos^2(ka)\!\int_a^\infty e^{-2\alpha(x-a)}\,dx\right],\\
   \|\psi_{\text{odd}}\|^2
   &= 2\!\left[\int_0^a \sin^2(kx)\,dx
   + \sin^2(ka)\!\int_a^\infty e^{-2\alpha(x-a)}\,dx\right].
   \end{aligned}

Set :math:`A_{\text{even}}=1/\|\psi_{\text{even}}\|` and
:math:`A_{\text{odd}}=1/\|\psi_{\text{odd}}\|` to obtain closed-form
normalization constants.

**Step 5 — Symbolic-to-numeric interface.**  
Keep the derivation symbolic in SymPy; later substitute numerical parameters,
solve the transcendental equations, and evaluate normalized wavefunctions via
:func:`sympy.utilities.lambdify.lambdify`.

.. plot::
   :context: reset
   :include-source: True
   :nofigs:
   :format: python

   import sympy as sp

   m, hbar, a, V0 = sp.symbols("m hbar a V0", positive=True, real=True)
   eps = sp.symbols("eps", positive=True, real=True)
   E = -eps

   k     = sp.sqrt(2*m*(V0 - eps))/hbar
   alpha = sp.sqrt(2*m*eps)/hbar
   x = sp.symbols("x", real=True)

   even_eq = k*sp.sin(k*a) - alpha*sp.cos(k*a)
   odd_eq  = k*sp.cos(k*a) + alpha*sp.sin(k*a)

   A = sp.symbols("A", positive=True)

   psi_even_unnorm = sp.Piecewise(
       (A*sp.cos(k*x), sp.Abs(x) <= a),
       (A*sp.cos(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   psi_odd_unnorm = sp.Piecewise(
       (A*sp.sin(k*x), sp.Abs(x) <= a),
       (A*sp.sign(x)*sp.sin(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   int_even_inside = sp.integrate(sp.cos(k*x)**2, (x, 0, a))
   int_even_out    = sp.integrate(sp.exp(-2*alpha*(x - a)), (x, a, sp.oo)) * sp.cos(k*a)**2
   norm_sq_even    = 2*(int_even_inside + int_even_out) * A**2
   A_even_exact    = sp.simplify(sp.sqrt(1/norm_sq_even).subs({A:1}))

   int_odd_inside = sp.integrate(sp.sin(k*x)**2, (x, 0, a))
   int_odd_out    = sp.integrate(sp.exp(-2*alpha*(x - a)), (x, a, sp.oo)) * sp.sin(k*a)**2
   norm_sq_odd    = 2*(int_odd_inside + int_odd_out) * A**2
   A_odd_exact    = sp.simplify(sp.sqrt(1/norm_sq_odd).subs({A:1}))

   psi_even = sp.Piecewise(
       (A_even_exact*sp.cos(k*x), sp.Abs(x) <= a),
       (A_even_exact*sp.cos(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   psi_odd = sp.Piecewise(
       (A_odd_exact*sp.sin(k*x), sp.Abs(x) <= a),
       (A_odd_exact*sp.sign(x)*sp.sin(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

Numeric solving (nsolve) + plots
--------------------------------

We now assign numerical values to the physical parameters
(:math:`m=\hbar=a=1, V_0=50`) to compute and visualize bound states.

.. plot::
   :context:
   :include-source: True
   :nofigs:
   :format: python

   import numpy as np
   import matplotlib.pyplot as plt

   m_val, hbar_val, a_val, V0_val = 1.0, 1.0, 1.0, 50.0

   even_eps = sp.lambdify(
       eps, even_eq.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val}), "numpy"
   )
   odd_eps  = sp.lambdify(
       eps, odd_eq.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val}), "numpy"
   )

   def find_roots_in_eps(expr, lo, hi, ntry=200):
       roots = []
       expr_E = expr.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val})
       guesses = np.linspace(lo, hi, ntry)
       for g in guesses:
           try:
               r = float(sp.nsolve(expr_E, eps, g))
               if lo < r < hi and not any(abs(r - q) < 1e-6 for q in roots):
                   roots.append(r)
           except Exception:
               pass
       return sorted(roots)

   eps_even_roots = find_roots_in_eps(even_eq, 1e-6, V0_val-1e-6)
   eps_odd_roots  = find_roots_in_eps(odd_eq,  1e-6, V0_val-1e-6)

   E_even = [-r for r in eps_even_roots]
   E_odd  = [-r for r in eps_odd_roots]

We then plot the zero-crossing functions and mark the roots corresponding to
allowed energies.

.. plot::
   :context:
   :include-source: True
   :format: python

   eps_grid = np.linspace(1e-6, V0_val-1e-6, 3000)
   E_grid   = -eps_grid
   Ye = np.clip(even_eps(eps_grid), -10, 10)
   Yo = np.clip(odd_eps(eps_grid), -10, 10)

   plt.close('all')
   plt.figure(figsize=(7.5, 4))
   plt.axhline(0.0, lw=1)
   plt.plot(E_grid, Ye, label="even condition")
   plt.plot(E_grid, Yo, label="odd condition")
   if E_even: plt.scatter(E_even, [0]*len(E_even), marker='o', label="even roots")
   if E_odd:  plt.scatter(E_odd,  [0]*len(E_odd),  marker='x', label="odd roots")
   plt.xlabel("Energy E")
   plt.ylabel("Zero-crossing function")
   plt.title("Finite Square Well: Transcendental Equations and Roots")
   plt.legend()
   plt.tight_layout()
   plt.show()

Next, we define :math:`V(x)` and reconstruct normalized
:math:`\psi_{\text{even}}(x)` and :math:`\psi_{\text{odd}}(x)`.

.. plot::
   :context:
   :include-source: True
   :nofigs:
   :format: python

   Vx_sym = sp.Piecewise(
       (-V0_val, sp.And(x >= -a_val, x <= a_val)),
       (0.0, True)
   )
   Vx = sp.lambdify(x, Vx_sym, "numpy")
   xs = np.linspace(-2*a_val, 2*a_val, 2000)
   V_vals = Vx(xs)

   def build_psi_at_eps(eps_val):
       k_val = float(sp.sqrt(2*m_val*(V0_val - eps_val))/hbar_val)
       apha  = float(sp.sqrt(2*m_val*eps_val)/hbar_val)
       Aeven = float(sp.N(A_even_exact.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val, eps:eps_val})))
       Aodd  = float(sp.N(A_odd_exact.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val, eps:eps_val})))
       psi_even_num = sp.lambdify(
           x,
           sp.Piecewise(
               (Aeven*sp.cos(k_val*x), sp.Abs(x) <= a_val),
               (Aeven*sp.cos(k_val*a_val)*sp.exp(-apha*(sp.Abs(x) - a_val)), True)
           ),
           "numpy"
       )
       psi_odd_num = sp.lambdify(
           x,
           sp.Piecewise(
               (Aodd*sp.sin(k_val*x), sp.Abs(x) <= a_val),
               (Aodd*sp.sign(x)*sp.sin(k_val*a_val)*sp.exp(-apha*(sp.Abs(x) - a_val)), True)
           ),
           "numpy"
       )
       return psi_even_num, psi_odd_num

Finally, we plot the potential and the lowest even/odd bound states.

.. plot::
   :context:
   :include-source: True
   :format: python

   if eps_even_roots and eps_odd_roots:
       eps0, eps1 = eps_even_roots[0], eps_odd_roots[0]
       psi_e_num, psi_o_num = build_psi_at_eps(eps0)[0], build_psi_at_eps(eps1)[1]
       psi_e_vals = psi_e_num(xs)
       psi_o_vals = psi_o_num(xs)

       plt.close('all')
       plt.figure(figsize=(8, 5))
       plt.plot(xs, V_vals, "k-", lw=2, label="V(x)")
       plt.plot(xs, -eps0 + psi_e_vals, "b", label=fr"even, $E={-eps0:.3f}$")
       plt.plot(xs, -eps1 + psi_o_vals, "r", label=fr"odd,  $E={-eps1:.3f}$")
       plt.axhline(0, color="black", lw=1)
       plt.xlabel("x")
       plt.ylabel(r"Energy / $\psi(x)$")
       plt.title("Square Well: Potential and Normalized Bound States")
       plt.legend()
       plt.tight_layout()
       plt.show()
   else:
       plt.close('all')
       plt.figure(figsize=(8, 4))
       plt.plot(xs, V_vals, "k-", lw=2, label="V(x)")
       plt.axhline(0, color="black", lw=1)
       plt.xlabel("x")
       plt.ylabel("Energy")
       plt.title("Square Well Potential (no bound states found)")
       plt.legend()
       plt.tight_layout()
       plt.show()

Summary
-------

The finite square well supports a **finite** number of bound states determined
by parity-resolved transcendental equations. Compared to the infinite well,
eigenfunctions have **exponential tails** outside the well, and the number of
levels increases with both depth :math:`V_0` and width :math:`a`. Analytic
normalization ensures correctly scaled wavefunctions, while the symbolic-to-
numeric workflow cleanly separates derivation from evaluation and visualization.

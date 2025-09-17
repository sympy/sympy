.. -*- coding: utf-8 -*-
.. _finite_square_well_sympy:

Finite Square Well with SymPy (Analytic → Numeric)
==================================================

This tutorial shows how to use **SymPy** to:

1. Derive the **finite square well** equations **analytically**:
   - Bound-state conditions for even/odd parity.
   - **Exact normalization constants** for the wavefunctions (symbolic integrals).
2. Solve **numerically** for bound-state energies using :func:'sympy.nsolve'.
3. Plot the **zero-crossing functions**, the **square-well potential**, and the
   **analytically normalized** wavefunctions.

We consider a symmetric well:
:math:`V(x) = -V_0` for :math:`\lvert x\rvert\le a` and :math:`V(x)=0` otherwise, with :math:`V_0>0`.

Analytic derivation (all SymPy)
-------------------------------

.. plot::
   :context: reset
   :include-source: True
   :format: python

   import sympy as sp

   # --- Symbols (all positive by physics) ---
   m, hbar, a, V0 = sp.symbols("m hbar a V0", positive=True, real=True)
   eps = sp.symbols("eps", positive=True, real=True)  # binding energy ε > 0
   E = -eps                                           # physical energy E < 0

   # Wave numbers (REAL for 0 < eps < V0)
   k     = sp.sqrt(2*m*(V0 - eps))/hbar
   alpha = sp.sqrt(2*m*eps)/hbar

   x = sp.symbols("x", real=True)

   # --- Bounded transcendental equations (same roots as tan/cot forms) ---
   # Even: k*sin(ka) - alpha*cos(ka) = 0
   # Odd : k*cos(ka) + alpha*sin(ka) = 0
   even_eq = k*sp.sin(k*a) - alpha*sp.cos(k*a)
   odd_eq  = k*sp.cos(k*a) + alpha*sp.sin(k*a)

   # --- Unnormalized bound-state wavefunctions (symbolic Piecewise) ---
   A = sp.symbols("A", positive=True)

   psi_even_unnorm = sp.Piecewise(
       (A*sp.cos(k*x), sp.Abs(x) <= a),
       (A*sp.cos(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   psi_odd_unnorm = sp.Piecewise(
       (A*sp.sin(k*x), sp.Abs(x) <= a),
       (A*sp.sign(x)*sp.sin(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   # --- Analytic normalization: ∫ |ψ|^2 dx = 1  →  solve for A ---
   # Use symmetry (integrate x ≥ 0 and double).
   int_even_inside = sp.integrate(sp.cos(k*x)**2, (x, 0, a))
   int_even_out    = sp.integrate(sp.exp(-2*alpha*(x - a)), (x, a, sp.oo)) * sp.cos(k*a)**2
   norm_sq_even    = 2*(int_even_inside + int_even_out) * A**2
   A_even_exact    = sp.simplify(sp.sqrt(1/norm_sq_even).subs({A:1}))  # constant when A=1 inside

   int_odd_inside = sp.integrate(sp.sin(k*x)**2, (x, 0, a))
   int_odd_out    = sp.integrate(sp.exp(-2*alpha*(x - a)), (x, a, sp.oo)) * sp.sin(k*a)**2
   norm_sq_odd    = 2*(int_odd_inside + int_odd_out) * A**2
   A_odd_exact    = sp.simplify(sp.sqrt(1/norm_sq_odd).subs({A:1}))

   # Build analytically normalized symbolic ψ(x)
   psi_even = sp.Piecewise(
       (A_even_exact*sp.cos(k*x), sp.Abs(x) <= a),
       (A_even_exact*sp.cos(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   psi_odd = sp.Piecewise(
       (A_odd_exact*sp.sin(k*x), sp.Abs(x) <= a),
       (A_odd_exact*sp.sign(x)*sp.sin(k*a)*sp.exp(-alpha*(sp.Abs(x) - a)), True)
   )

   # For reference, SymPy is effectively proving the compact forms:
   # ||ψ_even||^2 = 2 [ a/2 + sin(2ka)/(4k) + cos^2(ka)/(2α) ]
   # ||ψ_odd ||^2 = 2 [ a/2 - sin(2ka)/(4k) + sin^2(ka)/(2α) ]
   # so A_even = 1/sqrt(||ψ_even||^2) and A_odd = 1/sqrt(||ψ_odd||^2).


Numeric solving (nsolve) + plots
--------------------------------

.. plot::
   :context: close-figs
   :include-source: True
   :format: python

   import numpy as np
   import matplotlib.pyplot as plt

   # --- Choose physical parameters (numbers only now) ---
   m_val, hbar_val, a_val, V0_val = 1.0, 1.0, 1.0, 50.0

   # Zero-crossing functions of ε (bounded forms), for plotting and root finding
   even_eps = sp.lambdify(
       eps, even_eq.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val}), "numpy"
   )
   odd_eps  = sp.lambdify(
       eps, odd_eq.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val}), "numpy"
   )

   # Find roots in ε (0 < ε < V0) using sympy.nsolve with a grid of initial guesses
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

   eps_even_roots = find_roots_in_eps(even_eq, 1e-6, V0_val-1e-6, ntry=200)
   eps_odd_roots  = find_roots_in_eps(odd_eq,  1e-6, V0_val-1e-6, ntry=200)

   # Map to physical energies E = -ε
   E_even = [-r for r in eps_even_roots]
   E_odd  = [-r for r in eps_odd_roots]

   # --- Plot the bounded zero-crossing equations vs E ---
   eps_grid = np.linspace(1e-6, V0_val-1e-6, 3000)
   E_grid   = -eps_grid
   Ye = even_eps(eps_grid)
   Yo = odd_eps(eps_grid)
   Lclip = 10.0
   Ye = np.clip(Ye, -Lclip, Lclip)
   Yo = np.clip(Yo, -Lclip, Lclip)

   plt.figure(figsize=(7.5, 4))
   plt.axhline(0.0, lw=1)
   plt.plot(E_grid, Ye, label=r"even: $k\sin(ka)-\alpha\cos(ka)$")
   plt.plot(E_grid, Yo, label=r"odd:  $k\cos(ka)+\alpha\sin(ka)$")
   if E_even: plt.scatter(E_even, [0]*len(E_even), s=35, marker='o', label="even roots (nsolve)")
   if E_odd:  plt.scatter(E_odd,  [0]*len(E_odd),  s=35, marker='x', label="odd roots (nsolve)")
   plt.xlim(E_grid[0], E_grid[-1])   # (-V0, 0)
   plt.xlabel("Energy E")
   plt.ylabel("Zero-crossing function")
   plt.title("Finite Square Well: bounded equations (SymPy) & roots (nsolve)")
   plt.legend()
   plt.tight_layout()
   plt.show()

   # --- Potential V(x) and analytically normalized ψ(x) at the lowest even/odd levels ---
   x = sp.symbols("x", real=True)
   Vx_sym = sp.Piecewise(
       (-V0_val, sp.And(x >= -a_val, x <= a_val)),
       (0.0, True)
   )
   Vx = sp.lambdify(x, Vx_sym, "numpy")
   xs = np.linspace(-2*a_val, 2*a_val, 2000)
   V_vals = Vx(xs)

   # Helper: build normalized ψ(x) (even/odd) at a given ε using analytic A_even_exact/A_odd_exact
   def build_psi_at_eps(eps_val):
       # numeric k, alpha
       k_val = float(sp.sqrt(2*m_val*(V0_val - eps_val))/hbar_val)
       apha  = float(sp.sqrt(2*m_val*eps_val)/hbar_val)
       # analytic normalization constants evaluated numerically
       Aeven = float(sp.N(A_even_exact.subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val, eps:eps_val})))
       Aodd  = float(sp.N(A_odd_exact .subs({m:m_val, hbar:hbar_val, a:a_val, V0:V0_val, eps:eps_val})))
       # piecewise ψ_even
       psi_even_num = sp.lambdify(
           x,
           sp.Piecewise(
               (Aeven*sp.cos(k_val*x), sp.Abs(x) <= a_val),
               (Aeven*sp.cos(k_val*a_val)*sp.exp(-apha*(sp.Abs(x) - a_val)), True)
           ),
           "numpy"
       )
       # piecewise ψ_odd
       psi_odd_num = sp.lambdify(
           x,
           sp.Piecewise(
               (Aodd*sp.sin(k_val*x), sp.Abs(x) <= a_val),
               (Aodd*sp.sign(x)*sp.sin(k_val*a_val)*sp.exp(-apha*(sp.Abs(x) - a_val)), True)
           ),
           "numpy"
       )
       return psi_even_num, psi_odd_num

   # Pick the lowest even/odd (if available)
   if eps_even_roots and eps_odd_roots:
       eps0, eps1 = eps_even_roots[0], eps_odd_roots[0]
       psi_e_num, psi_o_num = build_psi_at_eps(eps0)[0], build_psi_at_eps(eps1)[1]
       psi_e_vals = psi_e_num(xs)
       psi_o_vals = psi_o_num(xs)

       plt.figure(figsize=(8, 5))
       plt.plot(xs, V_vals, "k-", lw=2, label="V(x)")
       plt.plot(xs, -eps0 + psi_e_vals, "b", label=fr"even, $E={-eps0:.3f}$")
       plt.plot(xs, -eps1 + psi_o_vals, "r", label=fr"odd,  $E={-eps1:.3f}$")
       plt.axhline(0, color="black", lw=1)
       plt.xlabel("x")
       plt.ylabel(r"Energy / $\psi(x)$")
       plt.title("Square Well: Potential and Analytically Normalized Bound States")
       plt.legend()
       plt.tight_layout()
       plt.show()
   else:
       # Fallback: just show the potential if no roots were found (e.g., extremely shallow well)
       plt.figure(figsize=(8, 4))
       plt.plot(xs, V_vals, "k-", lw=2, label="V(x)")
       plt.axhline(0, color="black", lw=1)
       plt.xlabel("x")
       plt.ylabel("Energy")
       plt.title("Square Well Potential (no bound states found for given V0, a)")
       plt.legend()
       plt.tight_layout()
       plt.show()

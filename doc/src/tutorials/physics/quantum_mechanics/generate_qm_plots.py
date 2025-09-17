# generate_qm_plots.py
import os
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, sqrt, tan, cot, lambdify

# Assume this script is run from doc/src/tutorials/physics/quantum_mechanics/
STATIC_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "_static")
os.makedirs(STATIC_DIR, exist_ok=True)

# 1) Pauli heatmaps
mats = {'σx': [[0, 1],[1, 0]], 'σy': [[0, -1j],[1j, 0]], 'σz': [[1, 0],[0, -1]]}
fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.2))
for ax, (name, M) in zip(axes, mats.items()):
    A = np.real(np.array(M, dtype=complex))
    im = ax.imshow(A, vmin=-1, vmax=1, cmap="bwr")
    ax.set_title(name)
    ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label="value")
fig.suptitle("Pauli matrices (real part)")
fig.tight_layout()
plt.savefig(os.path.join(STATIC_DIR, "pauli_matrices.png"))
plt.close()

# 2) Creation/annihilation energy ladder
nmax = 6
ns = np.arange(0, nmax+1)
En = ns + 0.5
plt.figure()
plt.stem(ns, En)
plt.xlabel("n"); plt.ylabel("E_n (arb. units)")
plt.title("Harmonic-oscillator energy ladder")
plt.tight_layout()
plt.savefig(os.path.join(STATIC_DIR, "creation_annihilation_ladder.png"))
plt.close()

# 3) Particle in a box eigenfunctions
L = 1.0
xs = np.linspace(0, L, 400)
def psi(n, x):
    return (2.0/L)**0.5 * np.sin(n * np.pi * x / L)
plt.figure()
for n in (1, 2, 3):
    plt.plot(xs, psi(n, xs), label=fr"$\psi_{{{n}}}(x)$")
plt.title("Particle in a box: first three eigenfunctions")
plt.xlabel("x"); plt.ylabel(r"$\psi_n(x)$")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(STATIC_DIR, "particle_in_box.png"))
plt.close()

# 4) Finite well even/odd equations
E = symbols('E')
m, hbar, a, V0 = symbols('m hbar a V0', positive=True)
params = {m:1.0, hbar:1.0, a:1.0, V0:50.0}
k = sqrt(2*m*(E + V0))/hbar
alpha = sqrt(-2*m*E)/hbar
f_even = k*tan(k*a) - alpha
f_odd  = -k*cot(k*a) - alpha
fe = lambdify(E, f_even.subs(params), 'numpy')
fo = lambdify(E, f_odd.subs(params),  'numpy')
Es = np.linspace(-params[V0]+1e-6, -1e-3, 4000)

def safe_eval(func, xs):
    ys = np.empty_like(xs, dtype=float); ys[:] = np.nan
    for i, x in enumerate(xs):
        try:
            y = func(x)
            ys[i] = y if np.isfinite(y) and abs(y)<50 else np.nan
        except: ys[i] = np.nan
    return ys

Ye, Yo = safe_eval(fe, Es), safe_eval(fo, Es)
plt.figure()
plt.axhline(0, lw=1)
plt.plot(Es, Ye, label="even")
plt.plot(Es, Yo, label="odd")
plt.xlabel("Energy E"); plt.ylabel("Equation value")
plt.title("Finite square well transcendental equations")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(STATIC_DIR, "finite_well.png"))
plt.close()

print("Plots generated in:", STATIC_DIR)

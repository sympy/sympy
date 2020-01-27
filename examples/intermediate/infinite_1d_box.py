#!/usr/bin/env python

"""
Applying perturbation theory to calculate the ground state energy
of the infinite 1D box of width ``a`` with a perturbation
which is linear in ``x``, up to second order in perturbation
"""

from sympy.core import pi
from sympy import Integral, var, S
from sympy.functions import sin, sqrt


def X_n(n, a, x):
    """
    Returns the wavefunction X_{n} for an infinite 1D box

    ``n``
        the "principal" quantum number. Corresponds to the number of nodes in
        the wavefunction.  n >= 0
    ``a``
        width of the well. a > 0
    ``x``
        x coordinate.
    """
    n, a, x = map(S, [n, a, x])
    C = sqrt(2 / a)
    return C * sin(pi * n * x / a)


def E_n(n, a, mass):
    """
    Returns the Energy psi_{n} for a 1d potential hole with infinity borders

    ``n``
        the "principal" quantum number. Corresponds to the number of nodes in
        the wavefunction.  n >= 0
    ``a``
        width of the well. a > 0
    ``mass``
        mass.
    """
    return ((n * pi / a)**2) / mass


def energy_corrections(perturbation, n, a=10, mass=0.5):
    """
    Calculating first two order corrections due to perturbation theory and
    returns tuple where zero element is unperturbated energy, and two second
    is corrections

    ``n``
        the "nodal" quantum number. Corresponds to the number of nodes in the
        wavefunction.  n >= 0
    ``a``
        width of the well. a > 0
    ``mass``
        mass.

    """
    x, _a = var("x _a")

    Vnm = lambda n, m, a: Integral(X_n(n, a, x) * X_n(m, a, x)
        * perturbation.subs({_a: a}), (x, 0, a)).n()

    # As we know from theory for V0*r/a we will just V(n, n-1) and V(n, n+1)
    #   wouldn't equals zero

    return (E_n(n, a, mass).evalf(),

            Vnm(n, n, a).evalf(),

            (Vnm(n, n - 1, a)**2/(E_n(n, a, mass) - E_n(n - 1, a, mass))
           + Vnm(n, n + 1, a)**2/(E_n(n, a, mass) - E_n(n + 1, a, mass))).evalf())


def main():
    print()
    print("Applying perturbation theory to calculate the ground state energy")
    print("of the infinite 1D box of width ``a`` with a perturbation")
    print("which is linear in ``x``, up to second order in perturbation.")
    print()

    x, _a = var("x _a")
    perturbation = .1 * x / _a

    E1 = energy_corrections(perturbation, 1)
    print("Energy for first term (n=1):")
    print("E_1^{(0)} = ", E1[0])
    print("E_1^{(1)} = ", E1[1])
    print("E_1^{(2)} = ", E1[2])
    print()

    E2 = energy_corrections(perturbation, 2)
    print("Energy for second term (n=2):")
    print("E_2^{(0)} = ", E2[0])
    print("E_2^{(1)} = ", E2[1])
    print("E_2^{(2)} = ", E2[2])
    print()

    E3 = energy_corrections(perturbation, 3)
    print("Energy for third term (n=3):")
    print("E_3^{(0)} = ", E3[0])
    print("E_3^{(1)} = ", E3[1])
    print("E_3^{(2)} = ", E3[2])
    print()


if __name__ == "__main__":
    main()

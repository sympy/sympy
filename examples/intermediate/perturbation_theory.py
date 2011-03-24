#!/usr/bin/env python
"""
Applying perturbation theory to calculate the ground state energy
of the 1D infinite potential well of width $a$ with a perturbation
which is linear in $x$, up to second order in perturbation by
K. Anatolii.
"""

from sympy import Integral, var, S
from sympy.core import pi
from sympy.functions import sin, sqrt

def R_n(n, a, r):
    """
    Returns the wavefunction R_{n} for a 1d potential hole with infinity borders

    ``n``
        the "nodal" quantum number.  Corresponds to the number of nodes in the
        wavefunction.  n >= 0
    ``a``
        Size of potential hole. a > 0
    ``r``
        r coordinate.
    """
    n, a, r = map(S, [n, a, r])
    C = sqrt(2/a)
    return C*sin(pi*n*r/a);


def E_n(n, a, mass):
    """
    Returns the Energy psi_{n} for a 1d potential hole with infinity borders

    ``n``
        the "nodal" quantum number.  Corresponds to the number of nodes in the
        wavefunction.  n >= 0
    ``a``
        Size of potential hole. a > 0
    ``mass``
        mass.
    """
    return ((n*pi/a)**2)/mass;

def EnergyCorrections(_n, _a=10, _mass=0.5):
    """
    Calculating first two order corrections due to perturbation theory and returns
    tuple where zero element is unperturbated energy, and two second is corrections

    ``_n``
        the "nodal" quantum number.  Corresponds to the number of nodes in the
        wavefunction.  n >= 0
    ``_a``
        Size of potential hole. a > 0
    ``_mass``
        mass.

    """
    r, a = var("r a");
    perturbation = .1*r/a;

    Vnm = lambda _n, _m, _a: Integral(R_n(_n, _a, r)*R_n(_m, _a, r)
        *perturbation.subs({a: _a}), (r, 0, _a)).n();

    # As we know from theory for V0*r/a we will just V(n, n-1) and V(n, n+1)
    #   wouldn't equals zero

    return (E_n(_n, _a, _mass).evalf(),

            Vnm(_n, _n, _a).evalf(),

            (Vnm(_n, _n-1, _a)**2/(E_n(_n, _a, _mass)- E_n(_n-1, _a, _mass))
           + Vnm(_n, _n+1, _a)**2/(E_n(_n, _a, _mass) - E_n(_n+1, _a, _mass))).evalf());

def main():
    print;
    print "Applying perturbation theory to calculate the ground state energy";
    print "of the 1D infinite potential well of width $a$ with a perturbation";
    print "which is linear in $x$, up to second order in perturbation.";
    print "by K. Anatolii";
    print;

    E1 = EnergyCorrections(1);
    print "Energy for first term (n=1):";
    print "E_1^{(0)} = ", E1[0];
    print "E_1^{(1)} = ", E1[1];
    print "E_1^{(2)} = ", E1[2];
    print;

    E2 = EnergyCorrections(2);
    print "Energy for second term (n=2):";
    print "E_2^{(0)} = ", E2[0];
    print "E_2^{(1)} = ", E2[1];
    print "E_2^{(2)} = ", E2[2];
    print;

    E3 = EnergyCorrections(3);
    print "Energy for third term (n=3):";
    print "E_3^{(0)} = ", E3[0];
    print "E_3^{(1)} = ", E3[1];
    print "E_3^{(2)} = ", E3[2];
    print;


if __name__ == "__main__":
    main();

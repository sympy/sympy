#!/usr/bin/env python

"""
This example shows how to work with the Hydrogen radial wavefunctions.
"""

from sympy import Eq, Integral, oo, pprint, symbols
from sympy.physics.hydrogen import R_nl


def main():
    print("Hydrogen radial wavefunctions:")
    a, r = symbols("a r")
    print("R_{21}:")
    pprint(R_nl(2, 1, a, r))
    print("R_{60}:")
    pprint(R_nl(6, 0, a, r))

    print("Normalization:")
    i = Integral(R_nl(1, 0, 1, r)**2 * r**2, (r, 0, oo))
    pprint(Eq(i, i.doit()))
    i = Integral(R_nl(2, 0, 1, r)**2 * r**2, (r, 0, oo))
    pprint(Eq(i, i.doit()))
    i = Integral(R_nl(2, 1, 1, r)**2 * r**2, (r, 0, oo))
    pprint(Eq(i, i.doit()))

if __name__ == '__main__':
    main()

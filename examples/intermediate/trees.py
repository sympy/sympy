#!/usr/bin/env python

"""
Calculates the Sloane's A000055 integer sequence, i.e. the "Number of
trees with n unlabeled nodes."

You can also google for "The On-Line Encyclopedia of Integer Sequences"
and paste in the sequence returned by this script:

1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106

and it will shows you the A000055
"""

from sympy import Symbol, Poly


def T(x):
    return x + x**2 + 2*x**3 + 4*x**4 + 9*x**5 + 20*x**6 + 48 * x**7 + \
        115*x**8 + 286*x**9 + 719*x**10


def A(x):
    return 1 + T(x) - T(x)**2/2 + T(x**2)/2


def main():
    x = Symbol("x")
    s = Poly(A(x), x)
    num = list(reversed(s.coeffs()))[:11]

    print(s.as_expr())
    print(num)

if __name__ == "__main__":
    main()

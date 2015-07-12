#!/usr/bin/env python

"""Pretty print example

Demonstrates pretty printing.
"""

from __future__ import division, print_function

from sympy import Symbol, pprint, sin, cos, exp, sqrt


def main():
    x = Symbol("x")
    y = Symbol("y")

    pprint( x**x )
    print('\n')  # separate with two blank likes

    pprint(x**2 + y + x)
    print('\n')

    pprint(sin(x)**x)
    print('\n')

    pprint( sin(x)**cos(x) )
    print('\n')

    pprint( sin(x)/(cos(x)**2 * x**x + (2*y)) )
    print('\n')

    pprint( sin(x**2 + exp(x)) )
    print('\n')

    pprint( sqrt(exp(x)) )
    print('\n')

    pprint( sqrt(sqrt(exp(x))) )
    print('\n')

    pprint( (1/cos(x)).series(x, 0, 10) )
    print('\n')

if __name__ == "__main__":
    main()

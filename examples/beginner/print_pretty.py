#!/usr/bin/env python

"""Pretty print example

Demonstrates pretty printing.
"""

from sympy import Symbol, pprint, sin, cos, exp, sqrt, MatrixSymbol, KroneckerProduct


def main():
    x = Symbol("x")
    y = Symbol("y")

    a = MatrixSymbol("a", 1, 1)
    b = MatrixSymbol("b", 1, 1)
    c = MatrixSymbol("c", 1, 1)

    pprint( x**x )
    print('\n')  # separate with two blank lines

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

    pprint(a*(KroneckerProduct(b, c)))
    print('\n')

if __name__ == "__main__":
    main()

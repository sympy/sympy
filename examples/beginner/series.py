#!/usr/bin/env python

"""Series example

Demonstrates series.
"""

from __future__ import division, print_function

from sympy import Symbol, cos, sin, pprint


def main():
    x = Symbol('x')

    e = 1/cos(x)
    print()
    print("Series for sec(x):")
    print()
    pprint(e.series(x, 0, 10))
    print("\n")

    e = 1/sin(x)
    print("Series for csc(x):")
    print()
    pprint(e.series(x, 0, 4))
    print()

if __name__ == "__main__":
    main()

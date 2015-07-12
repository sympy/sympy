#!/usr/bin/env python

"""Basic example

Demonstrates how to create symbols and print some algebra operations.
"""

from __future__ import division, print_function

import sympy
from sympy import pprint


def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    c = sympy.Symbol('c')
    e = ( a*b*b + 2*b*a*b )**c

    print
    pprint(e)
    print

if __name__ == "__main__":
    main()

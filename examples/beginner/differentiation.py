#!/usr/bin/env python

"""Differentiation example

Demonstrates some differentiation operations.
"""

import sympy
from sympy import pprint

def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = (a + 2*b)**5

    print
    pprint(e)
    print
    pprint(e.diff(a))
    print
    pprint(e.diff(b))
    print
    pprint(e.diff(b).diff(a, 3))
    print
    pprint(e.expand().diff(b).diff(a, 3))
    print

if __name__ == "__main__":
    main()

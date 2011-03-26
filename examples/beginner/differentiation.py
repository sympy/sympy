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

    pprint(e)
    pprint(e.diff(a))
    pprint(e.diff(b))
    pprint(e.diff(b).diff(a, 3))
    pprint(e.expand().diff(b).diff(a, 3))

if __name__ == "__main__":
    main()

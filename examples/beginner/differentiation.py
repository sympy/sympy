#!/usr/bin/env python
"""Differentiation example

Demonstrates some differentiation operations.
"""

import sympy

def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = (a + 2*b)**5

    print e
    print e.diff(a)
    print e.diff(b)
    print e.diff(b).diff(a, 3)
    print e.expand().diff(b).diff(a, 3)

if __name__ == "__main__":
    main()

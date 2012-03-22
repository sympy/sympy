#!/usr/bin/env python

"""Precision Example

Demonstrates SymPy's arbitrary precision abilities
"""

import sympy
from sympy import pprint, Symbol

def main():
    x = Pow(2, 50, evaluate=False)
    y = Pow(10, -50, evaluate=False)
    print
    print Mul(x, y, evaluate=False)
    print
    e = sympy.Rational(2)**50/sympy.Rational(10)**50
    pprint(e)

if __name__ == "__main__":
    main()

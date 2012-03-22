#!/usr/bin/env python

"""Precision Example

Demonstrates SymPy's arbitrary precision abilities
"""

import sympy
from sympy import pprint, Symbol

def main():
    x = Symbol('2')
    y = Symbol('10')
    pprint(x**50)
    print '__'
    pprint(y**50)
    print
    e = sympy.Rational(2)**50/sympy.Rational(10)**50
    pprint(e)

if __name__ == "__main__":
    main()

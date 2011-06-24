#!/usr/bin/env python
"""Precision Example

Demonstrates SymPy's arbitrary precision abilities
"""

import sympy
from sympy import pprint

def main():
    e = sympy.Rational(2)**50/sympy.Rational(10)**50
    pprint(e)

if __name__ == "__main__":
    main()

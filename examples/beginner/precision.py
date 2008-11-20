#!/usr/bin/env python
"""Precision Example

Demonstrates SymPy's arbitrary precision abilities
"""

import sympy

def main():
    e = sympy.Rational(2)**50/sympy.Rational(10)**50
    print e

if __name__ == "__main__":
    main()

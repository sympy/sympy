#!/usr/bin/env python

"""Precision Example

Demonstrates SymPy's arbitrary integer precision abilities
"""

import sympy
from sympy import Mul, Pow, S


def main():
    x = Pow(2, 50, evaluate=False)
    y = Pow(10, -50, evaluate=False)
    # A large, unevaluated expression
    m = Mul(x, y, evaluate=False)
    # Evaluating the expression
    e = S(2)**50/S(10)**50
    print("%s == %s" % (m, e))

if __name__ == "__main__":
    main()

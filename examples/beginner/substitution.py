#!/usr/bin/env python
"""Substitution example

Demonstrates substitution.
"""

import sympy
from sympy import pprint

def main():
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    e = 1/sympy.cos(x)
    pprint(e)
    pprint(e.subs(sympy.cos(x), y))
    pprint(e.subs(sympy.cos(x), y).subs(y, x**2))

    e = 1/sympy.log(x)
    e = e.subs(x, sympy.Real("2.71828"))
    pprint(e)
    pprint(e.evalf())

if __name__ == "__main__":
    main()

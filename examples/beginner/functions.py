#!/usr/bin/env python
"""Functions example

Demonstrates sympy defined functions.
"""

import sympy
from sympy import pprint

def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = sympy.log((a + b)**5)
    pprint(e)

    e = sympy.exp(e)
    pprint(e)

    e = sympy.log(sympy.exp((a + b)**5))
    pprint(e)

if __name__ == "__main__":
    main()

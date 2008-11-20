#!/usr/bin/env python
"""Functions example

Demonstrates sympy defined functions.
"""

import sympy

def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = sympy.log((a + b)**5)
    print e

    e = sympy.exp(e)
    print e

    e = sympy.log(sympy.exp((a + b)**5))
    print e

if __name__ == "__main__":
    main()

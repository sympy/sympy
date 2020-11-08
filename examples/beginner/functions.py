#!/usr/bin/env python

"""Functions example

Demonstrates functions defined in SymPy.
"""

from sympy import pprint, Symbol, log, exp

def main():
    a = Symbol('a')
    b = Symbol('b')
    e = log((a + b)**5)
    print()
    pprint(e)
    print('\n')

    e = exp(e)
    pprint(e)
    print('\n')

    e = log(exp((a + b)**5))
    pprint(e)
    print()

if __name__ == "__main__":
    main()

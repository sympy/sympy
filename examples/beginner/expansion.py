#!/usr/bin/env python

"""Expansion Example

Demonstrates how to expand expressions.
"""

import sympy
from sympy import pprint

def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = (a + b)**5

    print
    pprint(e)
    print '\n'
    pprint(e.expand())
    print

if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""Expansion Example

Demonstrates how to expand expressions.
"""

from __future__ import division, print_function

import sympy
from sympy import pprint


def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = (a + b)**5

    print("\nExpression:")
    pprint(e)
    print('\nExpansion of the above expression:')
    pprint(e.expand())
    print()

if __name__ == "__main__":
    main()

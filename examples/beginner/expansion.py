#!/usr/bin/env python

"""Expansion Example

Demonstrates how to expand expressions.
"""

from sympy import pprint, Symbol

def main():
    a = Symbol('a')
    b = Symbol('b')
    e = (a + b)**5

    print("\nExpression:")
    pprint(e)
    print('\nExpansion of the above expression:')
    pprint(e.expand())
    print()

if __name__ == "__main__":
    main()

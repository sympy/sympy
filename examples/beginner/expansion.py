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

    pprint(e)
    pprint(e.expand())

if __name__ == "__main__":
    main()

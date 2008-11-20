#!/usr/bin/env python
"""Expansion Example

Demonstrates how to expand expressions.
"""

import sympy

def main():
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')
    e = (a + b)**5

    print e
    print e.expand()

if __name__ == "__main__":
    main()

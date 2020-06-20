#!/usr/bin/env python

"""Differentiation example

Demonstrates some differentiation operations.
"""

from sympy import pprint, Symbol

def main():
    a = Symbol('a')
    b = Symbol('b')
    e = (a + 2*b)**5

    print("\nExpression : ")
    print()
    pprint(e)
    print("\n\nDifferentiating w.r.t. a:")
    print()
    pprint(e.diff(a))
    print("\n\nDifferentiating w.r.t. b:")
    print()
    pprint(e.diff(b))
    print("\n\nSecond derivative of the above result w.r.t. a:")
    print()
    pprint(e.diff(b).diff(a, 2))
    print("\n\nExpanding the above result:")
    print()
    pprint(e.expand().diff(b).diff(a, 2))
    print()

if __name__ == "__main__":
    main()

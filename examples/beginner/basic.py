#!/usr/bin/env python

"""Basic example

Demonstrates how to create symbols and print some algebra operations.
"""

from sympy import Symbol, pprint

def normal_equation():
    "This example is for adding the example of normal equation"
    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    e = (a * b +b*c+c*a)

    print('')
    pprint(e)
    print('')




def squared_equation_example():
    "This Function is for having the example of Squared Equation"
    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    e = (a * b * b + 2 * b * a * b) ** c

    print('')
    pprint(e)
    print('')


if __name__ == "__main__":
    squared_equation_example()
    normal_equation()

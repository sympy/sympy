#!/usr/bin/env python

"""Basic example

Demonstrates how to create symbols and print some algebra operations.
"""

from sympy import Symbol, pprint



def Normal_Equation():
    a=Symbol('a')
    b=Symbol('b')
    c=Symbol('c')
    e=(a*b)+(a*c)+(b*c)
    print('')
    print(e)
    print('')

def for_Squared_Expression():
    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    e = ( a*b*b + 2*b*a*b )**c


    print('')
    pprint(e)
    print('')

if __name__ == "__main__":
    for_Squared_Expression()
    Normal_Equation()

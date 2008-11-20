#!/usr/bin/env python
"""print_gtk example

Demonstrates printing with gtkmathview using mathml
"""

from sympy import Symbol, integrate, exp
from sympy.printing import print_gtk

def main():
    x = Symbol('x')

    #l1 = limit(sin(x)/x, x, 0, evaluate=False)
    #print_gtk(l1)

    l2 = integrate(exp(x), (x,0,1), evaluate=False)
    print_gtk(l2)

if __name__ == "__main__":
    main()

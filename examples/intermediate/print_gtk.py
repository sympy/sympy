#!/usr/bin/env python

"""print_gtk example

Demonstrates printing with gtkmathview using mathml
"""

from sympy import Integral, Limit, print_gtk, sin, Symbol


def main():
    x = Symbol('x')

    example_limit = Limit(sin(x)/x, x, 0)
    print_gtk(example_limit)

    example_integral = Integral(x, (x, 0, 1))
    print_gtk(example_integral)

if __name__ == "__main__":
    main()

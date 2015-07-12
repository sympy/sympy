#!/usr/bin/env python

"""Limits Example

Demonstrates limits.
"""

from __future__ import division, print_function

from sympy import exp, log, Symbol, Rational, sin, limit, sqrt, oo


def sqrt3(x):
    return x**Rational(1, 3)


def show(computed, correct):
    print("computed:", computed, "correct:", correct)


def main():
    x = Symbol("x")
    a = Symbol("a")
    h = Symbol("h")

    show( limit(sqrt(x**2 - 5*x + 6) - x, x, oo), -Rational(5)/2 )

    show( limit(x*(sqrt(x**2 + 1) - x), x, oo), Rational(1)/2 )

    show( limit(x - sqrt3(x**3 - 1), x, oo), Rational(0) )

    show( limit(log(1 + exp(x))/x, x, -oo), Rational(0) )

    show( limit(log(1 + exp(x))/x, x, oo), Rational(1) )

    show( limit(sin(3*x)/x, x, 0), Rational(3) )

    show( limit(sin(5*x)/sin(2*x), x, 0), Rational(5)/2 )

    show( limit(((x - 1)/(x + 1))**x, x, oo), exp(-2))

if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""Differential equations example

Demonstrates solving 1st and 2nd degree linear ordinary differential
equations.
"""

from sympy import dsolve, Eq, Function, sin, Symbol


def main():
    x = Symbol("x")
    f = Function("f")

    eq = Eq(f(x).diff(x), f(x))
    print("Solution for ", eq, " : ", dsolve(eq, f(x)))

    eq = Eq(f(x).diff(x, 2), -f(x))
    print("Solution for ", eq, " : ", dsolve(eq, f(x)))

    eq = Eq(x**2*f(x).diff(x), -3*x*f(x) + sin(x)/x)
    print("Solution for ", eq, " : ", dsolve(eq, f(x)))


if __name__ == "__main__":
    main()

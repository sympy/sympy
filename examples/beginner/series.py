#!/usr/bin/env python
"""Series example

Demonstrates series.
"""

from sympy import Symbol, cos, sin

def main():
    x = Symbol('x')

    e = 1/cos(x)
    print "Series for sec(x):"
    print e.series(x, 0, 10)
    print ""

    e = 1/sin(x)
    print "Series for csc(x):"
    print e.series(x, 0, 4)

if __name__ == "__main__":
    main()

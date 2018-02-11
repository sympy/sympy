#!/usr/bin/env python

"""Plotting example

Demonstrates simple plotting.
"""

from sympy import Symbol, cos, sin, log, tan
from sympy.plotting import PygletPlot
from sympy.abc import x, y


def main():
    fun1 = cos(x)*sin(y)
    fun2 = sin(x)*sin(y)
    fun3 = cos(y) + log(tan(y/2)) + 0.2*x

    PygletPlot(fun1, fun2, fun3, [x, -0.00, 12.4, 40], [y, 0.1, 2, 40])

if __name__ == "__main__":
    main()

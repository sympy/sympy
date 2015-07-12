#!/usr/bin/env python

"""Matplotlib 2D plotting example

Demonstrates plotting with matplotlib.
"""

from __future__ import division, print_function

import sys

from sample import sample

from sympy import log, pi, sqrt, sin, Symbol
from sympy.core.compatibility import is_sequence
from sympy.external import import_module


def mplot2d(f, var, show=True):
    """
    Plot a 2d function using matplotlib/Tk.
    """

    import warnings
    warnings.filterwarnings("ignore", "Could not match \S")

    p = import_module('pylab')
    if not p:
        sys.exit("Matplotlib is required to use mplot2d.")

    if not is_sequence(f):
        f = [f, ]

    for f_i in f:
        x, y = sample(f_i, var)
        p.plot(x, y)

    p.draw()
    if show:
        p.show()


def main():
    x = Symbol('x')

    # mplot2d(log(x), (x, 0, 2, 100))
    # mplot2d([sin(x), -sin(x)], (x, float(-2*pi), float(2*pi), 50))
    mplot2d([sqrt(x), -sqrt(x), sqrt(-x), -sqrt(-x)], (x, -40.0, 40.0, 80))

if __name__ == "__main__":
    main()

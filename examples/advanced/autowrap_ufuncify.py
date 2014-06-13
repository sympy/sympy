#!/usr/bin/env python
"""
Setup ufuncs for the legendre polynomials
-----------------------------------------

This example demonstrates how you can use the ufuncify utility in SymPy
to create fast, customized universal functions for use with numpy
arrays. An autowrapped sympy expression can be significantly faster than
what you would get by applying a sequence of the ufuncs shipped with
numpy. [0]

You need to have numpy installed to run this example, as well as a
working fortran compiler.


[0]:
http://ojensen.wordpress.com/2010/08/10/fast-ufunc-ish-hydrogen-solutions/
"""

import sys

from sympy.external import import_module

np = import_module('numpy')
if not np:
    sys.exit("Cannot import numpy. Exiting.")

import mpmath
from sympy.utilities.autowrap import ufuncify
from sympy.utilities.lambdify import implemented_function
from sympy import symbols, legendre, Plot, pprint


def main():

    print(__doc__)

    x = symbols('x')

    # a numpy array we can apply the ufuncs to
    grid = np.linspace(-1, 1, 1000)

    # set mpmath precision to 20 significant numbers for verification
    mpmath.mp.dps = 20

    print("Compiling legendre ufuncs and checking results:")

    # Let's also plot the ufunc's we generate
    plot1 = Plot(visible=False)
    for n in range(6):

        # Setup the SymPy expression to ufuncify
        expr = legendre(n, x)
        print("The polynomial of degree %i is" % n)
        pprint(expr)

        # This is where the magic happens:
        binary_poly = ufuncify(x, expr)

        # It's now ready for use with numpy arrays
        polyvector = binary_poly(grid)

        # let's check the values against mpmath's legendre function
        maxdiff = 0
        for j in range(len(grid)):
            precise_val = mpmath.legendre(n, grid[j])
            diff = abs(polyvector[j] - precise_val)
            if diff > maxdiff:
                maxdiff = diff
        print("The largest error in applied ufunc was %e" % maxdiff)
        assert maxdiff < 1e-14

        # We can also attach the autowrapped legendre polynomial to a sympy
        # function and plot values as they are calculated by the binary function
        g = implemented_function('g', binary_poly)
        plot1[n] = g(x), [200]

    print("Here's a plot with values calculated by the wrapped binary functions")
    plot1.show()

if __name__ == '__main__':
    main()

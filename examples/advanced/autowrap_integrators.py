#!/usr/bin/env python
"""
Numerical integration with autowrap
-----------------------------------

This example demonstrates how you can use the autowrap module in SymPy
to create fast, numerical integration routines callable from python. See
in the code for detailed explanations of the various steps. An
autowrapped sympy expression can be significantly faster than what you
would get by applying a sequence of the ufuncs shipped with numpy. [0]

We will find the coefficients needed to approximate a quantum mechanical
Hydrogen wave function in terms of harmonic oscillator solutions. For
the sake of demonstration, this will be done by setting up a simple
numerical integration scheme as a SymPy expression, and obtain a binary
implementation with autowrap.

You need to have numpy installed to run this example, as well as a
working fortran compiler. If you have pylab installed, you will be
rewarded with a nice plot in the end.

[0]:
http://ojensen.wordpress.com/2010/08/10/fast-ufunc-ish-hydrogen-solutions/

----
"""

from __future__ import division, print_function

import sys
from sympy.external import import_module

np = import_module('numpy')
if not np:
    sys.exit("Cannot import numpy. Exiting.")
pylab = import_module('pylab', warn_not_installed=True)

from sympy.utilities.lambdify import implemented_function
from sympy.utilities.autowrap import autowrap, ufuncify
from sympy import Idx, IndexedBase, Lambda, pprint, Symbol, oo, Integral,\
    Function
from sympy.physics.sho import R_nl
from sympy.physics.hydrogen import R_nl as hydro_nl


# ***************************************************************************
# calculation parameters to play with
# ***************************************************************************

basis_dimension = 5         # Size of h.o. basis (n < basis_dimension)
omega2 = 0.1                # in atomic units: twice the oscillator frequency
orbital_momentum_l = 1      # the quantum number `l` for angular momentum
hydrogen_n = 2              # the nodal quantum number for the Hydrogen wave
rmax = 20                   # cut off in the radial direction
gridsize = 200              # number of points in the grid

# ***************************************************************************


def main():

    print(__doc__)

    # arrays are represented with IndexedBase, indices with Idx
    m = Symbol('m', integer=True)
    i = Idx('i', m)
    A = IndexedBase('A')
    B = IndexedBase('B')
    x = Symbol('x')

    print("Compiling ufuncs for radial harmonic oscillator solutions")

    # setup a basis of ho-solutions  (for l=0)
    basis_ho = {}
    for n in range(basis_dimension):

        # Setup the radial ho solution for this n
        expr = R_nl(n, orbital_momentum_l, omega2, x)

        # Reduce the number of operations in the expression by eval to float
        expr = expr.evalf(15)

        print("The h.o. wave function with l = %i and n = %i is" % (
            orbital_momentum_l, n))
        pprint(expr)

        # implement, compile and wrap it as a ufunc
        basis_ho[n] = ufuncify(x, expr)

    # now let's see if we can express a hydrogen radial wave in terms of
    # the ho basis.  Here's the solution we will approximate:
    H_ufunc = ufuncify(x, hydro_nl(hydrogen_n, orbital_momentum_l, 1, x))

    # The transformation to a different basis can be written like this,
    #
    #   psi(r) = sum_i c(i) phi_i(r)
    #
    # where psi(r) is the hydrogen solution, phi_i(r) are the H.O. solutions
    # and c(i) are scalar coefficients.
    #
    # So in order to express a hydrogen solution in terms of the H.O. basis, we
    # need to determine the coefficients c(i).  In position space, it means
    # that we need to evaluate an integral:
    #
    #  psi(r) = sum_i Integral(R**2*conj(phi(R))*psi(R), (R, 0, oo)) phi_i(r)
    #
    # To calculate the integral with autowrap, we notice that it contains an
    # element-wise sum over all vectors.  Using the Indexed class, it is
    # possible to generate autowrapped functions that perform summations in
    # the low-level code.  (In fact, summations are very easy to create, and as
    # we will see it is often necessary to take extra steps in order to avoid
    # them.)
    # we need one integration ufunc for each wave function in the h.o. basis
    binary_integrator = {}
    for n in range(basis_dimension):

        #
        # setup basis wave functions
        #
        # To get inline expressions in the low level code, we attach the
        # wave function expressions to a regular SymPy function using the
        # implemented_function utility.  This is an extra step needed to avoid
        # erronous summations in the wave function expressions.
        #
        # Such function objects carry around the expression they represent,
        # but the expression is not exposed unless explicit measures are taken.
        # The benefit is that the routines that searches for repeated indices
        # in order to make contractions will not search through the wave
        # function expression.
        psi_ho = implemented_function('psi_ho',
                Lambda(x, R_nl(n, orbital_momentum_l, omega2, x)))

        # We represent the hydrogen function by an array which will be an input
        # argument to the binary routine.  This will let the integrators find
        # h.o. basis coefficients for any wave function we throw at them.
        psi = IndexedBase('psi')

        #
        # setup expression for the integration
        #

        step = Symbol('step')  # use symbolic stepsize for flexibility

        # let i represent an index of the grid array, and let A represent the
        # grid array.  Then we can approximate the integral by a sum over the
        # following expression (simplified rectangular rule, ignoring end point
        # corrections):

        expr = A[i]**2*psi_ho(A[i])*psi[i]*step

        if n == 0:
            print("Setting up binary integrators for the integral:")
            pprint(Integral(x**2*psi_ho(x)*Function('psi')(x), (x, 0, oo)))

        # Autowrap it.  For functions that take more than one argument, it is
        # a good idea to use the 'args' keyword so that you know the signature
        # of the wrapped function.  (The dimension m will be an optional
        # argument, but it must be present in the args list.)
        binary_integrator[n] = autowrap(expr, args=[A.label, psi.label, step, m])

        # Lets see how it converges with the grid dimension
        print("Checking convergence of integrator for n = %i" % n)
        for g in range(3, 8):
            grid, step = np.linspace(0, rmax, 2**g, retstep=True)
            print("grid dimension %5i, integral = %e" % (2**g,
                    binary_integrator[n](grid, H_ufunc(grid), step)))

    print("A binary integrator has been set up for each basis state")
    print("We will now use them to reconstruct a hydrogen solution.")

    # Note: We didn't need to specify grid or use gridsize before now
    grid, stepsize = np.linspace(0, rmax, gridsize, retstep=True)

    print("Calculating coefficients with gridsize = %i and stepsize %f" % (
        len(grid), stepsize))

    coeffs = {}
    for n in range(basis_dimension):
        coeffs[n] = binary_integrator[n](grid, H_ufunc(grid), stepsize)
        print("c(%i) = %e" % (n, coeffs[n]))

    print("Constructing the approximate hydrogen wave")
    hydro_approx = 0
    all_steps = {}
    for n in range(basis_dimension):
        hydro_approx += basis_ho[n](grid)*coeffs[n]
        all_steps[n] = hydro_approx.copy()
        if pylab:
            line = pylab.plot(grid, all_steps[n], ':', label='max n = %i' % n)

    # check error numerically
    diff = np.max(np.abs(hydro_approx - H_ufunc(grid)))
    print("Error estimate: the element with largest deviation misses by %f" % diff)
    if diff > 0.01:
        print("This is much, try to increase the basis size or adjust omega")
    else:
        print("Ah, that's a pretty good approximation!")

    # Check visually
    if pylab:
        print("Here's a plot showing the contribution for each n")
        line[0].set_linestyle('-')
        pylab.plot(grid, H_ufunc(grid), 'r-', label='exact')
        pylab.legend()
        pylab.show()

    print("""Note:
    These binary integrators were specialized to find coefficients for a
    harmonic oscillator basis, but they can process any wave function as long
    as it is available as a vector and defined on a grid with equidistant
    points. That is, on any grid you get from numpy.linspace.

    To make the integrators even more flexible, you can setup the harmonic
    oscillator solutions with symbolic parameters omega and l.  Then the
    autowrapped binary routine will take these scalar variables as arguments,
    so that the integrators can find coefficients for *any* isotropic harmonic
    oscillator basis.

    """)


if __name__ == '__main__':
    main()

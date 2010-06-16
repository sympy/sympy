"""
Algorithms for solving the Risch differential equation.

Given a differential field K of characteristic 0 that is a simple monomial
extension of a base field k and f, g in K, the Risch Differential Equation
problem is to decide if there exist y in K such that Dy + f*y == g and to find
one if there are some.  If t is a monomial over k and the coefficients of f and
g are in k(t), then y is in k(t), and the out line of the algorithm here is
given as:

1. Compute the normal part n of the denominator of y.  The problem is then
reduced to finding y' in k<t>, where y == y'/n.
2. Compute the special part s of the denominator of y.   The problem is then
reduced to finding y'' in k[t], where y == y''/(n*s)
3. Bound the degree of y''.
4. Reduce the equation Dy + f*y == g to a similar equation with f, g in k[t].
5. Find the solutions in k[t] of bounded degree of the reduced equation.

See Chapter 6 of "Symbolic Integration I: Transcendental Functions" by Manuel
Bronstein.
"""
from sympy.polys import Poly
#    from pudb import set_trace; set_trace() # Debugging

"""Module with some routines for polynomials"""

from sympy.polynomials.base import PolynomialException, Polynomial
from sympy.polynomials.roots_ import count_real_roots, roots, solve_system
from sympy.polynomials.wrapper import div, factor, gcd, groebner, lcm, \
     resultant, sqf, sqf_part
#from sympy.polynomials.ideals import Ideal

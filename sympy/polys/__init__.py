
from monomial import monomials, monomial_count
from polynomial import Poly, PolynomialError

from algorithms import poly_div, poly_pdiv

def poly_quo(f, g, *symbols):
    return poly_div(f, g, *symbols)[0]

def poly_rem(f, g, *symbols):
    return poly_div(f, g, *symbols)[1]

def poly_pquo(f, g, *symbols):
    return poly_pdiv(f, g, *symbols)[0]

def poly_prem(f, g, *symbols):
    return poly_pdiv(f, g, *symbols)[1]

from algorithms import poly_gcdex, poly_half_gcdex

from algorithms import poly_lcm, poly_gcd

from algorithms import poly_subresultants
from algorithms import poly_resultant

from algorithms import poly_groebner

from algorithms import poly_sqf

from algorithms import poly_decompose
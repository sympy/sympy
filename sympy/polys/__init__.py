
from monomial import monomials, monomial_count
from polynomial import Poly, PolynomialError

from algorithms import poly_div, poly_pdiv, poly_groebner, poly_lcm, poly_gcd, \
    poly_half_gcdex, poly_gcdex, poly_sqf, poly_resultant, poly_subresultants, \
    poly_decompose

from rootfinding import RootOf, RootsOf, RootSum, roots, poly_factors, poly_sturm

def _map_basic(f, *args, **kwargs):
    return tuple(g.as_basic() for g in f(*args, **kwargs))

def poly_quo(f, g, *symbols):
    """Wrapper for poly_div() """
    return poly_div(f, g, *symbols)[0]

def poly_rem(f, g, *symbols):
    """Wrapper for poly_div() """
    return poly_div(f, g, *symbols)[1]

def poly_pquo(f, g, *symbols):
    """Wrapper for poly_pdiv() """
    return poly_pdiv(f, g, *symbols)[0]

def poly_prem(f, g, *symbols):
    """Wrapper for poly_pdiv() """
    return poly_pdiv(f, g, *symbols)[1]

def div(*args, **kwargs):
    """Wrapper for poly_div() """
    q, r = poly_div(*args, **kwargs)

    if type(q) is not tuple:
        q = q.as_basic()
    else:
        q = tuple(p.as_basic() for p in q)

    return q, r.as_basic()

def quo(*args, **kwargs):
    """Wrapper for poly_quo() """
    q = poly_div(*args, **kwargs)[0]

    if type(q) is not tuple:
        return q.as_basic()
    else:
        return tuple(p.as_basic() for p in q)

def rem(*args, **kwargs):
    """Wrapper for poly_rem() """
    return poly_div(*args, **kwargs)[1].as_basic()

def pdiv(*args, **kwargs):
    """Wrapper for poly_pdiv() """
    return _map_basic(poly_pdiv, *args, **kwargs)

def pquo(*args, **kwargs):
    """Wrapper for poly_pquo() """
    return poly_pdiv(*args, **kwargs)[0].as_basic()

def prem(*args, **kwargs):
    """Wrapper for poly_prem() """
    return poly_pdiv(*args, **kwargs)[1].as_basic()

def groebner(*args, **kwargs):
    """Wrapper for poly_groebner() """
    return _map_basic(poly_groebner, *args, **kwargs)

def lcm(*args, **kwargs):
    """Wrapper for poly_lcm() """
    return _map_basic(poly_lcm, *args, **kwargs)

def gcd(*args, **kwargs):
    """Wrapper for poly_gcd() """
    return _map_basic(poly_gcd, *args, **kwargs)

def gcdex(*args, **kwargs):
    """Wrapper for poly_gcdex() """
    return _map_basic(poly_gcdex, *args, **kwargs)

def half_gcdex(*args, **kwargs):
    """Wrapper for poly_half_gcdex() """
    return _map_basic(poly_half_gcdex, *args, **kwargs)

def subresultants(*args, **kwargs):
    """Wrapper for poly_subresultants() """
    return _map_basic(poly_subresultants, *args, **kwargs)

def sqf(*args, **kwargs):
    """Wrapper for poly_sqf() """
    return _map_basic(poly_sqf, *args, **kwargs)

def decompose(*args, **kwargs):
    """Wrapper for poly_decompose() """
    return _map_basic(poly_decompose, *args, **kwargs)

def factors(*args, **kwargs):
    """Wrapper for poly_factors() """
    return _map_basic(poly_factors, *args, **kwargs)

def sturm(*args, **kwargs):
    """Wrapper for poly_sturm() """
    return _map_basic(poly_sturm, *args, **kwargs)

resultant = poly_resultant

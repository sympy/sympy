"""

"""

from sympy.polynomials import roots

def filter_roots(poly, n, predicate):
    return [ r for r in set(roots(poly, n)) if predicate(r) ]

def nni_roots(poly, n):
    return filter_roots(poly, n, lambda r: r.is_integer and r.is_nonnegative)
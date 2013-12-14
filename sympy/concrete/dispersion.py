from __future__ import print_function, division

from sympy.polys.polytools import factor_list


def dispersionset(p, q):
    r"""
    """
    # Check for valid input
    if not p.is_univariate or not q.is_univariate:
        raise ValueError("Polynomials need to be univariate")
    # We define the dispersion of constant polynomials to be zero
    if p.degree() < 1 or q.degree() < 1:
        return set([0])
    # The generator
    x = p.gens[0]
    # Factor p and q over the rationals
    fp = factor_list(p)
    fq = factor_list(q)
    J = set([])
    # Iterate over all pairs of factors
    for s, unused in fp[1]:
        for t, unused in fq[1]:
            m = s.degree()
            n = t.degree()
            if not m == n:
                D = []
            else:
                a = s.coeff_monomial(x**n)
                b = s.coeff_monomial(x**(n-1))
                c = t.coeff_monomial(x**n)
                d = t.coeff_monomial(x**(n-1))
                j = (b*c-a*d) / (a*c*n)
                if j < 0:
                    D = []
                else:
                    if (c*s - a*t.shift(j)).is_zero:
                        D = [j]
                    else:
                        D = []
            J = J.union(D)
    return J


def dispersion(p, q=None):
    r"""
    """
    if q is None:
        q = p

    J = dispersionset(p, q)
    return max(J)

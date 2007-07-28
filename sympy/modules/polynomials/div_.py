"""Polynomial division algorithms for use with class Polynomial"""

from sympy.modules.polynomials.base import *

def mv(f, g):
    """Devides f by the polynomials in g, returns the quotients and the remainder.

    g can be a single or a list of polynomials, assumed to be objects
    of the class Polynomials with matching var, orders and
    coeffs.
    Behaviour might be a bit surprising, as this algorithm doesn't
    stop when the leading terms don't divide any more, but instead
    tries to reduce the polynomial, by eliminating smaller terms.
    
    """
    if not isinstance(g, list):
        g = [g]
    r = Polynomial(S.Zero, var=f.var, order=f.order)
    q = []
    for i in range(0,len(g)):
        q.append(r)

    while f.sympy_expr is not S.Zero:
        for g_i in g:
            if g_i.sympy_expr is S.Zero: # avoid division by zero with bad arguments
                continue
            # Check if leading term of f is divisible by that of g_i.
            # TODO: Don't repeat term_div!
            if term_is_mult(f.coeffs[0], g_i.coeffs[0]):
                quot = Polynomial(coeffs=(term_div(f.coeffs[0], g_i.coeffs[0]),),
                                  var=f.var, order=f.order)
                q[g.index(g_i)] += quot
                f -= quot*g_i
                break
        else: # No division occured, add the leading term to remainder.
            lt = f.leading_term()
            r += lt
            f -= lt
    return q, r

def mv_int(f, g):
    """Devides f by the polynomials in g, returns the quotients and the remainder.

    g can be a single or a list of polynomials, assumed to be objects
    of the class Polynomials with matching var, orders and
    coeffs. Coeffs are assumed integers, in addition, the
    result are polynomials with integer coeffs.
    Behaviour might be a bit surprising, as this algorithm doesn't
    stop when the leading terms don't divide any more, but instead
    tries to reduce the polynomial, by eliminating smaller terms.

    """
    if not isinstance(g, list):
        g = [g]
    r = Polynomial(S.Zero, var=f.var, order=f.order)
    q = []
    for i in range(0,len(g)):
        q.append(r)

    while f.sympy_expr is not S.Zero:    # f != 0
        for g_i in g:
            if g_i.sympy_expr is S.Zero: # avoid division by zero with bad arguments
                continue
            # Check if leading term of f is divisible by that of g_i.
            # TODO: Don't repeat term_div!
            if term_is_mult(f.coeffs[0], g_i.coeffs[0]) \
                   and (f.coeffs[0][0]/g_i.coeffs[0][0]).is_integer:
                quot = Polynomial(coeffs=(term_div(f.coeffs[0],g_i.coeffs[0]),),
                                  var=f.var, order=f.order)
                q[g.index(g_i)] += quot
                f -= quot*g_i
                break
        else: # No division occured, add the leading term to remainder.
            lt = f.leading_term()
            r += lt
            f -= lt
    return q, r

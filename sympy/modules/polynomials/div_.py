"""Polynomial division algorithms for use with class Polynomial"""

from sympy.modules.polynomials.base import *

def mv(f, g):
    """Devides f by the polynomials in g, returns the quotients and the remainder.

    g can be a single or a list of polynomials, assumed to be objects
    of the class Polynomials with matching variables, orders and
    coefficients.
    Behaviour might be a bit surprising, as this algorithm doesn't
    stop when the leading terms don't divide any more, but instead
    tries to reduce the polynomial, by eliminating smaller terms.
    
    """
    if not isinstance(g, list):
        g = [g]
    r = Polynomial([[S.Zero]+[0]*len(f.var)], f.var, f.order, f.coeff)
    q = []
    for i in range(0,len(g)):
        q.append(r.copy())

    while f.cl[0][0] != 0:    # f != 0
        for g_i in g:
            if g_i.cl[0][0] == 0: # avoid division by zero with bad arguments
                continue
            # Check if leading term of f is divisible by that of g_i.
            # TODO: Don't repeat term_div!
            if term_is_mult(f.cl[0], g_i.cl[0]):
                quot = Polynomial([term_div(f.cl[0], g_i.cl[0])],
                                  f.var, f.order, f.coeff)
                q[g.index(g_i)] += quot
                f -= quot*g_i
                break
        else: # No division occured, add the leading term to remainder.
            # TODO: Don't create Polynomial, act on lists directly. 
            lt = Polynomial([f.cl[0]], f.var, f.order, f.coeff)
            r += lt # Append to the end.
            f -= lt # Remove from the beginning.
    return q, r

def mv_int(f, g):
    """Devides f by the polynomials in g, returns the quotients and the remainder.

    g can be a single or a list of polynomials, assumed to be objects
    of the class Polynomials with matching variables, orders and
    coefficients. Coefficients are assumed integers, in addition, the
    result are polynomials with integer coefficients.
    Behaviour might be a bit surprising, as this algorithm doesn't
    stop when the leading terms don't divide any more, but instead
    tries to reduce the polynomial, by eliminating smaller terms.

    """
    if not isinstance(g, list):
        g = [g]
    r = Polynomial([[S.Zero]+[0]*len(f.var)], f.var, f.order, f.coeff)
    q = []
    for i in range(0,len(g)):
        q.append(r.copy())

    while f.cl[0][0] != 0:    # f != 0
        for g_i in g:
            if g_i.cl[0][0] == 0: # avoid division by zero with bad arguments
                continue
            # Check if leading term of f is divisible by that of g_i.
            # TODO: Don't repeat term_div!
            if term_is_mult(f.cl[0], g_i.cl[0]) \
               and (f.cl[0][0]/g_i.cl[0][0]).is_integer:
                quot = Polynomial([term_div(f.cl[0], g_i.cl[0])],
                                  f.var, f.order, f.coeff)
                q[g.index(g_i)] += quot
                f -= quot*g_i
                break
        else: # No division occured, add the leading term to remainder.
            # TODO: Don't create Polynomial, act on lists directly. 
            lt = Polynomial([f.cl[0]], f.var, f.order, f.coeff)
            r += lt # Append to the end.
            f -= lt # Remove from the beginning.
    return q, r

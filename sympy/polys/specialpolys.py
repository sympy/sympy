"""Methods for creating interesting polynomials, e.g. for benchmarking """

from sympy import Add, Integer, Symbol

from sympy.polys.integerpolys import (
    zzX_zero, zzX_const, zzX_lift, zzx_from_dict,
    zzX_add_term, zzX_mul, zzX_sqr, zzX_neg)

def fateman_poly_F_1(n):
    """Fateman's GCD benchmark: trivial GCD """
    Y = [ Symbol('y_' + str(i)) for i in xrange(0, n+1) ]

    y_0, y_1 = Y[0], Y[1]

    u = y_0    + Add(*[ y    for y in Y[1:] ])
    v = y_0**2 + Add(*[ y**2 for y in Y[1:] ])

    F = ((u + 1)*(u + 2)).as_poly(*Y)
    G = ((v + 1)*(-3*y_1*y_0**2 + y_1**2 - 1)).as_poly(*Y)

    H = (Integer(1)).as_poly(*Y)

    return F, G, H

def zzX_fateman_poly_F_1(n):
    """Fateman's GCD benchmark: trivial GCD """
    u = [1, 0]

    for i in xrange(1, n+1):
        u = [zzX_const(i, 1), u]

    v = [1, 0, 0]

    for i in xrange(1, n+1):
        v = [zzX_const(i, 1), zzX_zero(i), v]

    U = zzX_add_term(u, zzX_const(n, 1))
    V = zzX_add_term(u, zzX_const(n, 2))

    W = zzX_add_term(v, zzX_const(n, 1))
    Y = zzX_lift(n-1, [[-3, 0], [], [1, 0, -1]])

    F = zzX_mul(U, V)
    G = zzX_mul(W, Y)

    H = zzX_const(n+1, 1)

    return F, G, H

def fateman_poly_F_2(n):
    """Fateman's GCD benchmark: linearly dense quartic inputs """
    Y = [ Symbol('y_' + str(i)) for i in xrange(0, n+1) ]

    y_0 = Y[0]

    u = Add(*[ y for y in Y[1:] ])

    H = ((y_0 + u + 1)**2).as_poly(*Y)

    F = ((y_0 - u - 2)**2).as_poly(*Y)
    G = ((y_0 + u + 2)**2).as_poly(*Y)

    return H*F, H*G, H

def zzX_fateman_poly_F_2(n):
    """Fateman's GCD benchmark: linearly dense quartic inputs """
    u = [1, 0]

    for i in xrange(1, n):
        u = [zzX_const(i, 1), u]

    v = zzX_add_term(u, zzX_const(n-1, 2))

    f = zzX_sqr([zzX_const(n, 1), zzX_neg(v)])
    g = zzX_sqr([zzX_const(n, 1), v])

    v = zzX_add_term(u, zzX_const(n-1, 1))

    h = zzX_sqr([zzX_const(n, 1), v])

    return zzX_mul(f,h), zzX_mul(g,h), h

def fateman_poly_F_3(n):
    """Fateman's GCD benchmark: sparse inputs (deg f ~ vars f) """
    Y = [ Symbol('y_' + str(i)) for i in xrange(0, n+1) ]

    y_0 = Y[0]

    u = Add(*[ y**(n+1) for y in Y[1:] ])

    H = ((y_0**(n+1) + u + 1)**2).as_poly(*Y)

    F = ((y_0**(n+1) - u - 2)**2).as_poly(*Y)
    G = ((y_0**(n+1) + u + 2)**2).as_poly(*Y)

    return H*F, H*G, H

def zzX_fateman_poly_F_3(n):
    """Fateman's GCD benchmark: sparse inputs (deg f ~ vars f) """
    u = zzx_from_dict({n+1:1})

    for i in xrange(1, n):
        u = zzX_add_term([u], zzX_const(i, 1), n+1)

    v = zzX_add_term(u, zzX_const(n-1, 2))

    f = zzX_sqr(zzX_add_term([zzX_neg(v)], zzX_const(n, 1), n+1))
    g = zzX_sqr(zzX_add_term([v], zzX_const(n, 1), n+1))

    v = zzX_add_term(u, zzX_const(n-1, 1))

    h = zzX_sqr(zzX_add_term([v], zzX_const(n, 1), n+1))

    return zzX_mul(f,h), zzX_mul(g,h), h

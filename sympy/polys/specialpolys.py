"""Functions for generating interesting polynomials, e.g. for benchmarking. """

from sympy.core import S, Add, Mul, Symbol, Rational, sympify, symbols

from sympy.polys.polytools import Poly
from sympy.polys.polyutils import _analyze_gens

from sympy.polys.polyclasses import DMP

from sympy.polys.densebasic import (
    dmp_zero, dmp_one, dmp_ground, dmp_normal,
    dup_from_raw_dict, dmp_raise, dup_random
)

from sympy.polys.densearith import (
    dmp_add_term, dmp_neg, dmp_mul, dmp_sqr
)

from sympy.polys.factortools import (
    dup_zz_cyclotomic_poly
)

from sympy.polys.domains import ZZ

from sympy.ntheory import nextprime

from sympy.utilities import cythonized, subsets

@cythonized("n,i")
def swinnerton_dyer_poly(n, x=None, **args):
    """Generates n-th Swinnerton-Dyer polynomial in `x`.  """
    if n <= 0:
        raise ValueError("can't generate Swinnerton-Dyer polynomial of order %s" % n)

    if x is not None:
        x = sympify(x)
    else:
        x = Symbol('x', dummy=True)

    p, elts = 2, [[x, -2**Rational(1,2)],
                  [x,  2**Rational(1,2)]]

    for i in xrange(2, n+1):
        p, _elts = nextprime(p), []

        neg_sqrt = -p**Rational(1,2)
        pos_sqrt = +p**Rational(1,2)

        for elt in elts:
            _elts.append(elt + [neg_sqrt])
            _elts.append(elt + [pos_sqrt])

        elts = _elts

    poly = []

    for elt in elts:
        poly.append(Add(*elt))

    if not args.get('polys', False):
        return Mul(*poly).expand()
    else:
        return Poly(Mul(*poly))

def cyclotomic_poly(n, x=None, **args):
    """Generates cyclotomic polynomial of order `n` in `x`. """
    if n <= 0:
        raise ValueError("can't generate cyclotomic polynomial of order %s" % n)

    if x is not None:
        x = sympify(x)
    else:
        x = Symbol('x', dummy=True)

    poly = Poly.new(DMP(dup_zz_cyclotomic_poly(int(n), ZZ), ZZ), x)

    if not args.get('polys', False):
        return poly.as_basic()
    else:
        return poly

def symmetric_poly(n, *gens, **args):
    """Generates symmetric polynomial of order `n`. """
    gens = _analyze_gens(gens)

    if n < 0 or n > len(gens) or not gens:
        raise ValueError("can't generate symmetric polynomial of order %s for %s" % (n, gens))
    elif not n:
        poly = S.One
    else:
        poly = Add(*[ Mul(*s) for s in subsets(gens, int(n)) ])

    if not args.get('polys', False):
        return poly
    else:
        return Poly(poly, *gens)

def random_poly(x, n, inf, sup, domain=ZZ, polys=False):
    """Return a polynomial of degree ``n`` with coefficients in ``[inf, sup]``. """
    poly = Poly(dup_random(n, inf, sup, domain), x, domain=domain)

    if not polys:
        return poly.as_basic()
    else:
        return poly

@cythonized("n,i,j")
def interpolating_poly(n, x, X='x', Y='y'):
    """Construct Lagrange interpolating polynomial for ``n`` data points. """
    if isinstance(X, str):
        X = symbols("%s:%s" % (xn, n))

    if isinstance(Y, str):
        Y = symbols("%s:%s" % (yn, n))

    coeffs = []

    for i in xrange(0, n):
        numer = []
        denom = []

        for j in xrange(0, n):
            if i == j:
                continue

            numer.append(x    - X[j])
            denom.append(X[i] - X[j])

        numer = Mul(*numer)
        denom = Mul(*denom)

        coeffs.append(numer/denom)

    return Add(*[ coeff*y for coeff, y in zip(coeffs, Y) ])

@cythonized("n,i")
def fateman_poly_F_1(n):
    """Fateman's GCD benchmark: trivial GCD """
    Y = [ Symbol('y_' + str(i)) for i in xrange(0, n+1) ]

    y_0, y_1 = Y[0], Y[1]

    u = y_0    + Add(*[ y    for y in Y[1:] ])
    v = y_0**2 + Add(*[ y**2 for y in Y[1:] ])

    F = ((u + 1)*(u + 2)).as_poly(*Y)
    G = ((v + 1)*(-3*y_1*y_0**2 + y_1**2 - 1)).as_poly(*Y)

    H = Poly(1, *Y)

    return F, G, H

@cythonized("n,m,i")
def dmp_fateman_poly_F_1(n, K):
    """Fateman's GCD benchmark: trivial GCD """
    u = [K(1), K(0)]

    for i in xrange(0, n):
        u = [dmp_one(i, K), u]

    v = [K(1), K(0), K(0)]

    for i in xrange(0, n):
        v = [dmp_one(i, K), dmp_zero(i), v]

    m = n-1

    U = dmp_add_term(u, dmp_ground(K(1), m), 0, n, K)
    V = dmp_add_term(u, dmp_ground(K(2), m), 0, n, K)

    f = [[-K(3), K(0)], [], [K(1), K(0), -K(1)]]

    W = dmp_add_term(v, dmp_ground(K(1), m), 0, n, K)
    Y = dmp_raise(f, m, 1, K)

    F = dmp_mul(U, V, n, K)
    G = dmp_mul(W, Y, n, K)

    H = dmp_one(n, K)

    return F, G, H

@cythonized("n,i")
def fateman_poly_F_2(n):
    """Fateman's GCD benchmark: linearly dense quartic inputs """
    Y = [ Symbol('y_' + str(i)) for i in xrange(0, n+1) ]

    y_0 = Y[0]

    u = Add(*[ y for y in Y[1:] ])

    H = Poly((y_0 + u + 1)**2, *Y)

    F = Poly((y_0 - u - 2)**2, *Y)
    G = Poly((y_0 + u + 2)**2, *Y)

    return H*F, H*G, H

@cythonized("n,m,i")
def dmp_fateman_poly_F_2(n, K):
    """Fateman's GCD benchmark: linearly dense quartic inputs """
    u = [K(1), K(0)]

    for i in xrange(0, n-1):
        u = [dmp_one(i, K), u]

    m = n-1

    v = dmp_add_term(u, dmp_ground(K(2), m-1), 0, n, K)

    f = dmp_sqr([dmp_one(m, K), dmp_neg(v, m, K)], n, K)
    g = dmp_sqr([dmp_one(m, K), v], n, K)

    v = dmp_add_term(u, dmp_one(m-1, K), 0, n, K)

    h = dmp_sqr([dmp_one(m, K), v], n, K)

    return dmp_mul(f, h, n, K), dmp_mul(g, h, n, K), h

@cythonized("n,i")
def fateman_poly_F_3(n):
    """Fateman's GCD benchmark: sparse inputs (deg f ~ vars f) """
    Y = [ Symbol('y_' + str(i)) for i in xrange(0, n+1) ]

    y_0 = Y[0]

    u = Add(*[ y**(n+1) for y in Y[1:] ])

    H = Poly((y_0**(n+1) + u + 1)**2, *Y)

    F = Poly((y_0**(n+1) - u - 2)**2, *Y)
    G = Poly((y_0**(n+1) + u + 2)**2, *Y)

    return H*F, H*G, H

@cythonized("n,i")
def dmp_fateman_poly_F_3(n, K):
    """Fateman's GCD benchmark: sparse inputs (deg f ~ vars f) """
    u = dup_from_raw_dict({n+1: K.one}, K)

    for i in xrange(0, n-1):
        u = dmp_add_term([u], dmp_one(i, K), n+1, i+1, K)

    v = dmp_add_term(u, dmp_ground(K(2), n-2), 0, n, K)

    f = dmp_sqr(dmp_add_term([dmp_neg(v, n-1, K)], dmp_one(n-1, K), n+1, n, K), n, K)
    g = dmp_sqr(dmp_add_term([v], dmp_one(n-1, K), n+1, n, K), n, K)

    v = dmp_add_term(u, dmp_one(n-2, K), 0, n-1, K)

    h = dmp_sqr(dmp_add_term([v], dmp_one(n-1, K), n+1, n, K), n, K)

    return dmp_mul(f, h, n, K), dmp_mul(g, h, n, K), h

# A few useful polynomials from Wang's paper ('78).

f_0 = dmp_normal([
    [[1,2,3], [2]],
    [[3]],
    [[4,5,6], [1,2,1], [1]]
], 2, ZZ)

f_1 = dmp_normal([
    [[1, 0], []],
    [[1, 0, 1], [20, 30], [1, 10, 0]],
    [[1, 0], [30, 20], [1, 10, 1, 610], [20, 230, 300]],
    [[1, 10, 0], [30, 320, 200], [600, 6000]]
], 2, ZZ)

f_2 = dmp_normal([
    [[1], [1, 0], [1, 0, 0], [1, 0, 0, 0]],
    [[]],
    [[1], [1, 90], [90, 0]],
    [[1, -11], [], [1, -11, 0, 0]],
    [[]],
    [[1, -11], [90, -990]]
], 2, ZZ)

f_3 = dmp_normal([
    [[1], [], []],
    [[1, 0, 0, 0, 1]],
    [[1, 0], [], [], [1, 0]],
    [[1], [1, 0, 0, 0], [], [1, 0, 0, 0, 1, 0], []],
    [[1, 0, 0, 0, 1], [1, 0, 0, 0, 1, 1, 0, 0], []],
    [[1, 0], [1, 0, 0, 0, 0], []]
], 2, ZZ)

f_4 = dmp_normal([
    [[-1, 0], [], [], [], [], [], [], [], []],
    [[-1, 0, 0, 0], [], [], [], [], []],
    [[-1, 0, 0], [], [], [], [-5], [], [], [], [], [], [], [], []],
    [[-1, 0, 0, 0, 0], [], [1, 0, 3, 0], [], [-5, 0, 0], [-1, 0, 0, 0], [], [], [], []],
    [[1, 0, 3, 0, 0, 0], [], [], [-1, 0, 0, 0, 0, 0], []],
    [[1, 0, 3, 0, 0], [], [], [-1, 0, 0, 0, 0], [5, 0, 15], [], [], [-5, 0, 0], [], [], [], []],
    [[1, 0, 3, 0, 0, 0, 0], [], [], [-1, 0, 0, 0, 0, 0, 0], [5, 0, 15, 0, 0], [1, 0, 3, 0, 0, 0], [], [-5, 0, 0, 0, 0], []],
    [[1, 0, 3, 0, 0, 0, 0, 0]],
    [[1, 0, 3, 0, 0, 0, 0], [], [], [], [5, 0, 15, 0, 0], [], [], []],
    [[1, 0, 3, 0, 0, 0, 0, 0, 0], [], [], [], [5, 0, 15, 0, 0, 0, 0]]
], 2, ZZ)

f_5 = dmp_normal([
    [[-1]],
    [[-3], [3, 0]],
    [[-3], [6, 0], [-3, 0, 0]],
    [[-1], [3, 0], [-3, 0, 0], [1, 0, 0, 0]]
], 2, ZZ)

f_6 = dmp_normal([
    [[[2115]], [[]]],
    [[[45, 0, 0], [], [], [-45, 0, 0]]],
    [[[]]],
    [[[-423]], [[-47]], [[]], [[141], [], [94, 0], []], [[]]],
    [[[-9, 0, 0], [], [], [9, 0, 0]],
     [[-1, 0, 0], [], [], [1, 0, 0]],
     [[]],
     [[3, 0, 0], [], [2, 0, 0, 0], [-3, 0, 0], [], [-2, 0, 0, 0], []]
    ]
], 3, ZZ)


w_1 = dmp_normal([
    [[4, 0, 0], [4, 0, 0, 0], [-4, 0, 0, 0, 0], [-4, 0, 0, 0, 0, 0], []],
    [[1, 0, 0, 0], [12, 0], [-1, 0, 0, 12, 0, 0], [-12, 0, 0, 0], [-12, 0, 0, 0, 0]],
    [[8], [6, 8, 0], [-4, 4, -8, 0, 0], [-4, -2, -8, 0, 0, 0], []],
    [[2, 0], [1, 0, 0, 0], [-1, 0, -2 , 0, 9, 0], [-12, 12, 0, 0], [-12, 3, 0, 0, 0]],
    [[6], [-6, 8, 0], [-2, -8, 2, 0, 0], []],
    [[2, 0], [-2, 0, 0, 0], [-3, 0], [3, 0, 0, 0]],
    [[-2], [2, 0, 0], []]
], 2, ZZ)

w_2 = dmp_normal([
    [24, 48, 0, 0],
    [24, 0, 0, -72, 0, 0],
    [25, 2, 0, 4, 8],
    [1, 0, 0, 1, 0, 0, -12],
    [1, -1, -2, 292, 0, 0],
    [-1, 0, 0, 3, 0, 0, 0],
    [-1, 0, 12, 0, 0, 48],
    [],
    [-12, 0, 0, 0]
], 1, ZZ)


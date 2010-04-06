from sympy import Function, symbols, S, sqrt, rf, factorial
from sympy.solvers.recurr import rsolve, rsolve_poly, rsolve_ratio, rsolve_hyper

y = Function('y')
n, k = symbols('nk', integer=True)
C0, C1, C2 = symbols('C0', 'C1', 'C2')

def test_rsolve_poly():
    assert rsolve_poly([-1, -1, 1], 0, n) == 0
    assert rsolve_poly([-1, -1, 1], 1, n) == -1

    assert rsolve_poly([-1, n+1], n, n) == 1
    assert rsolve_poly([-1, 1], n, n) == C0 + (n**2 - n)/2
    assert rsolve_poly([-n-1, n], 1, n) == C1*n - 1
    assert rsolve_poly([-4*n-2, 1], 4*n+1, n) == -1

    assert rsolve_poly([-1, 1], n**5 + n**3, n) == C0 - n**3 / 2 - n**5 / 2 + n**2 / 6 + n**6 / 6 + 2*n**4 / 3

def test_rsolve_ratio():
    solution = rsolve_ratio([-2*n**3+n**2+2*n-1, 2*n**3+n**2-6*n,
        -2*n**3-11*n**2-18*n-9, 2*n**3+13*n**2+22*n+8], 0, n)

    assert solution in [
            (-3*C1 + 2*C1*n)/(-2 + 2*n**2),
            ( 3*C1 - 2*C1*n)/( 2 - 2*n**2),
            (-3*C2 + 2*C2*n)/(-2 + 2*n**2),
            ( 3*C2 - 2*C2*n)/( 2 - 2*n**2),
                         ]

def test_rsolve_hyper():
    assert rsolve_hyper([-1, -1, 1], 0, n) == C0*(S.Half - S.Half*sqrt(5))**n + C1*(S.Half + S.Half*sqrt(5))**n

    assert rsolve_hyper([n**2-2, -2*n-1, 1], 0, n) in [C0*rf(sqrt(2), n) + C1*rf(-sqrt(2), n),
                                                       C1*rf(sqrt(2), n) + C0*rf(-sqrt(2), n)]

    assert rsolve_hyper([n**2-k, -2*n-1, 1], 0, n) in [C0*rf(sqrt(k), n) + C1*rf(-sqrt(k), n),
                                                       C1*rf(sqrt(k), n) + C0*rf(-sqrt(k), n)]

    assert rsolve_hyper([2*n*(n+1), -n**2-3*n+2, n-1], 0, n) == C0*factorial(n) + C1*2**n

    assert rsolve_hyper([n + 2, -(2*n + 3)*(17*n**2 + 51*n + 39), n + 1], 0, n) == 0

    assert rsolve_hyper([-n-1, -1, 1], 0, n) == 0

    assert rsolve_hyper([-1, 1], n, n).expand() == C0 + n**2/2 - n/2

    assert rsolve_hyper([-1, 1], 1+n, n).expand() == C0 + n**2/2 + n/2

    assert rsolve_hyper([-1, 1], 3*(n+n**2), n).expand() == C0 + n**3 - n

def recurrence_term(c, f):
    """Compute RHS of recurrence in f(n) with coefficients in c."""
    return sum(c[i]*f.subs(n, n+i) for i in range(len(c)))

def rsolve_bulk_checker(solver, c, q, p):
    """Used by test_rsolve_bulk."""
    pp = solver(c, q, n)
    assert pp == p

def test_rsolve_bulk():
    """Some bulk-generated tests."""
    funcs = [ n, n+1, n**2, n**3, n**4, n+n**2, 27*n + 52*n**2 - 3*n**3 + 12*n**4 - 52*n**5 ]
    coeffs = [ [-2, 1], [-2, -1, 1], [-1, 1, 1, -1, 1], [-n, 1], [n**2-n+12, 1] ]
    for p in funcs:
        # compute difference
        for c in coeffs:
            q = recurrence_term(c, p)
            if p.is_polynomial(n):
                yield rsolve_bulk_checker, rsolve_poly, c, q, p
            #if p.is_hypergeometric(n):
            #    yield rsolve_bulk_checker, rsolve_hyper, c, q, p

def test_rsolve():
    f = y(n+2) - y(n+1) - y(n)
    g = C0*(S.Half - S.Half*sqrt(5))**n \
      + C1*(S.Half + S.Half*sqrt(5))**n
    h = sqrt(5)*(S.Half + S.Half*sqrt(5))**n \
      - sqrt(5)*(S.Half - S.Half*sqrt(5))**n

    assert rsolve(f, y(n)) == g

    assert rsolve(f, y(n), [      0,      5 ]) == h
    assert rsolve(f, y(n), {   0 :0,   1 :5 }) == h
    assert rsolve(f, y(n), { y(0):0, y(1):5 }) == h

    f = (n-1)*y(n+2) - (n**2+3*n-2)*y(n+1) + 2*n*(n+1)*y(n)
    g = C0*factorial(n) + C1*2**n
    h = -3*factorial(n) + 3*2**n

    assert rsolve(f, y(n)) == g

    assert rsolve(f, y(n), [      0,      3 ]) == h
    assert rsolve(f, y(n), {   0 :0,   1 :3 }) == h
    assert rsolve(f, y(n), { y(0):0, y(1):3 }) == h

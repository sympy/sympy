from sympy import Poly, QQ, cyclotomic_poly, factorint
from sympy.core.mul import prod
from sympy.ntheory.residue_ntheory import n_order
from sympy.polys.numberfields.basis import round_two
from sympy.polys.numberfields.primes import (
    prime_decomp, _two_elt_rep, prime_valuation
)

def test_two_elt_rep():
    ell = 7
    T = Poly(cyclotomic_poly(ell))
    ZK, dK = round_two(T)
    for p in [29, 13, 11, 5]:
        P = prime_decomp(p, T)
        for Pi in P:
            # We have Pi in two-element representation, and, because we are
            # looking at a cyclotomic field, this was computed by the "easy"
            # method that just factors T mod p. We will now convert this to
            # a set of Z-generators, then convert that back into a two-element
            # rep. The latter need not be identical to the two-elt rep we
            # already have, but it must have the same HNF.
            H2 = p*ZK + Pi.alpha*ZK
            gens = H2.standard_reps()
            # Note: we could supply f = Pi.f, but prefer to test behavior without it.
            b = _two_elt_rep(gens, T, p, ZK=ZK)
            if b != Pi.alpha:
                H3 = p*ZK + b*ZK
                assert H3 == H2

def test_valuation_at_prime_ideal():
    p = 7
    T = Poly(cyclotomic_poly(p))
    ZK, dK = round_two(T)
    P = prime_decomp(p, T, dK=dK, ZK=ZK)
    assert len(P) == 1
    P0 = P[0]
    # Work through the function, passing a triple as second arg:
    v = prime_valuation(p * ZK, (p, [p, P0.alpha], ZK))
    assert v == P0.e
    # Work through the PrimeIdeal:
    v = P0.valuation(p*ZK)
    assert v == P0.e

def test_decomp_1():
    # All prime decompositions in cyclotomic fields are in the "easy case,"
    # since the index is unity.
    # Here we check the ramified prime.
    T = Poly(cyclotomic_poly(7))
    P = prime_decomp(7, T)
    assert len(P) == 1
    P0 = P[0]
    assert P0.e == 6
    assert P0.f == 1

def test_decomp_2():
    # More easy cyclotomic cases, but here we check unramified primes.
    ell = 7
    T = Poly(cyclotomic_poly(ell))
    for p in [29, 13, 11, 5]:
        f_exp = n_order(p, ell)
        g_exp = (ell - 1) // f_exp
        P = prime_decomp(p, T)
        assert len(P) == g_exp
        for Pi in P:
            assert Pi.e == 1
            assert Pi.f == f_exp

def test_decomp_3():
    from sympy.abc import x
    T = Poly(x ** 2 - 35)
    rad = {}
    ZK, dK = round_two(T, radicals=rad)
    # 35 is 3 mod 4, so field disc is 4*5*7, and theory says each of the
    # rational primes 2, 5, 7 should be the square of a prime ideal.
    for p in [2, 5, 7]:
        P = prime_decomp(p, T, dK=dK, ZK=ZK, radical=rad.get(p))
        assert len(P) == 1
        assert P[0].e == 2
        assert P[0]**2 == p*ZK

def test_decomp_4():
    from sympy.abc import x
    T = Poly(x ** 2 - 21)
    rad = {}
    ZK, dK = round_two(T, radicals=rad)
    # 21 is 1 mod 4, so field disc is 3*7, and theory says the
    # rational primes 3, 7 should be the square of a prime ideal.
    for p in [3, 7]:
        P = prime_decomp(p, T, dK=dK, ZK=ZK, radical=rad.get(p))
        assert len(P) == 1
        assert P[0].e == 2
        assert P[0]**2 == p*ZK

def test_decomp_5():
    # Here is our first test of the "hard case" of prime decomposition.
    # We work in a quadratic extension Q(sqrt(d)) where d is 1 mod 4, and
    # we consider the factorization of the rational prime 2, which divides
    # the index.
    # Theory says the form of p's factorization depends on the residue of
    # d mod 8, so we consider both cases.
    from sympy.abc import x
    for d in [-7, -3]:
        T = Poly(x ** 2 - d)
        rad = {}
        ZK, dK = round_two(T, radicals=rad)
        p = 2
        P = prime_decomp(p, T, dK=dK, ZK=ZK, radical=rad.get(p))
        if d % 8 == 1:
            assert len(P) == 2
            assert all(P[i].e == 1 and P[i].f == 1 for i in range(2))
            assert prod(Pi**Pi.e for Pi in P) == p * ZK
        else:
            assert d % 8 == 5
            assert len(P) == 1
            assert P[0].e == 1
            assert P[0].f == 2
            assert P[0] == p * ZK

def test_decomp_6():
    # Another case where 2 divides the index. This is Dedekind's example of
    # an essential discriminant divisor. (See Cohen, Excercise 6.10.)
    from sympy.abc import x
    T = Poly(x ** 3 + x ** 2 - 2 * x + 8)
    rad = {}
    ZK, dK = round_two(T, radicals=rad)
    p = 2
    P = prime_decomp(p, T, dK=dK, ZK=ZK, radical=rad.get(p))
    assert len(P) == 3
    assert all(Pi.e == Pi.f == 1 for Pi in P)
    assert prod(Pi**Pi.e for Pi in P) == p*ZK

def test_decomp_7():
    # Try working through an AlgebraicField
    from sympy.abc import x, theta
    T = Poly(x ** 3 + x ** 2 - 2 * x + 8)
    K = QQ.algebraic_field((T, theta))
    p = 2
    P = K.primes_above(p)
    ZK = K.int_ring
    assert len(P) == 3
    assert all(Pi.e == Pi.f == 1 for Pi in P)
    assert prod(Pi**Pi.e for Pi in P) == p*ZK

def test_decomp_8(inspect=False):
    # This time we consider various cubics, and try factoring all primes
    # dividing the index.
    from sympy.abc import x, theta
    cases = (
        x ** 3 + 3 * x ** 2 - 4 * x + 4,
        x ** 3 + 3 * x ** 2 + 3 * x - 3,
        x ** 3 + 5 * x ** 2 - x + 3,
        x ** 3 + 5 * x ** 2 - 5 * x - 5,
        x ** 3 + 3 * x ** 2 + 5,
        x ** 3 + 6 * x ** 2 + 3 * x - 1,
        x ** 3 + 6 * x ** 2 + 4,
        x ** 3 + 7 * x ** 2 + 7 * x - 7,
        x ** 3 + 7 * x ** 2 - x + 5,
        x ** 3 + 7 * x ** 2 - 5 * x + 5,
        x ** 3 + 4 * x ** 2 - 3 * x + 7,
        x ** 3 + 8 * x ** 2 + 5 * x - 1,
        x ** 3 + 8 * x ** 2 - 2 * x + 6,
        x ** 3 + 6 * x ** 2 - 3 * x + 8,
        x ** 3 + 9 * x ** 2 + 6 * x - 8,
        x ** 3 + 15 * x ** 2 - 9 * x + 13,
    )
    def display(T, p, radical, P, I, J):
        """Useful for inspection, when running test manually."""
        print('=' * 20)
        print(T, p, radical)
        for Pi in P:
            print(f'  ({Pi.pretty()})')
        print("I: ", I)
        print("J: ", J)
        print(f'Equal: {I == J}')
    for g in cases:
        T = Poly(g)
        rad = {}
        ZK, dK = round_two(T, radicals=rad)
        dT = T.discriminant()
        f_squared = dT // dK
        F = factorint(f_squared)
        for p in F:
            radical = rad.get(p)
            P = prime_decomp(p, T, dK=dK, ZK=ZK, radical=radical)
            I = prod(Pi**Pi.e for Pi in P)
            J = p * ZK
            if inspect:
                display(T, p, radical, P, I, J)
            else:
                assert I == J

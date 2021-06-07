from sympy import (S, Symbol, Rational, Poly, Eq, ratsimp, checkodesol, sqrt,
    oo, I)
from random import randint
from sympy.solvers.ode.riccati import (riccati_normal, riccati_inverse_normal,
    inverse_transform_poly, find_poles, limit_at_inf, check_necessary_conds,
    val_at_inf)

x = Symbol('x')

# These are the functions used to generate the tests
# SHOULD NOT BE USED DIRECTLY IN TESTS

def rand_rational(maxint):
    return Rational(randint(-maxint, maxint), randint(1, maxint))

def rand_poly(x, degree, maxint):
    return Poly([rand_rational(maxint) for _ in range(degree+1)], x)

def rand_rational_function(x, degree, maxint):
    degnum = randint(1, degree)
    degden = randint(1, degree)
    num = rand_poly(x, degnum, maxint)
    den = rand_poly(x, degden, maxint)
    while den == Poly(0, x):
        den = rand_poly(x, degden, maxint)
    return num / den

def find_riccati_ode(ratfunc, x, yf):
    y = ratfunc
    yp = y.diff(x)
    q1 = rand_rational_function(x, 1, 3)
    q2 = rand_rational_function(x, 1, 3)
    while q2 == 0:
        q2 = rand_rational_function(x, 1, 3)
    q0 = ratsimp(yp - q1*y - q2*y**2)
    eq = Eq(yf.diff(), q0 + q1*yf + q2*yf**2)
    sol = Eq(yf, y)
    print(eq)
    print(sol)
    assert checkodesol(eq, sol) == (True, 0)
    return eq, q0, q1, q2


def test_riccati_transformation():
    """
    Test the transformation that converts a
    Riccati equation into its normal form.

    Code
    ====

    This is the code to generate random test
    cases and test the transformations.

    for i in range(10):
        w, b1, b2 = [rand_rational_function(x, 10, 5).cancel() for i in range(3)]
        assert w == riccati_inverse_normal(riccati_normal(w, x, b1, b2), x, b1, b2).cancel()
        print("    (")
        print("        ", w, end=",\n")
        print("        ", b1, end=",\n")
        print("        ", b2, end=",\n")
        print("    ),")
    """
    test_cases = [
    (
        (-50*x**4 - 12*x**3 - 30*x**2 - 30*x + 24)/(12*x**4 + 45*x**3 + 30*x**2 - 6*x),
        (-50*x**3 + 40*x**2 - 5*x - 60)/(30*x**3 - 12*x**2 + 25*x + 80),
        (-300*x**9 - 60*x**8 + 120*x**7 + 48*x**6 + 40*x**5 + 300*x**4 - 150*x**3 + 300*x**2 + 30*x \
        - 36)/(45*x**9 - 60*x**8 + 30*x**7 + 40*x**6 + 60*x**5 + 60*x**4 - 90*x**3 - 12*x - 100),
    ),
    (
        (-6*x**5 - 40*x**4 + 30*x**3 + 6*x**2 - 30*x - 12)/(15*x**5 + 30*x**4 + 15*x**3 - 50*x**2 + \
        6*x - 60),
        (18*x**3 + 16*x**2 - 12*x - 24)/(12*x**3 + 24*x**2 + 15*x - 24),
        (-6*x**4 - 6*x**3 + 15*x**2 + 12*x)/(6*x**4 + 15*x**3 - 6*x**2 - 6*x + 2),
    ),
    (
        (20*x**3 - 6*x**2 - 3)/(20*x**3 + 12*x**2 - 12*x + 12),
        (36*x**8 + 150*x**7 + 120*x**6 - 180*x**5 - 45*x**4 - 48*x**3 + 30*x**2 - 24*x + 60)/(240*x**8 \
        - 75*x**7 - 36*x**6 + 48*x**5 + 60*x**4 - 100*x**3 - 240*x**2 - 240*x - 12),
        (12*x**6 - 80*x**5 + 15*x**4 + 24*x**3 + 60*x**2 - 180)/(12*x**6 - 75*x**5 + 120*x**4 - 30*x**3 \
        - 24*x**2 + 20),
    ),
    (
        (60*x**5 + 60*x**4 - 48*x**3 - 60*x**2 + 20*x + 30)/(80*x**5 - 40*x**4 - 75*x**3 - 75*x**2 - \
        30*x + 100),
        (12*x**5 - 2*x**4 - 2*x - 1)/(15*x**5 + x**4 - 3*x**3 + 5*x**2 - 6*x - 4),
        (-30*x**6 - 80*x**5 + 48*x**4 - 90*x**3 - 60*x + 24)/(240*x**6 + 48*x**5 + 30*x**4 - 48*x**3 \
        - 40*x**2 + 180*x - 45),
    ),
    (
        (-180*x**10 - 120*x**9 + 36*x**8 + 120*x**7 + 60*x**6 - 120*x**4 + 30*x**3 - 90*x**2 + 60)/(\
        100*x**9 + 15*x**7 - 80*x**6 + 48*x**5 + 240*x**4 - 150*x**3 + 60*x + 60),
        (20*x + 4)/(3*x + 6),
        (36*x**4 + 20*x**3 - 30*x**2 + 60*x)/(45*x**4 + 60*x**3 + 48*x**2 - 60),
    ),
    (
        (4*x**3 + 4*x**2 + 16*x + 15)/(50*x**2 + 25*x + 12),
        (-100*x**8 + 36*x**7 + 120*x**5 + 80*x**4 - 30*x**3 + 40*x**2 - 20*x - 40)/(180*x**9 + 100*x**7 \
        - 15*x**6 + 100*x**5 + 90*x**4 + 40*x**3 + 60*x**2 + 60*x + 120),
        (75 - 50*x)/(20*x - 6),
    ),
    (
        (40*x**5 - 180*x**3 - 60*x**2 - 60*x + 300)/(180*x**6 + 150*x**5 - 100*x**4 - 48*x**3 + 45*x**2 \
        - 20*x + 300),
        (-150*x - 12)/(300*x**2 + 15*x - 20),
        (-20*x**10 - 75*x**9 + 120*x**8 + 120*x**7 - 60*x**6 + 120*x**5 - 40*x**4 - 48*x**3 - 80*x**2 + \
        60*x - 30)/(120*x**10 - 30*x**9 - 120*x**8 - 120*x**7 + 150*x**6 - 60*x**5 - 15*x**4 - 150*x**3 -\
        150*x**2 + 120*x + 45),
    ),
    (
        (4*x - 5)/(12*x + 4),
        (8*x**3 + 30*x**2 - 3*x + 60)/(15*x**3 - 12*x**2 + 60*x - 6),
        (25*x**6 - 20*x**5 + 10*x**4 + 15*x**3 + 16*x**2 + 12*x)/(40*x**6 + 80*x**5 - 40*x**4 - 10*x**3 - \
        25*x**2 - 30*x - 10),
    ),
    (
        (4*x**8 + 25*x**7 + 12*x**5 - 10*x**4 + 20*x**3 - 80*x**2 - 8*x)/(20*x**8 - 40*x**6 - 20*x**5 + \
        25*x**4 + 5*x**3 - 40*x**2 - 60*x - 5),
        (60*x**3 - 90*x**2 + 12)/(30*x**3 + 10*x**2 - 15*x),
        (16*x + 50)/(25*x + 25),
    ),
    (
        (-20*x**6 - 30*x**5 - 15*x**4 + 60*x**3 + 40*x**2 + 50*x - 20)/(45*x**6 + 12*x**5 + 60*x**4 - \
        24*x**3 + 18*x**2 + 10*x + 10),
        (30*x**10 + 75*x**9 - 120*x**8 - 120*x**7 - 30*x**6 - 300*x**5 - 240*x**3 - 45*x**2 - 15*x - 90)/\
        (40*x**9 - 120*x**8 - 80*x**7 - 40*x**6 - 24*x**5 - 24*x**4 + 24*x**3 - 30*x - 80),
        (-20*x**3 + 12*x**2 + 12*x - 12)/(9*x**3 + 4*x**2 + 24*x + 6),
    )]
    for w, b1, b2 in test_cases:
        assert w == riccati_inverse_normal(riccati_normal(w, x, b1, b2), x, b1, b2).cancel()


def test_poles():
    tests = [
    (
        Poly(3, x),
        Poly(3*x**2 - x, x),
        {0: 1, S(1)/3: 1}
    ),
    (
        Poly(6*x - x**3, x),
        Poly((2 - x)**2*(2*x - 1)**2, x),
        {2: 2, S(1)/2: 2}
    ),
    (
        Poly(x, x),
        Poly(x**2 - x - 1, x),
        {(1 + sqrt(5))/2: 1, (1 - sqrt(5))/2: 1}
    ),
    (
        Poly(x**3 + 1, x),
        Poly(x + 2, x),
        {-2: 1, oo: 2}
    ),
    (
        Poly(x**5 + x**3 + 4, x),
        Poly(x**4 + 3*x**3 - 3*x**2 + 3*x - 4, x),
        {I: 1, -I: 1, 1: 1, -4: 1, oo: 1}
    ),
    (
        Poly(1, x),
        Poly(x**5 - 4*x**4 + 4*x**3 + x**2 - 4*x + 4, x),
        {2: 2, -1: 1, (1 - sqrt(3)*I)/2: 1, (1 + sqrt(3)*I)/2: 1}
    ),
    (
        Poly(x**3 + 4*x**2 - x + 12, x),
        Poly(x**3 + 12*x**2 + 39*x + 28, x),
        {-1: 1, -4: 1, -7: 1}
    ),
    (
        Poly(5*x**5 + 4*x**4 + 4*x**3 + 2*x**2 + x + 7, x),
        Poly(x**4 + 6*x**3 + 11*x**2 + 6*x + 1, x),
        {-S(3)/2 - sqrt(5)/2: 2, -S(3)/2 + sqrt(5)/2: 2, oo: 1}
    ),
    (
        Poly(x**2 - x + 1, x),
        Poly(x**5 + 5*x**4 + 10*x**3 + 10*x**2 + 3*x + 1, x),
        {}
    ),
    (
        Poly(5*x**5 - 4*x**4 - 3*x**3 - 3*x**2 - 4, x),
        Poly(2*x**7 + 6*x**6 + 11*x**5 - 12*x**4 - 28*x**3 + 6*x**2 + 15*x, x),
        {1: 2, -1: 2, -S(3)/2 - sqrt(21)*I/2: 1, -S(3)/2 + sqrt(21)*I/2: 1, 0: 1}
    )]

    for num, den, poles in tests:
        assert find_poles(num, den, x) == poles


def test_val_at_inf():
    tests = [
    (
        Poly(12*x**8 - 12*x**7 - 11*x**6 + 8*x**5 + 3*x**4 - x**3 + x**2 - 11*x, x),
        Poly(-14*x**2 + x, x),
         -6
    ),
    (
        Poly(-8*x**4 + 14*x**3 + 10*x**2 + 3*x - 11, x),
        Poly(-9*x**4 + 3*x**3 + 15*x**2 - 6*x - 14, x),
         0
    ),
    (
        Poly(-6*x**3 - 8*x**2 + 8*x - 6, x),
        Poly(-5*x**3 + 12*x**2 - 6*x - 9, x),
         0
    ),
    (
        Poly(15*x**3 - 5*x**2 - 6*x - 2, x),
        Poly(9*x**5 + 12*x**4 + 7*x**3 + 11*x**2 - 2*x - 9, x),
         2
    ),
    (
        Poly(10*x**3 + 8*x**2 - 13*x + 6, x),
        Poly(-13*x**10 - x**9 + 5*x**8 + 7*x**7 + 10*x**6 + 6*x**5 - 7*x**4 + 11*x**3 - 8*x**2 + 5*x + 13, x),
         7
    ),
    (
        Poly(-15*x - 15, x),
        Poly(15*x**3 - 13*x**2 - 4*x + 2, x),
         2
    ),
    (
        Poly(-9*x**8 - 9*x**7 + 11*x**6 - 11*x**5 + 9*x**4 + 7*x**3 - 8*x**2 - 5*x - 12, x),
        Poly(9*x**10 + 4*x**9 + 12*x**8 - 8*x**7 - 3*x**6 - 14*x**5 + 8*x**4 - 15*x**3 + 4*x**2 + 9*x + 1, x),
         2
    ),
    (
        Poly(13*x - 9, x),
        Poly(-8*x**4 + 5*x**3 - 4*x**2 - 14*x + 4, x),
         3
    ),
    (
        Poly(12*x**7 - 7*x**6 + 9*x**5 + 6*x**4 - 14*x**3 + 15*x**2 - x + 8, x),
        Poly(11*x**10 + 2*x**9 + x**8 - 3*x**7 - 12*x**6 - 10*x**5 - 15*x**4 - 14*x**3 - 5*x**2 + 15*x + 15, x),
         3
    ),
    (
        Poly(5*x**6 + 9*x**5 - 11*x**4 - 9*x**3 + x**2 - 4*x + 4, x),
        Poly(15*x**4 + 3*x**3 - 8*x**2 + 15*x + 12, x),
         -2
    )]
    for num, den, val in tests:
        assert val_at_inf(num, den, x) == val


def test_necessary_conds():
    # Valuation at Infinity is an odd negative integer
    assert check_necessary_conds(-3, [1, 2, 4]) == False
    # Valuation at Infinity is a positive integer lesser than 2
    assert check_necessary_conds(1, [1, 2, 4]) == False
    # Multiplicity of a pole is an odd integer greater than 1
    assert check_necessary_conds(2, [3, 1, 6]) == False
    # All values are correct
    assert check_necessary_conds(-10, [1, 2, 8, 12]) == True


def test_inverse_transform_poly():
    fns = [
        (30*x**6 - 12*x**5 + 15*x**4 - 15*x**2 + 10*x + 60)/(3*x**10 - 45*x**9 + 15*x**5 + 15*x**4 - 5*x**3 \
        + 15*x**2 + 45*x - 15),
        (15*x**3 - 8*x**2 - 2*x - 6)/(18*x + 6),
        (-75*x**6 + 12*x**5 - 80*x**4 + 15*x**3 - 75*x**2 + 120*x + 120)/(80*x**9 - 48*x**8 - 40*x**6 - \
        80*x**5 - 100*x**4 - 60*x**2 + 60*x + 60),
        (-40*x**4 + 180*x**3 - 60*x**2 + 60*x + 150)/(30*x**6 - 40*x**5 - 20*x**4 - 120*x**3 - 45*x - 48),
        (30*x**9 - 60*x**8 - 40*x**7 - 60*x**6 + 15*x**5 - 240*x**4 + 300*x**3 + 180*x**2 + 30*x - 90)/(\
        60*x**10 - 120*x**9 - 60*x**8 - 240*x**7 - 48*x**6 + 60*x**5 + 30*x**4 + 150*x**3 + 30*x**2 - 60*x),
        (24*x**10 + 20*x**7 - 48*x**5 - 45*x**4 - 60*x**3 + 60*x**2 + 60*x + 20)/(30*x**5 - 45*x**4 - 48*x**3
        - 24*x**2 + 30*x),
        (180*x**5 + 40*x**4 + 80*x**3 + 30*x**2 - 60*x - 80)/(180*x**3 - 150*x**2 + 75*x + 12),
        (-180*x**9 - 40*x**8 + 120*x**7 - 20*x**6 - 75*x**5 + 60*x**4 + 15*x**3 - 20*x**2 - 90*x + 15)/(\
        60*x**7 - 60*x**6 - 120*x**4 - 48*x**3 - 120*x**2 + 12*x + 100),
        (-15*x**5 - 36*x**4 + 75*x**3 - 60*x**2 - 80*x - 60)/(80*x**4 + 60*x**3 + 60*x**2 + 60*x - 80),
        (60*x**7 + 24*x**6 - 15*x**5 - 20*x**4 + 30*x**2 + 100*x - 60)/(240*x**2 - 20*x - 30)
    ]
    for f in fns:
        num, den = [Poly(e, x) for e in f.as_numer_denom()]
        num, den = inverse_transform_poly(num, den, x)
        assert f.subs(x, 1/x).cancel() == num/den


def test_limit_at_inf():
    tests = [
    (
        Poly(-12*x**2 + 20*x + 32, x),
        Poly(32*x**3 + 72*x**2 + 3*x - 32, x),
        0
    ),
    (
        Poly(25200*x**3 - 252*x**2 - 2205*x - 12600, x),
        Poly(840*x**9 - 630*x**8 + 280*x**7 + 1080*x**6 - 4032*x**5 - 10080*x**4 - \
        1400*x**3 + 4200*x**2 - 1575*x + 4410, x),
        0
    ),
    (
        Poly(12*x**5 - 840*x**3 - 140*x**2 + 840*x - 1260, x),
        Poly(2940*x**6 - 700*x**5 - 1890*x**4 + 7560*x**3 + 720*x**2 - 1134*x + 2205, x),
        0
    ),
    (
        Poly(1260*x**4 - 1260*x**3 - 700*x**2 - 1260*x + 1400, x),
        Poly(6300*x**3 - 1575*x**2 + 756*x - 540, x),
        oo
    ),
    (
        Poly(20*x**2 - 45*x + 105, x),
        Poly(30*x + 24, x),
        oo
    ),
    (
        Poly(-735*x**8 - 1400*x**7 + 1680*x**6 - 315*x**5 - 600*x**4 + 840*x**3 - 525*x**2 \
        + 630*x + 3780, x),
        Poly(1008*x**7 - 2940*x**6 - 84*x**5 + 2940*x**4 - 420*x**3 + 1512*x**2 + 105*x + 168, x),
        -oo
    ),
    (
        Poly(105*x**7 - 960*x**6 + 60*x**5 + 60*x**4 - 80*x**3 + 45*x**2 + 120*x + 15, x),
        Poly(735*x**7 + 525*x**6 + 720*x**5 + 720*x**4 - 8400*x**3 - 2520*x**2 + 2800*x + 280, x),
        S(1)/7
    ),
    (
        Poly(288*x**4 - 450*x**3 + 280*x**2 - 900*x - 90, x),
        Poly(607*x**4 + 840*x**3 - 1050*x**2 + 420*x + 420, x),
        S(288)/607
    ),
    (
        Poly(630*x**7 + 945*x**6 - 1890*x**5 + 756*x**4 + 1120*x**3 - 1620*x**2 + 420*x - 180, x),
        Poly(509*x**7 - 252*x**6 + 1470*x**5 + 3150*x**4 + 1764*x**3 - 1680*x**2 - 1134*x + 12600, x),
        S(630)/509
    )]
    for num, den, lim in tests:
        assert limit_at_inf(num, den, x) == lim

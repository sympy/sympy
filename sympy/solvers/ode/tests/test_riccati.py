from random import randint
from sympy import (S, symbols, Function, Rational, Poly, Eq, ratsimp,
    checkodesol, sqrt, oo, I, pi, Mul, sin, exp, log, tanh)
from sympy.testing.pytest import slow
from sympy.solvers.ode.riccati import (riccati_normal, riccati_inverse_normal,
    riccati_reduced, match_riccati, inverse_transform_poly, find_poles,
    limit_at_inf, check_necessary_conds, val_at_inf, construct_c_case_1,
    construct_c_case_2, construct_c_case_3, construct_c, construct_d_case_4,
    construct_d_case_5, construct_d_case_6, construct_d, rational_laurent_series)

f = Function('f')
x = symbols('x')

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
    assert checkodesol(eq, sol) == (True, 0)
    return eq, q0, q1, q2


# Testing functions start

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
        (-180*x**10 - 120*x**9 + 36*x**8 + 120*x**7 + 60*x**6 - 120*x**4 + 30*x**3 - 90*x**2 + 60)/( \
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
        60*x - 30)/(120*x**10 - 30*x**9 - 120*x**8 - 120*x**7 + 150*x**6 - 60*x**5 - 15*x**4 - 150*x**3 - \
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
        (30*x**10 + 75*x**9 - 120*x**8 - 120*x**7 - 30*x**6 - 300*x**5 - 240*x**3 - 45*x**2 - 15*x - 90)/ \
        (40*x**9 - 120*x**8 - 80*x**7 - 40*x**6 - 24*x**5 - 24*x**4 + 24*x**3 - 30*x - 80),
        (-20*x**3 + 12*x**2 + 12*x - 12)/(9*x**3 + 4*x**2 + 24*x + 6),
    )]
    for w, b1, b2 in test_cases:
        assert w == riccati_inverse_normal(riccati_normal(w, x, b1, b2), x, b1, b2).cancel()


def test_riccati_reduced():
    tests = [
    (
        f(x).diff(x) - x**2 - x*f(x) - x*f(x)**2,

        f(x).diff(x) + f(x)**2 + x**3 - x**2/4 - 3/(4*x**2)
    ),
    (
        f(x).diff(x) - 1/x + (x**2 - x)*f(x)/3 + f(x)**2/(x**2 - 2),

        -3*x**2/(x**2 - 2)**2 - 2*x*(1 - x**2/2)*(-x**2/3 + x/3)/( \
        x**2 - 2)**2 - x/3 - Mul(2, 1 - x**2/2, evaluate=False)* \
        (4*x**2/(x**2 - 2) - 1)/(x**2 - 2)**2 - (-x**2/3 + x/3 \
        )**2/4 + f(x)**2 + f(x).diff(x) + S(1)/6 - 1/(x*(x**2 - 2))
    ),
    (
        6*x/(2*x + 9) + f(x).diff(x) - (x + 1)*f(x)**2/x,

        -3*x**2*(1/x + (-x - 1)/x**2)**2/(4*(-x - 1)**2) + Mul(6, \
        -x - 1, evaluate=False)/(2*x + 9) + f(x)**2 + f(x).diff(x) \
        - (-1 + (x + 1)/x)/(x*(-x - 1))
    ),
    (
        -(S(3)/2 - 3*x)*f(x)/(-x - 3) + f(x).diff(x) - (2 - 3*x)/( \
        6*x) - f(x)**2/x,

        -(3 - 6*x)**2/(4*(2*x + 6)**2) + (3 - 6*x)/(2*x + 6)**2 \
        + f(x)**2 + f(x).diff(x) + 3/(2*x + 6) - (S(1)/2 - 1/(3*x)) \
        /x - (3 - 6*x)/(2*x*(2*x + 6)) + 1/(4*x**2)
    ),
    (
        -x*f(x)**2 + f(x).diff(x) + 1 - (x - S(1)/2)*f(x)/(x - S(2)/3),

        -x - (3 - 6*x)**2/(4*(6*x - 4)**2) + Mul(3, 3 - 6*x, evaluate= \
        False)/(6*x - 4)**2 + f(x)**2 + f(x).diff(x) + 3/(6*x - 4) + \
        (3 - 6*x)/(2*x*(6*x - 4)) - 3/(4*x**2)
    ),
    (
        f(x)**2 + f(x).diff(x) - (x - 1)*f(x)/(-x - S(1)/2),

        -(2*x - 2)**2/(4*(2*x + 1)**2) + (2*x - 2)/(2*x + 1)**2 + \
        f(x)**2 + f(x).diff(x) - 1/(2*x + 1)
    ),
    (
        f(x).diff(x) - f(x)**2/x,

        f(x)**2 + f(x).diff(x) + 1/(4*x**2)
    ),
    (
        -(3 - 4*x)/(2*x + 9) - (-3*x/2 - S(1)/2)*f(x)/(x - 2) + \
        f(x).diff(x) - (x - 1)*f(x)**2/x,

        -3*x**2*(1/x + (1 - x)/x**2)**2/(4*(1 - x)**2) - x*(1/x + \
        (1 - x)/x**2)*(3*x + 1)/(Mul(2, 1 - x, evaluate=False)* \
        (2*x - 4)) + f(x)**2 + f(x).diff(x) - 3/(Mul(2, 2*x - 4, \
        evaluate=False)) - (3*x + 1)**2/(4*(2*x - 4)**2) + (3*x \
        + 1)/(2*x - 4)**2 - (-1 + (x - 1)/x)/(x*(1 - x)) + (1 - \
        x)*(4*x/(2*x + 9) - 3/(2*x + 9))/x
    ),
    (
        x*f(x)/(3*(x + 1)) - (4*x**2 - 9)/(6*(x**2 - x + 1)) + \
        f(x).diff(x) - 3*(4*x - 1)*f(x)**2/(4*(x + 3)),

        -x**2/(4*(3*x + 3)**2) + 3*x/(2*(3*x + 3)**2) - x*(4*x \
        + 12)*(Mul(4, 3 - 12*x, evaluate=False)/(4*x + 12)**2 + \
        12/(4*x + 12))/(Mul(2, 3 - 12*x, evaluate=False)*(3*x + \
        3)) - Mul(3, -2 + (4*x - 1)/(Mul(2, x + 3, evaluate=\
        False)), evaluate=False)*(4*x + 12)/(Mul(2, 3 - 12*x, \
        evaluate=False)*(x + 3)**2) + (3 - 12*x)*(-4*x**2/( \
        6*x**2 - 6*x + 6) + 9/(6*x**2 - 6*x + 6))/(4*x + 12) \
        + f(x)**2 + f(x).diff(x) - 1/(Mul(2, 3*x + 3, evaluate \
        =False)) - 3*(4*x + 12)**2*(Mul(4, 3 - 12*x, evaluate= \
        False)/(4*x + 12)**2 + 12/(4*x + 12))**2/(4*(3 - 12*x)**2)
    ),
    (
        -3*(-x**2 - x + 1)/(x**2 + 6*x + 1) + f(x).diff(x) + f(x)**2/x,

        f(x)**2 + f(x).diff(x) + (3*x**2/(x**2 + 6*x + 1) + 3*x/(x**2 \
        + 6*x + 1) - 3/(x**2 + 6*x + 1))/x + 1/(4*x**2)
    )
    ]
    for eq, normal_eq in tests:
        assert normal_eq == riccati_reduced(eq, f, x)


def test_match_riccati():
    tests = [
    # Test Rational Riccati ODEs
    (
        f(x).diff(x) - (405*x**3 - 882*x**2 - 78*x + 92)/(243*x**4 \
        - 945*x**3 + 846*x**2 + 180*x - 72) - 2 - f(x)**2/(3*x + 1) \
        - (S(1)/3 - x)*f(x)/(S(1)/3 - 3*x/2),

        True,

        405*x**3/(243*x**4 - 945*x**3 + 846*x**2 + 180*x - 72) - \
        882*x**2/(243*x**4 - 945*x**3 + 846*x**2 + 180*x - 72) - \
        78*x/(243*x**4 - 945*x**3 + 846*x**2 + 180*x - 72) + 2 + \
        92/(243*x**4 - 945*x**3 + 846*x**2 + 180*x - 72),

        Mul(-1, 2 - 6*x, evaluate=False)/(9*x - 2),

        1/(3*x + 1)
    ),
    (
        f(x).diff(x) + 4*x/27 - (x/3 - 1)*f(x)**2 - (2*x/3 + \
        1)*f(x)/(3*x + 2) - S(10)/27 - (265*x**2 + 423*x + 162) \
        /(324*x**3 + 216*x**2),

        True,

        265*x**2/(324*x**3 + 216*x**2) - 4*x/27 + 423*x/(324*x**3 \
        + 216*x**2) + S(10)/27 + 162/(324*x**3 + 216*x**2),

        Mul(-1, -2*x - 3, evaluate=False)/(9*x + 6),

        x/3 - 1
    ),
    (
        f(x).diff(x) - (304*x**5 - 745*x**4 + 631*x**3 - 876*x**2 \
        + 198*x - 108)/(36*x**6 - 216*x**5 + 477*x**4 - 567*x**3 + \
        360*x**2 - 108*x) - S(17)/9 - (x - S(3)/2)*f(x)/(x/2 - \
        S(3)/2) - (x/3 - 3)*f(x)**2/(3*x),

        True,

        304*x**5/(36*x**6 - 216*x**5 + 477*x**4 - 567*x**3 + 360*x**2 \
        - 108*x) - 745*x**4/(36*x**6 - 216*x**5 + 477*x**4 - 567*x**3 \
        + 360*x**2 - 108*x) + 631*x**3/(36*x**6 - 216*x**5 + 477*x**4 \
        - 567*x**3 + 360*x**2 - 108*x) - 876*x**2/(36*x**6 - 216*x**5 \
        + 477*x**4 - 567*x**3 + 360*x**2 - 108*x) + 198*x/(36*x**6 - \
        216*x**5 + 477*x**4 - 567*x**3 + 360*x**2 - 108*x) + S(17)/9 - 108 \
        /(36*x**6 - 216*x**5 + 477*x**4 - 567*x**3 + 360*x**2 - 108*x),

        Mul(-1, 3 - 2*x, evaluate=False)/(x - 3),

        Mul(-1, 9 - x, evaluate=False)/(9*x)
    ),
    # Test Non-Rational Riccati ODEs
    (
        f(x).diff(x) - x**(S(3)/2)/(x**(S(1)/2) - 2) + x**2*f(x) + \
        x*f(x)**2/(x**(S(3)/4)),
        False, 0, 0, 0
    ),
    (
        f(x).diff(x) - sin(x**2) + exp(x)*f(x) + log(x)*f(x)**2,
        False, 0, 0, 0
    ),
    (
        f(x).diff(x) - tanh(x + sqrt(x)) + f(x) + x**4*f(x)**2,
        False, 0, 0, 0
    ),
    # Test Non-Riccati ODEs
    (
        (1 - x**2)*f(x).diff(x, 2) - 2*x*f(x).diff(x) + 20*f(x),
        False, 0, 0, 0
    ),
    (
        f(x).diff(x) - x**2 + x**3*f(x) + (x**2/(x + 1))*f(x)**3,
        False, 0, 0, 0
    ),
    (
        f(x).diff(x)*f(x)**2 + (x**2 - 1)/(x**3 + 1)*f(x) + 1/(2*x \
        + 3) + f(x)**2,
        False, 0, 0, 0
    )
    ]
    for eq, res, b0, b1, b2 in tests:
        match, funcs = match_riccati(eq, f, x)
        assert match == res
        if res:
            assert [b0, b1, b2] == funcs


def test_poles():
    tests = [
    (
        Poly(3*x**2 - x, x),
        {0: 1, S(1)/3: 1}
    ),
    (
        Poly((2 - x)**2*(2*x - 1)**2, x),
        {2: 2, S(1)/2: 2}
    ),
    (
        Poly(x**2 - x - 1, x),
        {(1 + sqrt(5))/2: 1, (1 - sqrt(5))/2: 1}
    ),
    (
        Poly(x + 2, x),
        {-2: 1}
    ),
    (
        Poly(x**4 + 3*x**3 - 3*x**2 + 3*x - 4, x),
        {I: 1, -I: 1, 1: 1, -4: 1}
    ),
    (
        Poly(x**5 - 4*x**4 + 4*x**3 + x**2 - 4*x + 4, x),
        {2: 2, -1: 1, (1 - sqrt(3)*I)/2: 1, (1 + sqrt(3)*I)/2: 1}
    ),
    (
        Poly(x**3 + 12*x**2 + 39*x + 28, x),
        {-1: 1, -4: 1, -7: 1}
    ),
    (
        Poly(x**4 + 6*x**3 + 11*x**2 + 6*x + 1, x),
        {-S(3)/2 - sqrt(5)/2: 2, -S(3)/2 + sqrt(5)/2: 2}
    ),
    (
        Poly(x**5 + 5*x**4 + 10*x**3 + 10*x**2 + 3*x + 1, x),
        {}
    ),
    (
        Poly(2*x**7 + 6*x**6 + 11*x**5 - 12*x**4 - 28*x**3 + 6*x**2 + 15*x, x),
        {1: 2, -1: 2, -S(3)/2 - sqrt(21)*I/2: 1, -S(3)/2 + sqrt(21)*I/2: 1, 0: 1}
    )]

    for den, poles in tests:
        assert find_poles(den, x) == poles


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
        (30*x**9 - 60*x**8 - 40*x**7 - 60*x**6 + 15*x**5 - 240*x**4 + 300*x**3 + 180*x**2 + 30*x - 90)/( \
        60*x**10 - 120*x**9 - 60*x**8 - 240*x**7 - 48*x**6 + 60*x**5 + 30*x**4 + 150*x**3 + 30*x**2 - 60*x),
        (24*x**10 + 20*x**7 - 48*x**5 - 45*x**4 - 60*x**3 + 60*x**2 + 60*x + 20)/(30*x**5 - 45*x**4 - 48*x**3
        - 24*x**2 + 30*x),
        (180*x**5 + 40*x**4 + 80*x**3 + 30*x**2 - 60*x - 80)/(180*x**3 - 150*x**2 + 75*x + 12),
        (-180*x**9 - 40*x**8 + 120*x**7 - 20*x**6 - 75*x**5 + 60*x**4 + 15*x**3 - 20*x**2 - 90*x + 15)/( \
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


def test_construct_c_case_1():
    tests = [
    (
        Poly(-3*x**3 + 3*x**2 + 4*x - 5, x, extension=True),
        Poly(4*x**8 + 16*x**7 + 9*x**5 + 12*x**4 + 6*x**3 + 12*x**2, x, extension=True),
        S(0),
        [[S(1)/2 + sqrt(6)*I/6], [S(1)/2 - sqrt(6)*I/6]]
    ),
    (
        Poly(1200*x**3 + 1440*x**2 + 816*x + 560, x, extension=True),
        Poly(128*x**5 - 656*x**4 + 1264*x**3 - 1125*x**2 + 385*x + 49, x, extension=True),
        S(7)/4,
        [[S(1)/2 + sqrt(16367978)/634], [S(1)/2 - sqrt(16367978)/634]]
    ),
    (
        Poly(4*x + 2, x, extension=True),
        Poly(18*x**4 + (2 - 18*sqrt(3))*x**3 + (14 - 11*sqrt(3))*x**2 + (4 - 6*sqrt(3))*x \
            + 8*sqrt(3) + 16, x, domain='QQ<sqrt(3)>'),
        (S(1) + sqrt(3))/2,
        [[S(1)/2 + sqrt(Mul(4, 2*sqrt(3) + 4, evaluate=False)/(19*sqrt(3) + 44) + 1)/2], \
            [S(1)/2 - sqrt(Mul(4, 2*sqrt(3) + 4, evaluate=False)/(19*sqrt(3) + 44) + 1)/2]]
    )]
    for num, den, pole, c in tests:
        assert construct_c_case_1(num, den, x, pole) == c


def test_construct_c_case_2():
    tests = [
    # Testing poles with multiplicity 2
    (
        Poly(1, x, extension=True),
        Poly((x - 1)**2*(x - 2), x, extension=True),
        1, 2,
        [[-I*(-1 - I)/2], [I*(-1 + I*(-1 - I)/2)/2]]
    ),
    (
        Poly(3*x**5 - 12*x**4 - 7*x**3 + 1, x, extension=True),
        Poly((3*x - 1)**2*(x + 2)**2, x, extension=True),
        S(1)/3, 2,
        [[-S(89)/98], [-S(103)/28]]
    ),
    # Testing poles with multiplicity 4
    (
        Poly(x**3 - x**2 + 4*x, x, extension=True),
        Poly((x - 2)**4*(x + 5)**2, x, extension=True),
        2, 4,
        [[7*sqrt(3)*(S(60)/343 - 4*sqrt(3)/7)/12, 2*sqrt(3)/7], \
        [-7*sqrt(3)*(S(60)/343 - 4*sqrt(3)/7)/12, -2*sqrt(3)/7]]
    ),
    (
        Poly(3*x**5 + x**4 + 3, x, extension=True),
        Poly((4*x + 1)**4*(x + 2), x, extension=True),
        -S(1)/4, 4,
        [[128*sqrt(439)*(-sqrt(439)/128 - S(55)/14336)/439, sqrt(439)/256], \
        [-128*sqrt(439)*(-sqrt(439)/128 - S(55)/14336)/439, -sqrt(439)/256]]
    ),
    # Testing poles with multiplicity 6
    (
        Poly(x**3 + 2, x, extension=True),
        Poly((3*x - 1)**6*(x**2 + 1), x, extension=True),
        S(1)/3, 6,
        [[27*sqrt(66)*(-sqrt(66)/54 - S(131)/267300)/22, -2*sqrt(66)/1485, sqrt(66)/162], \
        [-27*sqrt(66)*(-sqrt(66)/54 - S(131)/267300)/22, 2*sqrt(66)/1485, -sqrt(66)/162]]
    ),
    (
        Poly(x**2 + 12, x, extension=True),
        Poly((x - sqrt(2))**6, x, extension=True),
        sqrt(2), 6,
        [[sqrt(14)*(S(6)/7 - 3*sqrt(14))/28, sqrt(7)/7, sqrt(14)], \
        [-sqrt(14)*(S(6)/7 - 3*sqrt(14))/28, -sqrt(7)/7, -sqrt(14)]]
    )]
    for num, den, pole, mul, c in tests:
        assert construct_c_case_2(num, den, x, pole, mul) == c


def test_construct_c_case_3():
    assert construct_c_case_3() == [[1], [1]]


def test_construct_d_case_4():
    tests = [
    # Tests with multiplicity at oo = 2
    (
        Poly(-x**5 - 2*x**4 + 4*x**3 + 2*x + 5, x, extension=True),
        Poly(9*x**3 - 2*x**2 + 10*x - 2, x, extension=True),
        2,
        [[1193*sqrt(374)/90882, sqrt(374)/27, 27*sqrt(374)*(-sqrt(374)/27 \
        -S(2843831)/7361442)/748], [-1193*sqrt(374)/90882, -sqrt(374)/27, \
        -27*sqrt(374)*(-S(2843831)/7361442 + sqrt(374)/27)/748]]
    ),
    (
        Poly(-x**6 + 9*x**5 + 5*x**4 + 6*x**3 + 5*x**2 + 6*x + 7, x, extension=True),
        Poly(x**4 + 3*x**3 + 12*x**2 - x + 7, x, extension=True),
        2,
        [[-249*sqrt(82)*I/82, sqrt(82)*I, -sqrt(82)*I*(S(12227)/82 - sqrt(82)*I)/164], \
        [249*sqrt(82)*I/82, -sqrt(82)*I, sqrt(82)*I*(S(12227)/82 + sqrt(82)*I)/164]]
    ),
    # Tests with multiplicity at oo = 4
    (
        Poly(-2*x**6 - x**5 - x**4 - 2*x**3 - x**2 - 3*x - 3, x, extension=True),
        Poly(3*x**2 + 10*x + 7, x, extension=True),
        4,
        [[15481*sqrt(17)/20808, -137*sqrt(17)/306, sqrt(17)/3, 3*sqrt(17)*( \
        -S(9345037)/561816 - 2*sqrt(17)/3)/34], [-15481*sqrt(17)/20808, 137*sqrt \
        (17)/306, -sqrt(17)/3, -3*sqrt(17)*(-S(9345037)/561816 + 2*sqrt(17)/3)/34]]
    ),
    (
        Poly(-3*x**5 - 3*x**4 - 3*x**3 - x**2 - 1, x, extension=True),
        Poly(12*x - 2, x, extension=True),
        4,
        [[41*I/192, 7*I/24, I/2, -I*(-S(59)/6912 - I)], \
        [-41*I/192, -7*I/24, -I/2, I*(-S(59)/6912 + I)]]
    ),
    # Tests with multiplicity at oo = 4
    (
        Poly(-x**7 - x**5 - x**4 - x**2 - x, x, extension=True),
        Poly(x + 2, x, extension=True),
        6,
        [[-5*I/2, 2*I, -I, I, -I*(-9 - 3*I)/2], [5*I/2, -2*I, I, -I, I*(-9 + 3*I)/2]]
    ),
    (
        Poly(-x**7 - x**6 - 2*x**5 - 2*x**4 - x**3 - x**2 + 2*x - 2, x, extension=True),
        Poly(2*x - 2, x, extension=True),
        6,
        [[3*sqrt(2)*I/4, 3*sqrt(2)*I/4, sqrt(2)*I/2, sqrt(2)*I/2, -sqrt(2)*I*(-S(7)/8 - \
        3*sqrt(2)*I/2)/2], [-3*sqrt(2)*I/4, -3*sqrt(2)*I/4, -sqrt(2)*I/2, -sqrt(2)*I/2, \
        sqrt(2)*I*(-S(7)/8 + 3*sqrt(2)*I/2)/2]]
    )]
    for num, den, mul, d in tests:
        ser = rational_laurent_series(num, den, x, oo, mul, 1)
        assert construct_d_case_4(ser, mul//2, mul) == d


def test_construct_d_case_5():
    tests = [
    (
        Poly(2*x**3 + x**2 + x - 2, x, extension=True),
        Poly(9*x**3 + 5*x**2 + 2*x - 1, x, extension=True),
        [[5*sqrt(2)/27, -382*sqrt(2)/1215], [-5*sqrt(2)/27, 382*sqrt(2)/1215]]
    ),
    (
        Poly(3*x**5 + x**4 - x**3 + x**2 - 2*x - 2, x, domain='ZZ'),
        Poly(9*x**5 + 7*x**4 + 3*x**3 + 2*x**2 + 5*x + 7, x, domain='ZZ'),
        [[sqrt(27798)*I/243, 2914*sqrt(27798)*I/10132371], [-sqrt(27798)*I/243, \
        -2914*sqrt(27798)*I/10132371]]
    ),
    (
        Poly(x**2 - x + 1, x, domain='ZZ'),
        Poly(3*x**2 + 7*x + 3, x, domain='ZZ'),
        [[sqrt(10)*I/3, -7*sqrt(10)*I/18], [-sqrt(10)*I/3, 7*sqrt(10)*I/18]]
    )]
    for num, den, d in tests:
        # Multiplicity of oo is 0
        ser = rational_laurent_series(num, den, x, oo, 0, 1)
        assert construct_d_case_5(ser) == d


def test_construct_d_case_6():
    tests = [
    (
        Poly(-2*x**2 - 5, x, domain='ZZ'),
        Poly(4*x**4 + 2*x**2 + 10*x + 2, x, domain='ZZ'),
        [[S(1)/2 + I/2], [S(1)/2 - I/2]]
    ),
    (
        Poly(-2*x**3 - 4*x**2 - 2*x - 5, x, domain='ZZ'),
        Poly(x**6 - x**5 + 2*x**4 - 4*x**3 - 5*x**2 - 5*x + 9, x, domain='ZZ'),
        [[1], [0]]
    ),
    (
        Poly(-5*x**3 + x**2 + 11*x + 12, x, domain='ZZ'),
        Poly(6*x**8 - 26*x**7 - 27*x**6 - 10*x**5 - 44*x**4 - 46*x**3 - 34*x**2 \
        - 27*x - 42, x, domain='ZZ'),
        [[1], [0]]
    )]
    for num, den, d in tests:
        assert construct_d_case_6(num, den, x) == d


def test_rational_laurent_series():
    tests = [
    # Laurent series about simple pole (Multiplicity = 1)
    (
        Poly(x**2 - 3*x + 9, x, extension=True),
        Poly(x**2 - x, x, extension=True),
        S(1), 1, 6,
        [7, -8, 9, -9, 9, -9, 9]
    ),
    # Laurent series about multiple pole (Multiplicty > 1)
    (
        Poly(64*x**3 - 1728*x + 1216, x, extension=True),
        Poly(64*x**4 - 80*x**3 - 831*x**2 + 1809*x - 972, x, extension=True),
        S(9)/8, 2, 3,
        [S(1019)/984, S(209149)/75645, S(32177152)/46521675, S(11947565056)/28610830125, \
        S(3751586234368)/17595660526875]
    ),
    (
        Poly(1, x, extension=True),
        Poly(x**5 + (-4*sqrt(2) - 1)*x**4 + (4*sqrt(2) + 12)*x**3 + (-12 - 8*sqrt(2))*x**2 \
        + (4 + 8*sqrt(2))*x - 4, x, extension=True),
        sqrt(2), 4, 6,
        [1 + sqrt(2), -3 - 2*sqrt(2), 7 + 5*sqrt(2), -17 - 12*sqrt(2), 41 + 29*sqrt(2), \
        -99 - 70*sqrt(2), 239 + 169*sqrt(2), -577 - 408*sqrt(2), 1393 + 985*sqrt(2), \
        -3363 - 2378*sqrt(2)]
    ),
    # Laurent series about oo
    (
        Poly(x**5 - 4*x**3 + 6*x**2 + 10*x - 13, x, extension=True),
        Poly(x**2 - 5, x, extension=True),
        oo, 3, 6,
        [375, 85, 75, 17, 15, 6, 1, 0, 1]
    ),
    # Laurent series at x0 where x0 is not a pole of the function
    # Using multiplicity as 0 (as x0 will not be a pole)
    (
        Poly(3*x**3 + 6*x**2 - 2*x + 5, x, extension=True),
        Poly(9*x**4 - x**3 - 3*x**2 + 4*x + 4, x, extension=True),
        pi, 0, 1,
        [(-2*pi + 5 + 6*pi**2 + 3*pi**3)/(-pi**3 - 3*pi**2 + 4 + 4*pi + 9*pi**4), \
        (-108*pi**5 - 27*pi**6 - 160*pi**3 - 28 + 78*pi + 69*pi**2 + 51*pi**4)/ \
        (-18*pi**7 - 53*pi**6 - 32*pi**3 - 8*pi**2 + 16 + 32*pi + 73*pi**4 + \
        78*pi**5 + 81*pi**8)]
    ),
    (
        Poly(-7*x**2 + 2*x - 4, x, extension=True),
        Poly(7*x**5 + 9*x**4 + 8*x**3 + 3*x**2 + 6*x + 9, x, extension=True),
        oo, 0, 6,
        [-S(71)/49, S(11)/7, -S(1), 0, 0, 0]
    )]
    for num, den, x0, mul, n, ser in tests:
        assert ser == rational_laurent_series(num, den, x, x0, mul, n)

# def test_solve_riccati():
#     a, b, c, k, C1 = symbols('a b c k C1')
#     C0 = Dummy('C0')
#     # Type: 1st Order Rational Riccati, dy/dx = a + b*y + c*y**2,
#     # a, b, c are rational functions of x

#     eq = x**2 - (2*x + 1/x)*f(x) + f(x)**2 + f(x).diff(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (x**2 + 2)/x), Eq(f(x), (C0*x + x**3 + 2*x)/(C0 + x**2)), \
#         Eq(f(x), x)]
#     assert all([s1.dummy_eq(s2, C0) for s1, s2 in zip(sols, solse)])

#     eq = f(x)**2 + f(x).diff(x) - (4*x**6 - 8*x**5 + 12*x**4 + 4*x**3 + \
#             7*x**2 - 20*x + 4)/(4*x**4)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (2*x**5 - 2*x**4 - x**3 + 4*x**2 + 3*x - 2)/(2*x**4 - 2*x**2))]
#     assert all([s1.dummy_eq(s2, C0) for s1, s2 in zip(sols, solse)])

#     eq = -x*f(x)**2 + f(x).diff(x) - 2*f(x)/x
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), -4/x**2), Eq(f(x), -4*x**2/(C0 + x**4)), Eq(f(x), 0)]
#     assert all([s1.dummy_eq(s2, C0) for s1, s2 in zip(sols, solse)])

#     eq = -f(x)**2 + f(x).diff(x) + (15*x**2 - 20*x + 7)/((x - 1)**2*(2*x \
#             - 1)**2)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-15*x**2 + 15*x - 4)/(6*x**3 - 11*x**2 + 6*x - 1)), \
#         Eq(f(x), (-60*x**3 + 180*x**2 - 181*x + 62)/(24*x**4 - 92*x**3 + 130*x**2 \
#         - 79*x + 17)), Eq(f(x), (9*C0*x - 6*C0 - 15*x**5 + 60*x**4 - 94*x**3 + \
#         72*x**2 - 30*x + 6)/(6*C0*x**2 - 9*C0*x + 3*C0 + 6*x**6 - 29*x**5 + \
#         57*x**4 - 58*x**3 + 30*x**2 - 6*x)), Eq(f(x), (3*x - 2)/(2*x**2 - 3*x + 1))]
#     assert all([s1.dummy_eq(s2, C0) for s1, s2 in zip(sols, solse)])

#     eq = 9*x**2/4 - f(x)**2 + f(x).diff(x) - S(21)/2
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 3*x/2 - (3*x**2 - 1)/(x*(x**2 - 1)))]

#     eq = f(x)**2 + f(x).diff(x) - 15/(4*x**2)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-3*C1 + 5*x**4)/(2*x*(C1 + x**4))), Eq(f(x), -3/(2*x))]

#     eq = 3*f(x)**2 + f(x).diff(x) - 2/x**2
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-2*C1 + 3*x**5)/(3*x*(C1 + x**5))), Eq(f(x), -2/(3*x))]

#     eq = -a*b*k**2*f(x)**2 + 2*a*b*k*f(x) - a*b + f(x).diff(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (C1*a*b*k + a*b*k*x - 1)/(a*b*k**2*(C1 + x)))]

#     eq = f(x).diff(x) - 2*I*(f(x)**2 + 1)/x
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-I*C1 + I*x**4)/(C1 + x**4)), Eq(f(x), -I)]

#     eq = f(x).diff(x) - f(x)**2/x + 1/x
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (C1 - x**2)/(C1 + x**2)), Eq(f(x), 1)]

#     eq = a*f(x)**2 + c + x*f(x).diff(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), sqrt(-a*c)*(-C1 + x**(2*sqrt(-a*c)))/(a*(C1 + x**(2*sqrt(-a*c)))))]

#     eq = a*f(x)**2/x - b*f(x)/x - c/x + f(x).diff(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (C1*b - C1*sqrt(4*a*c + b**2) + b*x**(sqrt(4*a*c + b**2)) + \
#             x**(sqrt(4*a*c + b**2))*sqrt(4*a*c + b**2))/(2*a*(C1 + x**(sqrt(4*a*c + b**2)))))]

#     eq = f(x)**2 + f(x).diff(x) - f(x)/x
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 2*x/(C1 + x**2)), Eq(f(x), 0)]

#     eq = -x**2 - (2*x + 1/x)*f(x) - f(x)**2 + f(x).diff(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), x*(-C1 - x**2 - 2)/(C1 + x**2)), Eq(f(x), -x)]

#     eq = f(x)**2 + f(x).diff(x) + 4*f(x)/x + 2/x**2
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-2*C1 - x)/(x*(C1 + x)))]

#     eq = a*f(x)**2 - b/x**2 + f(x).diff(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-C1*sqrt(4*a*b + 1) + C1 + x**(sqrt(4*a*b + 1))*sqrt(4*a*b + 1) + \
#           x**(sqrt(4*a*b + 1)))/(2*a*x*(C1 + x**(sqrt(4*a*b + 1)))))]

#     eq = a*f(x)**2 + f(x).diff(x) + (b + c)/x**2
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-C1*sqrt(-4*a*b - 4*a*c + 1) + C1 + x**(sqrt(-4*a*b - 4*a*c + 1))*sqrt(\
#         -4*a*b - 4*a*c + 1) + x**(sqrt(-4*a*b - 4*a*c + 1)))/(2*a*x*(C1 + x**(sqrt(-4*a*b - \
#         4*a*c + 1)))))]

#     eq = 2*x**2*f(x).diff(x) - x*(4*f(x) + f(x).diff(x) - 4) + (f(x) - 1)*f(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (C1 + 2*x**2)/(C1 + x))]

#     eq = x**4*f(x).diff(x) + x**2 - x*(2*f(x)**2 + f(x).diff(x)) + f(x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), x*(C1*x + 1)/(C1 + x**2)), Eq(f(x), x**2)]

#     eq = f(x).diff(x) - a*f(x)**2 - b/x**2
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (C1*sqrt(-4*a*b + 1) - C1 - x**(sqrt(-4*a*b + 1))*sqrt(-4*a*b + 1) \
#         - x**(sqrt(-4*a*b + 1)))/(2*a*x*(C1 + x**(sqrt(-4*a*b + 1)))))]

#     eq = f(x).diff(x) - f(x)**2 - b/x**2 - a*f(x)/x,
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-C1*a + C1*sqrt(a**2 + 2*a - 4*b + 1) - C1 - a*x**(sqrt(a**2 + \
#         2*a - 4*b + 1)) - x**(sqrt(a**2 + 2*a - 4*b + 1))*sqrt(a**2 + 2*a - 4*b + 1) - \
#         x**(sqrt(a**2 + 2*a - 4*b + 1)))/(2*x*(C1 + x**(sqrt(a**2 + 2*a - 4*b + 1)))))]

#     eq = f(x).diff(x) - x*f(x)/(S(3)/2 - 2*x) - (x/2 - S(1)/3)*f(x)**2/\
#         (2*x/3 - S(1)/2) - S(5)/4 - (281*x**2 - 1260*x + 756)/(16*x**3 - 12*x**2))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (9 - x)/x), Eq(f(x), (40*x**14 + 28*x**13 + 420*x**12 + 2940*x**11 + \
#         18480*x**10 + 103950*x**9 + 519750*x**8 + 2286900*x**7 + 8731800*x**6 + 28378350*\
#         x**5 + 76403250*x**4 + 163721250*x**3 + 261954000*x**2 + 278326125*x + 147349125)/\
#         (x*(24*x**13 + 140*x**12 + 840*x**11 + 4620*x**10 + 23100*x**9 + 103950*x**8 + \
#         415800*x**7 + 1455300*x**6 + 4365900*x**5 + 10914750*x**4 + 21829500*x**3 + 32744250\
#         *x**2 + 32744250*x + 16372125)))]

#     eq = Eq(f(x).diff(x), (-12*x**2 - 48*x - 15)/(24*x**3 - 40*x**2 + 8*x + 8) \
#         + 3*f(x)**2/(6*x + 2))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (2*x + 1)/(Mul(2, x - 1, evaluate=False)))]

#     eq = Eq(f(x).diff(x), 2*x*f(x)**2/(9*x - 6) + (32 - 28*x)/(36*x**3 + \
#         12*x**2 - 15*x - 6) - 4*f(x)/(9*x - 6))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 2/(2*x + 1))]

#     eq = Eq(f(x).diff(x), x*f(x) + 2*x + (3*x - 2)*f(x)**2/(4*x + 2) + \
#         (8*x**2 - 7*x + 26)/(16*x**3 - 24*x**2 + 8) - S(3)/2)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (1 - 4*x)/(Mul(2, x - 1, evaluate=False)))]

#     eq = Eq(f(x).diff(x) + f(x)**2 - 2, 0)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), sqrt(2)), Eq(f(x), -sqrt(2))]

#     eq = Eq(f(x).diff(x), 8*x**3/729 - 8*x**2/243 + 4*x/9 + (S(2)/3 - 2*x/9)*\
#             f(x)**2 - S(4)/9 + 2/(3*x + 3) + (2*x - 1)*f(x)/(-x - 1))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 2*x/9)]

#     eq = Eq(f(x).diff(x), 2*x**3 + 4*x**2 - 7*x + (S(3)/2 - x/2)*f(x)**2 + \
#             (3*x - S(3)/2)*f(x) - 11)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), -2*x - 2)]

#     eq = Eq(f(x).diff(x), -9*x**3/8 + 3*x**2 + x/2 + (x/2 - S(2)/3)*f(x)**2 \
#             - S(1)/3 - 5/(6*x - 2) + (-2*x - 1)*f(x)/(S(1)/3 - x))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 1 - 3*x/2)]

#     eq = Eq(f(x).diff(x), -x**3/24 - x**2/8 - x/2 + (2*x/3 - S(2)/3)*f(x)**2 \
#             - S(11)/6 - 15/(4*x - 4) + (2*x + 3)*f(x)/(x - 1))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), x/4 + S(1)/2)]

#     eq = Eq(f(x).diff(x), 2*x**3/81 + x/18 - x*f(x)/(1 - x) + (S(3)/2 - x/2)*\
#             f(x)**2 + S(1)/6 + 5/(9*x - 9))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), -2*x/9 - S(1)/3)]

#     q = Eq(f(x).diff(x), (-x**6 + 15*x**4 - 40*x**3 + 45*x**2 - 24*x + 4)/\
#        (x**12 - 12*x**11 + 66*x**10 - 220*x**9 + 495*x**8 - 792*x**7 + 924*x**6 - \
#        792*x**5 + 495*x**4 - 220*x**3 + 66*x**2 - 12*x + 1) + f(x)**2 + f(x))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 1/(x**6 - 6*x**5 + 15*x**4 - 20*x**3 + 15*x**2 - 6*x + 1))]

#     eq = Eq(f(x).diff(x), x**3/4 - x**2/12 - 5*x/9 + (-6*x - 2)*f(x)/(9*x \
#         - 3) + (-x - 1)*f(x)**2 - S(7)/18 + 2/(27*x - 9))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), S(1)/3 - x/2)]

#     eq = Eq(f(x).diff(x), -3*x**3/4 - 3*x**2/4 - 3*x/4 + 3*x*f(x)/(2*x + \
#         -) + (3*x - 1)*f(x)**2 + S(85)/36 - 21/(4*x + 12))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), x/2 + S(1)/3)]

#     eq = Eq(f(x).diff(x), 2*x**3/3 + 8*x**2/3 + 10*x/3 + (-3*x/2 - 3)*\
#         -(x)**2 - f(x)/x - 2/(3*x))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), -2*x/3 - S(2)/3), Eq(f(x), -2*x/3 - S(2)/3)]

#     eq = Eq(f(x).diff(x), 18*x**3 + 18*x**2 + (-x/2 - S(1)/2)*f(x)**2 + 6),
#     sol = [Eq(f(x), 6*x)]

#     eq = Eq(f(x).diff(x), -3*x**3/4 + 15*x/2 + (x/3 - S(4)/3)*f(x)**2 \
#             + 9 + (1 - x)*f(x)/x + 3/x)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), -3*x/2 - 3)]

#     eq = Eq(f(x).diff(x), x**3/4 - 9*x**2/8 + 7*x/4 + (2 - 4*x)*f(x)**2\
#             - S(3)/4 + 1/(2*x - 6) + (x - 1)*f(x)/(x - 3))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), S(1)/2 - x/4)]

#     eq = Eq(f(x).diff(x), x**3/27 - 4*x**2/27 + 37*x/54 + (S(2)/3 - x/3)\
#             *f(x)**2 + (-9*x - 2)*f(x)/(6*x + 12) - S(61)/54 + 8/(3*x + 6))
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), x/3 - S(1)/3)]

#     eq = Eq(f(x).diff(x), -f(x)**2 - 2/(x**3 - x**2)),
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), 1/(x*(x - 1))), Eq(f(x), 1/(x*(x - 1)))]

# @slow
# def test_solve_riccati_slow():
#     eq = f(x).diff(x) + (3*x**2 + 1)*f(x)**2/x + (6*x**2 - x + 3)*f(x)/(x*(x \
#         - 1)) + (3*x**2 - 2*x + 2)/(x*(x - 1)**2)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), (-C1 - x**3 + x**2 - 2*x)/(C1*x - C1 + x**4 - x**3 + x**2 - x)), \
#         Eq(f(x), -1/(x - 1))]

#     eq = Eq(f(x).diff(x), (1 - x)*f(x)/(x - 3) + (2 - 12*x)*f(x)**2/(2*x - 9) + \
#         (54924*x**3 - 405264*x**2 + 1084347*x - 1087533)/(8*x**4 - 132*x**3 + 810*x**2 - \
#         2187*x + 2187) + 495)
#     _, funcs = match_riccati(eq, f, x)
#     sols = solve_riccati(f(x), x, *funcs)
#     solse = [Eq(f(x), Mul(6, 3*x + 1, evaluate=False)/(2*x - 9)), Eq(f(x), Mul(6, 3*x + \
#         1, evaluate=False)/(2*x - 9))]
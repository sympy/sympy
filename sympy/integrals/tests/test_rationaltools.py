from sympy import symbols, S, I, atan, log, Poly

from sympy.integrals.rationaltools import ratint, \
    ratint_ratpart, ratint_logpart, log_to_atan, log_to_real

from sympy.abc import a, b, x, t

half = S(1)/2

def test_ratint():
    assert ratint(S(0), x) == 0
    assert ratint(S(7), x) == 7*x

    assert ratint(x, x) == x**2/2
    assert ratint(2*x, x) == x**2

    assert ratint(8*x**7+2*x+1, x) == x**8+x**2+x

    f = S(1)
    g = x + 1

    assert ratint(f / g, x) == log(x + 1)
    assert ratint((f,g), x) == log(x + 1)

    f = x**3 - x
    g = x - 1

    assert ratint(f/g, x) == x**3/3 + x**2/2

    f = x
    g = (x - a)*(x + a)

    assert ratint(f/g, x) == log(x**2 - a**2)/2

    f = S(1)
    g = x**2 + 1

    assert ratint(f/g, x, real=None) == atan(x)
    assert ratint(f/g, x, real=True) == atan(x)

    assert ratint(f/g, x, real=False) == I*log(x + I)/2 - I*log(x - I)/2

    f = S(36)
    g = x**5-2*x**4-2*x**3+4*x**2+x-2

    assert ratint(f/g, x) == \
        -4*log(1 + x) + 4*log(-2 + x) - (6 + 12*x)/(1 - x**2)

    f = x**4-3*x**2+6
    g = x**6-5*x**4+5*x**2+4

    assert ratint(f/g, x) == \
        atan(x) + atan(x**3) + atan(x/2 - 3*x**S(3)/2 + S(1)/2*x**5)

    f = x**7-24*x**4-4*x**2+8*x-8
    g = x**8+6*x**6+12*x**4+8*x**2

    assert ratint(f/g, x) == \
        (4 + 6*x + 8*x**2 + 3*x**3)/(4*x + 4*x**3 + x**5) + log(x)

    assert ratint((x**3*f)/(x*g), x) == \
        -(12 - 16*x + 6*x**2 - 14*x**3)/(4 + 4*x**2 + x**4) - \
        5*2**(S(1)/2)*atan(x*2**(S(1)/2)/2) + S(1)/2*x**2 - 3*log(2 + x**2)

    f = x**5-x**4+4*x**3+x**2-x+5
    g = x**4-2*x**3+5*x**2-4*x+4

    assert ratint(f/g, x) == \
        x + S(1)/2*x**2 + S(1)/2*log(2-x+x**2) + (9-4*x)/(14-7*x+7*x**2) + \
        13*7**(S(1)/2)*atan(-S(1)/7*7**(S(1)/2) + 2*x*7**(S(1)/2)/7)/49

    assert ratint(1/(x**2+x+1), x) == \
        2*3**(S(1)/2)*atan(3**(S(1)/2)/3 + 2*x*3**(S(1)/2)/3)/3

    assert ratint(1/(x**3+1), x) == \
        -log(1 - x + x**2)/6 + log(1 + x)/3 + 3**(S(1)/2)*atan(-3**(S(1)/2)/3 + 2*x*3**(S(1)/2)/3)/3

    assert ratint(1/(x**2+x+1), x, real=False) == \
        -I*3**half*log(half + x - half*I*3**half)/3 + \
        I*3**half*log(half + x + half*I*3**half)/3

    assert ratint(1/(x**3+1), x, real=False) == log(1 + x)/3 - \
        (S(1)/6 - I*3**half/6)*log(-half + x + I*3**half/2) - \
        (S(1)/6 + I*3**half/6)*log(-half + x - I*3**half/2)

    # Issue 1892
    assert ratint(1/(x*(a+b*x)**3), x) == \
        (a**(-6))**(S(1)/2)*log(x + (2*a*b - 2*b*a**4*(a**(-6))**(S(1)/2))/(4*b**2)) + \
        (3*a + 2*b*x)/(2*a**2*b**2*x**2 + 4*b*x*a**3 + 2*a**4) - \
        (a**(-6))**(S(1)/2)*log(x + (2*a*b + 2*b*a**4*(a**(-6))**(S(1)/2))/(4*b**2))

def test_ratint_logpart():
    assert ratint_logpart(x, x**2-9, x, t) == \
        [(Poly(x**2 - 9, x), Poly(-2*t + 1, t))]
    assert ratint_logpart(x**2, x**3-5, x, t) == \
        [(Poly(x**3 - 5, x), Poly(-3*t + 1, t))]


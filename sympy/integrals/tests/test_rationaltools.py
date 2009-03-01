
from sympy import symbols, S, I, atan, log

from sympy.integrals.rationaltools import ratint, \
    ratint_ratpart, ratint_logpart, log_to_atan, log_to_real

x, t = symbols('x t')

def test_ratint():
    f = S(1)
    g = x + 1

    assert ratint(f / g, x) == log(x + 1)
    assert ratint((f,g), x) == log(x + 1)

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
        x + S(1)/2*x**2 + S(1)/2*log(2-x+x**2) + (S(9)/7-4*x/7)/(2-x+x**2) + \
        13*7**(S(1)/2)*atan(-S(1)/7*7**(S(1)/2) + 2*x*7**(S(1)/2)/7)/49


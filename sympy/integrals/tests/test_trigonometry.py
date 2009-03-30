from sympy import Symbol, Rational, sin, cos
from sympy.integrals.trigonometry import trigintegrate

x = Symbol('x')
y = Symbol('y')


def test_trigintegrate_odd():
    assert trigintegrate(Rational(1), x)    == x
    assert trigintegrate(x, x)      is None
    assert trigintegrate(x**2, x)   is None

    assert trigintegrate(sin(x), x) == -cos(x)
    assert trigintegrate(cos(x), x) ==  sin(x)

    assert trigintegrate(sin(3*x), x) == -cos(3*x)/3
    assert trigintegrate(cos(3*x), x) ==  sin(3*x)/3

    assert trigintegrate(sin(y*x), x) == -cos(y*x)/y
    assert trigintegrate(cos(y*x), x) ==  sin(y*x)/y

    assert trigintegrate(sin(x)*cos(x), x)      ==  sin(x)**2/2
    assert trigintegrate(sin(x)*cos(x)**2, x)   == -cos(x)**3/3
    assert trigintegrate(sin(x)**2*cos(x), x)   ==  sin(x)**3/3

    # check if it selects right function to substitute,
    # so the result is kept simple
    assert trigintegrate(sin(x)**7 * cos(x), x) ==  sin(x)**8/8
    assert trigintegrate(sin(x) * cos(x)**7, x) == -cos(x)**8/8

    assert trigintegrate(sin(x)**7 * cos(x)**3, x)  == -sin(x)**10/10 + sin(x)**8/8
    assert trigintegrate(sin(x)**3 * cos(x)**7, x)  ==  cos(x)**10/10 - cos(x)**8/8


def test_trigintegrate_even():
    assert trigintegrate(sin(x)**2, x)  == x/2 - cos(x)*sin(x)/2
    assert trigintegrate(cos(x)**2, x)  == x/2 + cos(x)*sin(x)/2

    assert trigintegrate(sin(3*x)**2, x)== x/2 - cos(3*x)*sin(3*x)/6
    assert trigintegrate(cos(3*x)**2, x)== x/2 + cos(3*x)*sin(3*x)/6
    assert trigintegrate(sin(x)**2 * cos(x)**2, x) == x/8 - cos(2*x)*sin(2*x)/16

    assert trigintegrate(sin(x)**4 * cos(x)**2, x) == x/16- sin(x)   *cos(x)/16 \
                                                          - sin(x)**3*cos(x)/24 \
                                                          + sin(x)**5*cos(x)/6

    assert trigintegrate(sin(x)**2 * cos(x)**4, x) == x/16+ cos(x)   *sin(x)/16 \
                                                          + cos(x)**3*sin(x)/24 \
                                                          - cos(x)**5*sin(x)/6

    assert trigintegrate(sin(x)**(-4),x) == -2*cos(x)/(3*sin(x)) \
                                            - cos(x)/(3*sin(x)**3)

    assert trigintegrate(cos(x)**(-6),x) == sin(x)/(5*cos(x)**5)\
                                            + 4*sin(x)/(15*cos(x)**3)\
                                            + 8*sin(x)/(15*cos(x))



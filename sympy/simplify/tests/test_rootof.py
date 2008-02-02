from sympy import Rational, symbols, I, sin, cos, pi
from sympy.simplify.rootof import roots, RootOf

def test_roots():
    a,b,c,x = symbols('abcx')

    assert roots(x**2-3*x+2, x) == [1, 2]
    assert roots(x**2-3*x/2+Rational(1,2), x) == [1, Rational(1,2)]
    assert roots(2*x**2-3*x+1, x) == [1, Rational(1,2)]
    assert roots(x**2-1, x) == [-1, 1]
    assert roots(x**2+1, x) == [I, -I]
    assert roots(x**3-1, x) == [1,
                             Rational(-1,2) + I*Rational(1,2)*3**Rational(1,2),
                             Rational(-1,2) - I*Rational(1,2)*3**Rational(1,2)]
    assert roots(x**3, x) == [0, 0, 0]
    assert roots(x**3-x, x) == [-1, 0, 1]
    assert roots(Rational(2), x) == []
    assert roots(a*x**2 + b*x + c, x) == \
           [-b/(a*2)+((b**2-4*a*c)**Rational(1,2))/2/a,
            -b/(a*2)-((b**2-4*a*c)**Rational(1,2))/2/a]
    assert roots(x**3 + x**2 + x + 1, x) == [-1, I, -I]
    assert roots(x**4 - 1, x) == [-1, 1, I, -I]
    assert roots(x**4 + 1, x) == [((-1)**Rational(1,4)).expand(complex=True),
                               ((-1)**Rational(3,4)).expand(complex=True),
                               (-(-1)**Rational(1,4)).expand(complex=True),
                               (-(-1)**Rational(3,4)).expand(complex=True)]
    assert roots(x**8 - 1, x) == [-1, 1,
        2**Rational(1,2)/2 + I*2**Rational(1,2)/2,
        I,
        -2**Rational(1,2)/2 + I*2**Rational(1,2)/2,
        -2**Rational(1,2)/2 - I*2**Rational(1,2)/2,
        -I,
        2**Rational(1,2)/2 - I*2**Rational(1,2)/2]
    assert roots(x**5 - Rational(3,2), x) == \
           [Rational(1,2)**Rational(1,5)*3**Rational(1,5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*cos(2*pi/5)
            + I*Rational(1,2)**Rational(1,5)*3**Rational(1,5)*sin(2*pi/5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*cos(4*pi/5)
            + I*Rational(1,2)**Rational(1,5)*3**Rational(1,5)*sin(4*pi/5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*cos(6*pi/5)
            + I*Rational(1,2)**Rational(1,5)*3**Rational(1,5)*sin(6*pi/5),
            Rational(1,2)**Rational(1,5)*3**Rational(1,5)*cos(8*pi/5)
            + I*Rational(1,2)**Rational(1,5)*3**Rational(1,5)*sin(8*pi/5)]

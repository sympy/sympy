from sympy import sqrt, Rational, S, sqrtdenest, Integral, cos
from sympy.simplify.sqrtdenest import _denester
from sympy.utilities.pytest import XFAIL

r2, r3, r5, r7, r29 = sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(29)

def test_sqrtdenest():
    d = {sqrt(5 + 2 * sqrt(6)): sqrt(2) + sqrt(3),
        sqrt(sqrt(2)): sqrt(sqrt(2)),
        sqrt(5+sqrt(7)): sqrt(5+sqrt(7)),
        sqrt(3+sqrt(5+2*sqrt(7))):
            3*sqrt(2)*(5 + 2*sqrt(7))**(S(1)/4)/(2*sqrt(6 + 3*sqrt(7))) +
            sqrt(2)*sqrt(6 + 3*sqrt(7))/(2*(5 + 2*sqrt(7))**(S(1)/4)),
        sqrt(3+2*sqrt(3)): 3**(S(3)/4)*(sqrt(6)/2 + 3*sqrt(2)/2)/3}
    for i in d:
        assert sqrtdenest(i) == d[i]

def test_sqrtdenest2():
    assert sqrtdenest(sqrt(16-2*sqrt(29)+2*sqrt(55-10*sqrt(29)))) == \
            sqrt(5) + sqrt(11-2*sqrt(29))
    assert sqrtdenest(sqrt(-r5+sqrt(-2*r29+2*sqrt(-10*r29+55)+16))) == \
    (-2*sqrt(29) + 11)**Rational(1,4)
    assert sqrtdenest(sqrt(1+sqrt(1+sqrt(7)))) == sqrt(1+sqrt(1+sqrt(7)))
    assert sqrtdenest(sqrt(((1+sqrt(1+2*sqrt(3+sqrt(2)+sqrt(5))))**2).expand())) == \
        1 + sqrt(1 + 2*sqrt(sqrt(2) + sqrt(5) + 3))
    assert sqrtdenest(sqrt(5*sqrt(3) + 6*sqrt(2))) == \
        sqrt(2)*3**Rational(1,4) + 3**Rational(3,4)

    assert sqrtdenest(sqrt(((1+sqrt(5)+sqrt(1+sqrt(3)))**2).expand())) == \
         1+sqrt(5)+sqrt(1+sqrt(3))

    assert sqrtdenest(sqrt(((1+r5+r7+sqrt(1+r3))**2).expand())) == \
        1+sqrt(1+r3)+r5+r7

    assert sqrtdenest(sqrt(((1+cos(2)+cos(3)+sqrt(1+r3))**2).expand())) == \
                cos(3)+cos(2)+1+sqrt(1+r3)

    assert sqrtdenest(sqrt(-2*sqrt(10)+2*r2*sqrt(-2*sqrt(10)+11)+14)) == \
        sqrt(-2*sqrt(10)+2*r2*(-1+sqrt(10))+14)

    # currencty cannot denest this; check one does not get a wrong answer
    z = sqrt(8 - sqrt(2)*sqrt(5-sqrt(5)) - 3*(1+sqrt(5)))
    assert (sqrtdenest(z) - z).evalf() < 1.0e-100

    # check that the result is not more complicated than the input
    z= sqrt(-2*sqrt(29) + cos(2) + 2*sqrt(-10*sqrt(29) + 55) + 16)
    assert sqrtdenest(z) == z

    assert sqrtdenest(sqrt(sqrt(6) + sqrt(15))) == sqrt(sqrt(6) + sqrt(15))

    # no assertion error when the 'r's are not the same in _denester
    z = sqrt(15 - 2*sqrt(31) + 2*sqrt(55 - 10*sqrt(29)))
    assert sqrtdenest(z) == z

def test_sqrtdenest_four_terms():
    assert sqrtdenest(sqrt(-4*sqrt(14)-2*sqrt(6)+4*sqrt(21)+33)) == \
      -sqrt(2) + sqrt(3) + 2*sqrt(7)
    assert sqrtdenest(sqrt(-28*sqrt(7)-14*sqrt(5)+4*sqrt(35) + 82)) == \
      -7 + sqrt(5) + 2*sqrt(7)

@XFAIL
def test_sqrtdenest2a():
    assert sqrtdenest(sqrt(13-2*sqrt(10)+2*sqrt(2)*sqrt(-2*sqrt(10)+11))) == \
            -1 + sqrt(2) + sqrt(10)

def test_issue_2758():
    from sympy.abc import x, y
    z = sqrt(1/(4*sqrt(3) + 7) + 1)
    ans = (sqrt(2) + sqrt(6))/(sqrt(3) + 2)
    assert sqrtdenest(z) == ans
    assert sqrtdenest(1 + z) == 1 + ans
    assert sqrtdenest(Integral(z + 1, (x, 1, 2))) == Integral(1 + ans, (x, 1, 2))
    assert sqrtdenest(x + sqrt(y)) == x + sqrt(y)

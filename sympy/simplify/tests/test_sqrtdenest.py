from sympy import sqrt, Rational, sqrtdenest
from sympy.utilities.pytest import XFAIL

def test_sqrtdenest():
    d = {sqrt(5 + 2 * sqrt(6)): sqrt(2) + sqrt(3),
        sqrt(sqrt(2)): sqrt(sqrt(2)),
        sqrt(5+sqrt(7)): sqrt(5+sqrt(7)),
        sqrt(3+sqrt(5+2*sqrt(7))):
            sqrt(6+3*sqrt(7))/(sqrt(2)*(5+2*sqrt(7))**Rational(1,4)) +
            3*(5+2*sqrt(7))**Rational(1,4)/(sqrt(2)*sqrt(6+3*sqrt(7))),
        sqrt(3+2*sqrt(3)): 3**Rational(1,4)/sqrt(2)+3/(sqrt(2)*3**Rational(1,4))}
    for i in d:
        assert sqrtdenest(i) == d[i]

# more complex example:
@XFAIL # this fails on amd64
def test_sqrtdenest2():
    assert sqrtdenest(sqrt(16-2*sqrt(29)+2*sqrt(55-10*sqrt(29)))) == \
            sqrt(5) + sqrt(11-2*sqrt(29))

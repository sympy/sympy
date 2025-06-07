from pytest import raises

from sympy import QQ, RR
from sympy.abc import x
from sympy.polys.domains.fpsring import PowerSeriesElement, PowerSeriesRing
from sympy.polys.orderings import lex
from sympy.polys.polyclasses import DMP_Python
from sympy.polys.polyerrors import GeneratorsError

def test_PowerSeriesRing__init__():
    R = PowerSeriesRing(QQ, 'x', 10)

    assert R.is_PowerSeriesRing
    assert R.domain == QQ
    assert R.symbol == (x,)
    assert R.prec == 10
    assert R.order == lex
    assert type(R) is PowerSeriesRing
    assert str(R.gen) == 'x + O(x**10)'
    assert type(R.gen) is PowerSeriesElement
    assert R.ring == 'Power Series Ring in x over QQ with lex order'

    raises (ValueError, lambda: PowerSeriesRing(QQ, 'x', -1))
    raises (GeneratorsError, lambda: PowerSeriesRing(QQ, 'x y'))


def test_PowerSeriesRing__hash__():
    R = PowerSeriesRing(QQ, 'x')
    assert hash(PowerSeriesRing(QQ, 'x', 6)) == hash(R)
    assert hash(PowerSeriesRing(QQ, 'x', 7)) != hash(R)


def test_PowerSeriesRing__eq__():
    assert PowerSeriesRing(QQ, 'x', 10) == PowerSeriesRing(QQ, [x], 10)
    assert PowerSeriesRing(QQ, x) == PowerSeriesRing(QQ, (x,), 6)
    assert PowerSeriesRing(QQ, 'x', 10) != PowerSeriesRing(QQ, 'y', 10)
    assert PowerSeriesRing(QQ, 'x') != PowerSeriesRing(RR, 'x')

def test_PowerSeriesRing_new():
    R = PowerSeriesRing(QQ, 'x', 10)

    assert R.ground_new(2) == R.from_list([2])
    assert R.term_new(2, 4) == R.from_list([2, 0, 0, 0])


def test_PowerSeriesRing_from_list():
    R = PowerSeriesRing(QQ, 'x', 8)
    x = R.gen
    O = R.order_term()

    assert R.from_list([3,2,1]) == 1 + 2*x + 3*x**2 + O
    assert R.from_list([1,0,1,0,1,0]) == x + x**3 + x**5 + O


def test_PowerSeriesElement__init__():
    R = PowerSeriesRing(QQ, 'x', 10)
    x = R.gen
    O = R.order_term()

    assert type(x) is PowerSeriesElement
    assert type(O) is PowerSeriesElement

    assert type(x.poly) is DMP_Python
    assert x.poly == DMP_Python([1, 0], QQ)
    assert x.ring == R
    assert O.ring == R

def test_PowerSeriesElement_add():
    R = PowerSeriesRing(QQ, 'x', 8)
    x = R.gen
    O = R.order_term()

    f = 1 + 2*x + 3*x**2
    g = 1 + 0*x + 5*x**2

    assert f + g == 2 + 2*x + 8*x**2 + O
    assert f + 3 == 4 + 2*x + 3*x**2 + O
    assert 3 + f == 4 + 2*x + 3*x**2 + O


def test_PowerSeriesElement_sub():
    R = PowerSeriesRing(QQ, 'x', 8)
    x = R.gen
    O = R.order_term()

    f = 1 + 2*x + 3*x**2
    g = 1 + 5*x**2

    assert f - g == 0 + 2*x - 2*x**2 + O
    assert f - 3 == -2 + 2*x + 3*x**2 + O
    assert 3 - f == 2 - 2*x - 3*x**2 + O


def test_PowerSeriesElement_mul():
    R = PowerSeriesRing(QQ, 'x', 8)
    x = R.gen
    O = R.order_term()

    f = R.from_list([3, 2, 1])
    g = R.from_list([5, 0, 1])

    assert f * g == 1 + 2*x + 8*x**2 + 10*x**3 + 15*x**4 + O
    assert x * x == x**2 + O
    assert f * 2 == 2 + 4*x + 6*x**2 + O
    assert 2 * f == 2 + 4*x + 6*x**2 + O


def test_PowerSeriesElement_pow():
    R = PowerSeriesRing(QQ, 'x', 8)
    x = R.gen
    O = R.order_term()

    # Basic exponentiation
    assert x**0 == 1 + O
    assert x**1 == x + O
    assert x**2 == x**2 + O
    assert x**3 == x**3 + O

    # More complex expression
    f = R.from_list([0, 1, 1])
    assert f**2 == 1 + 2*x + x**2 + O
    assert f**3 == 1 + 3*x + 3*x**2 + x**3 + O

    # Zero exponent
    g = R.from_list([3, 2, 1])
    assert g**0 == 1 + O

    raises(ValueError, lambda: x**-1)
    raises(ValueError, lambda: x**0.5)


def test_basics():
    R = PowerSeriesRing(QQ, 'x', 8)
    x = R.gen
    O = R.order_term()

    f = 7 + x**2 + x**4 + x

    assert f.get_constant_term() == QQ(7, 1)
    assert f.truncate(3) == 7 + x + x**2 + O

from sympy.abc import x, y
from sympy.core.evaluate import evaluate
from sympy.core import Mul, Add, Pow, S
from sympy import sqrt

def test_add():
    with evaluate(False):
        expr = x + x
        assert isinstance(expr, Add)
        assert expr.args == (x, x)

        with evaluate(True):
            assert (x + x).args == (2, x)

        assert (x + x).args == (x, x)

    assert isinstance(x + x, Mul)

    with evaluate(False):
        assert S(1) + 1 == Add(1, 1)
        assert S(4) - 3 == Add(4, -3)
        assert S(2) * 4 == Mul(2, 4)
        assert S(6) / 3 == Mul(6, S(1)/ 3)
        assert S(2) ** 9 == Pow(2, 9)
        assert S(2) / 2 == 2*S(1) / 2

        assert S(2) / 3 + 1 == Add(S(2) / 3, 1)
        assert S(4) / 7 - 3 == Add(S(4) / 7, -3)
        assert S(2) / 4 * 4 == Mul(S(2) / 4, 4)
        assert S(6) / 3 == Mul(6, S(1) / 3)
        assert S(2) ** 9 == Pow(S(2), 9)
        assert S(2) / 2 == 2*S(1) / 2

        assert S(1) / 3 + sqrt(3) == Add(S(1) / 3, sqrt(3))
        assert S(1) / 2 * 10.333 == Mul(S(1) / 2, 10.333)
        assert sqrt(2) * sqrt(2) == Mul(sqrt(2), sqrt(2))

        assert S(1) / 2 + x == Add(S(1) / 2, x)
        assert S(1) / x * x == Mul(S(1) / x, x)

def test_nested():
    with evaluate(False):
        expr = (x + x) + (y + y)
        assert expr.args == ((x + x), (y + y))
        assert expr.args[0].args == (x, x)

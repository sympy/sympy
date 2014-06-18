from sympy import sin, exp, Derivative
from sympy.core import Function
from sympy.functions import airyai
from sympy.series.formal import FormalSeries, simpleDE, Stream
from sympy.abc import x

def test_simpleDE():
    f = Function('f')

    assert simpleDE(x, sin(x), f(x)) == f(x) + Derivative(f(x), x, x)
    assert simpleDE(x, exp(x), f(x)) == -f(x) + Derivative(f(x), x)
    assert simpleDE(x, airyai(x), f(x)) == -x*f(x) + Derivative(f(x), x, x)


def test_DEtoRE():


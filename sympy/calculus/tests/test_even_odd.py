from sympy.calculus.even_odd import is_even_function, is_odd_function
from sympy.core.symbol import Symbol
from sympy.core.singleton import S
from sympy.functions import Abs, sign, sin, cos, exp, sinh, cosh, tanh, erf, fresnelc
from sympy.core.function import expand_trig


def test_is_even_odd_wolfram():
    # Test problems from
    # https://resources.wolframcloud.com/FunctionRepository/resources/EvenFunctionQ/
    # https://resources.wolframcloud.com/FunctionRepository/resources/OddFunctionQ/
    # https://resources.wolframcloud.com/FunctionRepository/resources/FunctionParity/

    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert is_even_function(0, x) is True
    assert is_odd_function(0, x) is True
    assert is_even_function(x**2, x) is True
    assert is_odd_function(x**3, x) is True
    assert is_even_function(1, x) is True
    assert is_even_function(Abs(x), x) is True
    assert is_odd_function(sin(x), x) is True
    assert is_even_function(cos(x), x) is True
    assert is_even_function(exp(-x**2), x) is True
    assert is_odd_function(sign(x) * exp(-x**2), x) is True
    assert is_even_function(cosh(x + y), x, y) is True
    assert is_odd_function(sinh(x + y), x, y) is True
    assert is_odd_function(tanh(x), x) is True
    assert is_odd_function(erf(x), x) is True
    assert is_odd_function(fresnelc(x), x) is True
    assert is_odd_function(sin(x * y * z), x, y, z) is True
    assert is_even_function(sin(x * y * Abs(z)), x, y) is True
    assert is_even_function(x**2 + y**2, x, y) is True

    assert is_odd_function(x / (x**2 + 1), x) is True

    # XXX It's not possible to decide even function for the pattern
    # if not simplified
    assert is_even_function(sin(x + S.Pi/4) + cos(x + S.Pi/4), x) is None
    assert is_even_function(expand_trig(sin(x + S.Pi/4) + cos(x + S.Pi/4)), x) is True

    # XXX May return False instead of None
    assert is_even_function(fresnelc(x) + 1, x) is None
    assert is_odd_function(fresnelc(x) + 1, x) is None
    assert is_even_function(x + 1, x) is None
    assert is_odd_function(x + 1, x) is None

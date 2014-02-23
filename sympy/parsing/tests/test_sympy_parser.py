from sympy.core import Symbol, Float, Rational, Integer, I, Mul, Pow
from sympy.functions import exp, factorial, sin
from sympy.logic import And
from sympy.series import Limit
from sympy.utilities.pytest import raises

from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations, rationalize, TokenError
)


def test_sympy_parser():
    x = Symbol('x')
    inputs = {
        '2*x': 2 * x,
        '3.00': Float(3),
        '22/7': Rational(22, 7),
        '2+3j': 2 + 3*I,
        'exp(x)': exp(x),
        'x!': factorial(x),
        '3.[3]': Rational(10, 3),
        '10!': 3628800,
        '-(2)': -Integer(2),
        '[-1, -2, 3]': [Integer(-1), Integer(-2), Integer(3)],
        'Symbol("x").free_symbols': x.free_symbols,
        "S('S(3).n(n=3)')": 3.00,
        'factorint(12, visual=True)': Mul(
            Pow(2, 2, evaluate=False),
            Pow(3, 1, evaluate=False),
            evaluate=False),
        'Limit(sin(x), x, 0, dir="-")': Limit(sin(x), x, 0, dir='-'),

    }
    for text, result in inputs.items():
        assert parse_expr(text) == result


def test_rationalize():
    inputs = {
        '0.123': Rational(123, 1000)
    }
    transformations = standard_transformations + (rationalize,)
    for text, result in inputs.items():
        assert parse_expr(text, transformations=transformations) == result


def test_factorial_fail():
    inputs = ['x!!!', 'x!!!!', '(!)']

    for text in inputs:
        try:
            parse_expr(text)
            assert False
        except TokenError:
            assert True

def test_local_dict():
    local_dict = {
        'my_function': lambda x: x + 2
    }
    inputs = {
        'my_function(2)': Integer(4)
    }
    for text, result in inputs.items():
        assert parse_expr(text, local_dict=local_dict) == result

def test_global_dict():
    global_dict = {
        'Symbol': Symbol
    }
    inputs = {
        'Q & S': And(Symbol('Q'), Symbol('S'))
    }
    for text, result in inputs.items():
        assert parse_expr(text, global_dict=global_dict) == result

def test_issue_2515_github():
    raises(TokenError, lambda: parse_expr('(()'))
    raises(TokenError, lambda: parse_expr('"""'))

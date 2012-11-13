from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations, rationalize, TokenError
)
from sympy import *

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
        '(_kern)': Symbol('_kern'),
        '-(2)': -Integer(2),
        '[-1, -2, 3]': [Integer(-1), Integer(-2), Integer(3)]
    }
    for text, result in inputs.items():
        assert parse_expr(text) == result

def test_rationalize():
    inputs = {
        '0.123': Rational(123,1000)
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

# -*- coding: utf-8 -*-


import sys
from random import choice


from sympy.core.function import Function
from sympy.core.mul import Mul
from sympy.core.numbers import (Float, I, Integer, Rational)
from sympy.core.power import Pow
from sympy.core.relational import Eq
from sympy.core.symbol import (Symbol, symbols)
from sympy.core.compatibility import PY3
from sympy.functions import exp, factorial, factorial2, sin, log
from sympy.logic import And
from sympy.series import Limit
from sympy.utilities.iterables import subsets
from sympy.utilities.pytest import raises, skip

from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations, rationalize, TokenError,
    split_symbols, implicit_multiplication, convert_equals_signs,
    convert_xor, function_exponentiation,
    implicit_multiplication_application, valid_name,
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
        'x!!': factorial2(x),
        '(x + 1)! - 1': factorial(x + 1) - 1,
        '3.[3]': Rational(10, 3),
        '.0[3]': Rational(1, 30),
        '3.2[3]': Rational(97, 30),
        '1.3[12]': Rational(433, 330),
        '1 + 3.[3]': Rational(13, 3),
        '1 + .0[3]': Rational(31, 30),
        '1 + 3.2[3]': Rational(127, 30),
        '.[0011]': Rational(1, 909),
        '0.1[00102] + 1': Rational(366697, 333330),
        '1.[0191]': Rational(10190, 9999),
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
        assert parse_expr(text) == result, (text, result)

    raises(TypeError, lambda:
        parse_expr('x', standard_transformations))
    raises(TypeError, lambda:
        parse_expr('x', transformations=lambda x,y: 1))
    raises(TypeError, lambda:
        parse_expr('x', transformations=(lambda x,y: 1,)))
    raises(TypeError, lambda: parse_expr('x', transformations=((),)))
    raises(TypeError, lambda: parse_expr('x', {}, [], []))
    raises(TypeError, lambda: parse_expr('x', [], [], {}))
    raises(TypeError, lambda: parse_expr('x', [], [], {}))


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


def test_repeated_fail():
    # All are valid Python, so only raise TypeError for invalid indexing
    inputs = ['1[1]', '.1e1[1]', '0x1[1]', '1.1j[1]', '1.1[1 + 1]',
        '0.1[[1]]', '0x1.1[1]']
    for text in inputs:
        raises(TypeError, lambda: parse_expr(text, pre=''))
    inputs = ['0.1[]']
    for text in inputs:
        raises(SyntaxError, lambda: parse_expr(text, pre=''))

    # these are not valid Python
    inputs = ['0.1[', '0.1[1']
    for text in inputs:
        raises(TokenError, lambda: parse_expr(text, pre=''))


def test_repeated_dot_only():
    assert parse_expr('.[1]') == Rational(1, 9)
    assert parse_expr('1 + .[1]') == Rational(10, 9)


def test_local_dict():
    local_dict = {
        'my_function': lambda x: x + 2
    }
    inputs = {
        'my_function(2)': Integer(4)
    }
    for text, result in inputs.items():
        assert parse_expr(text, local_dict=local_dict) == result


def test_local_dict_split_implmult():
    t = standard_transformations + (split_symbols, implicit_multiplication,)
    w = Symbol('w', real=True)
    y = Symbol('y')
    assert parse_expr('yx', local_dict={'x':w}, transformations=t) == y*w


def test_local_dict_symbol_to_fcn():
    x = Symbol('x')
    d = {'foo': Function('bar')}
    assert parse_expr('foo(x)', local_dict=d) == d['foo'](x)
    d = {'foo': Symbol('baz')}
    raises(TypeError, lambda: parse_expr('foo(x)', local_dict=d))


def test_global_dict():
    global_dict = {
        'Symbol': Symbol
    }
    inputs = {
        'Q & S': And(Symbol('Q'), Symbol('S'))
    }
    for text, result in inputs.items():
        assert parse_expr(text, global_dict=global_dict) == result


def test_issue_2515():
    raises(TokenError, lambda: parse_expr('(()'))
    raises(TokenError, lambda: parse_expr('"""'))


def test_issue_7663():
    x = Symbol('x')
    e = '2*(x+1)'
    assert parse_expr(e, evaluate=0) == parse_expr(e, evaluate=False)


def test_issue_10560():
    inputs = {
        '4*-3' : '(-3)*4',
        '-4*3' : '(-4)*3',
    }
    for text, result in inputs.items():
        assert parse_expr(text, evaluate=False) == parse_expr(result, evaluate=False)


def test_issue_10773():
    inputs = {
    '-10/5': '(-10)/5',
    '-10/-5' : '(-10)/(-5)',
    }
    for text, result in inputs.items():
        assert parse_expr(text, evaluate=False) == parse_expr(result, evaluate=False)


def test_split_symbols():
    transformations = standard_transformations + \
                      (split_symbols, implicit_multiplication,)
    x = Symbol('x')
    y = Symbol('y')
    xy = Symbol('xy')


    assert parse_expr("xy") == xy
    assert parse_expr("xy", transformations=transformations) == x*y


def test_split_symbols_function():
    transformations = standard_transformations + \
                      (split_symbols, implicit_multiplication,)
    x = Symbol('x')
    y = Symbol('y')
    a = Symbol('a')
    f = Function('f')


    assert parse_expr("ay(x+1)", transformations=transformations) == a*y*(x+1)
    assert parse_expr("af(x+1)", transformations=transformations,
                      local_dict={'f':f}) == a*f(x+1)


def test_functional_exponent():
    t = standard_transformations + (convert_xor, function_exponentiation)
    x = Symbol('x')
    y = Symbol('y')
    a = Symbol('a')
    yfcn = Function('y')
    assert parse_expr("sin^2(x)", transformations=t) == (sin(x))**2
    assert parse_expr("sin^y(x)", transformations=t) == (sin(x))**y
    assert parse_expr("exp^y(x)", transformations=t) == (exp(x))**y
    assert parse_expr("E^y(x)", transformations=t) == exp(yfcn(x))
    assert parse_expr("a^y(x)", transformations=t) == a**(yfcn(x))


def test_match_parentheses_implicit_multiplication():
    transformations = standard_transformations + \
                      (implicit_multiplication,)
    raises(TokenError, lambda: parse_expr('(1,2),(3,4]',transformations=transformations))


def test_convert_equals_signs():
    transformations = standard_transformations + \
                        (convert_equals_signs, )
    x = Symbol('x')
    y = Symbol('y')
    assert parse_expr("1*2=x", transformations=transformations) == Eq(2, x)
    assert parse_expr("y = x", transformations=transformations) == Eq(y, x)
    assert parse_expr("(2*y = x) = False",
        transformations=transformations) == Eq(Eq(2*y, x), False)


def test_parse_function_issue_3539():
    x = Symbol('x')
    f = Function('f')
    assert parse_expr('f(x)') == f(x)


def test_split_symbols_numeric():
    SI = standard_transformations + (
        implicit_multiplication_application,)
    do = lambda x: parse_expr(x, transformations=SI, pre='D')

    n = Symbol('n')
    assert parse_expr('2**n * 3**n') == do('2**n3**n') == 2**n*3**n
    assert do('n12n34') == n*12*n*34

    # issue 16632
    x, y, z, x_, x_3, j, jy, e4 = symbols('x y z x_ x_3 j jy e4')
    assert do('x3y') == x*3*y
    assert do('x3y') == 3*x*y
    assert do('x3j') == 3*j*x
    assert do('x3.j') == x*3.0*j, do('x3.j')
    assert parse_expr('x3.j', j=I) == 3.0*I*x
    assert do('x3.jy') == 3.0*x*j*y
    assert parse_expr('x3.jy', transformations=SI, j=I) == 3.0*x*I*y
    assert do('x_3e4yz') == Symbol('x_3e4yz'), do('x_3e4yz')
    assert do('x_3.e4yz') == 30000.0*x_*y*z
    assert do('x_.3e4yz') == 3000.0*x_*y*z
    assert do('x_3.3e4yz') == 33000.0*x_*y*z
    assert do('x_3yz') == Symbol('x_3yz')
    assert do('n1.1n22') == n*1.1*n*22
    assert parse_expr('log10(3.2x)', dict(log10=log),
        transformations=SI) == log(3.2*x)
    assert parse_expr('x3y',
        transformations=standard_transformations) == Symbol('x3y')
    assert parse_expr('2x1e2y', pre='e') == 200.0*x*y
    raises(SyntaxError, lambda: parse_expr('2x', pre=''))

    # check pre values
    # single number always recogized
    assert parse_expr('2e4', pre='d') == 2e4
    assert parse_expr('2e4', pre='D') == 2e4
    assert parse_expr('2E4', pre='D') == 2e4
    assert parse_expr('2e4', pre='e') == 2e4
    assert parse_expr('2e4', pre='E') == 2e4
    assert parse_expr('2e4', pre='') == 2e4
    e4 = Symbol('e4')
    E4 = Symbol('E4')
    # check sign collapse
    s = ''.join([choice('- +') for i in range(3)]) + '2'
    ans = parse_expr(s)
    assert ans == 2*(-1)**s.count('-'), '%s -> %s' % (s, ans)
    raises(SyntaxError, lambda: do('--'))

    assert parse_expr('3e4') == 3e4
    assert parse_expr('3e1', pre='e') == 30.0
    e1x = Symbol('e1x')
    assert parse_expr('2x') == 2*x
    assert parse_expr('2_3x') == 23*x
    _x = Symbol('_x')
    assert parse_expr('2_x') == 2*_x
    assert parse_expr('2x') == 2*x
    assert parse_expr('2e1x') == 2*e1x
    assert parse_expr('2e1x', pre='e') == 20.0*x
    assert parse_expr('2.e1x') == 2.0*e1x
    assert parse_expr('2.e1x', pre='e') == 20.0*x
    raises(SyntaxError, lambda: parse_expr('.e1x'))
    raises(SyntaxError, lambda: parse_expr('.e1x', pre='e'))
    E1x = Symbol('E1x')
    assert parse_expr('2.E1x') == 2.0*E1x
    assert parse_expr('2.E1x', pre='e') == 20.0*x


def test_unicode_names():
    if not PY3:
        skip("test_unicode_names can only pass in Python 3")


    assert parse_expr(u'α') == Symbol(u'α')


def test_python3_features():
    # Make sure the tokenizer can handle Python 3-only features
    if sys.version_info < (3, 6):
        skip("test_python3_features requires Python 3.6 or newer")


    assert parse_expr("123_456") == 123456
    assert parse_expr("1.2[3_4]") == parse_expr("1.2[34]") == Rational(611, 495)
    assert parse_expr("1.2[012_012]") == parse_expr("1.2[012012]") == Rational(400, 333)
    assert parse_expr('.[3_4]') == parse_expr('.[34]') == Rational(34, 99)
    assert parse_expr('.1[3_4]') == parse_expr('.1[34]') == Rational(133, 990)
    assert parse_expr('123_123.123_123[3_4]') == parse_expr('123123.123123[34]') == Rational(12189189189211, 99000000)


def test_valid_name():
    assert not valid_name(3)
    assert not valid_name('3')
    assert not valid_name(' x')
    assert not valid_name('a=b')
    assert not valid_name('for')
    assert valid_name('_ok2')

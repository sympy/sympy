# -*- coding: utf-8 -*-
import inspect
import itertools

import pytest
from multiset import Multiset

from sympy.pattern.expressions.substitution import Substitution

from sympy.pattern.expressions import (
    Arity, Operation, Symbol, Wildcard, SymbolWildcard, make_dot_variable, make_plus_variable, make_star_variable,
    make_symbol_variable
)

from multiset import Multiset

class SpecialSymbol(Symbol):
    pass

f = Operation.new('f', Arity.variadic)
f2 = Operation.new('f2', Arity.variadic)
f_u = Operation.new('f_u', Arity.unary)
f_i = Operation.new('f_i', Arity.variadic, one_identity=True)
f_c = Operation.new('f_c', Arity.variadic, commutative=True)
f2_c = Operation.new('f2_c', Arity.variadic, commutative=True)
f_a = Operation.new('f_a', Arity.variadic, associative=True)
f_ac = Operation.new('f_ac', Arity.variadic, associative=True, commutative=True)
a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')
s = SpecialSymbol('s')
_ = Wildcard.dot()
_s = Wildcard.symbol()
_ss = Wildcard.symbol(SpecialSymbol)
x_ = make_dot_variable('x')
s_ = make_symbol_variable('s')
ss_ = make_symbol_variable('ss', SpecialSymbol)
y_ = make_dot_variable('y')
z_ = make_dot_variable('z')
__ = Wildcard.plus()
x__ = make_plus_variable('x')
y__ = make_plus_variable('y')
z__ = make_plus_variable('z')
___ = Wildcard.star()
x___ = make_star_variable('x')
y___ = make_star_variable('y')
z___ = make_star_variable('z')


def test_union_with_var():
    test = [
        ({},                            'x',        a,                      {'x': a}),
        ({'x': a},                      'x',        a,                      {'x': a}),
        ({'x': a},                      'x',        b,                      ValueError),
        ({'x': a},                      'x',        (a, b),                 ValueError),
        ({'x': (a, )},                  'x',        a,                      ValueError),
        ({'x': (a, b)},                 'x',        (a, b),                 {'x': (a, b)}),
        ({'x': (a, b)},                 'x',        (a, a),                 ValueError),
        ({'x': (a, b)},                 'x',        Multiset([a, b]),       {'x': (a, b)}),
        ({'x': (a, b)},                 'x',        Multiset([a]),          ValueError),
        ({'x': Multiset([a, b])},       'x',        Multiset([a, b]),       {'x': Multiset([a, b])}),
        ({'x': Multiset([a, b])},       'x',        Multiset([]),           ValueError),
        ({'x': Multiset([a, b])},       'x',        (a, b),                 {'x': (a, b)}),
        ({'x': Multiset([a, b])},       'x',        (a, a),                 ValueError),
        ({'x': Multiset([a])},          'x',        (a, ),                  {'x': (a, )}),
        ({'x': Multiset([a])},          'x',        (b, ),                  ValueError),
        ({'x': Multiset([a])},          'x',        a,                      ValueError),
        ({'x': Multiset([a])},          'x',        b,                      ValueError),
    ]
    for substitution, variable, value, expected_result in test:
        substitution = Substitution(substitution)
        if expected_result is ValueError:
            with pytest.raises(ValueError):
                _ = substitution.union_with_variable(variable, value)
        else:
            result = substitution.union_with_variable(variable, value)
            assert result == expected_result


def test_union():
    test = [
        ({},                            {},                             {}),
        ({'x': a},                      {},                             {'x': a}),
        ({'x': a},                      {'y': b},                       {'x': a, 'y': b}),
        ({'x': a},                      {'x': b},                       ValueError),
        ({'x': a},                      {'x': a},                       {'x': a}),
    ]

    for substitution1, substitution2, expected_result in test:
        substitution1 = Substitution(substitution1)
        substitution2 = Substitution(substitution2)
        if expected_result is ValueError:
            with pytest.raises(ValueError):
                _ = substitution1.union(substitution2)
            with pytest.raises(ValueError):
                _ = substitution2.union(substitution1)
        else:
            result = substitution1.union(substitution2)
            assert result == expected_result
            assert result is not substitution1
            assert result is not substitution2
            result = substitution2.union(substitution1)
            assert result == expected_result
            assert result is not substitution1
            assert result is not substitution2


def test_extract_substitution():
    test = [
        ({},                            a,          a,                      {}),
        ({},                            a,          x_,                     {'x': a}),
        ({'x': a},                      a,          x_,                     {'x': a}),
        ({'x': b},                      a,          x_,                     False),
        ({},                            f(a),       f(a),                   {}),
        ({},                            f(a),       f(x_),                  {'x': a}),
        ({'x': a},                      f(a),       f(x_),                  {'x': a}),
        ({'x': b},                      f(a),       f(x_),                  False),
        ({},                            f(a, a),    f(x_, x_),              {'x': a}),
        ({},                            f(a, b),    f(x_, x_),              False),
        ({},                            f(a, b),    f(x_, y_),              {'x': a, 'y': b}),
    ]
    for substitution, subject, pattern, expected_result in test:
        substitution = Substitution(substitution)
        if expected_result is False:
            assert substitution.extract_substitution(subject, pattern) is False
        else:
            assert substitution.extract_substitution(subject, pattern) is True
            assert substitution == expected_result


def test_rename():
    test = [
        ({},                            {},                       {}),
        ({'x': a},                      {},                       {'x': a}),
        ({'x': a},                      {'x': 'y'},               {'y': a}),
        ({'x': a},                      {'y': 'x'},               {'x': a}),
    ]

    for substitution, renaming, expected_result in test:
        assert Substitution(substitution).rename(renaming) == expected_result

def test_copy():
    substitution = Substitution({'x': a})

    copy = substitution.__copy__()

    assert copy == substitution
    assert copy is not substitution

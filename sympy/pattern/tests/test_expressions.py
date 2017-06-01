from sympy.pattern.expressions import (
    Arity, Operation, Symbol, Wildcard, SymbolWildcard, make_dot_variable,
    make_plus_variable, make_star_variable, make_symbol_variable
)

from sympy.utilities import pytest
from multiset import Multiset
import inspect
from itertools import product

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


SIMPLE_EXPRESSIONS = [
    a,
    b,
    f(a, b),
    x_,
    ___,
    f(_, variable_name='x'),
    s_,
    _s,
]

class SpecialF(f):
    name = 'special'


def test_operation_simplify():
    test = [
        (f_i(a),                                                            a),
        (f_i(a, b),                                                         f_i(a, b)),
        (f_i(_),                                                            _),
        (f_i(___),                                                          f_i(___)),
        (f_i(__),                                                           f_i(__)),
        (f_i(x_),                                                           x_),
        (f_i(x___),                                                         f_i(x___)),
        (f_i(x__),                                                          f_i(x__)),
        (f_a(f_a(a)),                                                       f_a(a)),
        (f_a(f_a(a, b)),                                                    f_a(a, b)),
        (f_a(a, f_a(b)),                                                    f_a(a, b)),
        (f_a(f_a(a), b),                                                    f_a(a, b)),
        (f_a(f(a)),                                                         f_a(f(a))),
        (f_c(a, b),                                                         f_c(a, b)),
        (f_c(b, a),                                                         f_c(a, b)),
    ]

    for expression, simplified in test:
        assert expression == simplified


def test_operation_errors():
    test = [
        (Operation.new('f', Arity.unary),                       [],             ValueError),
        (Operation.new('f', Arity.unary),                       [a, b],         ValueError),
        (Operation.new('f', Arity.variadic),                    [],             None),
        (Operation.new('f', Arity.variadic),                    [a],            None),
        (Operation.new('f', Arity.variadic),                    [a, b],         None),
        (Operation.new('f', Arity.binary, associative=True),    [a, a, b],      ValueError),
        (Operation.new('f', Arity.binary),                      [x_, x___],     None),
        (Operation.new('f', Arity.binary),                      [x_, x__],      None),
        (Operation.new('f', Arity.binary),                      [x_, x_, x__],  ValueError),
        (Operation.new('f', Arity.binary),                      [x_, x_, x___], None),
        (Operation.new('f', Arity.binary),                      [x_, x_],       None),
        (Operation.new('f', Arity.binary),                      [x_, x_, x_],   ValueError),
    ]

    for operation, operands, expected_error in test:
        if expected_error is not None:
            with pytest.raises(expected_error):
                operation(*operands)
        else:
            _ = operation(*operands)


def test_is_constant():
    test = [
        (a,             True),
        (x_,            False),
        (_,             False),
        (f(a),          True),
        (f(a, b),       True),
        (f(x_),         False),
    ]

    for expression, is_constant in test:
        assert expression.is_constant == is_constant


def test_is_syntactic():
    test = [
        (a,             True),
        (x_,            True),
        (_,             True),
        (x___,          False),
        (___,           False),
        (x__,           False),
        (__,            False),
        (f(a),          True),
        (f(a, b),       True),
        (f(x_),         True),
        (f(x__),        False),
        (f_a(a),        False),
        (f_a(a, b),     False),
        (f_a(x_),       False),
        (f_a(x__),      False),
        (f_c(a),        False),
        (f_c(a, b),     False),
        (f_c(x_),       False),
        (f_c(x__),      False),
        (f_ac(a),       False),
        (f_ac(a, b),    False),
        (f_ac(x_),      False),
        (f_ac(x__),     False),
    ]

    for expression, is_syntactic in test:
        assert expression.is_syntactic == is_syntactic

def test_symbols():
    test = [
        (a,                 ['a']),
        (x_,                []),
        (_,                 []),
        (f(a),              ['a', 'f']),
        (f(a, b),           ['a', 'b', 'f']),
        (f(x_),             ['f']),
        (f(a, a),           ['a', 'a', 'f']),
        (f(f(a), f(b, c)),  ['a', 'b', 'c', 'f', 'f', 'f']),
    ]

    for expression, symbols in test:
        assert expression.symbols == Multiset(symbols)


def test_variables():
    test = [
        (a,                         []),
        (x_,                        ['x']),
        (_,                         []),
        (f(a),                      []),
        (f(x_),                     ['x']),
        (f(x_, x_),                 ['x', 'x']),
        (f(x_, a),                  ['x']),
        (f(x_, a, y_),              ['x', 'y']),
        (f(f(x_), f(b, x_)),        ['x', 'x']),
        (f(a, variable_name='x'),        ['x']),
        (f(f(y_), variable_name='x'),    ['x', 'y']),
    ]

    for expression, variables in test:
        assert expression.variables == Multiset(variables)

def test_preorder_iter():
    test = [
            (f(a, x_),      None,                       [(f(a, x_),         ()),
                                                         (a,                (0, )),
                                                         (x_,               (1, ))]),
            (f(a, f(x_)),   lambda e: e.head is None,   [(x_,               (1, 0))]),
            (f(a, f(x_)),   lambda e: e.head == f,      [(f(a, f(x_)),      ()),
                                                         (f(x_),            (1, ))])
    ]

    for expression, predicate, preorder_list in test:
        result = list(expression.preorder_iter(predicate))
        assert result == preorder_list

GETITEM_TEST_EXPRESSION = f(a, f(x_, b), _)

def test_getitem():
    test = [
        ((),            GETITEM_TEST_EXPRESSION),
        ((0, ),         a),
        ((0, 0),        IndexError),
        ((1, ),         f(x_, b)),
        ((1, 0),        x_),
        ((1, 0, 0),     IndexError),
        ((1, 1),        b),
        ((1, 1, 0),     IndexError),
        ((1, 2),        IndexError),
        ((2, ),         _),
        ((3, ),         IndexError),
    ]

    for position, expected_result in test:
        if inspect.isclass(expected_result) and issubclass(expected_result, Exception):
            with pytest.raises(expected_result):
                result = GETITEM_TEST_EXPRESSION[position]
                print(result)
        else:
            result = GETITEM_TEST_EXPRESSION[position]
            assert result == expected_result


def test_getitem_slice():
    test = [
        ((),            (),     [GETITEM_TEST_EXPRESSION]),
        ((0, ),         (0, ),  [a]),
        ((0, ),         (1, ),  [a, f(x_, b)]),
        ((0, ),         (2, ),  [a, f(x_, b), _]),
        ((0, ),         (3, ),  [a, f(x_, b), _]),
        ((1, ),         (2, ),  [f(x_, b), _]),
        ((1, 0),        (1, 1), [x_, b]),
        ((1, 0),        (2, ),  IndexError),
        ((1, ),         (0, ),  IndexError),
        ((1, 0),        (2, 0), IndexError),
    ]

    for start, end, expected_result in test:
        if inspect.isclass(expected_result) and issubclass(expected_result, Exception):
            with pytest.raises(expected_result):
                result = GETITEM_TEST_EXPRESSION[start:end]
                print(result)
        else:
            result = GETITEM_TEST_EXPRESSION[start:end]
            assert result == expected_result

def test_getitem_slice_symbol():
    with pytest.raises(IndexError):
        print(a[(0, ):()])
    with pytest.raises(IndexError):
        print(a[(0, ):(1, )])
    assert a[():()] == [a]


def test_lt():
    test = [
        (a,                             b),
        (a,                             Symbol('a', variable_name='x')),
        (Symbol('a', variable_name='x'),     Symbol('a', variable_name='y')),
        (a,                             _),
        (a,                             _s),
        (a,                             x_),
        (_,                             x_),
        (_s,                            x_),
        (x_,                            y_),
        (x_,                            x__),
        (f(a),                          f(b)),
        (f(a),                          f2(a)),
        (f(a),                          f(a, a)),
        (f(b),                          f(a, a)),
        (f(a, a),                       f(a, b)),
        (f(a, a),                       f(a, a, a)),
        (a,                             f(a)),
        (x_,                            f(a)),
        (_,                             f(a)),
        (_s,                            f(a)),
        (_s,                            s_),
        (SymbolWildcard(variable_name='x'),  SymbolWildcard(variable_name='y')),
        (s_,                            ss_),
        (_s,                            __),
        (_,                             _s),
        (SymbolWildcard(SpecialSymbol), SymbolWildcard(Symbol)),
        (f(a),                          SpecialF(a)),
    ]

    for expression1, expression2 in test:
        assert expression1 < expression2, "{!s} < {!s} did not hold".format(expression1, expression2)
        assert not (expression2 < expression1
                   ), "Inconsistent order: Both {0} < {1} and {1} < {0}".format(expression2, expression1)


def test_lt_error():
    test = [a, f(a), x_, _]
    for expression in test:
        with pytest.raises(TypeError):
            expression < object()


def test_operation_new_error():
    with pytest.raises(ValueError):
        _ = Operation.new('if', Arity.variadic)

    with pytest.raises(ValueError):
        _ = Operation.new('+', Arity.variadic)


def test_wildcard_error():
    with pytest.raises(ValueError):
        _ = Wildcard(-1, False)

    with pytest.raises(ValueError):
        _ = Wildcard(0, True)

def test_symbol_wildcard_error():
    with pytest.raises(TypeError):
        _ = SymbolWildcard(object)


def test_with_renamed_vars():
    test = [
            (a,                                 {},             a),
            (a,                                 {'x': 'y'},     a),
            (x_,                                {},             x_),
            (x_,                                {'x': 'y'},     y_),
            (SymbolWildcard(),                  {},             SymbolWildcard()),
            (SymbolWildcard(),                  {'x': 'y'},     SymbolWildcard()),
            (f(x_),                             {},             f(x_)),
            (f(x_),                             {'x': 'y'},     f(y_)),
    ]

    for expression, renaming, expected_result in test:
        new_expr = expression.with_renamed_vars(renaming)
        assert new_expr == expected_result


def test_hash():
    for expression, other in product(SIMPLE_EXPRESSIONS, SIMPLE_EXPRESSIONS):
        if expression != other:
            assert hash(expression) != hash(other), "hash({!s}) == hash({!s})".format(expression, other)
        else:
            assert hash(expression) == hash(other), "hash({!s}) != hash({!s})".format(expression, other)



def test_copy():
    for expression in SIMPLE_EXPRESSIONS:
        other = expression.__copy__()
        assert other == expression
        assert other is not expression



def test_contains():
    test = [
        (a,             a,              True),
        (a,             b,              False),
        (f(a),          a,              True),
        (f(a),          b,              False),
        (f(a),          f(a),           True),
        (f(a, b),       f(a),           False),
        (f(a),          f(a, b),        False),
        (f(x_, y_),     x_,             True),
        (f(x_, y_),     y_,             True),
        (f(x_, y_),     a,              False),
    ]

    for expression, subexpression, contains in test:
        if contains:
            assert subexpression in expression, "{!s} should be contained in {!s}".format(subexpression, expression)
        else:
            assert subexpression not in expression, "{!s} should not be contained in {!s}".format(subexpression, expression)

def test_one_identity_error():
    with pytest.raises(TypeError):
        Operation.new('Invalid', Arity.unary, one_identity=True)
    with pytest.raises(TypeError):
        Operation.new('Invalid', Arity.binary, one_identity=True)

def test_infix_error():
    with pytest.raises(TypeError):
        Operation.new('Invalid', Arity.unary, infix=True)

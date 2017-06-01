from sympy.pattern.expressions import (
    Arity, Operation, Symbol, Wildcard, SymbolWildcard, make_dot_variable,
    make_plus_variable, make_star_variable, make_symbol_variable
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

def test_simplified():
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

def test_syntactic():
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
    test = [                                               # expression        position
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
        #(f(a),                          SpecialF(a)),
    ]

    for expression1, expression2 in test:
        assert expression1 < expression2, "{!s} < {!s} did not hold".format(expression1, expression2)
        assert not (expression2 < expression1
                   ), "Inconsistent order: Both {0} < {1} and {1} < {0}".format(expression2, expression1)


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

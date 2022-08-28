from typing import Callable

import sympy
from sympy import Add, Implies, sqrt
from sympy.core import Mul, Pow
from sympy.core import (S, pi, symbols, Function, Rational, Integer,
                        Symbol, Eq, Ne, Le, Lt, Gt, Ge)
from sympy.functions import Piecewise, exp, sin, cos
from sympy.printing.smtlib import smtlib_code
from sympy.testing.pytest import raises

x, y, z = symbols('x,y,z')


def test_Integer():
    assert smtlib_code(Integer(67)) == "67"
    assert smtlib_code(Integer(-1)) == "-1"


def test_Rational():
    assert smtlib_code(Rational(3, 7)) == "(/ 3 7)"
    assert smtlib_code(Rational(18, 9)) == "2"
    assert smtlib_code(Rational(3, -7)) == "(/ -3 7)"
    assert smtlib_code(Rational(-3, -7)) == "(/ 3 7)"
    assert smtlib_code(x + Rational(3, 7), auto_declare=False) == "(+ (/ 3 7) x)"
    assert smtlib_code(Rational(3, 7) * x) == "(declare-const x Real)\n" \
                                              "(* (/ 3 7) x)"


def test_Relational():
    assert smtlib_code(Eq(x, y), auto_declare=False) == "(assert (= x y))"
    assert smtlib_code(Ne(x, y), auto_declare=False) == "(assert (not (= x y)))"
    assert smtlib_code(Le(x, y), auto_declare=False) == "(assert (<= x y))"
    assert smtlib_code(Lt(x, y), auto_declare=False) == "(assert (< x y))"
    assert smtlib_code(Gt(x, y), auto_declare=False) == "(assert (> x y))"
    assert smtlib_code(Ge(x, y), auto_declare=False) == "(assert (>= x y))"


def test_Function():
    assert smtlib_code(sin(x) ** cos(x), auto_declare=False) == "(pow (sin x) (cos x))"

    assert smtlib_code(
        abs(x),
        symbol_table={x: int, y: bool},
        known_types={int: "INTEGER_TYPE"},
        known_functions={sympy.Abs: "ABSOLUTE_VALUE_OF"}
    ) == "(declare-const x INTEGER_TYPE)\n" \
         "(ABSOLUTE_VALUE_OF x)"

    my_fun1 = Function('f1')
    assert smtlib_code(
        my_fun1(x),
        symbol_table={my_fun1: Callable[[bool], float]},
    ) == "(declare-const x Bool)\n" \
         "(declare-fun f1 (Bool) Real)\n" \
         "(f1 x)"

    assert smtlib_code(
        my_fun1(x),
        symbol_table={my_fun1: Callable[[bool], bool]},
    ) == "(declare-const x Bool)\n" \
         "(declare-fun f1 (Bool) Bool)\n" \
         "(assert (f1 x))"

    assert smtlib_code(
        Eq(my_fun1(x, z), y),
        symbol_table={my_fun1: Callable[[int, bool], bool]},
    ) == "(declare-const x Int)\n" \
         "(declare-const y Bool)\n" \
         "(declare-const z Bool)\n" \
         "(declare-fun f1 (Int Bool) Bool)\n" \
         "(assert (= (f1 x z) y))"

    assert smtlib_code(
        Eq(my_fun1(x, z), y),
        symbol_table={my_fun1: Callable[[int, bool], bool]},
        known_functions={my_fun1: "MY_KNOWN_FUN", Eq: '=='}
    ) == "(declare-const x Int)\n" \
         "(declare-const y Bool)\n" \
         "(declare-const z Bool)\n" \
         "(assert (== (MY_KNOWN_FUN x z) y))"

    assert smtlib_code(
        Eq(my_fun1(x, z), y),
        known_functions={my_fun1: "MY_KNOWN_FUN", Eq: '=='}
    ) == "(declare-const x Real)\n" \
         "(declare-const y Real)\n" \
         "(declare-const z Real)\n" \
         "(assert (== (MY_KNOWN_FUN x z) y))"


def test_Pow():
    assert smtlib_code(x ** 3, auto_declare=False) == "(pow x 3)"
    assert smtlib_code(x ** (y ** 3), auto_declare=False) == "(pow x (pow y 3))"
    assert smtlib_code(x ** Rational(2, 3), auto_declare=False) == '(pow x (/ 2 3))'

    a = Symbol('a', integer=True)
    b = Symbol('b', real=True)
    c = Symbol('c')

    def g(x): return 2 * x

    # if x=1, y=2, then expr=2.333...
    expr = 1 / (g(a) * 3.5) ** (a - b ** a) / (a ** 2 + b)

    assert smtlib_code([
        Eq(a < 2, c),
        Eq(b > a, c),
        c & True,
        Eq(expr, 2 + Rational(1, 3))
    ]) == '(declare-const a Int)\n' \
          '(declare-const b Real)\n' \
          '(declare-const c Bool)\n' \
          '(assert (= (< a 2) c))\n' \
          '(assert (= (> b a) c))\n' \
          '(assert c)\n' \
          '(assert (= ' \
          '(* (pow (* 7. a) (+ (pow b a) (* -1 a))) (pow (+ b (pow a 2)) -1)) ' \
          '(/ 7 3)' \
          '))'

    assert smtlib_code(
        Mul(-2, c, Pow(Mul(b, b, evaluate=False), -1, evaluate=False), evaluate=False)
    ) == '(declare-const b Real)\n' \
         '(declare-const c Real)\n' \
         '(* -2 c (pow (* b b) -1))'


def test_basic_ops():
    assert smtlib_code(x * y, auto_declare=False) == "(* x y)"
    assert smtlib_code(x + y, auto_declare=False) == "(+ x y)"
    # todo: implement re-write, currently does '(+ x (* -1 y))' instead
    # assert smtlib_code(x - y, auto_declare=False) == "(- x y)"
    assert smtlib_code(-x, auto_declare=False) == "(* -1 x)"


def test_quantifier_extensions():
    from sympy.logic.boolalg import Boolean
    from sympy import Interval, Tuple, sympify

    # start For-all quantifier class example
    class ForAll(Boolean):
        def _smtlib(self, printer):
            bound_symbol_declarations = [
                printer._s_expr(sym.name, [
                    printer._known_types[printer.symbol_table[sym]],
                    Interval(start, end)
                ]) for sym, start, end in self.limits
            ]
            return printer._s_expr('forall', [
                printer._s_expr('', bound_symbol_declarations),
                self.function
            ])

        @property
        def bound_symbols(self):
            return {s for s, _, _ in self.limits}

        @property
        def free_symbols(self):
            bound_symbol_names = {s.name for s in self.bound_symbols}
            return {
                s for s in self.function.free_symbols
                if s.name not in bound_symbol_names
            }

        def __new__(cls, *args):
            limits = [sympify(a) for a in args if isinstance(a, tuple) or isinstance(a, Tuple)]
            function = [sympify(a) for a in args if isinstance(a, Boolean)]
            assert len(limits) + len(function) == len(args)
            assert len(function) == 1
            function = function[0]

            if isinstance(function, ForAll): return ForAll.__new__(
                ForAll, *(limits + function.limits), function.function
            )
            inst = Boolean.__new__(cls)
            inst._args = tuple(limits + [function])
            inst.limits = limits
            inst.function = function
            return inst

    # end For-All Quantifier class example

    f = Function('f')
    assert smtlib_code(
        ForAll((x, -42, +21), Eq(f(x), f(x))),
        symbol_table={f: Callable[[float], float]}
    ) == '(assert (forall ( (x Real [-42, 21])) true))'

    assert smtlib_code(
        ForAll(
            (x, -42, +21), (y, -100, 3),
            Implies(Eq(x, y), Eq(f(x), f(y)))
        ),
        symbol_table={f: Callable[[float], float]}
    ) == '(declare-fun f (Real) Real)\n' \
         '(assert (' \
         'forall ( (x Real [-42, 21]) (y Real [-100, 3])) ' \
         '(=> (= x y) (= (f x) (f y)))' \
         '))'

    a = Symbol('a', integer=True)
    b = Symbol('b', real=True)
    c = Symbol('c')

    assert smtlib_code(
        ForAll(
            (a, 2, 100), ForAll(
                (b, 2, 100),
                Implies(a < b, sqrt(a) < b) | c
            ))
    ) == '(declare-const c Bool)\n' \
         '(assert (forall ( (a Int [2, 100]) (b Real [2, 100])) ' \
         '(or c (=> (< a b) (< (pow a (/ 1 2)) b)))' \
         '))'


def test_mix_number_mult_symbols():
    assert smtlib_code(
        1 / pi,
        known_constants={pi: "MY_PI"}
    ) == '(pow MY_PI -1)'

    assert smtlib_code(
        [
            Eq(pi, 3.14, evaluate=False),
            1 / pi,
        ],
        known_constants={pi: "MY_PI"}
    ) == '(assert (= MY_PI 3.14))\n' \
         '(pow MY_PI -1)'

    assert smtlib_code(
        Add(S.Zero, S.One, S.NegativeOne, S.Half,
            S.Exp1, S.Pi, S.GoldenRatio, evaluate=False),
        known_constants={
            S.Pi: 'p', S.GoldenRatio: 'g',
            S.Exp1: 'e'
        },
        known_functions={
            Add: 'plus',
            exp: 'exp'
        },
        precision=3
    ) == '(plus 0 1 -1 (/ 1 2) (exp 1) p g)'

    assert smtlib_code(
        Add(S.Zero, S.One, S.NegativeOne, S.Half,
            S.Exp1, S.Pi, S.GoldenRatio, evaluate=False),
        known_constants={
            S.Pi: 'p'
        },
        known_functions={
            Add: 'plus',
            exp: 'exp'
        },
        precision=3
    ) == '(plus 0 1 -1 (/ 1 2) (exp 1) p 1.62)'

    assert smtlib_code(
        Add(S.Zero, S.One, S.NegativeOne, S.Half,
            S.Exp1, S.Pi, S.GoldenRatio, evaluate=False),
        known_functions={Add: 'plus'},
        precision=3
    ) == '(plus 0 1 -1 (/ 1 2) 2.72 3.14 1.62)'

    assert smtlib_code(
        Add(S.Zero, S.One, S.NegativeOne, S.Half,
            S.Exp1, S.Pi, S.GoldenRatio, evaluate=False),
        known_constants={S.Exp1: 'e'},
        known_functions={Add: 'plus'},
        precision=3
    ) == '(plus 0 1 -1 (/ 1 2) e 3.14 1.62)'


def test_boolean():
    assert smtlib_code(x & y) == '(declare-const x Bool)\n' \
                                 '(declare-const y Bool)\n' \
                                 '(assert (and x y))'
    assert smtlib_code(x | y) == '(declare-const x Bool)\n' \
                                 '(declare-const y Bool)\n' \
                                 '(assert (or x y))'
    assert smtlib_code(~x) == '(declare-const x Bool)\n' \
                              '(assert (not x))'

    assert smtlib_code(x & y & z) == '(declare-const x Bool)\n' \
                                     '(declare-const y Bool)\n' \
                                     '(declare-const z Bool)\n' \
                                     '(assert (and x y z))'
    assert smtlib_code((x & ~y) | (z > 3)) == '(declare-const x Bool)\n' \
                                              '(declare-const y Bool)\n' \
                                              '(declare-const z Real)\n' \
                                              '(assert (or (> z 3) (and x (not y))))'

    f = Function('f')
    g = Function('g')
    h = Function('h')
    assert smtlib_code([
        Gt(f(x), y),
        Lt(y, g(z))
    ], symbol_table={
        f: Callable[[bool], int], g: Callable[[bool], int],
    }) == '(declare-const x Bool)\n' \
          '(declare-const y Real)\n' \
          '(declare-const z Bool)\n' \
          '(declare-fun f (Bool) Int)\n' \
          '(declare-fun g (Bool) Int)\n' \
          '(assert (> (f x) y))\n' \
          '(assert (< y (g z)))'

    assert smtlib_code([
        Eq(f(x), y),
        Lt(y, g(z))
    ], symbol_table={
        f: Callable[[bool], int], g: Callable[[bool], int],
    }) == '(declare-const x Bool)\n' \
          '(declare-const y Int)\n' \
          '(declare-const z Bool)\n' \
          '(declare-fun f (Bool) Int)\n' \
          '(declare-fun g (Bool) Int)\n' \
          '(assert (= (f x) y))\n' \
          '(assert (< y (g z)))'

    assert smtlib_code(
        [Eq(f(x), y),
         Eq(g(f(x)), z),
         Eq(h(g(f(x))), x)],
        symbol_table={
            f: Callable[[float], int],
            g: Callable[[int], bool],
            h: Callable[[bool], float]
        }
    ) == '(declare-const x Real)\n' \
         '(declare-const y Int)\n' \
         '(declare-const z Bool)\n' \
         '(declare-fun f (Real) Int)\n' \
         '(declare-fun g (Int) Bool)\n' \
         '(declare-fun h (Bool) Real)\n' \
         '(assert (= (f x) y))\n' \
         '(assert (= (g (f x)) z))\n' \
         '(assert (= (h (g (f x))) x))'


# todo: make smtlib_code support arrays
# def test_containers():
#     assert julia_code([1, 2, 3, [4, 5, [6, 7]], 8, [9, 10], 11]) == \
#            "Any[1, 2, 3, Any[4, 5, Any[6, 7]], 8, Any[9, 10], 11]"
#     assert julia_code((1, 2, (3, 4))) == "(1, 2, (3, 4))"
#     assert julia_code([1]) == "Any[1]"
#     assert julia_code((1,)) == "(1,)"
#     assert julia_code(Tuple(*[1, 2, 3])) == "(1, 2, 3)"
#     assert julia_code((1, x * y, (3, x ** 2))) == "(1, x .* y, (3, x .^ 2))"
#     # scalar, matrix, empty matrix and empty list
#     assert julia_code((1, eye(3), Matrix(0, 0, []), [])) == "(1, [1 0 0;\n0 1 0;\n0 0 1], zeros(0, 0), Any[])"


def test_smtlib_piecewise():
    assert smtlib_code(
        Piecewise((x, x < 1),
                  (x ** 2, True)),
        auto_declare=False
    ) == '(ite (< x 1) x (pow x 2))'

    assert smtlib_code(
        Piecewise((x ** 2, x < 1),
                  (x ** 3, x < 2),
                  (x ** 4, x < 3),
                  (x ** 5, True)),
        auto_declare=False
    ) == '(ite (< x 1) (pow x 2) ' \
         '(ite (< x 2) (pow x 3) ' \
         '(ite (< x 3) (pow x 4) ' \
         '(pow x 5))))'

    # Check that Piecewise without a True (default) condition error
    expr = Piecewise((x, x < 1), (x ** 2, x > 1), (sin(x), x > 0))
    raises(AssertionError, lambda: smtlib_code(expr))


def test_smtlib_piecewise_times_const():
    pw = Piecewise((x, x < 1), (x ** 2, True))
    assert smtlib_code(2 * pw) == '(declare-const x Real)\n(* 2 (ite (< x 1) x (pow x 2)))'
    assert smtlib_code(pw / x) == '(declare-const x Real)\n(* (pow x -1) (ite (< x 1) x (pow x 2)))'
    assert smtlib_code(pw / (x * y)) == '(declare-const x Real)\n(declare-const y Real)\n(* (pow x -1) (pow y -1) (ite (< x 1) x (pow x 2)))'
    assert smtlib_code(pw / 3) == '(declare-const x Real)\n(* (/ 1 3) (ite (< x 1) x (pow x 2)))'


# todo: make smtlib_code support arrays / matrices ?
# def test_smtlib_matrix_assign_to():
#     A = Matrix([[1, 2, 3]])
#     assert smtlib_code(A, assign_to='a') == "a = [1 2 3]"
#     A = Matrix([[1, 2], [3, 4]])
#     assert smtlib_code(A, assign_to='A') == "A = [1 2;\n3 4]"

# def test_julia_matrix_1x1():
#     A = Matrix([[3]])
#     B = MatrixSymbol('B', 1, 1)
#     C = MatrixSymbol('C', 1, 2)
#     assert julia_code(A, assign_to=B) == "B = [3]"
#     raises(ValueError, lambda: julia_code(A, assign_to=C))


# def test_julia_matrix_elements():
#     A = Matrix([[x, 2, x * y]])
#     assert julia_code(A[0, 0] ** 2 + A[0, 1] + A[0, 2]) == "x .^ 2 + x .* y + 2"
#     A = MatrixSymbol('AA', 1, 3)
#     assert julia_code(A) == "AA"
#     assert julia_code(A[0, 0] ** 2 + sin(A[0, 1]) + A[0, 2]) == \
#            "sin(AA[1,2]) + AA[1,1] .^ 2 + AA[1,3]"
#     assert julia_code(sum(A)) == "AA[1,1] + AA[1,2] + AA[1,3]"


def test_smtlib_boolean():
    assert smtlib_code(True, auto_assert=False) == 'true'
    assert smtlib_code(True) == '(assert true)'
    assert smtlib_code(S.true) == '(assert true)'
    assert smtlib_code(S.false) == '(assert false)'
    assert smtlib_code(False) == '(assert false)'
    assert smtlib_code(False, auto_assert=False) == 'false'


def test_not_supported():
    f = Function('f')
    raises(KeyError, lambda: smtlib_code(f(x).diff(x), symbol_table={f: Callable[[float], float]}))
    raises(KeyError, lambda: smtlib_code(S.ComplexInfinity))

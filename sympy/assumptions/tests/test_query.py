from __future__ import annotations
from sympy.abc import t, w, x, y, z, n, k, m, p, i
from sympy.assumptions import (ask, AssumptionsContext, Q)
from sympy.assumptions.assume import (assuming, global_assumptions, Predicate)
from sympy.assumptions.ask import _ask_recursive
from sympy.assumptions.cnf import CNF, Literal
from sympy.assumptions.facts import (single_fact_lookup,
    get_known_facts, generate_known_facts_dict, get_known_facts_keys)
from sympy.assumptions.ask_generated import (get_all_known_facts,
    get_known_facts_dict)
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.numbers import (I, Integer, Rational, oo, zoo, pi)
from sympy.core.singleton import S
from sympy.core.power import Pow
from sympy.core.symbol import Str, symbols, Symbol
from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.elementary.complexes import (Abs, im, re, sign)
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (
    acos, acot, asin, atan, cos, cot, sin, tan)
from sympy.logic.boolalg import Equivalent, Implies, Xor, And, to_cnf
from sympy.matrices import Matrix, SparseMatrix
from sympy.testing.pytest import XFAIL, slow, raises,  _both_exp_pow
import math


def test_int_1():
    z = 1
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is True
    assert _ask_recursive(Q.rational(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is True
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_int_11():
    z = 11
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is True
    assert _ask_recursive(Q.rational(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is True
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is True
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_int_12():
    z = 12
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is True
    assert _ask_recursive(Q.rational(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is True
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is True
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_float_1():
    z = 1.0
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is None
    assert _ask_recursive(Q.rational(z)) is None
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is None
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is None
    assert _ask_recursive(Q.odd(z)) is None
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is None
    assert _ask_recursive(Q.composite(z)) is None
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is None

    z = 7.2123
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is None
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is None
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is None

    # test for issue #12168
    assert _ask_recursive(Q.rational(math.pi)) is None


def test_zero_0():
    z = Integer(0)
    assert _ask_recursive(Q.nonzero(z)) is False
    assert _ask_recursive(Q.zero(z)) is True
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is True
    assert _ask_recursive(Q.rational(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is False
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is True
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is True
    assert _ask_recursive(Q.transcendental(z)) is False


def test_negativeone():
    z = Integer(-1)
    assert _ask_recursive(Q.nonzero(z)) is True
    assert _ask_recursive(Q.zero(z)) is False
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is True
    assert _ask_recursive(Q.rational(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is False
    assert _ask_recursive(Q.negative(z)) is True
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is True
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_infinity():
    assert _ask_recursive(Q.commutative(oo)) is True
    assert _ask_recursive(Q.integer(oo)) is False
    assert _ask_recursive(Q.rational(oo)) is False
    assert _ask_recursive(Q.algebraic(oo)) is False
    assert _ask_recursive(Q.real(oo)) is False
    assert _ask_recursive(Q.extended_real(oo)) is True
    assert _ask_recursive(Q.complex(oo)) is False
    assert _ask_recursive(Q.irrational(oo)) is False
    assert _ask_recursive(Q.imaginary(oo)) is False
    assert _ask_recursive(Q.positive(oo)) is False
    assert _ask_recursive(Q.extended_positive(oo)) is True
    assert _ask_recursive(Q.negative(oo)) is False
    assert _ask_recursive(Q.even(oo)) is False
    assert _ask_recursive(Q.odd(oo)) is False
    assert _ask_recursive(Q.finite(oo)) is False
    assert _ask_recursive(Q.infinite(oo)) is True
    assert _ask_recursive(Q.prime(oo)) is False
    assert _ask_recursive(Q.composite(oo)) is False
    assert _ask_recursive(Q.hermitian(oo)) is False
    assert _ask_recursive(Q.antihermitian(oo)) is False
    assert _ask_recursive(Q.positive_infinite(oo)) is True
    assert _ask_recursive(Q.negative_infinite(oo)) is False
    assert _ask_recursive(Q.transcendental(oo)) is False


def test_neg_infinity():
    mm = S.NegativeInfinity
    assert _ask_recursive(Q.commutative(mm)) is True
    assert _ask_recursive(Q.integer(mm)) is False
    assert _ask_recursive(Q.rational(mm)) is False
    assert _ask_recursive(Q.algebraic(mm)) is False
    assert _ask_recursive(Q.real(mm)) is False
    assert _ask_recursive(Q.extended_real(mm)) is True
    assert _ask_recursive(Q.complex(mm)) is False
    assert _ask_recursive(Q.irrational(mm)) is False
    assert _ask_recursive(Q.imaginary(mm)) is False
    assert _ask_recursive(Q.positive(mm)) is False
    assert _ask_recursive(Q.negative(mm)) is False
    assert _ask_recursive(Q.extended_negative(mm)) is True
    assert _ask_recursive(Q.even(mm)) is False
    assert _ask_recursive(Q.odd(mm)) is False
    assert _ask_recursive(Q.finite(mm)) is False
    assert _ask_recursive(Q.infinite(oo)) is True
    assert _ask_recursive(Q.prime(mm)) is False
    assert _ask_recursive(Q.composite(mm)) is False
    assert _ask_recursive(Q.hermitian(mm)) is False
    assert _ask_recursive(Q.antihermitian(mm)) is False
    assert _ask_recursive(Q.positive_infinite(-oo)) is False
    assert _ask_recursive(Q.negative_infinite(-oo)) is True
    assert _ask_recursive(Q.transcendental(-oo)) is False


def test_complex_infinity():
    assert _ask_recursive(Q.commutative(zoo)) is True
    assert _ask_recursive(Q.integer(zoo)) is False
    assert _ask_recursive(Q.rational(zoo)) is False
    assert _ask_recursive(Q.algebraic(zoo)) is False
    assert _ask_recursive(Q.real(zoo)) is False
    assert _ask_recursive(Q.extended_real(zoo)) is False
    assert _ask_recursive(Q.complex(zoo)) is False
    assert _ask_recursive(Q.irrational(zoo)) is False
    assert _ask_recursive(Q.imaginary(zoo)) is False
    assert _ask_recursive(Q.positive(zoo)) is False
    assert _ask_recursive(Q.negative(zoo)) is False
    assert _ask_recursive(Q.zero(zoo)) is False
    assert _ask_recursive(Q.nonzero(zoo)) is False
    assert _ask_recursive(Q.even(zoo)) is False
    assert _ask_recursive(Q.odd(zoo)) is False
    assert _ask_recursive(Q.finite(zoo)) is False
    assert _ask_recursive(Q.infinite(zoo)) is True
    assert _ask_recursive(Q.prime(zoo)) is False
    assert _ask_recursive(Q.composite(zoo)) is False
    assert _ask_recursive(Q.hermitian(zoo)) is False
    assert _ask_recursive(Q.antihermitian(zoo)) is False
    assert _ask_recursive(Q.positive_infinite(zoo)) is False
    assert _ask_recursive(Q.negative_infinite(zoo)) is False
    assert _ask_recursive(Q.transcendental(zoo)) is False


def test_positive_infinite_symbolic():
    assert _ask_recursive(Q.positive_infinite(x), Q.positive(x)) is False
    assert _ask_recursive(Q.positive_infinite(x), Q.finite(x)) is False
    assert _ask_recursive(Q.positive_infinite(x), Q.complex(x)) is False
    assert _ask_recursive(Q.positive_infinite(x + I), Q.real(x)) is False
    assert _ask_recursive(Q.positive_infinite(x*I), Q.real(x)) is False
    assert _ask_recursive(Q.positive_infinite(x + y), Q.real(x) & Q.imaginary(y)) is False


def test_negative_infinite_symbolic():
    assert _ask_recursive(Q.negative_infinite(x), Q.positive(x)) is False
    assert _ask_recursive(Q.negative_infinite(x), Q.finite(x)) is False
    assert _ask_recursive(Q.negative_infinite(x), Q.complex(x)) is False
    assert _ask_recursive(Q.negative_infinite(x + I), Q.real(x)) is False
    assert _ask_recursive(Q.negative_infinite(x*I), Q.real(x)) is False
    assert _ask_recursive(Q.negative_infinite(x + y), Q.real(x) & Q.imaginary(y)) is False


def test_nan():
    nan = S.NaN
    assert _ask_recursive(Q.commutative(nan)) is True
    assert _ask_recursive(Q.integer(nan)) is None
    assert _ask_recursive(Q.rational(nan)) is None
    assert _ask_recursive(Q.algebraic(nan)) is None
    assert _ask_recursive(Q.real(nan)) is None
    assert _ask_recursive(Q.extended_real(nan)) is None
    assert _ask_recursive(Q.complex(nan)) is None
    assert _ask_recursive(Q.irrational(nan)) is None
    assert _ask_recursive(Q.imaginary(nan)) is None
    assert _ask_recursive(Q.positive(nan)) is None
    assert _ask_recursive(Q.nonzero(nan)) is None
    assert _ask_recursive(Q.zero(nan)) is None
    assert _ask_recursive(Q.even(nan)) is None
    assert _ask_recursive(Q.odd(nan)) is None
    assert _ask_recursive(Q.finite(nan)) is None
    assert _ask_recursive(Q.infinite(nan)) is None
    assert _ask_recursive(Q.prime(nan)) is None
    assert _ask_recursive(Q.composite(nan)) is None
    assert _ask_recursive(Q.hermitian(nan)) is None
    assert _ask_recursive(Q.antihermitian(nan)) is None
    assert _ask_recursive(Q.transcendental(nan)) is None


def test_Rational_number():
    r = Rational(3, 4)
    assert _ask_recursive(Q.commutative(r)) is True
    assert _ask_recursive(Q.integer(r)) is False
    assert _ask_recursive(Q.rational(r)) is True
    assert _ask_recursive(Q.real(r)) is True
    assert _ask_recursive(Q.complex(r)) is True
    assert _ask_recursive(Q.irrational(r)) is False
    assert _ask_recursive(Q.imaginary(r)) is False
    assert _ask_recursive(Q.positive(r)) is True
    assert _ask_recursive(Q.negative(r)) is False
    assert _ask_recursive(Q.even(r)) is False
    assert _ask_recursive(Q.odd(r)) is False
    assert _ask_recursive(Q.finite(r)) is True
    assert _ask_recursive(Q.prime(r)) is False
    assert _ask_recursive(Q.composite(r)) is False
    assert _ask_recursive(Q.hermitian(r)) is True
    assert _ask_recursive(Q.antihermitian(r)) is False
    assert _ask_recursive(Q.transcendental(r)) is False

    r = Rational(1, 4)
    assert _ask_recursive(Q.positive(r)) is True
    assert _ask_recursive(Q.negative(r)) is False

    r = Rational(5, 4)
    assert _ask_recursive(Q.negative(r)) is False
    assert _ask_recursive(Q.positive(r)) is True

    r = Rational(5, 3)
    assert _ask_recursive(Q.positive(r)) is True
    assert _ask_recursive(Q.negative(r)) is False

    r = Rational(-3, 4)
    assert _ask_recursive(Q.positive(r)) is False
    assert _ask_recursive(Q.negative(r)) is True

    r = Rational(-1, 4)
    assert _ask_recursive(Q.positive(r)) is False
    assert _ask_recursive(Q.negative(r)) is True

    r = Rational(-5, 4)
    assert _ask_recursive(Q.negative(r)) is True
    assert _ask_recursive(Q.positive(r)) is False

    r = Rational(-5, 3)
    assert _ask_recursive(Q.positive(r)) is False
    assert _ask_recursive(Q.negative(r)) is True


def test_sqrt_2():
    z = sqrt(2)
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_pi():
    z = S.Pi
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is True

    z = S.Pi + 1
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is True

    z = 2*S.Pi
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is True

    z = S.Pi ** 2
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is True

    z = (1 + S.Pi) ** 2
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is True


def test_E():
    z = S.Exp1
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is False
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is True


def test_GoldenRatio():
    z = S.GoldenRatio
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_TribonacciConstant():
    z = S.TribonacciConstant
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is True
    assert _ask_recursive(Q.real(z)) is True
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is True
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is True
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is True
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False


def test_I():
    z = I
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is True
    assert _ask_recursive(Q.real(z)) is False
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is True
    assert _ask_recursive(Q.positive(z)) is False
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is False
    assert _ask_recursive(Q.antihermitian(z)) is True
    assert _ask_recursive(Q.transcendental(z)) is False

    z = 1 + I
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is True
    assert _ask_recursive(Q.real(z)) is False
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is False
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is False
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False

    z = I*(1 + I)
    assert _ask_recursive(Q.commutative(z)) is True
    assert _ask_recursive(Q.integer(z)) is False
    assert _ask_recursive(Q.rational(z)) is False
    assert _ask_recursive(Q.algebraic(z)) is True
    assert _ask_recursive(Q.real(z)) is False
    assert _ask_recursive(Q.complex(z)) is True
    assert _ask_recursive(Q.irrational(z)) is False
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.positive(z)) is False
    assert _ask_recursive(Q.negative(z)) is False
    assert _ask_recursive(Q.even(z)) is False
    assert _ask_recursive(Q.odd(z)) is False
    assert _ask_recursive(Q.finite(z)) is True
    assert _ask_recursive(Q.prime(z)) is False
    assert _ask_recursive(Q.composite(z)) is False
    assert _ask_recursive(Q.hermitian(z)) is False
    assert _ask_recursive(Q.antihermitian(z)) is False
    assert _ask_recursive(Q.transcendental(z)) is False

    z = I**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

    z = (-I)**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

    z = (3*I)**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is False

    z = (1)**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

    z = (-1)**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

    z = (1+I)**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is False

    z = (I)**(I+3)
    assert _ask_recursive(Q.imaginary(z)) is True
    assert _ask_recursive(Q.real(z)) is False

    z = (I)**(I+2)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

    z = (I)**(2)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

    z = (I)**(3)
    assert _ask_recursive(Q.imaginary(z)) is True
    assert _ask_recursive(Q.real(z)) is False

    z = (3)**(I)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is False

    z = (I)**(0)
    assert _ask_recursive(Q.imaginary(z)) is False
    assert _ask_recursive(Q.real(z)) is True

def test_bounded():
    x, y, z = symbols('x,y,z')
    a = x + y
    x, y = a.args
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(x)) is None
    assert _ask_recursive(Q.finite(x), Q.finite(x)) is True
    assert _ask_recursive(Q.finite(x), Q.finite(y)) is None
    assert _ask_recursive(Q.finite(x), Q.complex(x)) is True
    assert _ask_recursive(Q.finite(x), Q.extended_real(x)) is None

    assert _ask_recursive(Q.finite(x + 1)) is None
    assert _ask_recursive(Q.finite(x + 1), Q.finite(x)) is True
    a = x + y
    x, y = a.args
    # B + B
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)) is True
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.finite(y)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)) is True
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)) is True
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.finite(y)
        & ~Q.positive(y)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.positive(x)
        & Q.positive(y)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y) & ~Q.positive(x)
        & ~Q.positive(y)) is True
    # B + U
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x)
        & Q.positive_infinite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x)
        & Q.positive_infinite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)
        & ~Q.positive(y)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.positive(x)
        & Q.positive_infinite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.positive(x) & ~Q.finite(y)
        & ~Q.positive(y)) is False
    # B + ?
    assert _ask_recursive(Q.finite(a), Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x)
        & Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x)
        & Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.positive(x)
        & Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.positive(x)
        & ~Q.positive(y)) is None
    # U + U
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.positive_infinite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x) & ~Q.finite(y)
        & ~Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.extended_positive(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & ~Q.extended_positive(x) & ~Q.extended_positive(y)) is False
    # U + ?
    assert _ask_recursive(Q.finite(a), ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)
        & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)
        & Q.positive_infinite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)
        & ~Q.finite(y) & ~Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.extended_positive(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.extended_positive(x) & ~Q.finite(y)
        & ~Q.extended_positive(y)) is False
    # ? + ?
    assert _ask_recursive(Q.finite(a)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)
        & Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)
        & ~Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.extended_positive(x)
        & Q.extended_positive(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.extended_positive(x)
        & ~Q.extended_positive(y)) is None

    x, y, z = symbols('x,y,z')
    a = x + y + z
    x, y, z = a.args
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.negative(z)) is True
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.finite(z)) is True
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.positive(z)) is True
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & Q.finite(z)) is True
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & Q.positive(z)) is True
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.finite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive(y)
        & Q.positive(z)) is True
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.extended_positive(y)
        & Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative_infinite(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative_infinite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative_infinite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative_infinite(y)
        & Q.extended_negative(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x)
        & Q.negative_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.negative_infinite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & ~Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & ~Q.finite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & ~Q.finite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & ~Q.finite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive_infinite(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive_infinite(y)
        & Q.negative_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) &
         Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.positive_infinite(y)
        & Q.extended_positive(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.extended_negative(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.extended_negative(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative(x) & Q.extended_positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.finite(z)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.positive(z)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)
        & Q.positive(z)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.negative_infinite(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.negative_infinite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.negative_infinite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.negative_infinite(y)
        & Q.extended_negative(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x)
        & Q.negative_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.negative_infinite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive_infinite(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive_infinite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.positive_infinite(y)
        & Q.extended_positive(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.extended_negative(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.extended_negative(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.extended_positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)
        & Q.positive(z)) is True
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.negative_infinite(y)
        & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.negative_infinite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.negative_infinite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.negative_infinite(y)
        & Q.extended_negative(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x)
        & Q.negative_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.negative_infinite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & ~Q.finite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive_infinite(y)
        & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive_infinite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.positive_infinite(y)
        & Q.extended_positive(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.extended_negative(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.extended_negative(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive(x) & Q.extended_positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.negative_infinite(y) & Q.negative_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.negative_infinite(y) & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.negative_infinite(y)& Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.negative_infinite(y) & Q.extended_negative(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.negative_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.negative_infinite(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & ~Q.finite(y) & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & ~Q.finite(y) & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & ~Q.finite(y) & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & ~Q.finite(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.positive_infinite(y) & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.positive_infinite(y) & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.positive_infinite(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.extended_negative(y) & Q.extended_negative(z)) is False
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.extended_negative(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.negative_infinite(x)
        & Q.extended_positive(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.positive_infinite(z)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.positive_infinite(y)
        & Q.positive_infinite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.positive_infinite(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.positive_infinite(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.extended_negative(y)
        & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.extended_negative(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.extended_positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.positive_infinite(y) & Q.positive_infinite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.positive_infinite(y) & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.positive_infinite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.positive_infinite(y) & Q.extended_positive(z)) is False
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.extended_negative(y) & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.extended_negative(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.positive_infinite(x)
        & Q.extended_positive(y) & Q.extended_positive(z)) is False
    assert _ask_recursive(Q.finite(a), Q.extended_negative(x)
        & Q.extended_negative(y) & Q.extended_negative(z)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_negative(x)
        & Q.extended_negative(y)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_negative(x)
        & Q.extended_negative(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_negative(x)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_negative(x)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_negative(x)
        & Q.extended_positive(y) & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(y)
        & Q.extended_positive(z)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_positive(x)
        & Q.extended_positive(y) & Q.extended_positive(z)) is None

    assert _ask_recursive(Q.finite(2*x)) is None
    assert _ask_recursive(Q.finite(2*x), Q.finite(x)) is True

    x, y, z = symbols('x,y,z')
    a = x*y
    x, y = a.args
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.zero(x) & ~Q.finite(y)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.finite(y) &~Q.zero(y)) is False
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)) is False
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a)) is None
    a = x*y*z
    x, y, z = a.args
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & Q.finite(z)) is True
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.zero(x) & Q.finite(y)
        & ~Q.zero(y) & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.zero(x) & ~Q.finite(y)
        & Q.finite(z) & ~Q.zero(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.zero(x) & ~Q.finite(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.finite(y) & ~Q.zero(y)
        & Q.finite(z) & ~Q.zero(z)) is False
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.zero(x) & Q.finite(y)
        & ~Q.zero(y) & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & Q.finite(z) & ~Q.zero(z)) is False
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & ~Q.finite(z)) is False
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(y) & Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(y) & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(y) & Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(y) & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(z) & Q.extended_nonzero(x)
        & Q.extended_nonzero(y) & Q.extended_nonzero(z)) is None
    assert _ask_recursive(Q.finite(a), Q.extended_nonzero(x) & ~Q.finite(y)
        & Q.extended_nonzero(y) & ~Q.finite(z)
        & Q.extended_nonzero(z)) is False

    x, y, z = symbols('x,y,z')
    assert _ask_recursive(Q.finite(x**2)) is None
    assert _ask_recursive(Q.finite(2**x)) is None
    assert _ask_recursive(Q.finite(2**x), Q.finite(x)) is True
    assert _ask_recursive(Q.finite(x**x)) is None
    assert _ask_recursive(Q.finite(S.Half ** x)) is None
    assert _ask_recursive(Q.finite(S.Half ** x), Q.extended_positive(x)) is True
    assert _ask_recursive(Q.finite(S.Half ** x), Q.extended_negative(x)) is None
    assert _ask_recursive(Q.finite(2**x), Q.extended_negative(x)) is True
    assert _ask_recursive(Q.finite(sqrt(x))) is None
    assert _ask_recursive(Q.finite(2**x), ~Q.finite(x)) is False
    assert _ask_recursive(Q.finite(x**2), ~Q.finite(x)) is False

    # https://github.com/sympy/sympy/issues/27707
    assert _ask_recursive(Q.finite(x**y), Q.real(x) & Q.real(y)) is None
    assert _ask_recursive(Q.finite(x**y), Q.real(x) & Q.negative(y)) is None
    assert _ask_recursive(Q.finite(x**y), Q.zero(x) & Q.negative(y)) is False
    assert _ask_recursive(Q.finite(x**y), Q.real(x) & Q.positive(y)) is True
    assert _ask_recursive(Q.finite(x**y), Q.nonzero(x) & Q.real(y)) is True
    assert _ask_recursive(Q.finite(x**y), Q.nonzero(x) & Q.negative(y)) is True
    assert _ask_recursive(Q.finite(x**y), Q.zero(x) & Q.positive(y)) is True

    # sign function
    assert _ask_recursive(Q.finite(sign(x))) is True
    assert _ask_recursive(Q.finite(sign(x)), ~Q.finite(x)) is True

    # exponential functions
    assert _ask_recursive(Q.finite(log(x))) is None
    assert _ask_recursive(Q.finite(log(x)), Q.finite(x)) is None
    assert _ask_recursive(Q.finite(log(x)), ~Q.zero(x)) is True
    assert _ask_recursive(Q.finite(log(x)), Q.infinite(x)) is False
    assert _ask_recursive(Q.finite(log(x)), Q.zero(x)) is False
    assert _ask_recursive(Q.finite(exp(x))) is None
    assert _ask_recursive(Q.finite(exp(x)), Q.finite(x)) is True
    assert _ask_recursive(Q.finite(exp(2))) is True

    # trigonometric functions
    assert _ask_recursive(Q.finite(sin(x))) is True
    assert _ask_recursive(Q.finite(sin(x)), ~Q.finite(x)) is True
    assert _ask_recursive(Q.finite(cos(x))) is True
    assert _ask_recursive(Q.finite(cos(x)), ~Q.finite(x)) is True
    assert _ask_recursive(Q.finite(2*sin(x))) is True
    assert _ask_recursive(Q.finite(sin(x)**2)) is True
    assert _ask_recursive(Q.finite(cos(x)**2)) is True
    assert _ask_recursive(Q.finite(cos(x) + sin(x))) is True

@XFAIL
def test_unbounded_xfail():
    # TODO: Rewrite the logic for the zero and nonzero recursive
    # handlers in handlers/order.py so that these tests pass.
    assert _ask_recursive(Q.zero(1/y), Q.finite(y) & ~Q.zero(y)) is False
    # This is no longer handled as a side effect of the fix for issue #28129.
    assert _ask_recursive(Q.infinite(x / y), Q.infinite(x) & Q.finite(y) & ~Q.zero(y)) is True

def test_unbounded():
    assert _ask_recursive(Q.infinite(I * oo)) is True
    assert _ask_recursive(Q.infinite(1 + I*oo)) is True
    assert _ask_recursive(Q.infinite(3 * (I * oo))) is True
    assert _ask_recursive(Q.infinite(-I * oo)) is True
    assert _ask_recursive(Q.infinite(1 + zoo)) is True
    assert _ask_recursive(Q.infinite(I * zoo)) is True
    assert _ask_recursive(Q.infinite(I * oo - I * oo)) is None
    assert _ask_recursive(Q.infinite(x * I * oo)) is None
    assert _ask_recursive(Q.infinite(1 / x), Q.finite(x) & ~Q.zero(x)) is False
    assert _ask_recursive(Q.infinite(1 / (I * oo))) is False


def test_issue_27441():
    # https://github.com/sympy/sympy/issues/27441
    assert _ask_recursive(Q.composite(y), Q.integer(y) & Q.positive(y) & ~Q.prime(y)) is None


def test_issue_27447():
    x,y,z = symbols('x y z')
    a = x*y
    assert _ask_recursive(Q.finite(a), Q.finite(x)  & ~Q.finite(y)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x)  & Q.finite(y)) is None

    a = x*y*z
    assert _ask_recursive(Q.finite(a), Q.finite(x) & Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)
        & Q.finite(z) ) is None
    assert _ask_recursive(Q.finite(a), Q.finite(x) & ~Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.finite(y)
        & Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & Q.finite(y)
        & ~Q.finite(z)) is None
    assert _ask_recursive(Q.finite(a), ~Q.finite(x) & ~Q.finite(y)
        & Q.finite(z)) is None


@XFAIL
def test_issue_27662_xfail():
    assert ask(Q.finite(x*y), ~Q.finite(x)
        & Q.zero(y)) is None


@XFAIL
def test_bounded_xfail():
    """We need to support relations in ask for this to work"""
    assert _ask_recursive(Q.finite(sin(x)**x)) is True
    assert _ask_recursive(Q.finite(cos(x)**x)) is True


def test_commutative():
    """By default objects are Q.commutative that is why it returns True
    for both key=True and key=False"""
    assert _ask_recursive(Q.commutative(x)) is True
    assert _ask_recursive(Q.commutative(x), ~Q.commutative(x)) is False
    assert _ask_recursive(Q.commutative(x), Q.complex(x)) is True
    assert _ask_recursive(Q.commutative(x), Q.imaginary(x)) is True
    assert _ask_recursive(Q.commutative(x), Q.real(x)) is True
    assert _ask_recursive(Q.commutative(x), Q.positive(x)) is True
    assert _ask_recursive(Q.commutative(x), ~Q.commutative(y)) is True

    assert _ask_recursive(Q.commutative(2*x)) is True
    assert _ask_recursive(Q.commutative(2*x), ~Q.commutative(x)) is False

    assert _ask_recursive(Q.commutative(x + 1)) is True
    assert _ask_recursive(Q.commutative(x + 1), ~Q.commutative(x)) is False

    assert _ask_recursive(Q.commutative(x**2)) is True
    assert _ask_recursive(Q.commutative(x**2), ~Q.commutative(x)) is False

    assert _ask_recursive(Q.commutative(log(x))) is True


@_both_exp_pow
def test_complex():
    assert _ask_recursive(Q.complex(x)) is None
    assert _ask_recursive(Q.complex(x), Q.complex(x)) is True
    assert _ask_recursive(Q.complex(x), Q.complex(y)) is None
    assert _ask_recursive(Q.complex(x), ~Q.complex(x)) is False
    assert _ask_recursive(Q.complex(x), Q.real(x)) is True
    assert _ask_recursive(Q.complex(x), ~Q.real(x)) is None
    assert _ask_recursive(Q.complex(x), Q.rational(x)) is True
    assert _ask_recursive(Q.complex(x), Q.irrational(x)) is True
    assert _ask_recursive(Q.complex(x), Q.positive(x)) is True
    assert _ask_recursive(Q.complex(x), Q.imaginary(x)) is True
    assert _ask_recursive(Q.complex(x), Q.algebraic(x)) is True

    # a+b
    assert _ask_recursive(Q.complex(x + 1), Q.complex(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.real(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.rational(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.irrational(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.imaginary(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.integer(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.even(x)) is True
    assert _ask_recursive(Q.complex(x + 1), Q.odd(x)) is True
    assert _ask_recursive(Q.complex(x + y), Q.complex(x) & Q.complex(y)) is True
    assert _ask_recursive(Q.complex(x + y), Q.real(x) & Q.imaginary(y)) is True

    # a*x +b
    assert _ask_recursive(Q.complex(2*x + 1), Q.complex(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.real(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.positive(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.rational(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.irrational(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.imaginary(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.integer(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.even(x)) is True
    assert _ask_recursive(Q.complex(2*x + 1), Q.odd(x)) is True

    # x**2
    assert _ask_recursive(Q.complex(x**2), Q.complex(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.real(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.positive(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.rational(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.irrational(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.imaginary(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.integer(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.even(x)) is True
    assert _ask_recursive(Q.complex(x**2), Q.odd(x)) is True

    # 2**x
    assert _ask_recursive(Q.complex(2**x), Q.complex(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.real(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.positive(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.rational(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.irrational(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.imaginary(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.integer(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.even(x)) is True
    assert _ask_recursive(Q.complex(2**x), Q.odd(x)) is True
    assert _ask_recursive(Q.complex(x**y), Q.complex(x) & Q.complex(y)) is True

    # trigonometric expressions
    assert _ask_recursive(Q.complex(sin(x))) is True
    assert _ask_recursive(Q.complex(sin(2*x + 1))) is True
    assert _ask_recursive(Q.complex(cos(x))) is True
    assert _ask_recursive(Q.complex(cos(2*x + 1))) is True

    # exponential
    assert _ask_recursive(Q.complex(exp(x))) is True
    assert _ask_recursive(Q.complex(exp(x))) is True

    # Q.complexes
    assert _ask_recursive(Q.complex(Abs(x))) is True
    assert _ask_recursive(Q.complex(re(x))) is True
    assert _ask_recursive(Q.complex(im(x))) is True


def test_even_query():
    assert _ask_recursive(Q.even(x)) is None
    assert _ask_recursive(Q.even(x), Q.integer(x)) is None
    assert _ask_recursive(Q.even(x), ~Q.integer(x)) is False
    assert _ask_recursive(Q.even(x), Q.rational(x)) is None
    assert _ask_recursive(Q.even(x), Q.positive(x)) is None

    assert _ask_recursive(Q.even(2*x)) is None
    assert _ask_recursive(Q.even(2*x), Q.integer(x)) is True
    assert _ask_recursive(Q.even(2*x), Q.even(x)) is True
    assert _ask_recursive(Q.even(2*x), Q.irrational(x)) is False
    assert _ask_recursive(Q.even(2*x), Q.odd(x)) is True
    assert _ask_recursive(Q.even(2*x), ~Q.integer(x)) is None
    assert _ask_recursive(Q.even(3*x), Q.integer(x)) is None
    assert _ask_recursive(Q.even(3*x), Q.even(x)) is True
    assert _ask_recursive(Q.even(3*x), Q.odd(x)) is False

    assert _ask_recursive(Q.even(x + 1), Q.odd(x)) is True
    assert _ask_recursive(Q.even(x + 1), Q.even(x)) is False
    assert _ask_recursive(Q.even(x + 2), Q.odd(x)) is False
    assert _ask_recursive(Q.even(x + 2), Q.even(x)) is True
    assert _ask_recursive(Q.even(7 - x), Q.odd(x)) is True
    assert _ask_recursive(Q.even(7 + x), Q.odd(x)) is True
    assert _ask_recursive(Q.even(x + y), Q.odd(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.even(x + y), Q.odd(x) & Q.even(y)) is False
    assert _ask_recursive(Q.even(x + y), Q.even(x) & Q.even(y)) is True

    assert _ask_recursive(Q.even(2*x + 1), Q.integer(x)) is False
    assert _ask_recursive(Q.even(2*x*y), Q.rational(x) & Q.rational(x)) is None
    assert _ask_recursive(Q.even(2*x*y), Q.irrational(x) & Q.irrational(x)) is None

    assert _ask_recursive(Q.even(x + y + z), Q.odd(x) & Q.odd(y) & Q.even(z)) is True
    assert _ask_recursive(Q.even(x + y + z + t),
        Q.odd(x) & Q.odd(y) & Q.even(z) & Q.integer(t)) is None

    assert _ask_recursive(Q.even(Abs(x)), Q.even(x)) is True
    assert _ask_recursive(Q.even(Abs(x)), ~Q.even(x)) is None
    assert _ask_recursive(Q.even(re(x)), Q.even(x)) is True
    assert _ask_recursive(Q.even(re(x)), ~Q.even(x)) is None
    assert _ask_recursive(Q.even(im(x)), Q.even(x)) is True
    assert _ask_recursive(Q.even(im(x)), Q.real(x)) is True

    assert _ask_recursive(Q.even((-1)**n), Q.integer(n)) is False

    assert _ask_recursive(Q.even(k**2), Q.even(k)) is True
    assert _ask_recursive(Q.even(n**2), Q.odd(n)) is False
    assert _ask_recursive(Q.even(2**k), Q.even(k)) is None
    assert _ask_recursive(Q.even(x**2)) is None

    assert _ask_recursive(Q.even(k**m), Q.even(k) & Q.integer(m) & ~Q.negative(m)) is None
    assert _ask_recursive(Q.even(n**m), Q.odd(n) & Q.integer(m) & ~Q.negative(m)) is False

    assert _ask_recursive(Q.even(k**p), Q.even(k) & Q.integer(p) & Q.positive(p)) is True
    assert _ask_recursive(Q.even(n**p), Q.odd(n) & Q.integer(p) & Q.positive(p)) is False

    assert _ask_recursive(Q.even(m**k), Q.even(k) & Q.integer(m) & ~Q.negative(m)) is None
    assert _ask_recursive(Q.even(p**k), Q.even(k) & Q.integer(p) & Q.positive(p)) is None

    assert _ask_recursive(Q.even(m**n), Q.odd(n) & Q.integer(m) & ~Q.negative(m)) is None
    assert _ask_recursive(Q.even(p**n), Q.odd(n) & Q.integer(p) & Q.positive(p)) is None

    assert _ask_recursive(Q.even(k**x), Q.even(k)) is None
    assert _ask_recursive(Q.even(n**x), Q.odd(n)) is None

    assert _ask_recursive(Q.even(x*y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.even(x*x), Q.integer(x)) is None
    assert _ask_recursive(Q.even(x*(x + y)), Q.integer(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.even(x*(x + y)), Q.integer(x) & Q.even(y)) is None


@XFAIL
def test_evenness_in_ternary_integer_product_with_odd():
    # Tests that oddness inference is independent of term ordering.
    # Term ordering at the point of testing depends on SymPy's symbol order, so
    # we try to force a different order by modifying symbol names.
    assert _ask_recursive(Q.even(x*y*(y + z)), Q.integer(x) & Q.integer(y) & Q.odd(z)) is True
    assert _ask_recursive(Q.even(y*x*(x + z)), Q.integer(x) & Q.integer(y) & Q.odd(z)) is True


def test_evenness_in_ternary_integer_product_with_even():
    assert _ask_recursive(Q.even(x*y*(y + z)), Q.integer(x) & Q.integer(y) & Q.even(z)) is None


def test_extended_real():
    assert _ask_recursive(Q.extended_real(x), Q.positive_infinite(x)) is True
    assert _ask_recursive(Q.extended_real(x), Q.positive(x)) is True
    assert _ask_recursive(Q.extended_real(x), Q.zero(x)) is True
    assert _ask_recursive(Q.extended_real(x), Q.negative(x)) is True
    assert _ask_recursive(Q.extended_real(x), Q.negative_infinite(x)) is True

    assert _ask_recursive(Q.extended_real(-x), Q.positive(x)) is True
    assert _ask_recursive(Q.extended_real(-x), Q.negative(x)) is True

    assert _ask_recursive(Q.extended_real(x + S.Infinity), Q.real(x)) is True

    assert _ask_recursive(Q.extended_real(x), Q.infinite(x)) is None


@_both_exp_pow
def test_rational():
    assert _ask_recursive(Q.rational(x), Q.integer(x)) is True
    assert _ask_recursive(Q.rational(x), Q.irrational(x)) is False
    assert _ask_recursive(Q.rational(x), Q.real(x)) is None
    assert _ask_recursive(Q.rational(x), Q.positive(x)) is None
    assert _ask_recursive(Q.rational(x), Q.negative(x)) is None
    assert _ask_recursive(Q.rational(x), Q.nonzero(x)) is None
    assert _ask_recursive(Q.rational(x), ~Q.algebraic(x)) is False

    assert _ask_recursive(Q.rational(2*x), Q.rational(x)) is True
    assert _ask_recursive(Q.rational(2*x), Q.integer(x)) is True
    assert _ask_recursive(Q.rational(2*x), Q.even(x)) is True
    assert _ask_recursive(Q.rational(2*x), Q.odd(x)) is True
    assert _ask_recursive(Q.rational(2*x), Q.irrational(x)) is False

    assert _ask_recursive(Q.rational(x/2), Q.rational(x)) is True
    assert _ask_recursive(Q.rational(x/2), Q.integer(x)) is True
    assert _ask_recursive(Q.rational(x/2), Q.even(x)) is True
    assert _ask_recursive(Q.rational(x/2), Q.odd(x)) is True
    assert _ask_recursive(Q.rational(x/2), Q.irrational(x)) is False

    assert _ask_recursive(Q.rational(1/x), Q.rational(x) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(1/x), Q.integer(x) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(1/x), Q.even(x) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(1/x), Q.odd(x)) is True
    assert _ask_recursive(Q.rational(1/x), Q.irrational(x)) is False

    assert _ask_recursive(Q.rational(2/x), Q.rational(x) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(2/x), Q.integer(x) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(2/x), Q.even(x) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(2/x), Q.odd(x)) is True
    assert _ask_recursive(Q.rational(2/x), Q.irrational(x)) is False

    assert _ask_recursive(Q.rational(x), ~Q.algebraic(x)) is False

    # with multiple symbols
    assert _ask_recursive(Q.rational(x*y), Q.irrational(x) & Q.irrational(y)) is None
    assert _ask_recursive(Q.rational(y/x), Q.rational(x) & Q.rational(y) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(y/x), Q.integer(x) & Q.rational(y) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(y/x), Q.even(x) & Q.rational(y) & Q.nonzero(x)) is True
    assert _ask_recursive(Q.rational(y/x), Q.odd(x) & Q.rational(y)) is True
    assert _ask_recursive(Q.rational(y/x), Q.irrational(x) & Q.rational(y) & Q.nonzero(y)) is False

    for f in [exp, sin, tan, asin, atan, cos]:
        assert _ask_recursive(Q.rational(f(7))) is False
        assert _ask_recursive(Q.rational(f(7, evaluate=False))) is False
        assert _ask_recursive(Q.rational(f(0, evaluate=False))) is True
        assert _ask_recursive(Q.rational(f(x)), Q.rational(x)) is None
        assert _ask_recursive(Q.rational(f(x)), Q.rational(x) & Q.nonzero(x)) is False

    for g in [log, acos]:
        assert _ask_recursive(Q.rational(g(7))) is False
        assert _ask_recursive(Q.rational(g(7, evaluate=False))) is False
        assert _ask_recursive(Q.rational(g(1, evaluate=False))) is True
        assert _ask_recursive(Q.rational(g(x)), Q.rational(x)) is None
        assert _ask_recursive(Q.rational(g(x)), Q.rational(x) & Q.nonzero(x - 1)) is False

    for h in [cot, acot]:
        assert _ask_recursive(Q.rational(h(7))) is False
        assert _ask_recursive(Q.rational(h(7, evaluate=False))) is False
        assert _ask_recursive(Q.rational(h(x)), Q.rational(x)) is False

    # https://github.com/sympy/sympy/issues/27442
    assert _ask_recursive(Q.rational(x**y),Q.irrational(x) & Q.rational(y)) is None
    assert _ask_recursive(Q.rational(x**y),Q.integer(x) & Q.prime(x) & Q.rational(y)) is None
    assert _ask_recursive(Q.rational(x**y),Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.rational(x**y),Q.eq(x,1) & Q.rational(y)) is True
    assert _ask_recursive(Q.rational(x**y),Q.eq(x,-1) & Q.rational(y)) is None
    assert _ask_recursive(Q.rational(x**y), Q.prime(x) & Q.rational(y)) is None
    assert _ask_recursive(Q.rational(x**y), ~Q.rational(x) & Q.integer(y) ) is None
    assert _ask_recursive(Q.rational(Pow(-1, x, evaluate=False), Q.rational(x))) is None
    assert _ask_recursive(Q.rational(x**y), Q.integer(y) & ~Q. algebraic(x)) is None
    assert _ask_recursive(Q.rational(x**y), Q.integer(y) & ~Q. algebraic(x) & ~Q.zero(x)) is None
    assert _ask_recursive(Q.rational(x**y), Q.integer(y) & ~Q.algebraic(x) & Q.complex(x) & ~Q.real(x)) is None
    assert _ask_recursive(Q.rational(x**y), Q.integer(y) & ~Q.algebraic(x) & Q.complex(x)) is None


def test_hermitian():
    assert _ask_recursive(Q.hermitian(x)) is None
    assert _ask_recursive(Q.hermitian(x), Q.antihermitian(x)) is None
    assert _ask_recursive(Q.hermitian(x), Q.imaginary(x)) is False
    assert _ask_recursive(Q.hermitian(x), Q.prime(x)) is True
    assert _ask_recursive(Q.hermitian(x), Q.real(x)) is True
    assert _ask_recursive(Q.hermitian(x), Q.zero(x)) is True

    assert _ask_recursive(Q.hermitian(x + 1), Q.antihermitian(x)) is None
    assert _ask_recursive(Q.hermitian(x + 1), Q.complex(x)) is None
    assert _ask_recursive(Q.hermitian(x + 1), Q.hermitian(x)) is True
    assert _ask_recursive(Q.hermitian(x + 1), Q.imaginary(x)) is False
    assert _ask_recursive(Q.hermitian(x + 1), Q.real(x)) is True
    assert _ask_recursive(Q.hermitian(x + I), Q.antihermitian(x)) is None
    assert _ask_recursive(Q.hermitian(x + I), Q.complex(x)) is None
    assert _ask_recursive(Q.hermitian(x + I), Q.hermitian(x)) is False
    assert _ask_recursive(Q.hermitian(x + I), Q.imaginary(x)) is None
    assert _ask_recursive(Q.hermitian(x + I), Q.real(x)) is False
    assert _ask_recursive(
        Q.hermitian(x + y), Q.antihermitian(x) & Q.antihermitian(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.antihermitian(x) & Q.complex(y)) is None
    assert _ask_recursive(
        Q.hermitian(x + y), Q.antihermitian(x) & Q.hermitian(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.antihermitian(x) & Q.imaginary(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.antihermitian(x) & Q.real(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.hermitian(x) & Q.complex(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.hermitian(x) & Q.hermitian(y)) is True
    assert _ask_recursive(Q.hermitian(x + y), Q.hermitian(x) & Q.imaginary(y)) is False
    assert _ask_recursive(Q.hermitian(x + y), Q.hermitian(x) & Q.real(y)) is True
    assert _ask_recursive(Q.hermitian(x + y), Q.imaginary(x) & Q.complex(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.imaginary(x) & Q.imaginary(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.imaginary(x) & Q.real(y)) is False
    assert _ask_recursive(Q.hermitian(x + y), Q.real(x) & Q.complex(y)) is None
    assert _ask_recursive(Q.hermitian(x + y), Q.real(x) & Q.real(y)) is True

    assert _ask_recursive(Q.hermitian(I*x), Q.antihermitian(x)) is True
    assert _ask_recursive(Q.hermitian(I*x), Q.complex(x)) is None
    assert _ask_recursive(Q.hermitian(I*x), Q.hermitian(x)) is False
    assert _ask_recursive(Q.hermitian(I*x), Q.imaginary(x)) is True
    assert _ask_recursive(Q.hermitian(I*x), Q.real(x)) is False
    assert _ask_recursive(Q.hermitian(x*y), Q.hermitian(x) & Q.real(y)) is True

    assert _ask_recursive(
        Q.hermitian(x + y + z), Q.real(x) & Q.real(y) & Q.real(z)) is True
    assert _ask_recursive(Q.hermitian(x + y + z),
        Q.real(x) & Q.real(y) & Q.imaginary(z)) is False
    assert _ask_recursive(Q.hermitian(x + y + z),
        Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) is None
    assert _ask_recursive(Q.hermitian(x + y + z),
        Q.imaginary(x) & Q.imaginary(y) & Q.imaginary(z)) is None

    assert _ask_recursive(Q.antihermitian(x)) is None
    assert _ask_recursive(Q.antihermitian(x), Q.real(x)) is False
    assert _ask_recursive(Q.antihermitian(x), Q.prime(x)) is False

    assert _ask_recursive(Q.antihermitian(x + 1), Q.antihermitian(x)) is False
    assert _ask_recursive(Q.antihermitian(x + 1), Q.complex(x)) is None
    assert _ask_recursive(Q.antihermitian(x + 1), Q.hermitian(x)) is None
    assert _ask_recursive(Q.antihermitian(x + 1), Q.imaginary(x)) is False
    assert _ask_recursive(Q.antihermitian(x + 1), Q.real(x)) is None
    assert _ask_recursive(Q.antihermitian(x + I), Q.antihermitian(x)) is True
    assert _ask_recursive(Q.antihermitian(x + I), Q.complex(x)) is None
    assert _ask_recursive(Q.antihermitian(x + I), Q.hermitian(x)) is None
    assert _ask_recursive(Q.antihermitian(x + I), Q.imaginary(x)) is True
    assert _ask_recursive(Q.antihermitian(x + I), Q.real(x)) is False
    assert _ask_recursive(Q.antihermitian(x), Q.zero(x)) is True

    assert _ask_recursive(
        Q.antihermitian(x + y), Q.antihermitian(x) & Q.antihermitian(y)
    ) is True
    assert _ask_recursive(
        Q.antihermitian(x + y), Q.antihermitian(x) & Q.complex(y)) is None
    assert _ask_recursive(
        Q.antihermitian(x + y), Q.antihermitian(x) & Q.hermitian(y)) is None
    assert _ask_recursive(
        Q.antihermitian(x + y), Q.antihermitian(x) & Q.imaginary(y)) is True
    assert _ask_recursive(Q.antihermitian(x + y), Q.antihermitian(x) & Q.real(y)
        ) is False
    assert _ask_recursive(Q.antihermitian(x + y), Q.hermitian(x) & Q.complex(y)) is None
    assert _ask_recursive(Q.antihermitian(x + y), Q.hermitian(x) & Q.hermitian(y)
        ) is None
    assert _ask_recursive(
        Q.antihermitian(x + y), Q.hermitian(x) & Q.imaginary(y)) is None
    assert _ask_recursive(Q.antihermitian(x + y), Q.hermitian(x) & Q.real(y)) is None
    assert _ask_recursive(Q.antihermitian(x + y), Q.imaginary(x) & Q.complex(y)) is None
    assert _ask_recursive(Q.antihermitian(x + y), Q.imaginary(x) & Q.imaginary(y)) is True
    assert _ask_recursive(Q.antihermitian(x + y), Q.imaginary(x) & Q.real(y)) is False
    assert _ask_recursive(Q.antihermitian(x + y), Q.real(x) & Q.complex(y)) is None
    assert _ask_recursive(Q.antihermitian(x + y), Q.real(x) & Q.real(y)) is None

    assert _ask_recursive(Q.antihermitian(I*x), Q.real(x)) is True
    assert _ask_recursive(Q.antihermitian(I*x), Q.antihermitian(x)) is False
    assert _ask_recursive(Q.antihermitian(I*x), Q.complex(x)) is None
    assert _ask_recursive(Q.antihermitian(x*y), Q.antihermitian(x) & Q.real(y)) is True

    assert _ask_recursive(Q.antihermitian(x + y + z),
        Q.real(x) & Q.real(y) & Q.real(z)) is None
    assert _ask_recursive(Q.antihermitian(x + y + z),
        Q.real(x) & Q.real(y) & Q.imaginary(z)) is None
    assert _ask_recursive(Q.antihermitian(x + y + z),
        Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) is False
    assert _ask_recursive(Q.antihermitian(x + y + z),
        Q.imaginary(x) & Q.imaginary(y) & Q.imaginary(z)) is True


@_both_exp_pow
def test_imaginary():
    assert _ask_recursive(Q.imaginary(x)) is None
    assert _ask_recursive(Q.imaginary(x), Q.real(x)) is False
    assert _ask_recursive(Q.imaginary(x), Q.prime(x)) is False

    assert _ask_recursive(Q.imaginary(x + 1), Q.real(x)) is False
    assert _ask_recursive(Q.imaginary(x + 1), Q.imaginary(x)) is False
    assert _ask_recursive(Q.imaginary(x + I), Q.real(x)) is False
    assert _ask_recursive(Q.imaginary(x + I), Q.imaginary(x)) is True
    assert _ask_recursive(Q.imaginary(x + y), Q.imaginary(x) & Q.imaginary(y)) is True
    assert _ask_recursive(Q.imaginary(x + y), Q.real(x) & Q.real(y)) is False
    assert _ask_recursive(Q.imaginary(x + y), Q.imaginary(x) & Q.real(y)) is False
    assert _ask_recursive(Q.imaginary(x + y), Q.complex(x) & Q.real(y)) is None
    assert _ask_recursive(
        Q.imaginary(x + y + z), Q.real(x) & Q.real(y) & Q.real(z)) is False
    assert _ask_recursive(Q.imaginary(x + y + z),
        Q.real(x) & Q.real(y) & Q.imaginary(z)) is None
    assert _ask_recursive(Q.imaginary(x + y + z),
        Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) is False

    assert _ask_recursive(Q.imaginary(I*x), Q.real(x)) is True
    assert _ask_recursive(Q.imaginary(I*x), Q.imaginary(x)) is False
    assert _ask_recursive(Q.imaginary(I*x), Q.complex(x)) is None
    assert _ask_recursive(Q.imaginary(x*y), Q.imaginary(x) & Q.real(y)) is True
    assert _ask_recursive(Q.imaginary(x*y), Q.real(x) & Q.real(y)) is False

    assert _ask_recursive(Q.imaginary(I**x), Q.negative(x)) is None
    assert _ask_recursive(Q.imaginary(I**x), Q.positive(x)) is None
    assert _ask_recursive(Q.imaginary(I**x), Q.even(x)) is False
    assert _ask_recursive(Q.imaginary(I**x), Q.odd(x)) is True
    assert _ask_recursive(Q.imaginary(I**x), Q.imaginary(x)) is False
    assert _ask_recursive(Q.imaginary((2*I)**x), Q.imaginary(x)) is False
    assert _ask_recursive(Q.imaginary(x**0), Q.imaginary(x)) is False
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(x) & Q.imaginary(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(x) & Q.real(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.real(x) & Q.imaginary(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.real(x) & Q.real(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(y) & Q.integer(x)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(x) & Q.rational(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.imaginary(x) & Q.even(y)) is False

    assert _ask_recursive(Q.imaginary(x**y), Q.real(x) & Q.integer(y)) is False
    assert _ask_recursive(Q.imaginary(x**y), Q.positive(x) & Q.real(y)) is False
    assert _ask_recursive(Q.imaginary(x**y), Q.negative(x) & Q.real(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.negative(x) & Q.real(y) & ~Q.rational(y)) is False
    assert _ask_recursive(Q.imaginary(x**y), Q.integer(x) & Q.imaginary(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.negative(x) & Q.rational(y) & Q.integer(2*y)) is True
    assert _ask_recursive(Q.imaginary(x**y), Q.negative(x) & Q.rational(y) & ~Q.integer(2*y)) is False
    assert _ask_recursive(Q.imaginary(x**y), Q.negative(x) & Q.rational(y)) is None
    assert _ask_recursive(Q.imaginary(x**y), Q.real(x) & Q.rational(y) & ~Q.integer(2*y)) is False
    assert _ask_recursive(Q.imaginary(x**y), Q.real(x) & Q.rational(y) & Q.integer(2*y)) is None

    # logarithm
    assert _ask_recursive(Q.imaginary(log(I))) is True
    assert _ask_recursive(Q.imaginary(log(2*I))) is False
    assert _ask_recursive(Q.imaginary(log(I + 1))) is False
    assert _ask_recursive(Q.imaginary(log(x)), Q.complex(x)) is None
    assert _ask_recursive(Q.imaginary(log(x)), Q.imaginary(x)) is None
    assert _ask_recursive(Q.imaginary(log(x)), Q.positive(x)) is False
    assert _ask_recursive(Q.imaginary(log(exp(x))), Q.complex(x)) is None
    assert _ask_recursive(Q.imaginary(log(exp(x))), Q.imaginary(x)) is None  # zoo/I/a+I*b
    assert _ask_recursive(Q.imaginary(log(exp(I)))) is True

    # exponential
    assert _ask_recursive(Q.imaginary(exp(x)**x), Q.imaginary(x)) is False
    eq = Pow(exp(pi*I*x, evaluate=False), x, evaluate=False)
    assert _ask_recursive(Q.imaginary(eq), Q.even(x)) is False
    eq = Pow(exp(pi*I*x/2, evaluate=False), x, evaluate=False)
    assert _ask_recursive(Q.imaginary(eq), Q.odd(x)) is True
    assert _ask_recursive(Q.imaginary(exp(3*I*pi*x)**x), Q.integer(x)) is False
    assert _ask_recursive(Q.imaginary(exp(2*pi*I, evaluate=False))) is False
    assert _ask_recursive(Q.imaginary(exp(pi*I/2, evaluate=False))) is True

    # issue 7886
    assert _ask_recursive(Q.imaginary(Pow(x, Rational(1, 4))), Q.real(x) & Q.negative(x)) is False


def test_integer():
    assert _ask_recursive(Q.integer(x)) is None
    assert _ask_recursive(Q.integer(x), Q.integer(x)) is True
    assert _ask_recursive(Q.integer(x), ~Q.integer(x)) is False
    assert _ask_recursive(Q.integer(x), ~Q.real(x)) is False
    assert _ask_recursive(Q.integer(x), ~Q.positive(x)) is None
    assert ask(Q.integer(x), Q.even(x) | Q.odd(x)) is True

    assert _ask_recursive(Q.integer(2*x), Q.integer(x)) is True
    assert _ask_recursive(Q.integer(2*x), Q.even(x)) is True
    assert _ask_recursive(Q.integer(2*x), Q.prime(x)) is True
    assert _ask_recursive(Q.integer(2*x), Q.rational(x)) is None
    assert _ask_recursive(Q.integer(2*x), Q.real(x)) is None
    assert _ask_recursive(Q.integer(sqrt(2)*x), Q.integer(x)) is False
    assert _ask_recursive(Q.integer(sqrt(2)*x), Q.irrational(x)) is None

    assert _ask_recursive(Q.integer(x/2), Q.odd(x)) is False
    assert _ask_recursive(Q.integer(x/2), Q.even(x)) is True
    assert _ask_recursive(Q.integer(x/3), Q.odd(x)) is None
    assert _ask_recursive(Q.integer(x/3), Q.even(x)) is None

    # https://github.com/sympy/sympy/issues/7286
    assert _ask_recursive(Q.integer(Abs(x)),Q.integer(x)) is True
    assert _ask_recursive(Q.integer(Abs(-x)),Q.integer(x)) is True
    assert _ask_recursive(Q.integer(Abs(x)), ~Q.integer(x)) is None
    assert _ask_recursive(Q.integer(Abs(x)),Q.complex(x)) is None
    assert _ask_recursive(Q.integer(Abs(x+I*y)),Q.real(x) & Q.real(y)) is None

    # https://github.com/sympy/sympy/issues/27739
    assert _ask_recursive(Q.integer(x/y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.integer(1/x), Q.integer(x)) is None
    assert _ask_recursive(Q.integer(x**y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.integer(sqrt(5))) is False
    assert _ask_recursive(Q.integer(x**y), Q.nonzero(x) & Q.zero(y)) is True
    assert _ask_recursive(Q.integer(x**y), Q.integer(x) & Q.integer(y) & Q.positive(y)) is True
    assert _ask_recursive(Q.integer(-1**x), Q.integer(x)) is True
    assert _ask_recursive(Q.integer(x**y), Q.integer(x) & Q.integer(y) & Q.positive(y)) is True
    assert _ask_recursive(Q.integer(x**y), Q.zero(x) & Q.integer(y) & Q.positive(y)) is True
    assert _ask_recursive(Q.integer(pi**x), Q.zero(x)) is True
    assert _ask_recursive(Q.integer(x**y), Q.imaginary(x) & Q.zero(y)) is True
    assert _ask_recursive(Q.integer(x**y), Q.integer(x) & Q.integer(y) & ~Q.negative(y)) is True
    assert _ask_recursive(Q.integer(x**y), Q.integer(x) & Q.integer(y) & Q.negative(y) & ~Q.zero(x - 1) & ~Q.zero(x + 1)) is False


def test_negative():
    assert _ask_recursive(Q.negative(x), Q.negative(x)) is True
    assert _ask_recursive(Q.negative(x), Q.positive(x)) is False
    assert _ask_recursive(Q.negative(x), ~Q.real(x)) is False
    assert _ask_recursive(Q.negative(x), Q.prime(x)) is False
    assert _ask_recursive(Q.negative(x), ~Q.prime(x)) is None

    assert _ask_recursive(Q.negative(-x), Q.positive(x)) is True
    assert _ask_recursive(Q.negative(-x), ~Q.positive(x)) is None
    assert _ask_recursive(Q.negative(-x), Q.negative(x)) is False
    assert _ask_recursive(Q.negative(-x), Q.positive(x)) is True

    assert _ask_recursive(Q.negative(x - 1), Q.negative(x)) is True
    assert _ask_recursive(Q.negative(x + y)) is None
    assert _ask_recursive(Q.negative(x + y), Q.negative(x)) is None
    assert _ask_recursive(Q.negative(x + y), Q.negative(x) & Q.negative(y)) is True
    assert _ask_recursive(Q.negative(x + y), Q.negative(x) & Q.nonpositive(y)) is True
    assert _ask_recursive(Q.negative(2 + I)) is False
    # although this could be False, it is representative of expressions
    # that don't evaluate to a zero with precision
    assert _ask_recursive(Q.negative(cos(I)**2 + sin(I)**2 - 1)) is None
    assert _ask_recursive(Q.negative(-I + I*(cos(2)**2 + sin(2)**2))) is None

    assert _ask_recursive(Q.negative(x**2)) is None
    assert _ask_recursive(Q.negative(x**2), Q.real(x)) is False
    assert _ask_recursive(Q.negative(x**1.4), Q.real(x)) is None

    assert _ask_recursive(Q.negative(x**I), Q.positive(x)) is None

    assert _ask_recursive(Q.negative(x*y)) is None
    assert _ask_recursive(Q.negative(x*y), Q.positive(x) & Q.positive(y)) is False
    assert _ask_recursive(Q.negative(x*y), Q.positive(x) & Q.negative(y)) is True
    assert _ask_recursive(Q.negative(x*y), Q.complex(x) & Q.complex(y)) is None

    assert _ask_recursive(Q.negative(x**y)) is None
    assert _ask_recursive(Q.negative(x**y), Q.negative(x) & Q.even(y)) is False
    assert _ask_recursive(Q.negative(x**y), Q.negative(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.negative(x**y), Q.positive(x) & Q.integer(y)) is False

    assert _ask_recursive(Q.negative(Abs(x))) is False


def test_nonzero():
    assert _ask_recursive(Q.nonzero(x)) is None
    assert _ask_recursive(Q.nonzero(x), Q.real(x)) is None
    assert _ask_recursive(Q.nonzero(x), Q.positive(x)) is True
    assert _ask_recursive(Q.nonzero(x), Q.negative(x)) is True
    assert ask(Q.nonzero(x), Q.negative(x) | Q.positive(x)) is True

    assert _ask_recursive(Q.nonzero(x + y)) is None
    assert _ask_recursive(Q.nonzero(x + y), Q.positive(x) & Q.positive(y)) is True
    assert _ask_recursive(Q.nonzero(x + y), Q.positive(x) & Q.negative(y)) is None
    assert _ask_recursive(Q.nonzero(x + y), Q.negative(x) & Q.negative(y)) is True

    assert _ask_recursive(Q.nonzero(2*x)) is None
    assert _ask_recursive(Q.nonzero(2*x), Q.positive(x)) is True
    assert _ask_recursive(Q.nonzero(2*x), Q.negative(x)) is True
    assert _ask_recursive(Q.nonzero(x*y), Q.nonzero(x)) is None
    assert _ask_recursive(Q.nonzero(x*y), Q.nonzero(x) & Q.nonzero(y)) is True

    # https://github.com/sympy/sympy/pull/29225
    assert _ask_recursive(Q.nonzero(Pow(x, y), Q.nonzero(x))) is None
    assert _ask_recursive(Q.nonzero(Pow(x, y)), Q.positive(x) & Q.real(y)) is True
    assert _ask_recursive(Q.nonzero(Pow(5, 2*I*n*pi)), Q.integer(n)) is False
    assert _ask_recursive(Q.nonzero(Pow(0, x)), Q.positive(x)) is False
    assert _ask_recursive(Q.nonzero(Pow(5, I*n), Q.integer(n))) is None
    assert _ask_recursive(Q.nonzero(Pow(-1, x)), Q.real(x)) is None
    assert _ask_recursive(Q.nonzero(Pow(x, x)), Q.zero(x)) is True
    assert _ask_recursive(Q.nonzero(Pow(I, x)), Q.zero(x)) is True

    assert _ask_recursive(Q.nonzero(Abs(x))) is None
    assert _ask_recursive(Q.nonzero(Abs(x)), Q.nonzero(x)) is True

    assert _ask_recursive(Q.nonzero(log(exp(2*I)))) is False
    # although this could be False, it is representative of expressions
    # that don't evaluate to a zero with precision
    assert _ask_recursive(Q.nonzero(cos(1)**2 + sin(1)**2 - 1)) is None


def test_zero():
    assert _ask_recursive(Q.zero(x)) is None
    assert _ask_recursive(Q.zero(x), Q.real(x)) is None
    assert _ask_recursive(Q.zero(x), Q.positive(x)) is False
    assert _ask_recursive(Q.zero(x), Q.negative(x)) is False
    assert ask(Q.zero(x), Q.negative(x) | Q.positive(x)) is False

    assert ask(Q.zero(x), Q.nonnegative(x) & Q.nonpositive(x)) is True

    assert _ask_recursive(Q.zero(x + y)) is None
    assert _ask_recursive(Q.zero(x + y), Q.positive(x) & Q.positive(y)) is False
    assert _ask_recursive(Q.zero(x + y), Q.positive(x) & Q.negative(y)) is None
    assert _ask_recursive(Q.zero(x + y), Q.negative(x) & Q.negative(y)) is False

    assert _ask_recursive(Q.zero(2*x)) is None
    assert _ask_recursive(Q.zero(2*x), Q.positive(x)) is False
    assert _ask_recursive(Q.zero(2*x), Q.negative(x)) is False
    assert _ask_recursive(Q.zero(x*y), Q.nonzero(x)) is None

    assert _ask_recursive(Q.zero(Abs(x))) is None
    assert _ask_recursive(Q.zero(Abs(x)), Q.zero(x)) is True

    assert _ask_recursive(Q.integer(x), Q.zero(x)) is True
    assert _ask_recursive(Q.even(x), Q.zero(x)) is True
    assert _ask_recursive(Q.odd(x), Q.zero(x)) is False
    assert _ask_recursive(Q.zero(x), Q.even(x)) is None
    assert _ask_recursive(Q.zero(x), Q.odd(x)) is False
    assert ask(Q.zero(x) | Q.zero(y), Q.zero(x*y)) is True


def test_odd_query():
    assert _ask_recursive(Q.odd(x)) is None
    assert _ask_recursive(Q.odd(x), Q.odd(x)) is True
    assert _ask_recursive(Q.odd(x), Q.integer(x)) is None
    assert _ask_recursive(Q.odd(x), ~Q.integer(x)) is False
    assert _ask_recursive(Q.odd(x), Q.rational(x)) is None
    assert _ask_recursive(Q.odd(x), Q.positive(x)) is None

    assert _ask_recursive(Q.odd(-x), Q.odd(x)) is True

    assert _ask_recursive(Q.odd(2*x)) is None
    assert _ask_recursive(Q.odd(2*x), Q.integer(x)) is False
    assert _ask_recursive(Q.odd(2*x), Q.odd(x)) is False
    assert _ask_recursive(Q.odd(2*x), Q.irrational(x)) is False
    assert _ask_recursive(Q.odd(2*x), ~Q.integer(x)) is None
    assert _ask_recursive(Q.odd(3*x), Q.integer(x)) is None

    assert _ask_recursive(Q.odd(x/3), Q.odd(x)) is None
    assert _ask_recursive(Q.odd(x/3), Q.even(x)) is None

    assert _ask_recursive(Q.odd(x + 1), Q.even(x)) is True
    assert _ask_recursive(Q.odd(x + 2), Q.even(x)) is False
    assert _ask_recursive(Q.odd(x + 2), Q.odd(x)) is True
    assert _ask_recursive(Q.odd(3 - x), Q.odd(x)) is False
    assert _ask_recursive(Q.odd(3 - x), Q.even(x)) is True
    assert _ask_recursive(Q.odd(3 + x), Q.odd(x)) is False
    assert _ask_recursive(Q.odd(3 + x), Q.even(x)) is True
    assert _ask_recursive(Q.odd(x + y), Q.odd(x) & Q.odd(y)) is False
    assert _ask_recursive(Q.odd(x + y), Q.odd(x) & Q.even(y)) is True
    assert _ask_recursive(Q.odd(x - y), Q.even(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.odd(x - y), Q.odd(x) & Q.odd(y)) is False

    assert _ask_recursive(Q.odd(x + y + z), Q.odd(x) & Q.odd(y) & Q.even(z)) is False
    assert _ask_recursive(Q.odd(x + y + z + t),
        Q.odd(x) & Q.odd(y) & Q.even(z) & Q.integer(t)) is None

    assert _ask_recursive(Q.odd(2*x + 1), Q.integer(x)) is True
    assert _ask_recursive(Q.odd(2*x + y), Q.integer(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.odd(2*x + y), Q.integer(x) & Q.even(y)) is False
    assert _ask_recursive(Q.odd(2*x + y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.odd(x*y), Q.odd(x) & Q.even(y)) is False
    assert _ask_recursive(Q.odd(x*y), Q.odd(x) & Q.odd(y)) is True
    assert _ask_recursive(Q.odd(2*x*y), Q.rational(x) & Q.rational(x)) is None
    assert _ask_recursive(Q.odd(2*x*y), Q.irrational(x) & Q.irrational(x)) is None

    assert _ask_recursive(Q.odd(Abs(x)), Q.odd(x)) is True

    assert _ask_recursive(Q.odd((-1)**n), Q.integer(n)) is True

    assert _ask_recursive(Q.odd(k**2), Q.even(k)) is False
    assert _ask_recursive(Q.odd(n**2), Q.odd(n)) is True
    assert _ask_recursive(Q.odd(3**k), Q.even(k)) is None

    assert _ask_recursive(Q.odd(k**m), Q.even(k) & Q.integer(m) & ~Q.negative(m)) is None
    assert _ask_recursive(Q.odd(n**m), Q.odd(n) & Q.integer(m) & ~Q.negative(m)) is True

    assert _ask_recursive(Q.odd(k**p), Q.even(k) & Q.integer(p) & Q.positive(p)) is False
    assert _ask_recursive(Q.odd(n**p), Q.odd(n) & Q.integer(p) & Q.positive(p)) is True

    assert _ask_recursive(Q.odd(m**k), Q.even(k) & Q.integer(m) & ~Q.negative(m)) is None
    assert _ask_recursive(Q.odd(p**k), Q.even(k) & Q.integer(p) & Q.positive(p)) is None

    assert _ask_recursive(Q.odd(m**n), Q.odd(n) & Q.integer(m) & ~Q.negative(m)) is None
    assert _ask_recursive(Q.odd(p**n), Q.odd(n) & Q.integer(p) & Q.positive(p)) is None

    assert _ask_recursive(Q.odd(k**x), Q.even(k)) is None
    assert _ask_recursive(Q.odd(n**x), Q.odd(n)) is None

    assert _ask_recursive(Q.odd(x*y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.odd(x*x), Q.integer(x)) is None
    assert _ask_recursive(Q.odd(x*(x + y)), Q.integer(x) & Q.odd(y)) is False
    assert _ask_recursive(Q.odd(x*(x + y)), Q.integer(x) & Q.even(y)) is None


@XFAIL
def test_oddness_in_ternary_integer_product_with_odd():
    # Tests that oddness inference is independent of term ordering.
    # Term ordering at the point of testing depends on SymPy's symbol order, so
    # we try to force a different order by modifying symbol names.
    assert _ask_recursive(Q.odd(x*y*(y + z)), Q.integer(x) & Q.integer(y) & Q.odd(z)) is False
    assert _ask_recursive(Q.odd(y*x*(x + z)), Q.integer(x) & Q.integer(y) & Q.odd(z)) is False


def test_oddness_in_ternary_integer_product_with_even():
    assert _ask_recursive(Q.odd(x*y*(y + z)), Q.integer(x) & Q.integer(y) & Q.even(z)) is None


def test_prime():
    assert _ask_recursive(Q.prime(x), Q.prime(x)) is True
    assert _ask_recursive(Q.prime(x), ~Q.prime(x)) is False
    assert _ask_recursive(Q.prime(x), Q.integer(x)) is None
    assert _ask_recursive(Q.prime(x), ~Q.integer(x)) is False

    assert _ask_recursive(Q.prime(2*x), Q.integer(x)) is None
    assert _ask_recursive(Q.prime(x*y)) is None
    assert _ask_recursive(Q.prime(x*y), Q.prime(x)) is None
    assert _ask_recursive(Q.prime(x*y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.prime(4*x), Q.integer(x)) is False
    assert _ask_recursive(Q.prime(4*x)) is None

    assert _ask_recursive(Q.prime(x**2), Q.integer(x)) is False
    assert _ask_recursive(Q.prime(x**2), Q.prime(x)) is False

    # https://github.com/sympy/sympy/issues/27446
    assert _ask_recursive(Q.prime(4**x), Q.integer(x)) is False
    assert _ask_recursive(Q.prime(p**x), Q.prime(p) & Q.integer(x) & Q.ne(x, 1)) is False
    assert _ask_recursive(Q.prime(n**x), Q.integer(x) & Q.composite(n)) is False
    assert _ask_recursive(Q.prime(x**y), Q.integer(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.prime(2**x), Q.integer(x)) is None
    assert _ask_recursive(Q.prime(p**x), Q.prime(p) & Q.integer(x)) is None

    # Ideally, these should return True since the base is prime and the exponent is one,
    # but currently, they return None.
    assert _ask_recursive(Q.prime(x**y), Q.prime(x) & Q.eq(y,1)) is None
    assert _ask_recursive(Q.prime(x**y), Q.prime(x) & Q.integer(y) & Q.gt(y,0) & Q.lt(y,2)) is None

    assert _ask_recursive(Q.prime(Pow(x,1, evaluate=False)), Q.prime(x)) is True


@_both_exp_pow
def test_positive():
    assert _ask_recursive(Q.positive(cos(I) ** 2 + sin(I) ** 2 - 1)) is None
    assert _ask_recursive(Q.positive(x), Q.positive(x)) is True
    assert _ask_recursive(Q.positive(x), Q.negative(x)) is False
    assert _ask_recursive(Q.positive(x), Q.nonzero(x)) is None

    assert _ask_recursive(Q.positive(-x), Q.positive(x)) is False
    assert _ask_recursive(Q.positive(-x), Q.negative(x)) is True

    assert _ask_recursive(Q.positive(x + y), Q.positive(x) & Q.positive(y)) is True
    assert _ask_recursive(Q.positive(x + y), Q.positive(x) & Q.nonnegative(y)) is True
    assert _ask_recursive(Q.positive(x + y), Q.positive(x) & Q.negative(y)) is None
    assert _ask_recursive(Q.positive(x + y), Q.positive(x) & Q.imaginary(y)) is False

    assert _ask_recursive(Q.positive(2*x), Q.positive(x)) is True
    assumptions = Q.positive(x) & Q.negative(y) & Q.negative(z) & Q.positive(w)
    assert _ask_recursive(Q.positive(x*y*z)) is None
    assert _ask_recursive(Q.positive(x*y*z), assumptions) is True
    assert _ask_recursive(Q.positive(-x*y*z), assumptions) is False

    assert _ask_recursive(Q.positive(x**I), Q.positive(x)) is None

    assert _ask_recursive(Q.positive(x**2), Q.positive(x)) is True
    assert _ask_recursive(Q.positive(x**2), Q.negative(x)) is True
    assert _ask_recursive(Q.positive(x**3), Q.negative(x)) is False
    assert _ask_recursive(Q.positive(1/(1 + x**2)), Q.real(x)) is True
    assert _ask_recursive(Q.positive(2**I)) is False
    assert _ask_recursive(Q.positive(2 + I)) is False
    # although this could be False, it is representative of expressions
    # that don't evaluate to a zero with precision
    assert _ask_recursive(Q.positive(cos(I)**2 + sin(I)**2 - 1)) is None
    assert _ask_recursive(Q.positive(-I + I*(cos(2)**2 + sin(2)**2))) is None

    #exponential
    assert _ask_recursive(Q.positive(exp(x)), Q.real(x)) is True
    assert _ask_recursive(Q.positive(x**y), Q.nonzero(x) & Q.even(y)) is True
    assert _ask_recursive(~Q.negative(exp(x)), Q.real(x)) is True
    assert _ask_recursive(Q.positive(x + exp(x)), Q.real(x)) is None
    assert _ask_recursive(Q.positive(exp(x)), Q.imaginary(x)) is None
    assert _ask_recursive(Q.positive(exp(2*pi*I, evaluate=False)), Q.imaginary(x)) is True
    assert _ask_recursive(Q.negative(exp(pi*I, evaluate=False)), Q.imaginary(x)) is True
    assert _ask_recursive(Q.positive(exp(x*pi*I)), Q.even(x)) is True
    assert _ask_recursive(Q.positive(exp(x*pi*I)), Q.odd(x)) is False
    assert _ask_recursive(Q.positive(x**y), Q.zero(x) & Q.positive(y)) is None
    assert _ask_recursive(Q.positive(exp(x*pi*I)), Q.real(x)) is None
    assert _ask_recursive(Q.positive(x**2), Q.real(x)) is None
    assert _ask_recursive(Q.positive(x**y), Q.even(y) & Q.real(x)) is None
    assert _ask_recursive(Q.positive(x**y),Q.zero(x)) is None
    assert _ask_recursive(Q.positive(x**y), Q.zero(x) & Q.negative(y)) is None
    assert _ask_recursive(Q.positive(x**y), Q.zero(x) & Q.even(y)) is None

    # logarithm
    assert _ask_recursive(Q.positive(log(x)), Q.imaginary(x)) is False
    assert _ask_recursive(Q.positive(log(x)), Q.negative(x)) is False
    assert _ask_recursive(Q.positive(log(x)), Q.positive(x)) is None
    assert _ask_recursive(Q.positive(log(x + 2)), Q.positive(x)) is True

    # factorial
    assert _ask_recursive(Q.positive(factorial(x)), Q.integer(x) & Q.positive(x)) is True
    assert _ask_recursive(Q.positive(factorial(x)), Q.integer(x)) is None

    #absolute value
    assert _ask_recursive(Q.positive(Abs(x))) is None  # Abs(0) = 0
    assert _ask_recursive(Q.positive(Abs(x)), Q.positive(x)) is True


def test_nonpositive():
    assert _ask_recursive(Q.nonpositive(-1)) is True
    assert _ask_recursive(Q.nonpositive(0)) is True
    assert _ask_recursive(Q.nonpositive(1)) is False
    assert _ask_recursive(~Q.positive(x), Q.nonpositive(x)) is True
    assert _ask_recursive(Q.nonpositive(x), Q.positive(x)) is False
    assert _ask_recursive(Q.nonpositive(sqrt(-1))) is False
    assert _ask_recursive(Q.nonpositive(x), Q.imaginary(x)) is False


def test_nonnegative():
    assert _ask_recursive(Q.nonnegative(-1)) is False
    assert _ask_recursive(Q.nonnegative(0)) is True
    assert _ask_recursive(Q.nonnegative(1)) is True
    assert _ask_recursive(~Q.negative(x), Q.nonnegative(x)) is True
    assert _ask_recursive(Q.nonnegative(x), Q.negative(x)) is False
    assert _ask_recursive(Q.nonnegative(sqrt(-1))) is False
    assert _ask_recursive(Q.nonnegative(x), Q.imaginary(x)) is False

def test_real_basic():
    assert _ask_recursive(Q.real(x)) is None
    assert _ask_recursive(Q.real(x), Q.real(x)) is True
    assert _ask_recursive(Q.real(x), Q.nonzero(x)) is True
    assert _ask_recursive(Q.real(x), Q.positive(x)) is True
    assert _ask_recursive(Q.real(x), Q.negative(x)) is True
    assert _ask_recursive(Q.real(x), Q.integer(x)) is True
    assert _ask_recursive(Q.real(x), Q.even(x)) is True
    assert _ask_recursive(Q.real(x), Q.prime(x)) is True

    assert _ask_recursive(Q.real(x/sqrt(2)), Q.real(x)) is True
    assert _ask_recursive(Q.real(x/sqrt(-2)), Q.real(x)) is None # x can be zero

    assert _ask_recursive(Q.real(x + 1), Q.real(x)) is True
    assert _ask_recursive(Q.real(x + I), Q.real(x)) is False
    assert _ask_recursive(Q.real(x + I), Q.complex(x)) is None

    assert _ask_recursive(Q.real(2*x), Q.real(x)) is True
    assert _ask_recursive(Q.real(I*x), Q.real(x)) is None
    assert _ask_recursive(Q.real(I*x), Q.imaginary(x)) is True
    assert _ask_recursive(Q.real(I*x), Q.complex(x)) is None


def test_real_pow():
    assert _ask_recursive(Q.real(x**2), Q.real(x)) is True
    assert _ask_recursive(Q.real(sqrt(x)), Q.negative(x)) is False
    assert _ask_recursive(Q.real(x**y), Q.real(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.real(x**y), Q.real(x) & Q.real(y)) is None
    assert _ask_recursive(Q.real(x**y), Q.positive(x) & Q.real(y)) is True
    assert _ask_recursive(Q.real(x**y), Q.imaginary(x) & Q.imaginary(y)) is None  # I**I or (2*I)**I
    assert _ask_recursive(Q.real(x**y), Q.imaginary(x) & Q.real(y)) is None  # I**1 or I**0
    assert _ask_recursive(Q.real(x**y), Q.real(x) & Q.imaginary(y)) is None  # could be exp(2*pi*I) or 2**I
    assert _ask_recursive(Q.real(x**0), Q.imaginary(x)) is True
    assert _ask_recursive(Q.real(x**y), Q.positive(x) & Q.real(y)) is True
    assert _ask_recursive(Q.real(x**y), Q.real(x) & Q.rational(y)) is None
    assert _ask_recursive(Q.real(x**y), Q.imaginary(x) & Q.integer(y)) is None
    assert _ask_recursive(Q.real(x**y), Q.imaginary(x) & Q.odd(y)) is False
    assert _ask_recursive(Q.real(x**y), Q.imaginary(x) & Q.even(y)) is True
    assert _ask_recursive(Q.real(x**(y/z)), Q.real(x) & Q.real(y/z) & Q.rational(y/z) & Q.even(z) & Q.positive(x)) is True
    assert _ask_recursive(Q.real(x**(y/z)), Q.real(x) & Q.rational(y/z) & Q.even(z) & Q.negative(x)) is None
    assert _ask_recursive(Q.real(x**(y/z)), Q.real(x) & Q.integer(y/z)) is None
    assert _ask_recursive(Q.real(x**(y/z)), Q.real(x) & Q.real(y/z) & Q.positive(x)) is True
    assert _ask_recursive(Q.real(x**(y/z)), Q.real(x) & Q.real(y/z) & Q.negative(x)) is None
    assert _ask_recursive(Q.real((-I)**i), Q.imaginary(i)) is True
    assert _ask_recursive(Q.real(I**i), Q.imaginary(i)) is True
    assert _ask_recursive(Q.real(i**i), Q.imaginary(i)) is None  # i might be 2*I
    assert _ask_recursive(Q.real(x**i), Q.imaginary(i)) is None  # x could be 0
    assert _ask_recursive(Q.real(x**(I*pi/log(x))), Q.real(x)) is True

    # https://github.com/sympy/sympy/issues/27485
    assert _ask_recursive(Q.real(n**p), Q.negative(n) & Q.positive(p)) is None

    # https://github.com/sympy/sympy/issues/16530
    assert _ask_recursive(Q.real(1/Abs(x))) is None
    assert _ask_recursive(Q.real(x**y), Q.zero(x) & Q.real(y)) is None
    assert _ask_recursive(Q.real(x**y), Q.zero(x) & Q.positive(y)) is True

    # https://github.com/sympy/sympy/issues/28142
    assert _ask_recursive(Q.real(sqrt(x)), Q.real(x)) is None
    one_half = Mul(S.Half, 1, evaluate=False)
    assert _ask_recursive(Q.real(x ** one_half), Q.real(x)) is None
    assert _ask_recursive(Q.real(Pow(x, 1, evaluate=False)), Q.zero(x)) is True
    assert _ask_recursive(Q.real(x**x), Q.zero(x)) is True
    assert _ask_recursive(Q.real(Pow(x, -1, evaluate=False)), Q.zero(x)) is False
    two_thirds = Rational(2, 3)
    assert _ask_recursive(Q.real(x**two_thirds), Q.zero(x)) is True
    assert _ask_recursive(Q.real(x**pi), Q.zero(x)) is True
    assert _ask_recursive(Q.real(x**sqrt(2)), Q.zero(x)) is True
    assert _ask_recursive(Q.real(x**Rational(1/2)), Q.zero(x)) is True


@_both_exp_pow
def test_real_functions():
    # trigonometric functions
    assert _ask_recursive(Q.real(sin(x))) is None
    assert _ask_recursive(Q.real(cos(x))) is None
    assert _ask_recursive(Q.real(sin(x)), Q.real(x)) is True
    assert _ask_recursive(Q.real(cos(x)), Q.real(x)) is True

    # exponential function
    assert _ask_recursive(Q.real(exp(x))) is None
    assert _ask_recursive(Q.real(exp(x)), Q.real(x)) is True
    assert _ask_recursive(Q.real(x + exp(x)), Q.real(x)) is True
    assert _ask_recursive(Q.real(exp(2*pi*I, evaluate=False))) is True
    assert _ask_recursive(Q.real(exp(pi*I, evaluate=False))) is True
    assert _ask_recursive(Q.real(exp(pi*I/2, evaluate=False))) is False

    # logarithm
    assert _ask_recursive(Q.real(log(I))) is False
    assert _ask_recursive(Q.real(log(2*I))) is False
    assert _ask_recursive(Q.real(log(I + 1))) is False
    assert _ask_recursive(Q.real(log(x)), Q.complex(x)) is None
    assert _ask_recursive(Q.real(log(x)), Q.imaginary(x)) is False
    assert _ask_recursive(Q.real(log(exp(x))), Q.imaginary(x)) is None  # exp(2*pi*I) is 1, log(exp(pi*I)) is pi*I (disregarding periodicity)
    assert _ask_recursive(Q.real(log(exp(x))), Q.complex(x)) is None
    eq = Pow(exp(2*pi*I*x, evaluate=False), x, evaluate=False)
    assert _ask_recursive(Q.real(eq), Q.integer(x)) is True
    assert _ask_recursive(Q.real(exp(x)**x), Q.imaginary(x)) is True
    assert _ask_recursive(Q.real(exp(x)**x), Q.complex(x)) is None

    # Q.complexes
    assert _ask_recursive(Q.real(re(x))) is True
    assert _ask_recursive(Q.real(im(x))) is True


def test_matrix():

    # hermitian
    assert _ask_recursive(Q.hermitian(Matrix([[2, 2 + I, 4], [2 - I, 3, I], [4, -I, 1]]))) == True
    assert _ask_recursive(Q.hermitian(Matrix([[2, 2 + I, 4], [2 + I, 3, I], [4, -I, 1]]))) == False
    z = symbols('z', complex=True)
    assert _ask_recursive(Q.hermitian(Matrix([[2, 2 + I, z], [2 - I, 3, I], [4, -I, 1]]))) is None
    assert _ask_recursive(Q.hermitian(SparseMatrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11))))) == True
    assert _ask_recursive(Q.hermitian(SparseMatrix(((25, 15, -5), (15, I, 0), (-5, 0, 11))))) == False
    assert _ask_recursive(Q.hermitian(SparseMatrix(((25, 15, -5), (15, z, 0), (-5, 0, 11))))) is None

    # antihermitian
    A = Matrix([[0, -2 - I, 0], [2 - I, 0, -I], [0, -I, 0]])
    B = Matrix([[-I, 2 + I, 0], [-2 + I, 0, 2 + I], [0, -2 + I, -I]])
    assert _ask_recursive(Q.antihermitian(A)) is True
    assert _ask_recursive(Q.antihermitian(B)) is True
    assert _ask_recursive(Q.antihermitian(A**2)) is False
    C = (B**3)
    C.simplify()
    assert _ask_recursive(Q.antihermitian(C)) is True
    _A = Matrix([[0, -2 - I, 0], [z, 0, -I], [0, -I, 0]])
    assert _ask_recursive(Q.antihermitian(_A)) is None


@_both_exp_pow
def test_algebraic():
    assert _ask_recursive(Q.algebraic(x)) is None

    assert _ask_recursive(Q.algebraic(I)) is True
    assert _ask_recursive(Q.algebraic(2*I)) is True
    assert _ask_recursive(Q.algebraic(I/3)) is True

    assert _ask_recursive(Q.algebraic(sqrt(7))) is True
    assert _ask_recursive(Q.algebraic(2*sqrt(7))) is True
    assert _ask_recursive(Q.algebraic(sqrt(7)/3)) is True

    assert _ask_recursive(Q.algebraic(I*sqrt(3))) is True
    assert _ask_recursive(Q.algebraic(sqrt(1 + I*sqrt(3)))) is True

    assert _ask_recursive(Q.algebraic(1 + I*sqrt(3)**Rational(17, 31))) is True
    assert _ask_recursive(Q.algebraic(1 + I*sqrt(3)**(17/pi))) is None

    assert _ask_recursive(Q.algebraic(log(3, 2) ** 2)) is None

    for f in [exp, sin, tan, asin, atan, cos]:
        assert _ask_recursive(Q.algebraic(f(7))) is False
        assert _ask_recursive(Q.algebraic(f(7, evaluate=False))) is False
        assert _ask_recursive(Q.algebraic(f(0, evaluate=False))) is True
        assert _ask_recursive(Q.algebraic(f(x)), Q.algebraic(x)) is None
        assert _ask_recursive(Q.algebraic(f(x)), Q.algebraic(x) & Q.nonzero(x)) is False

    for g in [log, acos]:
        assert _ask_recursive(Q.algebraic(g(7))) is False
        assert _ask_recursive(Q.algebraic(g(7, evaluate=False))) is False
        assert _ask_recursive(Q.algebraic(g(1, evaluate=False))) is True
        assert _ask_recursive(Q.algebraic(g(x)), Q.algebraic(x)) is None
        assert _ask_recursive(Q.algebraic(g(x)), Q.algebraic(x) & Q.nonzero(x - 1)) is False

    for h in [cot, acot]:
        assert _ask_recursive(Q.algebraic(h(7))) is False
        assert _ask_recursive(Q.algebraic(h(7, evaluate=False))) is False
        assert _ask_recursive(Q.algebraic(h(x)), Q.algebraic(x)) is False

    assert _ask_recursive(Q.algebraic(sqrt(sin(7)))) is None
    assert _ask_recursive(Q.algebraic(sqrt(y + I*sqrt(7)))) is None

    assert _ask_recursive(Q.algebraic(2.47)) is None

    assert _ask_recursive(Q.algebraic(x), Q.transcendental(x)) is False
    assert _ask_recursive(Q.transcendental(x), Q.algebraic(x)) is False

    #https://github.com/sympy/sympy/issues/27445
    assert _ask_recursive(Q.algebraic(Pow(1, x, evaluate=False)), Q.algebraic(x)) is None
    assert _ask_recursive(Q.algebraic(Pow(x, y))) is None
    assert _ask_recursive(Q.algebraic(Pow(1, x, evaluate=False))) is None
    assert _ask_recursive(Q.algebraic(x**(pi*I))) is None
    assert _ask_recursive(Q.algebraic(pi**n),Q.integer(n) & Q.positive(n)) is False
    assert _ask_recursive(Q.algebraic(x**y),Q.algebraic(x) & Q.rational(y)) is True



def test_global():
    """Test ask with global assumptions"""
    assert ask(Q.integer(x)) is None
    global_assumptions.add(Q.integer(x))
    assert ask(Q.integer(x)) is True
    global_assumptions.clear()
    assert ask(Q.integer(x)) is None


def test_custom_context():
    """Test ask with custom assumptions context"""
    assert ask(Q.integer(x)) is None
    local_context = AssumptionsContext()
    local_context.add(Q.integer(x))
    assert ask(Q.integer(x), context=local_context) is True
    assert ask(Q.integer(x)) is None


def test_functions_in_assumptions():
    assert ask(Q.negative(x), Q.real(x) >> Q.positive(x)) is False
    assert ask(Q.negative(x), Equivalent(Q.real(x), Q.positive(x))) is False
    assert ask(Q.negative(x), Xor(Q.real(x), Q.negative(x))) is False


def test_composite_ask():
    assert ask(Q.negative(x) & Q.integer(x),
        assumptions=Q.real(x) >> Q.positive(x)) is False


def test_composite_proposition():
    assert _ask_recursive(True) is True
    assert _ask_recursive(False) is False
    assert _ask_recursive(~Q.negative(x), Q.positive(x)) is True
    assert _ask_recursive(~Q.real(x), Q.commutative(x)) is None
    assert _ask_recursive(Q.negative(x) & Q.integer(x), Q.positive(x)) is False
    assert _ask_recursive(Q.negative(x) & Q.integer(x)) is None
    assert _ask_recursive(Q.real(x) | Q.integer(x), Q.positive(x)) is True
    assert _ask_recursive(Q.real(x) | Q.integer(x)) is None
    assert _ask_recursive(Q.real(x) >> Q.positive(x), Q.negative(x)) is False
    assert _ask_recursive(Implies(
        Q.real(x), Q.positive(x), evaluate=False), Q.negative(x)) is False
    assert _ask_recursive(Implies(Q.real(x), Q.positive(x), evaluate=False)) is None
    assert _ask_recursive(Equivalent(Q.integer(x), Q.even(x)), Q.even(x)) is True
    assert _ask_recursive(Equivalent(Q.integer(x), Q.even(x))) is None
    assert _ask_recursive(Equivalent(Q.positive(x), Q.integer(x)), Q.integer(x)) is None
    assert ask(Q.real(x) | Q.integer(x), Q.real(x) | Q.integer(x)) is True

def test_tautology():
    assert ask(Q.real(x) | ~Q.real(x)) is True
    assert ask(Q.real(x) & ~Q.real(x)) is False

def test_composite_assumptions():
    assert _ask_recursive(Q.real(x), Q.real(x) & Q.real(y)) is True
    assert _ask_recursive(Q.positive(x), Q.positive(x) | Q.positive(y)) is None
    assert _ask_recursive(Q.positive(x), Q.real(x) >> Q.positive(y)) is None
    assert _ask_recursive(Q.real(x), ~(Q.real(x) >> Q.real(y))) is True

def test_key_extensibility():
    """test that you can add keys to the ask system at runtime"""
    # make sure the key is not defined
    raises(AttributeError, lambda: _ask_recursive(Q.my_key(x)))

    # New handler system
    class MyPredicate(Predicate):
        pass
    try:
        Q.my_key = MyPredicate()
        @Q.my_key.register(Symbol)
        def _(expr, assumptions):
            return True
        assert _ask_recursive(Q.my_key(x)) is True
        assert _ask_recursive(Q.my_key(x+1)) is None
    finally:
        del Q.my_key
    raises(AttributeError, lambda: _ask_recursive(Q.my_key(x)))


def test_type_extensibility():
    """test that new types can be added to the ask system at runtime
    """
    from sympy.core import Basic

    class MyType(Basic):
        pass

    @Q.prime.register(MyType)
    def _(expr, assumptions):
        return True

    assert _ask_recursive(Q.prime(MyType())) is True


def test_single_fact_lookup():
    known_facts = And(Implies(Q.integer, Q.rational),
                      Implies(Q.rational, Q.real),
                      Implies(Q.real, Q.complex))
    known_facts_keys = {Q.integer, Q.rational, Q.real, Q.complex}

    known_facts_cnf = to_cnf(known_facts)
    mapping = single_fact_lookup(known_facts_keys, known_facts_cnf)

    assert mapping[Q.rational] == {Q.real, Q.rational, Q.complex}


def test_generate_known_facts_dict():
    known_facts = And(Implies(Q.integer(x), Q.rational(x)),
                      Implies(Q.rational(x), Q.real(x)),
                      Implies(Q.real(x), Q.complex(x)))
    known_facts_keys = {Q.integer(x), Q.rational(x), Q.real(x), Q.complex(x)}

    assert generate_known_facts_dict(known_facts_keys, known_facts) == \
        {Q.complex: ({Q.complex}, set()),
         Q.integer: ({Q.complex, Q.integer, Q.rational, Q.real}, set()),
         Q.rational: ({Q.complex, Q.rational, Q.real}, set()),
         Q.real: ({Q.complex, Q.real}, set())}


@slow
def test_known_facts_consistent():
    """"Test that ask_generated.py is up-to-date"""
    x = Symbol('x')
    fact = get_known_facts(x)
    # test cnf clauses of fact between unary predicates
    cnf = CNF.to_CNF(fact)
    clauses = set()
    clauses.update(frozenset(Literal(lit.arg.function, lit.is_Not) for lit in sorted(cl, key=str)) for cl in cnf.clauses)
    assert get_all_known_facts() == clauses
    # test dictionary of fact between unary predicates
    keys = [pred(x) for pred in get_known_facts_keys()]
    mapping = generate_known_facts_dict(keys, fact)
    assert get_known_facts_dict() == mapping


def test_Add_queries():
    assert _ask_recursive(Q.prime(12345678901234567890 + (cos(1)**2 + sin(1)**2))) is True
    assert _ask_recursive(Q.even(Add(S(2), S(2), evaluate=False))) is True
    assert _ask_recursive(Q.prime(Add(S(2), S(2), evaluate=False))) is False
    assert _ask_recursive(Q.integer(Add(S(2), S(2), evaluate=False))) is True


def test_positive_assuming():
    with assuming(Q.positive(x + 1)):
        assert not ask(Q.positive(x))


def test_issue_5421():
    raises(TypeError, lambda: ask(pi/log(x), Q.real))


def test_issue_3906():
    raises(TypeError, lambda: ask(Q.positive))


def test_issue_5833():
    assert _ask_recursive(Q.positive(log(x)**2), Q.positive(x)) is None
    assert _ask_recursive(~Q.negative(log(x)**2), Q.positive(x)) is True


def test_issue_6732():
    raises(ValueError, lambda: ask(Q.positive(x), Q.positive(x) & Q.negative(x)))
    raises(ValueError, lambda: ask(Q.negative(x), Q.positive(x) & Q.negative(x)))


def test_issue_7246():
    assert _ask_recursive(Q.positive(atan(p)), Q.positive(p)) is True
    assert _ask_recursive(Q.positive(atan(p)), Q.negative(p)) is False
    assert _ask_recursive(Q.positive(atan(p)), Q.zero(p)) is False
    assert _ask_recursive(Q.positive(atan(x))) is None

    assert _ask_recursive(Q.positive(asin(p)), Q.positive(p)) is None
    assert _ask_recursive(Q.positive(asin(p)), Q.zero(p)) is None
    assert _ask_recursive(Q.positive(asin(Rational(1, 7)))) is True
    assert _ask_recursive(Q.positive(asin(x)), Q.positive(x) & Q.nonpositive(x - 1)) is True
    assert _ask_recursive(Q.positive(asin(x)), Q.negative(x) & Q.nonnegative(x + 1)) is False

    assert _ask_recursive(Q.positive(acos(p)), Q.positive(p)) is None
    assert _ask_recursive(Q.positive(acos(Rational(1, 7)))) is True
    assert _ask_recursive(Q.positive(acos(x)), Q.nonnegative(x + 1) & Q.nonpositive(x - 1)) is True
    assert _ask_recursive(Q.positive(acos(x)), Q.nonnegative(x - 1)) is None

    assert _ask_recursive(Q.positive(acot(x)), Q.positive(x)) is True
    assert _ask_recursive(Q.positive(acot(x)), Q.real(x)) is True
    assert _ask_recursive(Q.positive(acot(x)), Q.imaginary(x)) is False
    assert _ask_recursive(Q.positive(acot(x))) is None


@XFAIL
def test_issue_7246_failing():
    #Move this test to test_issue_7246 once
    #the new assumptions module is improved.
    assert _ask_recursive(Q.positive(acos(x)), Q.zero(x)) is True


def test_check_old_assumption():
    x = symbols('x', real=True)
    assert _ask_recursive(Q.real(x)) is True
    assert _ask_recursive(Q.imaginary(x)) is False
    assert _ask_recursive(Q.complex(x)) is True

    x = symbols('x', imaginary=True)
    assert _ask_recursive(Q.real(x)) is False
    assert _ask_recursive(Q.imaginary(x)) is True
    assert _ask_recursive(Q.complex(x)) is True

    x = symbols('x', complex=True)
    assert _ask_recursive(Q.real(x)) is None
    assert _ask_recursive(Q.complex(x)) is True

    x = symbols('x', positive=True)
    assert _ask_recursive(Q.positive(x)) is True
    assert _ask_recursive(Q.negative(x)) is False
    assert _ask_recursive(Q.real(x)) is True

    x = symbols('x', commutative=False)
    assert _ask_recursive(Q.commutative(x)) is False

    x = symbols('x', negative=True)
    assert _ask_recursive(Q.positive(x)) is False
    assert _ask_recursive(Q.negative(x)) is True

    x = symbols('x', nonnegative=True)
    assert _ask_recursive(Q.negative(x)) is False
    assert _ask_recursive(Q.positive(x)) is None
    assert _ask_recursive(Q.zero(x)) is None

    x = symbols('x', finite=True)
    assert _ask_recursive(Q.finite(x)) is True

    x = symbols('x', prime=True)
    assert _ask_recursive(Q.prime(x)) is True
    assert _ask_recursive(Q.composite(x)) is False

    x = symbols('x', composite=True)
    assert _ask_recursive(Q.prime(x)) is False
    assert _ask_recursive(Q.composite(x)) is True

    x = symbols('x', even=True)
    assert _ask_recursive(Q.even(x)) is True
    assert _ask_recursive(Q.odd(x)) is False

    x = symbols('x', odd=True)
    assert _ask_recursive(Q.even(x)) is False
    assert _ask_recursive(Q.odd(x)) is True

    x = symbols('x', nonzero=True)
    assert _ask_recursive(Q.nonzero(x)) is True
    assert _ask_recursive(Q.zero(x)) is False

    x = symbols('x', zero=True)
    assert _ask_recursive(Q.zero(x)) is True

    x = symbols('x', integer=True)
    assert _ask_recursive(Q.integer(x)) is True

    x = symbols('x', rational=True)
    assert _ask_recursive(Q.rational(x)) is True
    assert _ask_recursive(Q.irrational(x)) is False

    x = symbols('x', irrational=True)
    assert _ask_recursive(Q.irrational(x)) is True
    assert _ask_recursive(Q.rational(x)) is False

    x = symbols('x', transcendental=True)
    assert _ask_recursive(Q.transcendental(x)) is True


def test_issue_9636():
    assert _ask_recursive(Q.integer(1.0)) is None
    assert _ask_recursive(Q.prime(3.0)) is None
    assert _ask_recursive(Q.composite(4.0)) is None
    assert _ask_recursive(Q.even(2.0)) is None
    assert _ask_recursive(Q.odd(3.0)) is None


def test_autosimp_used_to_fail():
    # See issue #9807
    assert _ask_recursive(Q.imaginary(0**I)) is None
    assert _ask_recursive(Q.imaginary(0**(-I))) is None
    assert _ask_recursive(Q.real(0**I)) is None
    assert _ask_recursive(Q.real(0**(-I))) is None


def test_custom_AskHandler():
    from sympy.logic.boolalg import conjuncts

    class MersennePredicate(Predicate):
        pass
    try:
        Q.mersenne = MersennePredicate()
        @Q.mersenne.register(Integer)
        def _(expr, assumptions):
            if _ask_recursive(Q.integer(log(expr + 1, 2))):
                return True
        @Q.mersenne.register(Symbol)
        def _(expr, assumptions):
            if expr in conjuncts(assumptions):
                return True
        assert _ask_recursive(Q.mersenne(7)) is True
        assert _ask_recursive(Q.mersenne(n), Q.mersenne(n)) is None
    finally:
        del Q.mersenne


def test_polyadic_predicate():

    class SexyPredicate(Predicate):
        pass
    try:
        Q.sexyprime = SexyPredicate()

        @Q.sexyprime.register(Integer, Integer)
        def _(int1, int2, assumptions):
            args = sorted([int1, int2])
            if not all(_ask_recursive(Q.prime(a), assumptions) for a in args):
                return False
            return args[1] - args[0] == 6

        @Q.sexyprime.register(Integer, Integer, Integer)
        def _(int1, int2, int3, assumptions):
            args = sorted([int1, int2, int3])
            if not all(_ask_recursive(Q.prime(a), assumptions) for a in args):
                return False
            return args[2] - args[1] == 6 and args[1] - args[0] == 6

        assert _ask_recursive(Q.sexyprime(5, 11)) is True
        assert _ask_recursive(Q.sexyprime(7, 13, 19)) is True
    finally:
        del Q.sexyprime


def test_Predicate_handler_is_unique():

    # Undefined predicate does not have a handler
    assert Predicate('mypredicate').handler is None

    # Handler of defined predicate is unique to the class
    class MyPredicate(Predicate):
        pass
    mp1 = MyPredicate(Str('mp1'))
    mp2 = MyPredicate(Str('mp2'))
    assert mp1.handler is mp2.handler


def test_relational():
    assert _ask_recursive(Q.eq(x, 0), Q.zero(x)) is True
    assert _ask_recursive(Q.eq(x, 0), Q.nonzero(x)) is False
    assert _ask_recursive(Q.ne(x, 0), Q.zero(x)) is False
    assert _ask_recursive(Q.ne(x, 0), Q.nonzero(x)) is True


def test_issue_25221():
    assert _ask_recursive(Q.transcendental(x), Q.algebraic(x) | Q.positive(y,y)) is None
    assert _ask_recursive(Q.transcendental(x), Q.algebraic(x) | (0 > y)) is None
    assert _ask_recursive(Q.transcendental(x), Q.algebraic(x) | Q.gt(0,y)) is None


def test_issue_27440():
    nan = S.NaN
    assert _ask_recursive(Q.negative(nan)) is None

def test_issue_28127():
    assert _ask_recursive(Q.le(x,y), Q.gt(x,y)) is False
    assert _ask_recursive(Q.ge(x,y), Q.lt(x,y)) is False
    assert _ask_recursive(Q.gt(y,x), Q.lt(x,y)) is True
    assert _ask_recursive(Q.lt(y,x), Q.gt(x,y)) is True
    assert _ask_recursive(Q.le(y,x), Q.ge(x,y)) is True


# https://github.com/sympy/sympy/pull/30175
@XFAIL
def test_queries_require_satask():
    """
    Queries that `ask` only answers through `satask`. `_ask_recursive`
    returns None for every assertion below.
    """
    assert _ask_recursive(Q.integer(x), Q.even(x) | Q.odd(x)) is True
    assert _ask_recursive(Q.nonzero(x), Q.negative(x) | Q.positive(x)) is True
    assert _ask_recursive(Q.zero(x), Q.negative(x) | Q.positive(x)) is False
    assert _ask_recursive(Q.zero(x), Q.nonnegative(x) & Q.nonpositive(x)) is True
    assert _ask_recursive(Q.zero(x) | Q.zero(y), Q.zero(x*y)) is True
    assert _ask_recursive(Q.negative(x), Q.real(x) >> Q.positive(x)) is False
    assert _ask_recursive(Q.negative(x), Equivalent(Q.real(x), Q.positive(x))) is False
    assert _ask_recursive(Q.negative(x), Xor(Q.real(x), Q.negative(x))) is False
    assert _ask_recursive(Q.negative(x) & Q.integer(x),
        assumptions=Q.real(x) >> Q.positive(x)) is False
    assert _ask_recursive(Q.real(x) | Q.integer(x), Q.real(x) | Q.integer(x)) is True
    assert _ask_recursive(Q.real(x) | ~Q.real(x)) is True
    assert _ask_recursive(Q.real(x) & ~Q.real(x)) is False

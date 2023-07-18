from sympy.parsing.tests.test_latex import GOOD_PAIRS
from sympy.testing.pytest import raises, XFAIL
from sympy.external import import_module

from sympy.concrete.products import Product
from sympy.concrete.summations import Sum
from sympy.core.add import Add
from sympy.core.function import (Derivative, Function)
from sympy.core.mul import Mul
from sympy.core.numbers import (E, oo)
from sympy.core.power import Pow
from sympy.core.relational import (GreaterThan, LessThan, StrictGreaterThan, StrictLessThan, Unequality)
from sympy.core.symbol import Symbol
from sympy.functions.combinatorial.factorials import (binomial, factorial)
from sympy.functions.elementary.complexes import (Abs, conjugate)
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.integers import (ceiling, floor)
from sympy.functions.elementary.miscellaneous import (root, sqrt)
from sympy.functions.elementary.trigonometric import (asin, cos, csc, sec, sin, tan)
from sympy.integrals.integrals import Integral
from sympy.series.limits import Limit

from sympy.core.relational import Eq, Ne, Lt, Le, Gt, Ge
from sympy.physics.quantum.state import Bra, Ket
from sympy.abc import x, y, z, a, b, c, t, k, n
lark = import_module("lark")

if not lark:
    disabled = True

theta = Symbol('theta')
f = Function('f')


# shorthand definitions
def _Add(a, b):
    return Add(a, b, evaluate=False)


def _Mul(a, b):
    return Mul(a, b, evaluate=False)


def _Pow(a, b):
    return Pow(a, b, evaluate=False)


def _Sqrt(a):
    return sqrt(a, evaluate=False)


def _Conjugate(a):
    return conjugate(a, evaluate=False)


def _Abs(a):
    return Abs(a, evaluate=False)


def _factorial(a):
    return factorial(a, evaluate=False)


def _exp(a):
    return exp(a, evaluate=False)


def _log(a, b):
    return log(a, b, evaluate=False)


def _binomial(n, k):
    return binomial(n, k, evaluate=False)


# Temporary testing code to allow for quick prototyping
def determine_parseable_lark():
    from sympy.parsing.latex.lark import parse_latex_lark

    failure_list = []
    for i, (latex_str, sympy_expr) in enumerate(GOOD_PAIRS):
        try:
            parse_latex_lark(latex_str)
        except Exception as e:
            failure_list.append(i)

    return failure_list

# A list of tests we usually use to test the edge cases and stress test the WIP implementation
TRICKY_TESTS = [
    (r"\mathit{HELLO world}", Symbol('HELLO world')),
    (r"-3.14", -3.14),
    (r"(-7.13)(1.5)", _Mul(-7.13, 1.5)),
    (r"3x - 1", _Add(_Mul(3, x), -1)),
    (r"-c", -c),
    (r"a + b - a", _Add(a + b, -a)),
    (r"(x + y) z", _Mul(_Add(x, y), z)),
    (r"\frac12", _Pow(2, -1)),
    (r"\frac12y", _Mul(_Pow(2, -1), y)),
    (r"\frac1234", _Mul(_Pow(2, -1), 34)),
    (r"\frac2{3}", _Mul(2, _Pow(3, -1))),
    (r"\sin a \cos b", _Mul(sin(a), cos(b))),
    (r"100!", _factorial(100)),
    (r"a \negthinspace b", _Mul(a, b))
]


def tricky_tests():
    from sympy.parsing.latex.lark import parse_latex_lark

    failure_list = []
    for i, (latex_str, sympy_expr) in enumerate(TRICKY_TESTS):
        try:
            parse_latex_lark(latex_str)
        except Exception as e:
            failure_list.append(i)

    return failure_list

if __name__ == "__main__":
    fail_list = determine_parseable_lark()
    print(fail_list)
    print("Failures =", len(fail_list))
    print("Passes =", len(GOOD_PAIRS) - len(fail_list))

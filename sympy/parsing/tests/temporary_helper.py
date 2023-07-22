from sympy.parsing.latex.lark.latex_parser import pretty_print_lark_trees
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

disabled = lark is None

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


# # Temporary testing code to allow for quick prototyping
# def determine_parseable_lark():
#     from sympy.parsing.latex.lark import parse_latex_lark

#     failure_list = []
#     for i, (latex_str, sympy_expr) in enumerate(GOOD_PAIRS):
#         try:
#             parse_latex_lark(latex_str)
#         except Exception as e:
#             failure_list.append(i)

#     return failure_list

def determine_parseable_lark():
    from sympy.parsing.latex.lark import parse_latex_lark

    # If the expression doesn't raise an error, it means that the expression is parsed into _something_, so it's not a
    # complete failure for that test case. If even that doesn't happen, we add it into `complete_failure_list`
    complete_failure_list = []
    # If the expression parses into _something_, but the output doesn't equal the string given in the test case, then
    # it's a partial failure. If that happens, we add it here.
    partial_failure_list = []
    # If the expression is correctly handled, we add it to this list.
    success_list = []
    # checks outputted sympy expression for equality with the test case.
    for i, (latex_str, sympy_expr) in enumerate(GOOD_PAIRS):
        try:
            result = parse_latex_lark(latex_str)

            if result != sympy_expr:
                partial_failure_list.append(i)
                parse_latex_lark(latex_str, print_debug_output=True)  # for debugging purposes
            else:
                success_list.append(i)
        except Exception as e:
            complete_failure_list.append(i)

    return success_list, partial_failure_list, complete_failure_list

# A list of tests we usually use to test the edge cases and stress test the WIP implementation
TRICKY_TESTS = [
    (r"\mathit{HELLOworld}", Symbol('HELLOworld')),
    # (r"\mathit{HELLO world}", Symbol('HELLO world')),
    (r"-3.14", -3.14),
    (r"0", 0),
    (r".0", .0),
    (r"2.0", 2.0),
    (r"14", 14),
    (r"(-7.13)(1.5)", _Mul(-7.13, 1.5)),
    (r"3x - 1", _Add(_Mul(3, x), -1)),
    (r"xy", _Mul(x, y)),
    (r"-c", -c),
    (r"a + b - a", _Add(a + b, -a)),
    (r"a - b - c", None),
    (r"a(b)", None),
    (r"(x + y) z", _Mul(_Add(x, y), z)),
    (r"\frac12", _Pow(2, -1)),
    (r"\frac12y", _Mul(_Pow(2, -1), y)),
    (r"\frac1234", _Mul(_Pow(2, -1), 34)),
    (r"\frac2{3}", _Mul(2, _Pow(3, -1))),
    (r"\sin a \cos b", _Mul(sin(a), cos(b))),
    (r"100!", _factorial(100)),
    (r"a \negthinspace b", _Mul(a, b)),
]


def test_tricky_tests():
    from sympy.parsing.latex.lark import parse_latex_lark

    failure_list = []
    for i, (latex_str, sympy_expr) in enumerate(TRICKY_TESTS):
        try:
            out = parse_latex_lark(latex_str)
            pretty = pretty_print_lark_trees(out)
            print(pretty)
            assert "_ambig" not in pretty
        except Exception as e:
            failure_list.append(i)

    return failure_list

# if __name__ == "__main__":
#     fail_list = determine_parseable_lark()
#     print(fail_list)
#     print("Failures =", len(fail_list))
#     print("Passes =", len(GOOD_PAIRS) - len(fail_list))

if __name__ == "__main__":
    successes, partial_failures, complete_failures = determine_parseable_lark()

    print("List of failures =", complete_failures)
    print("Not even parsed =", len(complete_failures))
    print("Parsing but not (fully) transformed =", len(partial_failures))
    print("Passes =", len(successes))

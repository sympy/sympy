from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.parsing.ast_parser import parse_expr
from sympy.testing.pytest import raises
from sympy.core.sympify import SympifyError
import warnings

def test_parse_expr():
    a, b = symbols('a, b')
    # tests issue_16393
    assert parse_expr('a + b', {}) == a + b
    raises(SympifyError, lambda: parse_expr('a + ', {}))

    # tests Transform.visit_Num
    assert parse_expr('1 + 2', {}) == S(3)
    assert parse_expr('1 + 2.0', {}) == S(3.0)

    # tests Transform.visit_Name
    assert parse_expr('Rational(1, 2)', {}) == S(1)/2
    assert parse_expr('a', {'a': a}) == a

    # tests issue_23092
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert parse_expr('6 * 7', {}) == S(42)
        # Note: The test below can be removed when support for Python 3.7 is
        # dropped. This test exists to ensure that the visit_Num function
        # exists for Python 3.7, because in 3.7, Python didn't use the
        # visit_Constant to create AST Nodes yet.
        test_expr = parse_expr('1 / 3', {})
        assert test_expr == S(1)/3 # sanity check
        assert isinstance(test_expr, Rational)

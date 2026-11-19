from __future__ import annotations
from sympy.printing.codeprinter import CodePrinter, PrintMethodNotImplementedError
from sympy.core import symbols
from sympy.core.symbol import Dummy
from sympy.core.relational import Eq, Lt
from sympy.logic.boolalg import And
from sympy.printing import (ccode, fcode, glsl_code, jscode, julia_code,
    mathematica_code, octave_code, pycode, rcode, rust_code)
from sympy.testing.pytest import raises
from sympy import cos
from sympy.utilities.lambdify import lambdify
from math import cos as math_cos
from sympy.printing.lambdarepr import LambdaPrinter


def setup_test_printer(**kwargs):
    p = CodePrinter(settings=kwargs)
    p._not_supported = set()
    p._number_symbols = set()
    return p


def test_print_Dummy():
    d = Dummy('d')
    p = setup_test_printer()
    assert p._print_Dummy(d) == "d_%i" % d.dummy_index

def test_print_Symbol():

    x, y = symbols('x, if')

    p = setup_test_printer()
    assert p._print(x) == 'x'
    assert p._print(y) == 'if'

    p.reserved_words.update(['if'])
    assert p._print(y) == 'if_'

    p = setup_test_printer(error_on_reserved=True)
    p.reserved_words.update(['if'])
    with raises(ValueError):
        p._print(y)

    p = setup_test_printer(reserved_word_suffix='_He_Man')
    p.reserved_words.update(['if'])
    assert p._print(y) == 'if_He_Man'


def test_lambdify_LaTeX_symbols_issue_23374():
    # Create symbols with Latex style names
    x1, x2 = symbols("x_{1} x_2")

    # Lambdify the function
    f1 = lambdify([x1, x2], cos(x1 ** 2 + x2 ** 2))

    # Test that the function works correctly (numerically)
    assert f1(1, 2) == math_cos(1 ** 2 + 2 ** 2)

    # Explicitly generate a custom printer to verify the naming convention
    p = LambdaPrinter()
    expr_str = p.doprint(cos(x1 ** 2 + x2 ** 2))
    assert 'x_1' in expr_str
    assert 'x_2' in expr_str


def test_issue_15791():
    class CrashingCodePrinter(CodePrinter):
        def emptyPrinter(self, obj):
            raise NotImplementedError

    from sympy.matrices import (
        MutableSparseMatrix,
        ImmutableSparseMatrix,
    )

    c = CrashingCodePrinter()

    # these should not silently succeed
    with raises(PrintMethodNotImplementedError):
        c.doprint(ImmutableSparseMatrix(2, 2, {}))
    with raises(PrintMethodNotImplementedError):
        c.doprint(MutableSparseMatrix(2, 2, {}))


def test_issue_11366():
    x, y, z = symbols('x y z')
    nested = Eq(x, Eq(y, z, evaluate=False), evaluate=False)
    compared = Eq(x > 1, y > 1, evaluate=False)
    chained = Lt(x, Lt(y, z, evaluate=False), evaluate=False)
    boolean = Eq(And(x > 0, y > 0), z > 0, evaluate=False)
    for printer in (ccode, glsl_code, jscode, julia_code, mathematica_code,
                    octave_code, rcode, rust_code):
        assert printer(nested) == 'x == (y == z)'
        assert printer(chained) == 'x < (y < z)'
    for printer in (ccode, glsl_code, jscode, julia_code, mathematica_code):
        assert printer(compared) == '(x > 1) == (y > 1)'
        assert printer(boolean) == '(x > 0 && y > 0) == (z > 0)'
    assert octave_code(boolean) == '(x > 0 & y > 0) == (z > 0)'
    assert rcode(boolean) == '(x > 0 & y > 0) == (z > 0)'
    assert rust_code(compared) == '(x > 1.0) == (y > 1.0)'
    assert fcode(chained, source_format='free') == 'x < (y < z)'
    assert pycode(compared) == '((x > 1) == (y > 1))'
    assert pycode(boolean) == '(((x > 0) and (y > 0)) == (z > 0))'
    assert ccode(Eq(x + 1, y)) == 'x + 1 == y'

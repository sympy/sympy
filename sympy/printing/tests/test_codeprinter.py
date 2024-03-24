from sympy.printing.codeprinter import CodePrinter, PrintMethodNotImplementedError, InvalidVariableNameError
from sympy.core import symbols
from sympy.core.symbol import Symbol, Dummy
from sympy.testing.pytest import raises


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


    p2 = setup_test_printer(strict_names=True)
    for invalid_name in ["f'(x)", "a b", "a+b", "#error", "", "$_"]:
        invalid_symbol = Symbol(invalid_name)
        assert p._print(invalid_symbol) == invalid_symbol.name
        with raises(InvalidVariableNameError):
            p2._print(invalid_symbol)


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

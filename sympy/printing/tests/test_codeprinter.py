from sympy.printing.codeprinter import CodePrinter
from sympy.core import C

def setup_test_printer(*args, **kwargs):
    p = CodePrinter(*args, **kwargs)
    p._not_supported = set()
    p._number_symbols = set()
    return p

def test_print_Dummy():
    d = C.Dummy('d')
    p = setup_test_printer()
    assert p._print_Dummy(d) == "d_%i" % d.dummy_index

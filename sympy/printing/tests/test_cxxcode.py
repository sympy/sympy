from sympy import symbols
from sympy.functions import beta, Ei, zeta
from sympy.printing.cfunctions import log1p
from sympy.printing.cxxcode import CXX11CodePrinter, _CXX17CodePrinter

x, y = symbols('x y')

def test_CXX11CodePrinter():
    assert CXX11CodePrinter().doprint(log1p(x)) == 'std::log1p(x)'

def test_subclass_print_method():
    class MyPrinter(CXX11CodePrinter):
        def _print_log1p(self, expr):
            return 'my_library::log1p(%s)' % ', '.join(map(self._print, expr.args))

    assert MyPrinter().doprint(log1p(x)) == 'my_library::log1p(x)'


def test_CXX17CodePrinter():
    assert _CXX17CodePrinter().doprint(beta(x, y)) == 'std::beta(x, y)'
    assert _CXX17CodePrinter().doprint(Ei(x)) == 'std::expint(x)'
    assert _CXX17CodePrinter().doprint(zeta(x)) == 'std::riemann_zeta(x)'

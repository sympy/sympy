from sympy import symbols
from sympy.tensor import Indexed, Idx
from sympy.utilities.autowrap import autowrap
from sympy.utilities.pytest import XFAIL

#
# test runners used by several language-backend combinations
#

def runtest_autowrap_twice(language, backend):
    a, b, c = symbols('a b c')
    f = autowrap((((a + b)/c)**5).expand(), language, backend)
    g = autowrap((((a + b)/c)**4).expand(), language, backend)

    # check that autowrap updates the module name.  Else, g gives the same as f
    assert f(1, -2, 1) == -1.0
    assert g(1, -2, 1) ==  1.0

def runtest_autowrap_trace(language, backend):
    A = IndexedBase('A')
    n = symbols('n', integer=True)
    x = symbols('x')
    i = Idx('i', n)

    expr = Eq(x, A[i, i])
    tr = binary_function(expr, language, backend)
    trace = tr(A[i, i])
    assert str(trace) == "tr(A[i, i])"
    tr = binary_function(expr, language, backend)

def runtest_autowrap_arrays(language, backend):
    A, x, y = map(IndexedBase, ['A', 'x', 'y'])
    n, m = symbols('n m', integer=True)
    i = Idx('i', m)
    j = Idx('j', n)

    expr = Eq(y[i], A[i, j]*x[j])
    mv = autowrap(expr, language, backend)

#
# tests of language-backend combinations
#

def test_wrap_twice_c_cython():
    runtest_autowrap_twice('C', 'cython')

def test_wrap_twice_f95_f2py():
    runtest_autowrap_twice('f95', 'f2py')

@XFAIL
def test_autowrap_arrays_cython_C():
    runtest_autowrap_arrays('C', 'cython')

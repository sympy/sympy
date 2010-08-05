from sympy import symbols
from sympy.utilities.autowrap import autowrap

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

#
# tests of language-backend combinations
#

def test_wrap_twice_c_cython():
    runtest_autowrap_twice('C', 'cython')

def test_wrap_twice_f95_f2py():
    runtest_autowrap_twice('f95', 'f2py')

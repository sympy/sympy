from sympy import symbols
from sympy.utilities.autowrap import autowrap

def test_autowrap_f2py_f95():

    a, b, c = symbols('a b c')

    f = autowrap((((a + b)/c)**5).expand(), language='f95', backend='f2py')
    assert f(1, -2, 1) == -1.0
    assert f(1.0, -2.0, 1.0) == -1.0
    assert f(3, -4, 1) == -1.0
    assert f(-3, 4, 1) ==  1.0

    # check that autowrap updates the module name.  If not, the next line
    # will calculate the above expression again and assertion will fail:
    g = autowrap((((a + b)/ c)**4).expand(), language='f95', backend='f2py')
    assert g(1, -2, 1) ==  1.0



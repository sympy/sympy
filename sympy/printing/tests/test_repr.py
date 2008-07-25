from sympy.utilities.pytest import XFAIL
from sympy import Symbol, Function, WildFunction
from sympy.polynomials import Polynomial
from sympy.polys.polynomial import Poly
from sympy.printing import srepr

x = Symbol('x')
y = Symbol('y')

def test_function_repr():
    f = Function('f')
    fx= f(x)
    w = WildFunction('w')
    assert srepr(fx) == "Function('f')(Symbol('x'))"
    assert srepr(w)  == "WildFunction('w')"

@XFAIL
def test_unapplied_function_repr():
    f = Function('f')
    assert srepr(f) == "Function('f')"   # this does not work

def test_Polynomial():
    f = Polynomial(x+2)
    assert srepr(f) == "Polynomial(Add(Integer(2), Symbol('x')), ((One(1), One(1)), (Integer(2), Zero(0))), [Symbol('x')], 'grevlex')"

def test_Poly():
    assert srepr(Poly(7, x)) == \
        "Poly([(Integer(7), (0,))], Symbol('x'), order='grlex')"

    assert srepr(Poly(2*x*y + 7, x, y)) == \
        "Poly([(Integer(2), (1, 1)), (Integer(7), (0, 0))]," \
        " Symbol('x'), Symbol('y'), order='grlex')"

    assert srepr(Poly(2*x*y - 7, x, y, order='grevlex')) == \
        "Poly([(Integer(2), (1, 1)), (Integer(-7), (0, 0))]," \
        " Symbol('x'), Symbol('y'), order='grevlex')"
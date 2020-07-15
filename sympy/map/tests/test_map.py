from sympy import symbols, S
from sympy.core.symbol import Str
from sympy.map import (
    Map, UndefinedMap, InverseMap, IdentityMap, AppliedMap
)
from sympy.testing.pytest import raises

x, y = symbols('x y')

class F(Map):
    def eval(self, x):
        return x+1
f = F()

def test_Map():

    # test UndefinedMap generation
    undf = Map('f')
    assert isinstance(undf, UndefinedMap)
    assert undf.name == Str('f')

    # test nargs
    assert Map('f', domain=S.Reals).nargs == 1
    assert Map('f', domain=S.Reals**1).nargs == 1
    assert Map('f', domain=S.Reals*S.Integers, codomain=S.Reals).nargs == 2
    assert Map('f', domain=S.Reals**4).nargs == 4

def test_AppliedMap():

    arg = x.diff(x, evaluate=False) # unevaluated argument

    # do not evaluate by default
    assert f(arg) == AppliedMap(f, arg)

    # evaluate by evaluate=True keyword
    assert f(arg, evaluate=True) == AppliedMap(f, arg, evaluate=True) == arg+1

    # evaluation by doit(deep=False)
    assert AppliedMap(f, arg).doit(deep=False) == arg+1

    # recursive evaluation by doit()
    assert AppliedMap(f, arg).doit() == 2

    # undefined map cannot be evaluated
    assert Map('f')(x).doit() == AppliedMap(Map('f'), x)

def test_InverseMap():

    # cannot be evaluated when _eval_inverse is not implemented
    assert f.inv() == f.inv(evaluate=True) == InverseMap(f)

    class G(Map):
        def _eval_inverse(self):
            return H()
    g = G()
    class H(Map):
        def _eval_inverse(self):
            return G()
    h = H()

    # do not evaluate by default
    assert g.inv() == InverseMap(g)
    assert h.inv() == InverseMap(h)

    # evaluate by evaluate=True keyword
    assert g.inv(evaluate=True) == InverseMap(g, evaluate=True) == h
    assert h.inv(evaluate=True) == InverseMap(h, evaluate=True) == g

    # evaluation by doit()
    assert InverseMap(g).doit() == InverseMap(g).doit(deep=False) == h
    assert InverseMap(h).doit() == InverseMap(h).doit(deep=False) == g

def test_IdentityMap():

    I = IdentityMap()

    # inverse is self
    assert I.inv(evaluate=True) == I

    # application returns argument
    assert I(x, evaluate=True) == x
    assert I((x, y), evaluate=True) == (x, y)
    # does not take multiple arguments
    raises(TypeError, lambda: I(x,y))

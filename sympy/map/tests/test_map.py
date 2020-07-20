from sympy import symbols, S, assuming, Q, ask
from sympy.core.symbol import Str
from sympy.map import (
    Map, UndefinedMap, InverseMap, IdentityMap, AppliedMap,
    RestrictedMap
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

    # test arity
    assert Map('f', domain=S.Reals**0).arity == 0
    assert Map('f', domain=S.Reals).arity == 1
    assert Map('f', domain=S.Reals**1).arity == 1
    assert Map('f', domain=S.Reals*S.Integers, codomain=S.Reals).arity == 2
    assert Map('f', domain=S.Reals**4).arity == 4

    # invertibility cannot be determined by default
    assert ask(Q.invertible(f)) is None

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

    # cannot invert if mapping is assumed to be non-invertible
    with assuming(~Q.invertible(f)):
        raises(TypeError, lambda: f.inv())
    # can invert if mapping is assumed to be invertible
    with assuming(Q.invertible(f)):
        f.inv()
    # can invert if invertibility annot be determined as well
    f.inv()

    # inverted map is always assumed to be invertible back
    assert ask(Q.invertible(f.inv())) is True

    # domain and codomain
    m = Map('m', domain=S.Reals, codomain=S.Naturals0)
    assert m.inv().domain is S.Naturals0, m.inv().domain is S.Reals

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

def test_RestrictedMap():
    assert f.restrict(S.Naturals0)(2, evaluate=True) == 3
    assert f.restrict(S.Naturals0).restrict(S.Naturals, evaluate=True) == f.restrict(S.Naturals)
    assert f.restrict(S.Naturals0).is_restriction(f)

def test_restriction():
    f1, f2 = Map('f', S.Reals, S.Reals), Map('f', S.Naturals0, S.Naturals0)

    assert f2.is_restriction(f1)
    assert f2.inv().is_restriction(f1.inv())
    assert IdentityMap(S.Naturals0).is_restriction(IdentityMap(S.Reals))
    assert f2.restrict(S.Naturals).is_restriction(f2)

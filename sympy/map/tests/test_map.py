from sympy import symbols, S, assuming, Q, ask
from sympy.core.symbol import Str
from sympy.map import (
    Map, UndefinedMap, InverseMap, IdentityMap, AppliedMap,
    ConstantMap, RestrictedMap, isappliedmap
)
from sympy.testing.pytest import raises, warns_deprecated_sympy

x, y = symbols('x y')

class F(Map):
    def eval(self, x):
        return x+1
    def _eval_inverse(self):
        return G()
    def _eval_restrict(self, domain):
        return H()
class G(Map):
    invertible = True
    commutative = True
    associative = False
    def eval(self, x):
        return x-1
    def _eval_inverse(self):
        return F()
class H(Map):
    invertible = False
    commutative = False
    associative = True
    def _map_content(self):
        return F()._map_content()

f, g, h = F(), G(), H()

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
    # invertibility can be set
    assert ask(Q.invertible(g)) is True
    assert ask(Q.invertible(h)) is False
    # invertibility can be assumed
    assert ask(Q.invertible(f), Q.invertible(f)) is True
    assert ask(Q.invertible(g), ~Q.invertible(g)) is False
    assert ask(Q.invertible(h), Q.invertible(h)) is True

    # commutativity cannot be determined by default
    assert ask(Q.commutative(f)) is False
    # commutativity can be set
    assert ask(Q.commutative(g)) is True
    assert ask(Q.commutative(h)) is False
    # commutativity can be assumed
    assert ask(Q.commutative(f), Q.commutative(f)) is True
    assert ask(Q.commutative(g), ~Q.commutative(g)) is False
    assert ask(Q.commutative(h), Q.commutative(h)) is True

    # associativity cannot be determined by default
    assert ask(Q.associative(f)) is None
    # associativity can be set
    assert ask(Q.associative(g)) is False
    assert ask(Q.associative(h)) is True
    # associativity can be assumed
    assert ask(Q.associative(f), Q.associative(f)) is True
    assert ask(Q.associative(g), Q.associative(g)) is True
    assert ask(Q.associative(h), ~Q.associative(h)) is False

def test_AppliedMap():

    arg = x.diff(x, evaluate=False) # unevaluated argument
    args = (arg,)

    # do not evaluate by default
    assert f(arg) == AppliedMap(f, args)

    # evaluate by evaluate=True keyword
    assert f(arg, evaluate=True) == AppliedMap(f, args, evaluate=True) == arg+1

    # evaluation by doit(deep=False)
    assert AppliedMap(f, args).doit(deep=False) == arg+1

    # recursive evaluation by doit()
    assert AppliedMap(f, args).doit() == 2

    # undefined map cannot be evaluated
    assert Map('f')(x).doit() == AppliedMap(Map('f'), (x,))

    # map's unevaluated codomain contains the result of evaluation
    fx = Map('f', codomain=S.Integers)(x)
    assert fx in S.Integers
    assert S.Integers.contains(fx)
    assert fx in S.Reals
    assert S.Reals.contains(fx)

    # error if evaluated result is not in codomain of mapping
    class Phi(Map):
        codomain = S.Naturals
        def eval(self, x):
            return x/2
    phi = Phi()
    raises(TypeError, lambda: phi(S.One, evaluate=True))

    # commutative map sorts the argument
    with assuming(Q.commutative(h)):
        assert h(x, y, evaluate=True) == h(y, x, evaluate=True)

    # differentiation
    a, b, c = Map('a', domain=S.Reals**2), Map('b', codomain=S.Reals), Map('c', codomain=S.Reals)
    assert a(x,y).diff(x) == a.diff(1)(x,y)
    assert a(y,x).diff(x) == a.diff(2)(y,x)
    assert b(x).diff(x) == b.diff(1)(x)
    assert b(c(x)).diff(x) == b.diff(1)(c(x)) * c.diff(1)(x)
    assert a(b(x), c(x)).diff(x) == \
        a.diff(1)(b(x), c(x))*b.diff(1)(x) + a.diff(2)(b(x), c(x))*c.diff(1)(x)

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
    assert m.inv().domain is S.Naturals0 and m.inv().codomain is S.Reals

    # cannot be evaluated when _eval_inverse is not implemented
    assert m.inv() == m.inv(evaluate=True) == InverseMap(m)

    # do not evaluate by default
    assert f.inv() == InverseMap(f)
    assert g.inv() == InverseMap(g)

    # evaluate by evaluate=True keyword
    assert f.inv(evaluate=True) == InverseMap(f, evaluate=True) == g
    assert g.inv(evaluate=True) == InverseMap(g, evaluate=True) == f

    # evaluation by doit()
    assert InverseMap(f).doit() == InverseMap(f).doit(deep=False) == g
    assert InverseMap(g).doit() == InverseMap(g).doit(deep=False) == f

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
    # restriction can be evaluated by _eval_restrict method
    assert f.restrict(S.Naturals, evaluate=True) == h

    # else, cannot be evaluated
    assert g.restrict(S.Naturals, evaluate=True) == RestrictedMap(g, S.Naturals)

    # restricted function pertains the original structure
    assert g.restrict(S.Naturals0)(2, evaluate=True) == 1

    # restriction of restriction can be evaluated.
    assert g.restrict(S.Naturals0).restrict(S.Naturals, evaluate=True) == g.restrict(S.Naturals)

def test_restriction():
    # RestrictedMap is surely restriction of base function.
    assert g.restrict(S.Naturals0).is_restriction(g)

    # restriction checked by domain relations
    f1, f2 = Map('f', S.Reals, S.Reals), Map('f', S.Naturals0, S.Naturals0)

    assert f2.is_restriction(f1)
    assert f2.inv().is_restriction(f1.inv())
    assert IdentityMap(S.Naturals0).is_restriction(IdentityMap(S.Reals))
    assert f2.restrict(S.Naturals).is_restriction(f2)

    # restriction checked by overriding _map_content
    assert h.is_restriction(f)

def test_ConstantMap():
    f = ConstantMap(1)
    assert f(x, evaluate=True) == 1

def test_isappliedmap():
    f1, f2 = Map('f', domain=S.Reals), Map('f', domain=S.Integers)
    g = Map('g')

    assert isappliedmap(f1(1), f1)
    assert not isappliedmap(f1(1), f2)
    assert isappliedmap(f2(1), f1)
    assert isappliedmap(f2(1), (f2, g))

def test_compatibility():
    from sympy.core.function import Function

    with warns_deprecated_sympy():
        assert f(x).args[0] == x
        assert isinstance(f(x), Function)

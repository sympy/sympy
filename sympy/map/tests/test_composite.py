from sympy import Symbol, S, FiniteSet, ask, Q, assuming
from sympy.map import Map, CompositeMap, IdentityMap, IteratedMap
from sympy.testing.pytest import raises

x = Symbol('x')
f = Map('f', domain=FiniteSet(1,2), codomain=S.Naturals0)
g = Map('g', domain=S.Reals, codomain=FiniteSet(3,4))
h, h1, h2, h3 = Map('h'), Map('h1'), Map('h2'), Map('h3')

def test_CompositeMap():

    # domain and codomain
    raises(TypeError, lambda: CompositeMap(f, g))
    gf = CompositeMap(g, f)
    assert gf.domain == FiniteSet(1,2)
    assert gf.codomain == FiniteSet(3,4)

    # flattening
    hh = CompositeMap(h, h)
    assert CompositeMap(hh, h, evaluate=True) == CompositeMap(h, h, h)

    # inverse cancellation
    assert CompositeMap(h, h.inv(), evaluate=True) == IdentityMap()
    assert CompositeMap(h1, h2, h3, h3.inv(), h2.inv(), evaluate=True) == h1

    # identity cancellation only when domain and codomain are strictly same
    assert CompositeMap(IdentityMap(S.Naturals0), f, evaluate=True) == f
    assert CompositeMap(IdentityMap(S.Reals), f, evaluate=True) != f
    assert CompositeMap(g, IdentityMap(S.Reals), evaluate=True) == g
    assert CompositeMap(g, IdentityMap(S.Naturals0), evaluate=True) != g

    # not invertible when domains and codomains are not suitable
    assert ask(Q.invertible(gf)) is False
    # not invertible when any of its argument is not invertible
    with assuming(~Q.invertible(h2)):
        assert ask(Q.invertible(CompositeMap(h1, h2, h3))) is False
    # invertible when all arguments are invertible and domains/codomains are suitable
    with assuming(Q.invertible(h1) & Q.invertible(h2) & Q.invertible(h3)):
        assert ask(Q.invertible(CompositeMap(h1, h2, h3))) is True
    # invertible if assumed to be invertible
    with assuming(Q.invertible(hh)):
        assert ask(Q.invertible(hh)) is True

    # inverse evaluation
    assert CompositeMap(h1, h2).inv(evaluate=True) == CompositeMap(h2.inv(), h1.inv())

    # evaluation
    class F1(Map):
        def eval(self, x):
            return x+2
    class F2(Map):
        def eval(self, x):
            return x*3
    f1, f2 = F1(), F2()
    assert CompositeMap(f1, f2)(x, evaluate=True) == 3*x + 2
    assert CompositeMap(f2, f1)(x, evaluate=True) == 3*x + 6

    # multivariate
    class F3(Map):
        def eval(self, tup):
            a,b = tup
            return a+b
    class F4(Map):
        def eval(self, tup):
            a,b = tup
            return (a+1, b+1)
    f3, f4 = F3(), F4()
    assert CompositeMap(f2, f3, f4)((1, 2), evaluate=True) == 15

    # composition can be defined
    class A(Map):
        def _eval_composite(self, other):
            if other == f1:
                return f2
    assert CompositeMap(A(), f1, evaluate=True) == f2

    # __matmul__ returns evaluated composition
    assert A()@f1 == f2

def test_IteratedMap():
    # domain and codomain
    raises(TypeError, lambda: IteratedMap(f, 2))

    # 0 and 1 exponents
    assert IteratedMap(h, 0, evaluate=True) == IdentityMap()
    assert IteratedMap(h, 1, evaluate=True) == h

    # negative exponents
    assert IteratedMap(h, -1, evaluate=True) == h.inv()
    assert IteratedMap(h, -2, evaluate=True) == IteratedMap(h.inv(), 2)

    # iteration can be defined
    class A(Map):
        def _eval_iteration(self, n):
            if n == S.One/2:
                return h1
    assert IteratedMap(A(), S.One/2, evaluate=True) == h1

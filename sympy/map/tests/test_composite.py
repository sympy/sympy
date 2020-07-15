from sympy import Symbol, S
from sympy.map import Map, CompositeMap, IdentityMap

x = Symbol('x')
f, g = Map('f', codomain=S.Reals), Map('g', domain=S.Naturals0)

def test_CompositeMap():

    # domain and codomain
    fg = CompositeMap(f, g)
    assert fg.domain == S.Naturals0
    assert fg.codomain == S.Reals

    # flattening
    ff = CompositeMap(f, f)
    assert CompositeMap(f, ff, evaluate=True) == CompositeMap(f, f, f)

    # inverse cancellation
    assert CompositeMap(g, f, f.inv(), f, evaluate=True) == CompositeMap(g, f)
    assert CompositeMap(g, f, f.inv(), g.inv(), evaluate=True) == IdentityMap(S.Naturals0)

    # inverse
    assert fg.inv(evaluate=True) == CompositeMap(g.inv(), f.inv())

    # evaluation
    class H1(Map):
        def eval(cls, x):
            return x+2
    class H2(Map):
        def eval(cls, x):
            return x*3
    h1, h2 = H1(), H2()
    assert CompositeMap(h1, h2)(x, evaluate=True) == 3*x + 2
    assert CompositeMap(h2, h1)(x, evaluate=True) == 3*x + 6

    # multivariate
    class H3(Map):
        def eval(cls, tup):
            a,b = tup
            return a+b
    class H4(Map):
        def eval(cls, tup):
            a,b = tup
            return (a+1, b+1)
    h3, h4 = H3(), H4()
    assert CompositeMap(h2, h3, h4)((1, 2), evaluate=True) == 15

from sympy import symbols
from sympy.map import Map, AppliedMap, CompositeMap, IdentityMap

x, y = symbols('x y')
f, g = Map(name='f'), Map(name='g')

def test_CompositeMap():
    from sympy.core.compatibility import iterable

    # inverse cancellation
    assert CompositeMap(g, f, f.inv(), f, evaluate=True) == CompositeMap(g, f)
    assert CompositeMap(g, f, f.inv(), g.inv(), evaluate=True) == IdentityMap()

    # inverse
    assert CompositeMap(f, g).inv(evaluate=True) == CompositeMap(g.inv(), f.inv())

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


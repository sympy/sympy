from sympy.diffgeom import Manifold, Patch, CoordSystem, Point
from sympy import symbols, Function

m = Manifold('m', 2)
p = Patch('p', m)
cs = CoordSystem('cs', p, ['a', 'b'])
cs_noname = CoordSystem('cs', p)
x, y = symbols('x y')
f = Function('f')
s1, s2 = cs.coord_functions()
v1, v2 = cs.base_vectors()
f1, f2 = cs.base_oneforms()


def test_point():
    point = Point(cs, [x, y])
    assert point == point.func(*point.args)
    assert point != Point(cs, [2, y])
    #TODO assert point.subs(x, 2) == Point(cs, [2, y])
    #TODO assert point.free_symbols == set([x, y])


def test_atomicclass_args():
    assert m.args == ()
    assert p.args == ()
    assert cs.args == ()
    assert cs_noname.args == ()
    assert s1.args == ()
    assert v1.args == ()


def test_subs():
    assert s1.subs(s1, s2) == s2
    assert v1.subs(v1, v2) == v2
    assert f1.subs(f1, f2) == f2
    assert (x*f(s1) + y).subs(s1, s2) == x*f(s2) + y
    assert (f(s1)*v1).subs(v1, v2) == f(s1)*v2
    assert (y*f(s1)*f1).subs(f1, f2) == y*f(s1)*f2

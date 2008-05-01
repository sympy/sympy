from sympy import Symbol, Function, Number, diff, Derivative, expand, sin, cos, exp, log

def test_simplemulti():
    f = Function("f")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    a = Number(1)
    b = Number("1.1")
    c = Number(1)/4
    res = [f(x),f(y),f(z)]
    assert f([x,y,z])==res
    assert f((x,y,z))==res
    assert f([x, [y, [a, 1], [b, z]], [[c]]])==[f(x), [f(y), [f(a), f(1)], [f(b), f(z)]], [[f(c)]]]

def test_somefuncs():
    x = Symbol("x")
    y = Symbol("y")
    for f in (sin, cos, exp, log):
        f([x,y])==[f(x), f(y)]

def test_diffmulti():
    f = Function("f")
    g = Function("g")
    h = Function("h")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    r = [x, y, z]
    assert diff(f(r), x)==[Derivative(f(x),x), 0, 0]
    assert diff(f(*r), r)==[Derivative(f(*r),x), Derivative(f(*r),y), Derivative(f(*r),z)]
    assert diff([f(*r), g(*r), h(*r)], r)==[
            [Derivative(f(*r), x), Derivative(f(*r), y), Derivative(f(*r), z)],
            [Derivative(g(*r), x), Derivative(g(*r), y), Derivative(g(*r), z)],
            [Derivative(h(*r), x), Derivative(h(*r), y), Derivative(h(*r), z)]]
    assert diff(f(x),x,[0,1,2,3])==[diff(f(x),x,0), diff(f(x),x,1), diff(f(x),x,2), diff(f(x),x,3)]

def test_expandmulti():
    x = Symbol("x")
    y = Symbol("y")
    assert expand([(x+y)**2, (x-1)**2])==[expand((x+y)**2), expand((x-1)**2)]
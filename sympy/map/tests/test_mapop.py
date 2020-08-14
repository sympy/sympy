from sympy import Map, MapAdd, MapMul, MapPow, ConstantMap
from sympy.abc import x

f, g = Map('f'), Map('g')

def test_MapAdd():
    assert MapAdd(f, g, evaluate=True) == f+g
    assert MapAdd(f, f, evaluate=True) == MapMul(2, f) == f+f == 2*f
    assert f+g+f+g == 2*f + 2*g
    assert ConstantMap(0) + f == f

    assert ConstantMap(2) + ConstantMap(3) == ConstantMap(5)
    assert ConstantMap(2) + ConstantMap(x) == ConstantMap(2+x)

    assert (f+g)(x, evaluate=True) == f(x) + g(x)

def test_MapMul():
    assert MapMul(f, g, evaluate=True) == f*g
    assert MapMul(f, f, evaluate=True) == MapPow(f, 2) == f*f == f**2
    assert f*g*f*g == f**2 * g**2
    assert 1*f == f*1 == f
    assert 2*f*3 == 6*f

    assert f*0 == 0*f == ConstantMap(0)
    assert ConstantMap(3)*4 == ConstantMap(3)*ConstantMap(4) == ConstantMap(12)

    assert (f*g)(x, evaluate=True) == f(x)*g(x)

def test_MapPow():
    assert MapPow(f, 2, evaluate=True) == f**2
    assert f**0 == ConstantMap(1)
    assert f**1 == f

    assert (f**2)(x, evaluate=True) == f(x)**2

"""Test modules.py code"""

from sympy.polys.agca.modules import FreeModule
from sympy.polys import CoercionFailed, QQ
from sympy.abc import x, y

def test_FreeModule():
    M1 = FreeModule(QQ[x], 2)
    assert M1 == FreeModule(QQ[x], 2)
    assert M1 != FreeModule(QQ[y], 2)
    assert M1 != FreeModule(QQ[x], 3)
    M2 = FreeModule(QQ.poly_ring(x, order="ilex"), 2)

    assert [x, 1] in M1
    assert [x] not in M1
    assert [2, y] not in M1
    assert [1/(x + 1), 2] not in M1

    e = M1.convert([x, x**2 + 1])
    X = QQ[x].convert(x)
    assert e == [X, X**2 + 1]
    assert e == [x, x**2 + 1]
    assert 2*e == [2*x, 2*x**2 + 2]
    assert e*2 == [2*x, 2*x**2 + 2]
    assert e/2 == [x/2, (x**2 + 1)/2]
    assert x*e == [x**2, x**3 + x]
    assert e*x == [x**2, x**3 + x]
    assert X*e == [x**2, x**3 + x]
    assert e*X == [x**2, x**3 + x]

    assert [x, 1] in M2
    assert [x] not in M2
    assert [2, y] not in M2
    assert [1/(x + 1), 2] in M2

    e = M2.convert([x, x**2 + 1])
    X = QQ.poly_ring(x, order="ilex").convert(x)
    assert e == [X, X**2 + 1]
    assert e == [x, x**2 + 1]
    assert 2*e == [2*x, 2*x**2 + 2]
    assert e*2 == [2*x, 2*x**2 + 2]
    assert e/2 == [x/2, (x**2 + 1)/2]
    assert x*e == [x**2, x**3 + x]
    assert e*x == [x**2, x**3 + x]
    assert e/(1 + x) == [x/(1 + x), (x**2 + 1)/(1 + x)]
    assert X*e == [x**2, x**3 + x]
    assert e*X == [x**2, x**3 + x]

    M3 = FreeModule(QQ[x, y], 2)
    assert M3.convert(e) == M3.convert([x, x**2 + 1])

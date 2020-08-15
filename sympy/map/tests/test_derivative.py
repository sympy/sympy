from sympy.map import Map

f = Map('f')

def test_DerivativeFunction():
    assert f.diff(1).diff(1) == f.diff((1,2))

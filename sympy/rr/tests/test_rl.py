from sympy.rr.rl import rmid, glom
from sympy import Basic

def test_rmid():
    rmzeros = rmid(lambda x: x == 0)
    assert rmzeros(Basic(0, 1)) == Basic(1)
    assert rmzeros(Basic(0, 0)) == Basic(0)
    assert rmzeros(Basic(2, 1)) == Basic(2, 1)

def test_glom():
    conglomerate = glom(lambda num, x: num * x)
    assert conglomerate(Basic(1, 2, 2)) == Basic(1, 4)
    conglomerate = glom(lambda num, x: x ** num)
    assert conglomerate(Basic(1, 3, 3)) == Basic(1, 9)

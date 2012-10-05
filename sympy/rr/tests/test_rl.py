from sympy.rr.rl import rmid, glom, flatten, unpack, sort
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

def test_flatten():
    assert flatten(Basic(1, 2, Basic(3, 4))) == Basic(1, 2, 3, 4)

def test_unpack():
    assert unpack(Basic(2)) == 2

def test_sort():
    assert sort(str)(Basic(3,1,2)) == Basic(1,2,3)


from sympy.rr.rl import rmid
from sympy import Basic

class Container(Basic):
    pass

def test_rmid():
    rmzeros = rmid(lambda x: x == 0)
    x = Basic(0, 0, 0)
    assert rmzeros(Basic(0, 1)) == Basic(1)
    assert rmzeros(Basic(0, 0)) == Basic(0)
    assert rmzeros(Basic(2, 1)) == Basic(2, 1)


from sympy.rr import typed, rmid
from sympy import Basic

def test_typed():
    class A(Basic): pass
    class B(Basic): pass
    rmzeros = rmid(lambda x: x == 0)
    rmones  = rmid(lambda x: x == 1)
    remove_something = typed({A: rmzeros, B: rmones})

    assert remove_something(A(0, 1)) == A(1)
    assert remove_something(B(0, 1)) == B(0)


from sympy.rules import typed, rm_id, do_one
from sympy import Basic

def test_typed():
    class A(Basic):
        pass
    class B(Basic):
        pass
    rmzeros = rm_id(lambda x: x == 0)
    rmones  = rm_id(lambda x: x == 1)
    remove_something = typed({A: rmzeros, B: rmones})

    assert remove_something(A(0, 1)) == A(1)
    assert remove_something(B(0, 1)) == B(0)

def test_do_one():
    rl1 = lambda x: 2 if x == 1 else x
    rl2 = lambda x: 3 if x == 2 else x

    rule = do_one(rl1, rl2)
    assert rule(1) == 2
    assert rule(rule(1)) == 3

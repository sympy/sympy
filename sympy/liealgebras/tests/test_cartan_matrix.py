from sympy.liealgebras.cartan_matrix import CartanMatrix
from sympy.matrices import Matrix

def test_CartanMatrix_TypeA():
    c = CartanMatrix("A3")
    m = Matrix(3, 3, [2, -1, 0, -1, 2, -1, 0, -1, 2])
    assert c == m

def test_CartanMatrix_TypeB():
    c = CartanMatrix("B4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2])
    assert c == m

def test_CartanMatrix_TypeC():
    c = CartanMatrix("C4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -2, 2])
    assert c == m

def test_CartanMatrix_TypeD():
    c = CartanMatrix("D4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, -1, 0, -1, 2, 0, 0, -1, 0, 2])
    assert c == m

def test_CartanMatrix_TypeE6():
    e6 = CartanMatrix(["E", 6])
    e6t = Matrix(6,6,[
        2,-1,0,0,0,0,
        -1,2,-1,0,0,0,
        0,-1,2,-1,0,-1,
        0,0,-1,2,-1,0,
        0,0,0,-1,2,0,
        0,0,-1,0,0,2
    ])
    assert e6t == e6

def test_CartanMatrix_TypeE7():
    e7 = CartanMatrix(["E", 7])
    e7t = Matrix([
        [2, -1, 0, 0, 0, 0, 0],
        [-1, 2, -1, 0, 0, 0, 0],
        [0, -1, 2, -1, 0, 0, 0],
        [0, 0, -1, 2, -1, 0, -1],
        [0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, -1, 2, 0],
        [0, 0, 0, -1, 0, 0, 2]
    ])
    assert e7t == e7

def test_CartanMatrix_TypeE8():
    e8 = CartanMatrix(["E", 8])
    e8t = Matrix([
        [2, -1, 0, 0, 0, 0, 0, 0],
        [-1, 2, -1, 0, 0, 0, 0, 0],
        [0, -1, 2, -1, 0, 0, 0, 0],
        [0, 0, -1, 2, -1, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, -1],
        [0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, -1, 2, 0],
        [0, 0, 0, 0, -1, 0, 0, 2]
    ])
    assert e8t == e8

def test_CartanMatrix_TypeF4():
    a = CartanMatrix(["F",4])
    mt = Matrix(4, 4, [
        2, -1, 0, 0,
        -1,2,-2,0,
        0,-1,2,-1,
        0,0,-1,2
    ])
    assert a == mt

def test_CartanMatrix_TypeG2():
    a = CartanMatrix(["G",2])
    mt = Matrix(2, 2, [2, -1, -3, 2])
    assert a == mt

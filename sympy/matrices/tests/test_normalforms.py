from sympy import Symbol, Poly
from sympy.polys.solvers import RawMatrix as Matrix
from sympy.matrices.normalforms import invariant_factors, smith_normal_form, hermite_normal_form
from sympy.polys.domains import ZZ, QQ

def test_smith_normal():
    m = Matrix([[12, 6, 4,8],[3,9,6,12],[2,16,14,28],[20,10,10,20]])
    setattr(m, 'ring', ZZ)
    smf = Matrix([[1, 0, 0, 0], [0, 10, 0, 0], [0, 0, -30, 0], [0, 0, 0, 0]])
    assert smith_normal_form(m) == smf

    x = Symbol('x')
    m = Matrix([[Poly(x-1), Poly(1, x),Poly(-1,x)],
                [0, Poly(x), Poly(-1,x)],
                [Poly(0,x),Poly(-1,x),Poly(x)]])
    setattr(m, 'ring', QQ[x])
    invs = (Poly(1, x), Poly(x - 1), Poly(x**2 - 1))
    assert invariant_factors(m) == invs

    m = Matrix([[2, 4]])
    setattr(m, 'ring', ZZ)
    smf = Matrix([[2, 0]])
    assert smith_normal_form(m) == smf

def test_hermite_normal():
    m = Matrix([[3, 3, 1, 4], [0, 1, 0, 0], [0, 0, 19, 16], [0, 0, 0, 3]])
    setattr(m, 'ring', ZZ)
    H, U = hermite_normal_form(m)
    assert H == Matrix([[3, 0, 1, 1], [0, 1, 0, 0], [0, 0, 19, 1], [0, 0, 0, 3]])
    assert U == Matrix([[1, -3, 0, -1], [0, 1, 0, 0], [0, 0, 1, -5], [0, 0, 0, 1]])

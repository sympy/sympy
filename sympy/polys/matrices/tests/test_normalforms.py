from sympy.testing.pytest import raises

from sympy import Symbol, sympify
from sympy.polys.matrices.normalforms import invariant_factors, smith_normal_form
from sympy.polys.domains import ZZ, QQ
from sympy.polys.matrices import DomainMatrix


def test_smith_normal():

    def DM(elems, domain):
        conv = lambda e: domain.from_sympy(sympify(e))
        elems = [[conv(e) for e in row] for row in elems]
        return DomainMatrix(elems, (len(elems), len(elems[0])), domain)

    m = DM([[12, 6, 4, 8], [3, 9, 6, 12], [2, 16, 14, 28], [20, 10, 10, 20]], ZZ)
    smf = DM([[1, 0, 0, 0], [0, 10, 0, 0], [0, 0, -30, 0], [0, 0, 0, 0]], ZZ)
    assert smith_normal_form(m).to_dense() == smf

    x = Symbol('x')
    m = DM([[x-1,  1, -1],
            [  0,  x, -1],
            [  0, -1,  x]], QQ[x])
    dx = m.domain.gens[0]
    assert invariant_factors(m) == (1, dx-1, dx**2-1)

    zr = DomainMatrix([], (0, 2), ZZ)
    zc = DomainMatrix([[], []], (2, 0), ZZ)
    assert smith_normal_form(zr).to_dense() == zr
    assert smith_normal_form(zc).to_dense() == zc

    assert smith_normal_form(DM([[2, 4]], ZZ)).to_dense() == DM([[2, 0]], ZZ)
    assert smith_normal_form(DM([[0, -2]], ZZ)).to_dense() == DM([[-2, 0]], ZZ)
    assert smith_normal_form(DM([[0], [-2]], ZZ)).to_dense() == DM([[-2], [0]], ZZ)

    m =   DM([[3, 0, 0, 0], [0, 0, 0, 0], [0, 0, 2, 0]], ZZ)
    snf = DM([[1, 0, 0, 0], [0, 6, 0, 0], [0, 0, 0, 0]], ZZ)
    assert smith_normal_form(m).to_dense() == snf

    raises(ValueError, lambda: smith_normal_form(DM([[1]], ZZ[x])))

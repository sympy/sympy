from sympy.polys.domains import GF, QQ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.exceptions import DDMRankError
from sympy.polys.matrices.noninvertible import invertible_supplement
from sympy.testing.pytest import raises


def test_invertible_supplement_1():
    M = DomainMatrix([[QQ(1), QQ(7), QQ(0)], [QQ(2), QQ(3), QQ(4)]], (2, 3), QQ).transpose()

    # First supplement over QQ:
    B = invertible_supplement(M)
    assert B[:, :2] == M
    assert B[:, 2] == DomainMatrix.eye(3, QQ).to_dense()[:, 0]

    # Now supplement over GF(7):
    M = M.convert_to(GF(7))
    B = invertible_supplement(M)
    assert B[:, :2] == M
    # When we work mod 7, first col of M goes to [1, 0, 0],
    # so the supplementary vector cannot equal this, as it did
    # when we worked over QQ. Instead, we get the second std basis vector:
    assert B[:, 2] == DomainMatrix.eye(3, GF(7)).to_dense()[:, 1]


def test_invertible_supplement_2():
    M = DomainMatrix([[QQ(1), QQ(0), QQ(0)], [QQ(2), QQ(0), QQ(0)]], (2, 3), QQ).transpose()
    with raises(DDMRankError):
        invertible_supplement(M)

from sympy import QQ
from sympy.matrices import Matrix
from sympy.polys import Poly, cyclotomic_poly


def test_submodule():
    from sympy.abc import zeta
    T = Poly(cyclotomic_poly(5))
    K = QQ.algebraic_field((T, zeta))
    ZK = K.int_ring
    C = 2 * Matrix.eye(4)
    M = ZK.submodule(C)
    assert M.matrix == C
    for u in range(4):
        for v in range(u, 4):
            for k in range(4):
                assert M.mult_tab[u][v][k] == 2 * ZK.mult_tab[u][v][k]

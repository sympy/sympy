from sympy import symbols

from sympy.physics.quantum.kronecker import KroneckerDelta

def test_kronecker_delta():
    i, j, k = symbols('i j k')
    D = KroneckerDelta
    assert D(i, i) == 1
    assert D(i, i + 1) == 0
    assert D(0, 0) == 1
    assert D(0, 1) == 0
    # assert D(i, i + k) == D(0, k)
    assert D(i + k, i + k) == 1
    assert D(i + k, i + 1 + k) == 0
    assert D(i, j).subs(dict(i=1, j=0)) == 0
    assert D(i, j).subs(dict(i=3, j=3)) == 1

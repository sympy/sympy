from sympy.tensor.indexed_integrals import _IndexedIntegral
from sympy import IndexedBase, Idx, symbols, sin, cos


def test_indexed_integrals():
    A = IndexedBase('A')
    i, j = symbols('i j', integer=True)
    a1, a2 = symbols('a1:3', cls=Idx)
    assert isinstance(a1, Idx)

    assert _IndexedIntegral(1, A[i]).doit() == A[i]
    assert _IndexedIntegral(A[i], A[i]).doit() == A[i]**2/2
    assert _IndexedIntegral(A[j], A[i]).doit() == A[i]*A[j]
    assert _IndexedIntegral(A[i]*A[j], A[i]).doit() == A[i]**2*A[j]/2
    assert _IndexedIntegral(sin(A[i]), A[i]).doit() == -cos(A[i])
    assert _IndexedIntegral(sin(A[j]), A[i]).doit() == sin(A[j])*A[i]

    assert _IndexedIntegral(1, A[a1]).doit() == A[a1]
    assert _IndexedIntegral(A[a1], A[a1]).doit() == A[a1]**2/2
    assert _IndexedIntegral(A[a2], A[a1]).doit() == A[a1]*A[a2]
    assert _IndexedIntegral(A[a1]*A[a2], A[a1]).doit() == A[a1]**2*A[a2]/2
    assert _IndexedIntegral(sin(A[a1]), A[a1]).doit() == -cos(A[a1])
    assert _IndexedIntegral(sin(A[a2]), A[a1]).doit() == sin(A[a2])*A[a1]

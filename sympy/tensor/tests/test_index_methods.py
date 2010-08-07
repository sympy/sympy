from sympy.core import symbols, S
from sympy.utilities.pytest import raises
from sympy.tensor.indexed import Idx, IndexedBase, Indexed
from sympy.tensor.index_methods import IndexConformanceException

# import test
from sympy import get_contraction_structure, get_indices

def test_trivial_indices():
    x, y = symbols('x y')
    assert get_indices(x) == (set([]), {})
    assert get_indices(x*y) == (set([]), {})
    assert get_indices(x + y) == (set([]), {})
    assert get_indices(x**y) == (set([]), {})

def test_get_indices_Indexed():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    assert get_indices(x[i, j]) == (set([i, j]), {})
    assert get_indices(x[j, i]) == (set([j, i]), {})

def test_get_indices_mul():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    assert get_indices(x[j]*y[i]) == (set([i, j]), {})
    assert get_indices(x[i]*y[j]) == (set([i, j]), {})

def test_get_indices_exceptions():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    raises(IndexConformanceException, 'get_indices(x[i] + y[j])')

def test_scalar_broadcast():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    assert get_indices(x[i] + y[i, i]) == (set([i]), {})

def test_get_indices_add():
    x = IndexedBase('x')
    y = IndexedBase('y')
    A = IndexedBase('A')
    i, j, k = Idx('i'), Idx('j'), Idx('k')
    assert get_indices(x[i] + 2*y[i]) == (set([i,]), {})
    assert get_indices(y[i] + 2*A[i, j]*x[j]) == (set([i,]), {})
    assert get_indices(y[i] + 2*(x[i] + A[i, j]*x[j])) == (set([i,]), {})
    assert get_indices(y[i] + x[i]*(A[j, j] + 1)) == (set([i,]), {})
    assert get_indices(y[i] + x[i]*x[j]*(y[j] + A[j, k]*x[k])) == (set([i,]), {})

def test_get_contraction_structure_basic():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    assert get_contraction_structure(x[i]*y[j]) == {None: set([x[i]*y[j]])}
    assert get_contraction_structure(x[i] + y[j]) == {None: set([x[i], y[j]])}
    assert get_contraction_structure(x[i]*y[i]) == {(i,): set([x[i]*y[i]])}
    assert get_contraction_structure(1 + x[i]*y[i]) == {None: set([S.One]), (i,): set([x[i]*y[i]])}

def test_get_contraction_structure_complex():
    x = IndexedBase('x')
    y = IndexedBase('y')
    A = IndexedBase('A')
    i, j, k = Idx('i'), Idx('j'), Idx('k')
    expr1 = y[i] + A[i, j]*x[j]
    d1 = {None: set([y[i]]), (j,): set([A[i, j]*x[j]])}
    assert get_contraction_structure(expr1) == d1
    expr2 = expr1*A[k, i] + x[k]
    d2 = {None: set([x[k]]), (i,): set([expr1*A[k, i]]), expr1*A[k, i]: [d1]}
    assert get_contraction_structure(expr2) == d2

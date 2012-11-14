from sympy.computations.core import Computation, unique

def test_Computation():
    C = Computation()
    assert hasattr(C, 'inputs')
    assert hasattr(C, 'outputs')
    assert hasattr(C, 'edges')

def test_unique():
    assert tuple(unique((1, 3, 1, 2))) == (1, 3, 2)

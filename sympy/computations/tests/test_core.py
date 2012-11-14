from sympy.computations.core import Computation

def test_Computation():
    C = Computation()
    assert hasattr(C, 'inputs')
    assert hasattr(C, 'outputs')

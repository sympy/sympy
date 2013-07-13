from sympy.physics.quantum.circuitplot import labler
from sympy.physics.quantum.gate import CNOT, H, X, Z, SWAP, CGate, S, T

def test_labler():
    """Test the labler utility"""
    assert labler(2) == ['q_1', 'q_0']
    assert labler(3,'j') == ['j_2', 'j_1', 'j_0']
    return

if __name__ == '__main__':
    test_ex4()

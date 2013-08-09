from sympy.physics.quantum.qasm import Qasm
from sympy.physics.quantum.gate import CNOT,H,SWAP,CGate,T,S,Z,X
from sympy.physics.quantum.circuitplot import Mz

def test_qasm_ex1():
    q = Qasm('qubit q0','qubit q1','h q0','cnot q0,q1')
    assert q.get_circuit() == CNOT(1,0)*H(1)

def test_qasm_ex1_methodcalls():
    q = Qasm()
    q.qubit('q_0')
    q.qubit('q_1')
    q.h('q_0')
    q.cnot('q_0','q_1')    
    assert q.get_circuit() == CNOT(1,0)*H(1)

def test_qasm_swap():
    q = Qasm('qubit q0','qubit q1','cnot q0,q1','cnot q1,q0','cnot q0,q1')
    assert q.get_circuit() == CNOT(1,0)*CNOT(0,1)*CNOT(1,0)


def test_qasm_ex2():
    q = Qasm('qubit q_0','qubit q_1','qubit q_2','h  q_1',
             'cnot q_1,q_2','cnot q_0,q_1','h q_0',
             'measure q_1','measure q_0',
             'c-x q_1,q_2','c-z q_0,q_2')
    assert q.get_circuit() == CGate(2,Z(0))*CGate(1,X(0))*Mz(2)*Mz(1)*H(2)*CNOT(2,1)*CNOT(1,0)*H(1)

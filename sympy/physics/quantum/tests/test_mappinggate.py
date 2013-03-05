"""Tests for mappinggate.py"""

from sympy import I, Integer, Mul, Add
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.mappinggate import MappingGate

M = MappingGate(('00', I, '11'), ('01', '10'), ('10', '01'), ('11', -I, '00'))
M_rep = represent(M, format='sympy')

def test_MappingGate():
    assert M.get_final_state('00') == I*Qubit('11')
    assert M.mapping['00'] == I*Qubit('11')
    assert qapply(M*Qubit('01')) == Qubit('10')
    for i in range(2**M.nqubits):
        terms = []
        temp = bin(i)[2:]
        string = (M.nqubits - len(temp))*'0' + temp
        initial = Qubit(string)
        fin = M.get_final_state(initial)
        terms.append(fin*Dagger(initial))
        assert M.rewrite('op').doit() == Add(*terms)
    for i, f in M.mapping.items():
        col = IntQubit(i).as_int()
        terms = split_qexpr_parts(f)
        if len(terms[1]) == 0:
            row = IntQubit(*terms[0]).as_int()
            scalar = Integer(1)
        else:
            row = IntQubit(*terms[1]).as_int()
            scalar = Mul(*terms[0])
        assert M_rep[row, col] == scalar

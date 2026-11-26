"""Tests for mappinggate.py"""

from sympy import I, Integer, Mul, Add
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.qexpr import split_qexpr_parts
from sympy.physics.quantum.hilbert import ComplexSpace
from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.mappinggate import MappingGate
from sympy.external import import_module
from sympy.utilities.pytest import skip

np = import_module('numpy', min_python_version=(2, 6))
scipy = import_module('scipy', __import__kwargs={'fromlist': ['sparse']})

# All 3 ways produce same Qubit Mappings
M = MappingGate(('00', I, '11'), ('01', '10'), ('10', '01'), ('11', -I, '00'))
M_half = MappingGate(('00', I, '11'), ('01', '10'))
d = dict({Qubit('00'):I*Qubit('11'), Qubit('01'):Qubit('10')})
M_dict = MappingGate(d)

M_rep = represent(M, format='sympy')

def test_MappingGate():
    assert M.get_final_state('00') == I*Qubit('11')
    assert M.mapping[Qubit('00')] == I*Qubit('11')
    assert qapply(M*Qubit('01')) == Qubit('10')
    assert M.hilbert_space == ComplexSpace(2)**M.nqubits
    # Shows same qubit mappings
    assert M.args == M_half.args
    assert M.args == M_dict.args

    terms = []
    for i in range(2**M.nqubits):
        initial = Qubit(IntQubit(i, M.nqubits))
        fin = M.get_final_state(initial)
        terms.append(fin*Dagger(initial))
    result = Add(*terms)
    assert M.rewrite('op') == result

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

    if not np:
        skip("numpy not installed or Python too old.")

    M_rep_numpy = represent(M, format='numpy')
    for i, f in M.mapping.items():
        col = IntQubit(i).as_int()
        terms = split_qexpr_parts(f)
        if len(terms[1]) == 0:
            row = IntQubit(*terms[0]).as_int()
            scalar = Integer(1)
        else:
            row = IntQubit(*terms[1]).as_int()
            scalar = Mul(*terms[0])
        assert M_rep_numpy[row, col] == complex(scalar)

    if not np:
        skip("numpy not installed or Python too old.")
    if not scipy:
        skip("scipy not installed.")
    else:
        sparse = scipy.sparse

    M_rep_scipy = represent(M, format='scipy.sparse')
    for i, f in M.mapping.items():
        col = IntQubit(i).as_int()
        terms = split_qexpr_parts(f)
        if len(terms[1]) == 0:
            row = IntQubit(*terms[0]).as_int()
            scalar = Integer(1)
        else:
            row = IntQubit(*terms[1]).as_int()
            scalar = Mul(*terms[0])
        col = float(col)
        row = float(row)
        scalar = complex(scalar)
        assert M_rep_scipy[row, col] == scalar

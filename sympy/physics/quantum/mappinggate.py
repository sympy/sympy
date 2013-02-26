"""Mapping Quantum Gates"""

from sympy import Integer, conjugate, Add
from sympy.core.containers import Dict
from sympy.core.basic import Basic
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.gate import Gate
from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.matrixutils import matrix_eye
from sympy.physics.quantum.qexpr import split_qexpr_parts

#-----------------------------------------------------------------------------------

class MappingGate(Gate):
    """

    Parameters
    ==========

    args : tuple
        
        The list of initial state and final state pairs. There must be an even
        number of arguments so every initial state has a final state. The initial
        and final states can be separated by a scalar. If no scalar is given, 1 is
        assumed.

    """

    @classmethod
    def _eval_args(cls, args):
        temp = {}
        for arg in args:
            i = Qubit(arg[0])
            if len(arg) == 2:
                scalar = Integer(1)
                f = Qubit(arg[1])
            elif len(arg) == 3:
                scalar = arg[1]
                f = Qubit(arg[2])
            else:
                raise ValueError('Too many scalar arguments')
            if i.nqubits != f.nqubits:
                raise ValueError('Number of qubits for each state do not match')
            temp[i] = scalar*f
            #temp[Dagger(i)] = conjugate(scalar)*Dagger(f)
        new_args = Dict(temp)
        return (new_args,)

    @classmethod
    def _eval_hilbert_space(cls, args):
        pass

    @property
    def size(self):
        """Total number of initial and final state pairs"""
        return int(len(self.args[0])/2)

    @property
    def mapping(self):
        return self.args[0]

    @property
    def nqubits(self):
        """Gives the dimension of the matrix representation"""
        return self.args[0].args[0].args[0].nqubits

    def get_final_state(self, qubit):
        i = Qubit(qubit)
        f = self.mapping.get(i, None)
        if f is None:
            return i
        else:
            return f

    def _apply_operator_Qubit(self, qubit):
        return self.get_final_state(qubit)

    def _eval_rewrite_as_op(self, *args):
        terms = []
        for i in range(2**self.nqubits):
            temp = bin(i)[2:]
            string = (self.nqubits - len(temp))*'0' + temp
            initial = Qubit(string)
            fin = self.get_final_state(initial)
            terms.append(fin*Dagger(initial))
        return Add(*terms)

    def _represent_default_basis(self, **options):
        return self._represent_ZGate(None, **options)

    def _represent_ZGate(self, basis, **options):
        format = options.get('format','sympy')
        spmatrix = options.get('spmatrix', 'csr')
        matrix = matrix_eye(2**self.nqubits, **options)
        for i, f in self.mapping.items():
            col = IntQubit(i).as_int()
            terms = split_qexpr_parts(f)
            if len(terms[0]) == 1:
                row = IntQubit(*terms[1]).as_int()
                scalar = terms[0]
            else:
                row = IntQubit(*terms[0]).as_int()
                scalar = Integer(1)
            matrix[col, col] = Integer(0)
            matrix[row, col] = scalar
        return matrix

"""Mapping Quantum Gates"""

from sympy import Integer
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
        number of arguments so every initial state has a final state.

    """

    @classmethod
    def _eval_args(cls, args):
        temp = {}
        for arg in args:
            i = Qubit(arg[0])
            f = Qubit(arg[1])
            #parts = split_qexpr_parts(f)
            #f = Qubit(parts[1])
            #scalar = parts[0]
            if i.nqubits != f.nqubits:
                raise ValueError('Number of qubits for each state do not match')
            temp[i] = f
            #temp[Dagger(i)] = conjugate(scalar)*Dagger(f)
        new_args = Dict(temp)
        print new_args
        return new_args

    @classmethod
    def _eval_hilbert_space(cls, args):
        pass

    @property
    def size(self):
        """Total number of initial and final state pairs"""
        return int(len(self.args)/2)

    @property
    def nqubits(self):
        """Gives the dimension of the matrix representation"""
        return self.args[0].nqubits

    def get_final_state(self, qubit):
        i = Qubit(qubit)
        f = self.args.get(i, None)
        if f is None:
            return i
        else:
            return f

    def _apply_operator_Qubit(self, qubit):
        return self.get_final_state(qubit)

    def _eval_rewrite_as_op(self, *args):
        terms = []
        for i in range(2**self.nqubits):
            init = Qubit(IntQubit(i))
            fin = self.get_final_state(init)
            term[i] = fin*Dagger(init)
        return Add(*terms)

    def _represent_default_basis(self, **options):
        return self._represent_ZGate(None, **options)

    def _represent_ZGate(self, basis, *args, **options):
        format = options.get('format','sympy')
        spmatrix = options.get('spmatrix', 'csr')
        matrix = matrix_eye(2**self.nqubits, **options)
        for i, f in self.args.items():
            col = IntQubit(i).as_int()
            row = IntQubit(f).as_int()
            matrix[col, col] = Integer(0)
            matrix[row, col] = f
        if format == 'scipy.sparse':
            matrix = matrix.tocsr()
        return matrix

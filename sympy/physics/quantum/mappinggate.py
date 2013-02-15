"""Mapping Quantum Gates"""

from sympy import Integer
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.gate import Gate
from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.matrixutils import matrix_zeros

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
        if len(args) % 2 == 1:
            raise ValueError('Need even number of arguments') 
        values = matrix_zeros(1, len(args))
        for i in range(len(args)):
            values[i] = Qubit(str(args[i]))
        args = values
        for i in range(len(args) - 1):
            if args[i].nqubits != args[i + 1].nqubits:
                raise ValueError('Number of qubits for each state do not match')
        return args

    @classmethod
    def _eval_hilbert_space(cls, args):
        pass

    @property
    def size(self, *args):
        """Total number of initial and final state pairs"""
        return int(len(self.args)/2)

    @property
    def mat(self, *args):
        """Gives the dimension of the matrix representation"""
        return 2**self.args[0].nqubits

    def _apply_operator_Qubit(self, qubits, *args):
        for i in range(self.size):
            if qubits == self.args[2*i]:
                return self.args[2*i + 1]
            elif qubits.nqubits != self.args[2*i].nqubits:
                raise ValueError('Initial state not specified')

    def _eval_rewrite_as_op(self, *args):
        product = matrix_zeros(1, len(self.args))
        for i in range(self.size):
            value = self.args[2*i + 1]*Dagger(self.args[2*i])
            product[i] = value
        return sum(product)

    def _represent_default_basis(self, **options):
        return self._represent_ZGate(None, **options)

    def _represent_ZGate(self, basis, *args, **options):
        format = options.get('format','sympy')
        spmatrix = options.get('spmatrix', 'csr')
        matrix = matrix_zeros(self.mat, self.mat, **options)
        for i in range(self.size):
            for j in range(self.size):
                for k in range(self.size):
                    if IntQubit(j).as_int() == IntQubit(self.args[2*i + 1]).as_int():
                        if IntQubit(k).as_int() == IntQubit(self.args[2*i]).as_int():
                            if format == 'scipy.sparse':
                                matrix[j, k] = matrix[j, k] + float(Integer(1))
                            else:
                                matrix[j, k] = matrix[j, k] + Integer(1)
                        else:
                            matrix[j, k] = matrix[j, k]
                    else:
                        matrix[j, k] = matrix[j, k]
        if format == 'scipy.sparse':
            matrix = matrix.tocsr()
        return matrix

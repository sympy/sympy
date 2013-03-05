"""Mapping Quantum Gates"""

from sympy import Integer, conjugate, Add, Mul
from sympy.core.containers import Dict
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

        arg[0] = initial state
        arg[1] = scalar or final state
        arg[2] = None or final state

        The list of initial state and final state pairs. The final states are to be
        multiplied by some scalar. If no scalar is given, 1 is assumed.

    Examples
    ========

    Creating a Mapping Gate and checking its arguments and properties. Getting final state from
    initial state.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy import I
        >>> M = MappingGate(('00', I, '11'), ('01', '10'), ('10', '01'), ('11', -I,
        '00'))
        >>> M.args
        ({|00>: I*|11>, |01>: |10>, |10>: |01>, |11>: -I*|00>},)
        >>> M.nqubits
        2
        >>> M.mapping[Qubit('00')]
        I*|11>
        >>> M.get_final_state('00')
        I*|11>

    Using qapply on initial states returns the final states.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy import I
        >>> from sympy.physics.quantum.qapply import qapply
        >>> from sympy.physics.quantum.qubit import Qubit
        >>> M = MappingGate(('00', I, '11'), ('01', '10'), ('10', '01'), ('11', -I,
        '00'))
        >>> q = Qubit('00') + Qubit('01')
        >>> qapply(M*q)
        |10> + I*|11>

    The MappingGate can be rewritten as an outer product of states. We will show two
    examples: one where all four states are given and one where only one state is
    given. If not all initial states are specified they return themselves as final
    states.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy import I
        >>> M = MappingGate(('00', I, '11'), ('01', '10'), ('10', '01'), ('11', -I,
        '00'))
        >>> M.rewrite('op')
        |01><10| + |10><01| - I*|00>*<11| + I*|11>*<00|
        >>> M = MappingGate(('00', -1, '00'))
        >>> M.rewrite('op')
        |01><01| + |10><10| + |11><11| - |00>*<00|

    The MappingGate is also expressed as a matrix where the rows and columns
    represent the Qubits.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy.physics.quantum.represent import represent
        >>> from sympy import I
        >>> M = MappingGate(('00', I, '11'), ('01', '10'), ('10', '01'), ('11', -I,
        '00'))
        >>> represent(M)
        [0, 0, 0, -I]
        [0, 0, 1,  0]
        [0, 1, 0,  0]
        [I, 0, 0,  0]

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
    def mapping(self):
        return self.args[0]

    @property
    def nqubits(self):
        """Gives the dimension of the matrix representation"""
        return self.args[0].args[0].args[0].nqubits

    def get_final_state(self, qubit):
        """Returns the final state for a given initial state, if initial state is
        not mapped to a final state the initial state is returned."""
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
            if len(terms[1]) == 0:
                row = IntQubit(*terms[0]).as_int()
                scalar = Integer(1)
            else:
                row = IntQubit(*terms[1]).as_int()
                scalar = Mul(*terms[0])
            matrix[col, col] = Integer(0)
            matrix[row, col] = scalar
        return matrix

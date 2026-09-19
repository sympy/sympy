"""Mapping Quantum Gates

TODO:
- Enable sparse mappings for large numbers of qubits
"""

from sympy import Integer, conjugate, Add, Mul
from sympy.core.containers import Dict
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.gate import Gate
from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.matrixutils import matrix_eye
from sympy.physics.quantum.qexpr import split_qexpr_parts
from sympy.physics.quantum.hilbert import ComplexSpace

#-----------------------------------------------------------------------------------


class MappingGate(Gate):
    """

    Parameters
    ==========

    args : tuple, dict

        arg[0] = initial state
        arg[1] = scalar or final state
        arg[2] = None or final state

        The list of initial state and final state pairs. The final states are to be
        multiplied by some scalar. If no scalar is given, 1 is assumed. Only need to
        supply the qubit mapping for half of the matrix representation of the gate.
        Since a quantum gate is required to be unitary, the other half is created
        to ensure it is unitary. Can pass either a python or sympy dictionary that
        already has the qubit mappings.

    Examples
    ========

    Creating a Mapping Gate and checking its arguments and properties. Getting final
    state from initial state.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy.physics.quantum.qubit import Qubit
        >>> from sympy import I
        >>> M = MappingGate(('00',I,'11'), ('01','10'), ('10','01'), ('11',-I,'00'))
        >>> M.args
        ({|00>: I*|11>, |01>: |10>, |10>: |01>, |11>: -I*|00>},)
        >>> M.nqubits
        2
        >>> M.mapping[Qubit('00')]
        I*|11>
        >>> M.get_final_state('00')
        I*|11>

    Create a Mapping Gate by only giving half of the initial and final state pairs,
    the resulting arguments are the same as the example above. Also passing a python
    or sympy dictionary to MappingGate can have the same result.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy import I
        >>> M = MappingGate(('00', I, '11'), ('01', '10'))
        >>> M.args
        ({|00>: I*|11>, |01>: |10>, |10>: |01>, |11>: -I*|00>},)
        >>> d = dict({Qubit('00'):I*Qubit('11'), Qubit('01'):Qubit('10')})
        >>> M_dict = MappingGate(d)
        >>> M.args
        ({|00>: I*|11>, |01>: |10>, |10>: |01>, |11>: -I*|00>},)

    Using qapply on initial states returns the final states.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy import I
        >>> from sympy.physics.quantum.qapply import qapply
        >>> from sympy.physics.quantum.qubit import Qubit
        >>> M = MappingGate(('00', I, '11'), ('01', '10'))
        >>> q = Qubit('00') + Qubit('01')
        >>> qapply(M*q)
        |10> + I*|11>

    The MappingGate can be rewritten as an outer product of states. We will show two
    examples: one where all four states are given and one where only one state is
    given. If not all initial states are specified they return themselves as final
    states.

        >>> from sympy.physics.quantum.mappinggate import MappingGate
        >>> from sympy import I
        >>> M = MappingGate(('00', I, '11'), ('01', '10'))
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
        >>> M = MappingGate(('00', I, '11'), ('01', '10'))
        >>> represent(M)
        [0, 0, 0, -I]
        [0, 0, 1,  0]
        [0, 1, 0,  0]
        [I, 0, 0,  0]

    """

    @classmethod
    def _eval_args(cls, args):
        if len(args) == 1 and isinstance(args[0], (dict, Dict)):
            temp = {}
            for i, f in args[0].items():
                terms = split_qexpr_parts(f)
                if len(terms[1]) == 0:
                    temp[f] = i
                else:
                    temp[terms[1][0]] = conjugate(Mul(*terms[0]))*i
                temp[i] = f
            new_args = Dict(temp)
        else:
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
                temp[f] = conjugate(scalar)*i
                temp[i] = scalar*f
            new_args = Dict(temp)
        return (new_args,)

    @classmethod
    def _eval_hilbert_space(cls, args):
        return ComplexSpace(2)**args[0].keys()[0].nqubits

    @property
    def mapping(self):
        return self.args[0]

    @property
    def nqubits(self):
        """Gives the dimension of the matrix representation"""
        return self.args[0].keys()[0].nqubits

    def get_final_state(self, qubit):
        """Returns the final state for a given initial state, if initial state is
        not mapped to a final state the initial state is returned."""
        i = Qubit(qubit)
        return self.mapping.get(i, i)

    def _apply_operator_Qubit(self, qubit):
        return self.get_final_state(qubit)

    def _eval_rewrite_as_op(self, *args):
        terms = []
        for i in range(2**self.nqubits):
            initial = Qubit(IntQubit(i, self.nqubits))
            fin = self.get_final_state(initial)
            terms.append(fin*Dagger(initial))
        return Add(*terms)

    def _represent_default_basis(self, **options):
        return self._represent_ZGate(None, **options)

    def _represent_ZGate(self, basis, **options):
        format = options.get('format','sympy')
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
            if format == 'scipy.sparse':
                matrix = matrix.tolil()
                col = float(col)
                row = float(row)
                scalar = complex(scalar)
                matrix[col, col] = 0.0
                matrix[row, col] = scalar
            elif format == 'numpy':
                scalar = complex(scalar)
                matrix[col, col] = 0.0
                matrix[row, col] = scalar
            else:
                matrix[col, col] = Integer(0)
                matrix[row, col] = scalar
        if format == 'scipy.sparse':
            matrix = matrix.tocsr()
        return matrix

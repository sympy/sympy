"""Grover's algorithm and helper functions.

Todo:
* W gate construction (or perhaps -W gate based on Mermin's book)
* Perform the Grover iteration (application of WV on a state) many times
  in order to achieve the desired outcome with a high probability.
  Look back at the description to understand what is meant by this.
* Generalize the algorithm for an unknown function that returns 1 on
  multiple qubit states, not just one.  
"""

from sympy import sqrt  
from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import ComplexSpace
from sympy.physics.quantum.operator import UnitaryOperator
from sympy.physics.quantum.gate import Gate, HadamardGate
from sympy.physics.quantum.qubit import Qubit

__all__ = [
    'OracleGate',
    'WGate'
]

"""Provides a string representation of 0 of some length
Parameters: 
numzeroes - int
    The number of zeroes in the string

Return: string
    Returns a string of numzeroes zeroes
"""
def _create_zeroes(numzeroes):
    zeroes = ['0' for i in range(numzeroes)]
    return reduce((lambda x, y: x + y), zeroes, '')

"""Creates a superposition of basis states for n-qubits
Parameters:
-----------
nqubits - int
    The number of qubits for the state

Return: qexpr
    An equal superposition of nqubit basis states
"""
def create_basis_states(nqubits):
    basis = [i for i in range(2**nqubits)]
    # Return the string starting from index 2 to remove the substring
    # '0b' from the binary string.  ex. ('0b1111')[2:] -> '1111'
    bin_convert = (lambda num: 
                      _create_zeroes(nqubits-len(bin(num)[2:])) + bin(num)[2:]
                  )
    qubit_states = (lambda bin_rep:
                       (1/sqrt(2**nqubits))*Qubit(bin_rep)
                   )
    basis_st = map(qubit_states, map(bin_convert, basis))
    return reduce(lambda cur_st, st: cur_st + st, basis_st[1:], basis_st[0])

class OracleGate(Gate):
    """A black box gate

    The gate flips the sign of the basis state if the unknown function
    applied to the basis state is 1.

    Parameters
    ----------
    label : tuple (int, callable)
        Number of qubits
        A callable function that returns 1 or 0 for a basis state.

    Examples
    --------
    """

    gate_name = u'V'
    gate_name_latex = u'V'
    
    #-------------------------------------------------------------------------
    # Initialization/creation
    #-------------------------------------------------------------------------

    # args : tuple
    # Returns: tuple, (int tuple - target qubits, callable)
    @classmethod
    def _eval_args(cls, args):
        if len(args) != 2:
            raise QuantumError(
                'Insufficient/excessive arguments to Oracle.  Please ' +
                    'supply the number of qubits and an unknown function.'
            )
        numqubits = args[0], 
        numqubits = UnitaryOperator._eval_args(numqubits)
        if not numqubits[0].is_Integer:
           raise TypeError('Integer expected, got: %r' % numqubits[0])
        if not callable(args[1]):
           raise TypeError('Callable expected, got: %r' % args[1])
        numqubits = UnitaryOperator._eval_args(tuple(range(args[0])))
        # Ask about because returning (int tuple, callable) b/c args is
        # used to create Hilbert space
        return (numqubits, args[1])

    @classmethod
    def _eval_hilbert_space(cls, args):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(args[0])+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def search_function(self):
        """The unknown function that helps find the sought after state"""
        return self.label[1]

    @property
    def targets(self):
        """A tuple of target qubits."""
        return self.label[0]

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    # qubits : a set of qubits (Qubit)
    # Return: quantum object (quantum expression - QExpr)
    def _apply_operator_Qubit(self, qubits, **options):
        #print 'Oracle _apply_operator_Qubit'
        if qubits.nqubits != self.nqubits:
            raise QuantumError(
                'OracleGate operates on %r qubits, got: %r' 
                    (self.nqubits, qubits.nqubits)
            )
        # If function returns 1 on qubits
            # return the negative of the qubits (flip the sign)
        return -qubits if (self.search_function)(qubits) == 1 else qubits

    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    # Still To Do
    def _represent_ZGate(self, basis, **options):
        return basis
 

class WGate(Gate):
    """General n qubit W Gate in Grover's algorithm.

    The gate performs the operation 2|phi><phi| - 1 on a basis state.
    |phi> = (tensor product of n Hadamards)*(|0> with n qubits)

    Parameters
    ----------
    nqubits : int
        The number of qubits to operate on
    """

    gate_name = u'W'
    gate_name_latex = u'W'

    @classmethod
    def _eval_args(cls, args):
        if len(args) != 1:
            raise QuantumError(
                'Insufficient/excessive arguments to W gate.  Please ' +
                    'supply the number of qubits to operate on.'
            )
        args = UnitaryOperator._eval_args(args)
        if not args[0].is_Integer:
           raise TypeError('Integer expected, got: %r' % args[0])
        return tuple([args[0]-1-i for i in range(args[0])])

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    # qubits : a set of qubits (Qubit)
    # Return: quantum object (quantum expression - QExpr)
    def _apply_operator_Qubit(self, qubits, **options):
        #print 'WGate _apply_operator_Qubit'
        if qubits.nqubits != self.nqubits:
            raise QuantumError(
                'WGate operates on %r qubits, got: %r' 
                    (self.nqubits, qubits.nqubits)
            )

        # See 'Quantum Computer Science' by David Mermin p.92 -> W|a> result
        # Return (2/(sqrt(2^n)))|phi> - |a> where |a> is the current basis
        # state and phi is the superposition of basis states (see function
        # create_basis_states above)
        basis_states = create_basis_states(self.nqubits)
        change_to_basis = (2/sqrt(2**self.nqubits))*basis_states
        return change_to_basis - qubits

"""Applies one application of the Oracle and W Gate, WV
Parameters:
qstate : QExpr (quantum object)
    A quantum state
oracle : OracleGate
    The black box operator that flips the sign of the desired basis state

Returns:
QExpr : The quantum state after applying the oracle and W gate.
"""
def grover_iteration(qstate, oracle):
    wgate = WGate(oracle.nqubits)
    return wgate*oracle*qstate

"""Applies grover's algorithm
Parameters:
bbox : callable 
    The unknown callable function that returns 1 applied to the 
    desired basis state and 0 otherwise.  

Returns:
Qubit : The sought after qubits for the unknown callable function
"""
def grover(bbox):
    # Ask about the number of iterations and measurement
    return 0

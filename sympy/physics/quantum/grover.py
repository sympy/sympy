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
from sympy.physics.quantum.qubit import IntQubit

__all__ = [
    'OracleGate',
    'WGate'
]

def _create_computational_basis(nqubits):
    """Creates an equal superposition of the computational basis.

    Parameters:
    -----------
    nqubits : int
        The number of qubits

    Return
    ------
    Qubit : An equal superposition of the computational basis with nqubits
    """
    amp = 1/sqrt(2**nqubits)
    return sum(amp*IntQubit(n, nqubits) for n in range(2**nqubits))

class OracleGate(Gate):
    """A black box gate

    The gate marks the desired qubits of an unknown function by flipping
    the sign of the qubits.  The unknown function returns true when it
    finds its desired qubits and false otherwise.  

    Parameters
    ----------
    label : tuple (int, callable)
        Number of qubits
        A callable function that returns a boolean on a computational basis

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
        sub_args = args[0], 
        sub_args = UnitaryOperator._eval_args(sub_args)
        if not sub_args[0].is_Integer:
           raise TypeError('Integer expected, got: %r' % sub_args[0])
        if not callable(args[1]):
           raise TypeError('Callable expected, got: %r' % args[1])
        sub_args = UnitaryOperator._eval_args(tuple(range(args[0])))
        return (sub_args, args[1])

    @classmethod
    def _eval_hilbert_space(cls, args):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**(max(args[0])+1)

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def search_function(self):
        """The unknown function that helps find the sought after qubits"""
        return self.label[1]

    @property
    def targets(self):
        """A tuple of target qubits."""
        return self.label[0]

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    # qubits : a set of qubits (Qubit)
    # Return: Qubit
    def _apply_operator_Qubit(self, qubits, **options):
        if qubits.nqubits != self.nqubits:
            raise QuantumError(
                'OracleGate operates on %r qubits, got: %r' 
                    (self.nqubits, qubits.nqubits)
            )
        # If function returns 1 on qubits
            # return the negative of the qubits (flip the sign)
        return -qubits if (self.search_function)(qubits) else qubits

    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    # Still To Do
    def _represent_ZGate(self, basis, **options):
        raise NotImplementedError(
            "Represent for the Oracle has not been implemented yet"
        )
 

class WGate(Gate):
    """General n qubit W Gate in Grover's algorithm.

    The gate performs the operation 2|phi><phi| - 1 on some qubits.
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
        return tuple(reversed(range(args[0])))

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    # qubits : a set of qubits (Qubit)
    # Return: quantum object (quantum expression - QExpr)
    def _apply_operator_Qubit(self, qubits, **options):
        if qubits.nqubits != self.nqubits:
            raise QuantumError(
                'WGate operates on %r qubits, got: %r' 
                    (self.nqubits, qubits.nqubits)
            )

        # See 'Quantum Computer Science' by David Mermin p.92 -> W|a> result
        # Return (2/(sqrt(2^n)))|phi> - |a> where |a> is the current basis
        # state and phi is the superposition of basis states (see function
        # create_computational_basis above)
        basis_states = _create_computational_basis(self.nqubits)
        change_to_basis = (2/sqrt(2**self.nqubits))*basis_states
        return change_to_basis - qubits

def grover_iteration(qstate, oracle):
    """Applies one application of the Oracle and W Gate, WV

    Parameters
    ----------
    qstate : Qubit
        A superposition of qubits
    oracle : OracleGate
        The black box operator that flips the sign of the desired basis qubits

    Returns
    -------
    Qubit : The qubits after applying the Oracle and W gate.

    """
    wgate = WGate(oracle.nqubits)
    return wgate*oracle*qstate

def grover(bbox):
    """Applies grover's algorithm

    Parameters
    ----------
    bbox : callable 
        The unknown callable function that returns true when applied to the 
        desired qubits and false otherwise.  

    Returns
    -------
    Qubit : The sought after qubits for the unknown callable function

    """
    # Ask about the number of iterations and measurement
    return 0

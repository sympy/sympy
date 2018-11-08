"""Grover's algorithm and helper functions.

Todo:

* W gate construction (or perhaps -W gate based on Mermin's book)
* Implement _represent_ZGate in OracleGate
"""

from __future__ import print_function, division

from sympy import floor, pi, sqrt, sympify, eye
from sympy.core.compatibility import range
from sympy.core.numbers import NegativeOne
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.qexpr import QuantumError
from sympy.physics.quantum.hilbert import ComplexSpace
from sympy.physics.quantum.operator import UnitaryOperator
from sympy.physics.quantum.gate import Gate
from sympy.physics.quantum.qubit import IntQubit, Qubit

__all__ = [
    'OracleGate',
    'WGate',
    'superposition_basis',
    'grover_iteration',
    'apply_grover',
    'random_oracle',
    'g_bbht_search'
]


def superposition_basis(nqubits):
    """Creates an equal superposition of the computational basis.

    Parameters
    ==========

    nqubits : int
        The number of qubits.

    Returns
    =======

    state : Qubit
        An equal superposition of the computational basis with nqubits.

    Examples
    ========

    Create an equal superposition of 2 qubits::

        >>> from sympy.physics.quantum.grover import superposition_basis
        >>> superposition_basis(2)
        |0>/2 + |1>/2 + |2>/2 + |3>/2
    """

    amp = 1/sqrt(2**nqubits)
    return sum([amp*IntQubit(n, nqubits) for n in range(2**nqubits)])


class OracleGate(Gate):
    """A black box gate.

    The gate marks the desired qubits of an unknown function by flipping
    the sign of the qubits.  The unknown function returns true when it
    finds its desired qubits and false otherwise.

    Parameters
    ==========

    qubits : int
        Number of qubits.

    oracle : callable
        A callable function that returns a boolean on a computational basis.

    Examples
    ========

    Apply an Oracle gate that flips the sign of ``|2>`` on different qubits::

        >>> from sympy.physics.quantum.qubit import IntQubit
        >>> from sympy.physics.quantum.qapply import qapply
        >>> from sympy.physics.quantum.grover import OracleGate
        >>> f = lambda qubits: qubits == IntQubit(2)
        >>> v = OracleGate(2, f)
        >>> qapply(v*IntQubit(2))
        -|2>
        >>> qapply(v*IntQubit(3))
        |3>
    """

    gate_name = u'V'
    gate_name_latex = u'V'

    #-------------------------------------------------------------------------
    # Initialization/creation
    #-------------------------------------------------------------------------

    @classmethod
    def _eval_args(cls, args):
        # TODO: args[1] is not a subclass of Basic
        if len(args) != 2:
            raise QuantumError(
                'Insufficient/excessive arguments to Oracle.  Please ' +
                'supply the number of qubits and an unknown function.'
            )
        sub_args = (args[0],)
        sub_args = UnitaryOperator._eval_args(sub_args)
        if not sub_args[0].is_Integer:
            raise TypeError('Integer expected, got: %r' % sub_args[0])

        if not callable(args[1]):
            raise TypeError('Callable expected, got: %r' % args[1])
        return (sub_args[0], args[1])

    @classmethod
    def _eval_hilbert_space(cls, args):
        """This returns the smallest possible Hilbert space."""
        return ComplexSpace(2)**args[0]

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def search_function(self):
        """The unknown function that helps find the sought after qubits."""
        return self.label[1]

    @property
    def targets(self):
        """A tuple of target qubits."""
        return sympify(tuple(range(self.args[0])))

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    def _apply_operator_Qubit(self, qubits, **options):
        """Apply this operator to a Qubit subclass.

        Parameters
        ==========

        qubits : Qubit
            The qubit subclass to apply this operator to.

        Returns
        =======

        state : Expr
            The resulting quantum state.
        """
        if qubits.nqubits != self.nqubits:
            raise QuantumError(
                'OracleGate operates on %r qubits, got: %r'
                % (self.nqubits, qubits.nqubits)
            )
        # If function returns 1 on qubits
            # return the negative of the qubits (flip the sign)
        if self.search_function(qubits):
            return -qubits
        else:
            return qubits

    #-------------------------------------------------------------------------
    # Represent
    #-------------------------------------------------------------------------

    def _represent_ZGate(self, basis, **options):
        """
        Represent the OracleGate in the computational basis.
        """
        nbasis = 2**self.nqubits  # compute it only once
        matrixOracle = eye(nbasis)
        # Flip the sign given the output of the oracle function
        for i in range(nbasis):
            if self.search_function(IntQubit(i, self.nqubits)):
                matrixOracle[i, i] = NegativeOne()
        return matrixOracle


class WGate(Gate):
    """General n qubit W Gate in Grover's algorithm.

    The gate performs the operation ``2|phi><phi| - 1`` on some qubits.
    ``|phi> = (tensor product of n Hadamards)*(|0> with n qubits)``

    Parameters
    ==========

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
        return args

    #-------------------------------------------------------------------------
    # Properties
    #-------------------------------------------------------------------------

    @property
    def targets(self):
        return sympify(tuple(reversed(range(self.args[0]))))

    #-------------------------------------------------------------------------
    # Apply
    #-------------------------------------------------------------------------

    def _apply_operator_Qubit(self, qubits, **options):
        """
        qubits: a set of qubits (Qubit)
        Returns: quantum object (quantum expression - QExpr)
        """
        if qubits.nqubits != self.nqubits:
            raise QuantumError(
                'WGate operates on %r qubits, got: %r'
                % (self.nqubits, qubits.nqubits)
            )

        # See 'Quantum Computer Science' by David Mermin p.92 -> W|a> result
        # Return (2/(sqrt(2^n)))|phi> - |a> where |a> is the current basis
        # state and phi is the superposition of basis states (see function
        # create_computational_basis above)
        basis_states = superposition_basis(self.nqubits)
        change_to_basis = (2/sqrt(2**self.nqubits))*basis_states
        return change_to_basis - qubits


def grover_iteration(qstate, oracle):
    """Applies one application of the Oracle and W Gate, WV.

    Parameters
    ==========

    qstate : Qubit
        A superposition of qubits.
    oracle : OracleGate
        The black box operator that flips the sign of the desired basis qubits.

    Returns
    =======

    Qubit : The qubits after applying the Oracle and W gate.

    Examples
    ========

    Perform one iteration of grover's algorithm to see a phase change::

        >>> from sympy.physics.quantum.qapply import qapply
        >>> from sympy.physics.quantum.qubit import IntQubit
        >>> from sympy.physics.quantum.grover import OracleGate
        >>> from sympy.physics.quantum.grover import superposition_basis
        >>> from sympy.physics.quantum.grover import grover_iteration
        >>> numqubits = 2
        >>> basis_states = superposition_basis(numqubits)
        >>> f = lambda qubits: qubits == IntQubit(2)
        >>> v = OracleGate(numqubits, f)
        >>> qapply(grover_iteration(basis_states, v))
        |2>

    """
    wgate = WGate(oracle.nqubits)
    return wgate*oracle*qstate


def apply_grover(oracle, nqubits, iterations=None):
    """Applies grover's algorithm.

    Parameters
    ==========

    oracle : callable
        The unknown callable function that returns true when applied to the
        desired qubits and false otherwise.

    Returns
    =======

    state : Expr
        The resulting state after Grover's algorithm has been iterated.

    Examples
    ========

    Apply grover's algorithm to an even superposition of 2 qubits::

        >>> from sympy.physics.quantum.qapply import qapply
        >>> from sympy.physics.quantum.qubit import IntQubit
        >>> from sympy.physics.quantum.grover import apply_grover
        >>> f = lambda qubits: qubits == IntQubit(2)
        >>> qapply(apply_grover(f, 2))
        |2>

    """
    if nqubits <= 0:
        raise QuantumError(
            'Grover\'s algorithm needs nqubits > 0, received %r qubits'
            % nqubits
        )
    if iterations is None:
        iterations = floor(sqrt(2**nqubits)*(pi/4))

    v = OracleGate(nqubits, oracle)
    iterated = superposition_basis(nqubits)
    for iter in range(iterations):
        iterated = grover_iteration(iterated, v)
        iterated = qapply(iterated)

    return iterated


def random_oracle(nqubits, min_img=1, max_img=1, q_type='bin'):
    """Create a random OracleGate under the given parameter

    Parameters
    ==========

    nqubits : int
        The number of qubits for OracleGate
    min_pic : int
        Minimum number of inverse images that are mapped to 1
    max_pic : int
        Maximum number of inverse images that are mapped to 1
    q_type : OracleGate
        Type of the Qubits that the oracle should be applied on.
        Can be 'bin' for binary (Qubit()) or 'int' for integer (IntQubit()).

    Returns
    =======

    OracleGate : random OracleGate under the given parameter

    Examples
    ========

    Generate random OracleGate that outputs 1 for 2-4 inputs::

        >>> from sympy.physics.quantum.grover import random_oracle
        >>> oracle = random_oracle(4, min_img=2, max_img=4, q_type="bin")

    """
    if q_type != 'bin' and q_type != 'int':
        raise QuantumError("q_type must be 'int' or 'bin'")

    if min_img < 1 or max_img < 1:
        raise QuantumError("min_pic, max_pic must be > 0")

    if min_img > max_img:
        raise QuantumError("max_pic must be >= min_pic")

    if min_img >= 2 ** nqubits or max_img > 2 ** nqubits:
        raise QuantumError("min_pic must be < 2**nqubits and max_pic must be <= 2**nqubits")

    import random
    pics = random.randint(min_img, max_img)
    integers = random.sample(range(2 ** nqubits), pics)
    if q_type == "int":
        items = [IntQubit(i) for i in integers]
    else:
        items = [Qubit(IntQubit(i)) for i in integers]

    return OracleGate(nqubits, lambda qubits: qubits in items)


def g_bbht_search(qstate, oracle):
    """G-BBHT-Search

    G-BBHT-Search is an algorithm based on Grover's algorithm,
    but designed for an unknown oracle function that returns 1
    for an unknown number of qubit states. It is based on the
    paper: Boyer, M., Brassard, G., Hoyer, P.,  Tapp, A. (1998).
    Tight bounds on quantum searching: Progress of Physics,
    46(4-5), 493-505.

    Parameters
    ==========

    qstate : Qubit
        State the G-BBHT-Search should be applied on
    oracle : OracleGate
        The black box operator that flips the sign of the desired basis qubits.

    Returns
    =======

    (Qubit, int) : (single Qubit that fulfills oracle, number of grover iterations applied)

    Examples
    ========

    G-BBHT-Search for an orcale that returns 1 for an unkown number of statesbetween 5 - 15::

        >>> from sympy.physics.quantum.grover import random_oracle, superposition_basis, g_bbht_search
        >>> basis_states = superposition_basis(4)
        >>> oracle = random_oracle(4, min_img=5, max_img=15, q_type="bin")
        >>> x = g_bbht_search(basis_states, oracle)[0]
        >>> oracle.search_function(x)
        True

    """
    import random
    from sympy.physics.quantum.qubit import measure_all_oneshot

    max_iterations = 1
    factor_iterations = 6.0 / 5.0
    count_grover_iterations = 0

    while True:
        for _ in range(1 if max_iterations == 1 else random.randint(0, int(max_iterations))):
            count_grover_iterations += 1
            qstate = qapply(grover_iteration(qstate, oracle))

        measure = measure_all_oneshot(qstate)
        max_iterations *= factor_iterations
        if oracle.search_function(measure) is True:
            return measure, count_grover_iterations

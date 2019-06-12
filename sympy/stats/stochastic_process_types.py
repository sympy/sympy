from sympy import (Symbol, Matrix, MatrixSymbol, S, Indexed, Basic,
                    Set, And, Tuple, Eq, FiniteSet, ImmutableMatrix,
                    nsimplify)
from sympy.stats.rv import (RandomIndexedSymbol, random_symbols, RandomSymbol,
                            _symbol_converter)
from sympy.core.compatibility import string_types
from sympy.core.relational import Relational
from sympy.stats.symbolic_probability import Probability
from sympy.stats.stochastic_process import StochasticPSpace
from sympy.logic.boolalg import Boolean

__all__ = [
    'StochasticProcess',
    'DiscreteTimeStochasticProcess',
    'DiscreteMarkovChain',
    'TransitionMatrixOf',
    'StochasticStateSpaceOf'
]

def _set_converter(itr):
    """
    Helper function for converting list/tuple/set to Set.
    If parameter is not an instance of list/tuple/set then
    no operation is performed.

    Returns
    =======

    Set
        The argument converted to Set.


    Raises
    ======

    TypeError
        If the argument is not an instance of list/tuple/set.
    """
    if isinstance(itr, (list, tuple, set)):
        itr = FiniteSet(*itr)
    if not isinstance(itr, Set):
        raise TypeError("%s is not an instance of list/tuple/set."%(itr))
    return itr

def _matrix_checks(matrix):
    if not isinstance(matrix, (Matrix, MatrixSymbol, ImmutableMatrix)):
        raise TypeError("Transition probabilities etiher should "
                            "be a Matrix or a MatrixSymbol.")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("%s is not a square matrix"%(matrix))
    if isinstance(matrix, Matrix):
        matrix = ImmutableMatrix(matrix.tolist())
    return matrix

class StochasticProcess(Basic):
    """
    Base class for all the stochastic processes whether
    discrete or continuous.

    Parameters
    ==========

    sym: Symbol or string_types
    state_space: Set
        The state space of the stochastic process, by default S.Reals.
        For discrete sets it is zero indexed.

    See Also
    ========

    DiscreteTimeStochasticProcess
    """

    def __new__(cls, sym, state_space=S.Reals):
        sym = _symbol_converter(sym)
        state_space = _set_converter(state_space)
        return Basic.__new__(cls, sym, state_space)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def state_space(self):
        return self.args[1]

    def __call__(self, time):
        """
        Overrided in ContinuousTimeStochasticProcess.
        """
        raise NotImplementedError("Use [] for indexing discrete time stochastic process.")

    def __getitem__(self, time):
        """
        Overrided in DiscreteTimeStochasticProcess.
        """
        raise NotImplementedError("Use () for indexing continuous time stochastic process.")

    def probability(self, condition):
        raise NotImplementedError()

class DiscreteTimeStochasticProcess(StochasticProcess):
    """
    Base class for all discrete stochastic processes.
    """
    def __getitem__(self, time):
        """
        For indexing discrete time stochastic processes.

        Returns
        =======

        RandomIndexedSymbol
        """
        if time not in self.index_set:
            raise IndexError("%s is not in the index set of %s"%(time, self.symbol))
        idx_obj = Indexed(self.symbol, time)
        pspace_obj = StochasticPSpace(self.symbol, self)
        return RandomIndexedSymbol(idx_obj, pspace_obj)

class TransitionMatrixOf(Boolean):
    """
    Assumes that the matrix is the transition matrix
    of the process.
    """

    def __new__(cls, process, matrix):
        if not isinstance(process, DiscreteMarkovChain):
            raise ValueError("Currently only DiscreteMarkovChain "
                                "support TransitionMatrixOf.")
        matrix = _matrix_checks(matrix)
        return Basic.__new__(cls, process, matrix)

    process = property(lambda self: self.args[0])
    matrix = property(lambda self: self.args[1])

class StochasticStateSpaceOf(Boolean):

    def __new__(cls, process, state_space):
        if not isinstance(process, DiscreteMarkovChain):
            raise ValueError("Currently only DiscreteMarkovChain "
                                "support StochasticStateSpaceOf.")
        state_space = _set_converter(state_space)
        return Basic.__new__(cls, process, state_space)

    process = property(lambda self: self.args[0])
    state_space = property(lambda self: self.args[1])

class DiscreteMarkovChain(DiscreteTimeStochasticProcess):
    """
    Represents discrete Markov chain.

    Parameters
    ==========

    sym: Symbol
    state_space: Set
        Optional, by default, S.Reals
    trans_probs: Matrix/ImmutableMatrix/MatrixSymbol
        Optional, by default, None

    Examples
    ========

    >>> from sympy.stats import DiscreteMarkovChain, TransitionMatrixOf
    >>> from sympy import Matrix, MatrixSymbol, Eq
    >>> from sympy.stats import P
    >>> T = Matrix([[0.5, 0.2, 0.3],[0.2, 0.5, 0.3],[0.2, 0.3, 0.5]])
    >>> Y = DiscreteMarkovChain("Y", [0, 1, 2], T)
    >>> YS = DiscreteMarkovChain("Y")
    >>> Y.state_space
    {0, 1, 2}
    >>> Y.transition_probabilities
    Matrix([
    [0.5, 0.2, 0.3],
    [0.2, 0.5, 0.3],
    [0.2, 0.3, 0.5]])
    >>> TS = MatrixSymbol('T', 3, 3)
    >>> P(Eq(YS[3], 2), Eq(YS[1], 1) & TransitionMatrixOf(YS, TS))
    T[0, 2]*T[1, 0] + T[1, 1]*T[1, 2] + T[1, 2]*T[2, 2]
    >>> P(Eq(Y[3], 2), Eq(Y[1], 1)).round(2)
    0.36
    """

    index_set = S.Naturals0

    def __new__(cls, sym, state_space=S.Reals, trans_probs=None):
        sym = _symbol_converter(sym)
        state_space = _set_converter(state_space)
        if trans_probs != None:
            trans_probs = _matrix_checks(trans_probs)
        return Basic.__new__(cls, sym, state_space, trans_probs)

    @property
    def transition_probabilities(self):
        """
        Transition probabilities of discrete Markov chain,
        either an instance of Matrix or MatrixSymbol.
        """
        return self.args[2]

    def probability(self, condition, given_condition, **kwargs):
        """
        Handles probability queries for discrete Markov chains.

        Parameters
        ==========

        condition: Relational
        given_condition: Relational/And

        Returns
        =======

        Probability
            If the transition probabilities are not available
        Expr
            If the transition probabilities is MatrixSymbol or Matrix

        Note
        ====

        Any information passed at the time of query overrides
        any information passed at the time of object creation like
        transition probabilities, state space.

        Pass the transition matrix using TransitionMatrixOf and state space
        using StochasticStateSpaceOf in given_condition using & or And.
        """

        # extracting transition matrix and state space
        trans_probs, state_space = self.transition_probabilities, self.state_space
        if isinstance(given_condition, And):
            gcs = given_condition.args
            for gc in gcs:
                if isinstance(gc, TransitionMatrixOf):
                    trans_probs = gc.matrix
                if isinstance(gc, StochasticStateSpaceOf):
                    state_space = gc.state_space
                if isinstance(gc, Eq):
                    given_condition = gc
        if isinstance(given_condition, TransitionMatrixOf):
            trans_probs = given_condition.matrix
        if isinstance(given_condition, StochasticStateSpaceOf):
            state_space = given_condition.state_space

        # given_condition does not have sufficient information
        # for computations
        if trans_probs == None or \
            given_condition == None:
            return Probability(condition, given_condition, **kwargs)

        # working out transition probabilities
        if not isinstance(trans_probs, MatrixSymbol):
            rows = trans_probs.tolist()
            for row in rows:
                if (sum(row) - 1) != 0:
                    raise ValueError("Probabilities in a row must sum to 1. "
                    "If you are using Float or floats then please use Rational.")

        # if given condition is None, then there is no need to work out
        # state_space from random variables
        if given_condition != None:
            rand_var = list(given_condition.atoms(RandomSymbol) -
                        given_condition.atoms(RandomIndexedSymbol))
            if len(rand_var) == 1:
                state_space = rand_var[0].pspace.set
        if not FiniteSet(*[i for i in range(trans_probs.shape[0])]).is_subset(state_space):
            raise ValueError("state space is not compatible with the transition probabilites.")

        if isinstance(condition, Eq) and \
            isinstance(given_condition, Eq) and \
            len(given_condition.atoms(RandomSymbol)) == 1:
            # handles simple queries like P(Eq(X[i], dest_state), Eq(X[i], init_state))
            lhsc, rhsc = condition.lhs, condition.rhs
            lhsg, rhsg = given_condition.lhs, given_condition.rhs
            if not isinstance(lhsc, RandomIndexedSymbol):
                lhsc, rhsc = (rhsc, lhsc)
            if not isinstance(lhsg, RandomIndexedSymbol):
                lhsg, rhsg = (rhsg, lhsg)
            keyc, statec, keyg, stateg = (lhsc.key, rhsc, lhsg.key, rhsg)
            if stateg >= trans_probs.shape[0] == False or statec >= trans_probs.shape[1]:
                raise IndexError("No information is avaliable for (%s, %s) in "
                    "transition probabilities of shape, (%s, %s). "
                    "State space is zero indexed."
                    %(stateg, statec, trans_probs.shape[0], trans_probs.shape[1]))
            if keyc < keyg:
                raise ValueError("Incorrect given condition is given, probability "
                      "of past state cannot be computed from future state.")
            nsteptp = trans_probs**(keyc - keyg)
            if hasattr(nsteptp, "__getitem__"):
                return nsteptp.__getitem__((stateg, statec))
            return Indexed(nsteptp, stateg, statec)

        if isinstance(condition, And):
            # handle queries like,
            # P(Eq(X[i+k], s1) & Eq(X[i+m], s2) . . . & Eq(X[i], sn), Eq(P(X[i]), prob))
            conds = condition.args
            i, result = -1, 1
            while i > -len(conds):
                result *= self.probability(conds[i], conds[i-1] & \
                            TransitionMatrixOf(self, trans_probs) & \
                            StochasticStateSpaceOf(self, state_space))
                i -= 1
            if isinstance(given_condition, (TransitionMatrixOf, StochasticStateSpaceOf)):
                return result * Probability(conds[i])
            if isinstance(given_condition, Eq):
                if not isinstance(given_condition.lhs, Probability) or \
                    given_condition.lhs.args[0] != conds[i]:
                    raise ValueError("Probability for %s needed", conds[i])
                return result * given_condition.rhs

        raise NotImplementedError("Mechanism for handling (%s, %s) queries hasn't been "
                                "implemented yet."%(condition, given_condition))

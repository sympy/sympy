from __future__ import print_function, division
import random

import itertools

from sympy import (Matrix, MatrixSymbol, S, Indexed, Basic,
                   Set, And, Eq, FiniteSet, ImmutableMatrix,
                   Lambda, Mul, Dummy, IndexedBase, Add, Interval, oo,
                   linsolve, Or, Not, Intersection, factorial, Contains,
                   Union, Expr, Function, exp, cacheit, sqrt, pi, gamma,
                   Ge, Piecewise, Symbol, NonSquareMatrixError, EmptySet,
                   Range, limit_seq, zeros, ones, Identity, FunctionMatrix)
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean
from sympy.stats.joint_rv import JointDistribution
from sympy.stats.joint_rv_types import JointDistributionHandmade
from sympy.stats.rv import (RandomIndexedSymbol, random_symbols, RandomSymbol,
                            _symbol_converter, _value_check, pspace, given,
                           dependent, is_random, sample_iter)
from sympy.stats.stochastic_process import StochasticPSpace
from sympy.stats.symbolic_probability import Probability, Expectation
from sympy.stats.frv_types import Bernoulli, BernoulliDistribution, FiniteRV
from sympy.stats.drv_types import Poisson, PoissonDistribution
from sympy.stats.crv_types import Normal, NormalDistribution, Gamma, GammaDistribution
from sympy.core.sympify import _sympify

__all__ = [
    'StochasticProcess',
    'DiscreteTimeStochasticProcess',
    'DiscreteMarkovChain',
    'TransitionMatrixOf',
    'StochasticStateSpaceOf',
    'GeneratorMatrixOf',
    'ContinuousMarkovChain',
    'BernoulliProcess',
    'PoissonProcess',
    'WienerProcess',
    'GammaProcess'
]


@is_random.register(Indexed)
def _(x):
    return is_random(x.base)

@is_random.register(RandomIndexedSymbol)
def _(x):
    return True

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

def _sym_sympify(arg):
    """
    Converts an arbitrary expression to a type that can be used inside SymPy.
    As generally strings are unwise to use in the expressions,
    it returns the Symbol of argument if the string type argument is passed.

    Parameters
    =========

    arg: The parameter to be converted to be used in Sympy.

    Returns
    =======

    The converted parameter.

    """
    if isinstance(arg, str):
        return Symbol(arg)
    else:
        return _sympify(arg)

def _matrix_checks(matrix):
    if not isinstance(matrix, (Matrix, MatrixSymbol, ImmutableMatrix, FunctionMatrix)):
        raise TypeError("Transition probabilities either should "
                            "be a Matrix or a MatrixSymbol.")
    if matrix.shape[0] != matrix.shape[1]:
        raise NonSquareMatrixError("%s is not a square matrix"%(matrix))
    if isinstance(matrix, Matrix):
        matrix = ImmutableMatrix(matrix.tolist())
    return matrix

class StochasticProcess(Basic):
    """
    Base class for all the stochastic processes whether
    discrete or continuous.

    Parameters
    ==========

    sym: Symbol or str
    state_space: Set
        The state space of the stochastic process, by default S.Reals.
        For discrete sets it is zero indexed.

    See Also
    ========

    DiscreteTimeStochasticProcess
    """

    index_set = S.Reals

    def __new__(cls, sym, state_space=S.Reals, **kwargs):
        sym = _symbol_converter(sym)
        state_space = _set_converter(state_space)
        return Basic.__new__(cls, sym, state_space)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def state_space(self):
        return self.args[1]

    @property
    def distribution(self):
        return None

    def __call__(self, time):
        """
        Overridden in ContinuousTimeStochasticProcess.
        """
        raise NotImplementedError("Use [] for indexing discrete time stochastic process.")

    def __getitem__(self, time):
        """
        Overridden in DiscreteTimeStochasticProcess.
        """
        raise NotImplementedError("Use () for indexing continuous time stochastic process.")

    def probability(self, condition):
        raise NotImplementedError()

    def joint_distribution(self, *args):
        """
        Computes the joint distribution of the random indexed variables.

        Parameters
        ==========

        args: iterable
            The finite list of random indexed variables/the key of a stochastic
            process whose joint distribution has to be computed.

        Returns
        =======

        JointDistribution
            The joint distribution of the list of random indexed variables.
            An unevaluated object is returned if it is not possible to
            compute the joint distribution.

        Raises
        ======

        ValueError: When the arguments passed are not of type RandomIndexSymbol
        or Number.
        """
        args = list(args)
        for i, arg in enumerate(args):
            if S(arg).is_Number:
                if self.index_set.is_subset(S.Integers):
                    args[i] = self.__getitem__(arg)
                else:
                    args[i] = self.__call__(arg)
            elif not isinstance(arg, RandomIndexedSymbol):
                raise ValueError("Expected a RandomIndexedSymbol or "
                                "key not  %s"%(type(arg)))

        if args[0].pspace.distribution == None: # checks if there is any distribution available
            return JointDistribution(*args)

        pdf = Lambda(tuple(args),
                expr=Mul.fromiter(arg.pspace.process.density(arg) for arg in args))
        return JointDistributionHandmade(pdf)

    def expectation(self, condition, given_condition):
        raise NotImplementedError("Abstract method for expectation queries.")

    def sample(self):
        raise NotImplementedError("Abstract method for sampling queries.")

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
        pspace_obj = StochasticPSpace(self.symbol, self, self.distribution)
        return RandomIndexedSymbol(idx_obj, pspace_obj)

class ContinuousTimeStochasticProcess(StochasticProcess):
    """
    Base class for all continuous time stochastic process.
    """
    def __call__(self, time):
        """
        For indexing continuous time stochastic processes.

        Returns
        =======

        RandomIndexedSymbol
        """
        if time not in self.index_set:
            raise IndexError("%s is not in the index set of %s"%(time, self.symbol))
        func_obj = Function(self.symbol)(time)
        pspace_obj = StochasticPSpace(self.symbol, self, self.distribution)
        return RandomIndexedSymbol(func_obj, pspace_obj)

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

class GeneratorMatrixOf(TransitionMatrixOf):
    """
    Assumes that the matrix is the generator matrix
    of the process.
    """

    def __new__(cls, process, matrix):
        if not isinstance(process, ContinuousMarkovChain):
            raise ValueError("Currently only ContinuousMarkovChain "
                                "support GeneratorMatrixOf.")
        matrix = _matrix_checks(matrix)
        return Basic.__new__(cls, process, matrix)

class StochasticStateSpaceOf(Boolean):

    def __new__(cls, process, state_space):
        if not isinstance(process, (DiscreteMarkovChain, ContinuousMarkovChain)):
            raise ValueError("Currently only DiscreteMarkovChain and ContinuousMarkovChain "
                                "support StochasticStateSpaceOf.")
        state_space = _set_converter(state_space)
        return Basic.__new__(cls, process, state_space)

    process = property(lambda self: self.args[0])
    state_space = property(lambda self: self.args[1])

class MarkovProcess(StochasticProcess):
    """
    Contains methods that handle queries
    common to Markov processes.
    """
    def _extract_information(self, given_condition):
        """
        Helper function to extract information, like,
        transition matrix/generator matrix, state space, etc.
        """
        if isinstance(self, DiscreteMarkovChain):
            trans_probs = self.transition_probabilities
        elif isinstance(self, ContinuousMarkovChain):
            trans_probs = self.generator_matrix
        state_space = self.state_space
        if isinstance(given_condition, And):
            gcs = given_condition.args
            given_condition = S.true
            for gc in gcs:
                if isinstance(gc, TransitionMatrixOf):
                    trans_probs = gc.matrix
                if isinstance(gc, StochasticStateSpaceOf):
                    state_space = gc.state_space
                if isinstance(gc, Relational):
                    given_condition = given_condition & gc
        if isinstance(given_condition, TransitionMatrixOf):
            trans_probs = given_condition.matrix
            given_condition = S.true
        if isinstance(given_condition, StochasticStateSpaceOf):
            state_space = given_condition.state_space
            given_condition = S.true
        return trans_probs, state_space, given_condition

    def _check_trans_probs(self, trans_probs, row_sum=1):
        """
        Helper function for checking the validity of transition
        probabilities.
        """
        if not isinstance(trans_probs, MatrixSymbol):
            rows = trans_probs.tolist()
            for row in rows:
                if (sum(row) - row_sum) != 0:
                    raise ValueError("Values in a row must sum to %s. "
                    "If you are using Float or floats then please use Rational."%(row_sum))

    def _work_out_state_space(self, state_space, given_condition, trans_probs):
        """
        Helper function to extract state space if there
        is a random symbol in the given condition.
        """
        # if given condition is None, then there is no need to work out
        # state_space from random variables
        if given_condition != None:
            rand_var = list(given_condition.atoms(RandomSymbol) -
                        given_condition.atoms(RandomIndexedSymbol))
            if len(rand_var) == 1:
                state_space = rand_var[0].pspace.set
        if not FiniteSet(*[i for i in range(trans_probs.shape[0])]).is_subset(state_space):
            raise ValueError("state space is not compatible with the transition probabilities.")
        state_space = FiniteSet(*[i for i in range(trans_probs.shape[0])])
        return state_space

    @cacheit
    def _preprocess(self, given_condition, evaluate):
        """
        Helper function for pre-processing the information.
        """
        is_insufficient = False

        if not evaluate: # avoid pre-processing if the result is not to be evaluated
            return (True, None, None, None)

        # extracting transition matrix and state space
        trans_probs, state_space, given_condition = self._extract_information(given_condition)

        # given_condition does not have sufficient information
        # for computations
        if trans_probs == None or \
            given_condition == None:
            is_insufficient = True
        else:
            # checking transition probabilities
            if isinstance(self, DiscreteMarkovChain):
                self._check_trans_probs(trans_probs, row_sum=1)
            elif isinstance(self, ContinuousMarkovChain):
                self._check_trans_probs(trans_probs, row_sum=0)

            # working out state space
            state_space = self._work_out_state_space(state_space, given_condition, trans_probs)

        return is_insufficient, trans_probs, state_space, given_condition

    def probability(self, condition, given_condition=None, evaluate=True, **kwargs):
        """
        Handles probability queries for Markov process.

        Parameters
        ==========

        condition: Relational
        given_condition: Relational/And

        Returns
        =======
        Probability
            If the information is not sufficient.
        Expr
            In all other cases.

        Note
        ====
        Any information passed at the time of query overrides
        any information passed at the time of object creation like
        transition probabilities, state space.
        Pass the transition matrix using TransitionMatrixOf,
        generator matrix using GeneratorMatrixOf and state space
        using StochasticStateSpaceOf in given_condition using & or And.
        """
        check, mat, state_space, new_given_condition = \
            self._preprocess(given_condition, evaluate)

        if check:
            return Probability(condition, new_given_condition)

        if isinstance(self, ContinuousMarkovChain):
            trans_probs = self.transition_probabilities(mat)
        elif isinstance(self, DiscreteMarkovChain):
            trans_probs = mat

        if isinstance(condition, Relational):
            rv, states = (list(condition.atoms(RandomIndexedSymbol))[0], condition.as_set())
            if isinstance(new_given_condition, And):
                gcs = new_given_condition.args
            else:
                gcs = (new_given_condition, )
            grvs = new_given_condition.atoms(RandomIndexedSymbol)

            min_key_rv = None
            for grv in grvs:
                if grv.key <= rv.key:  # TODO: remove. Knowledge of the future does give knowledge of the present
                    min_key_rv = grv
            if min_key_rv == None:
                return Probability(condition)

            prob, gstate = dict(), None
            for gc in gcs:
                if gc.has(min_key_rv):
                    if gc.has(Probability):
                        p, gp = (gc.rhs, gc.lhs) if isinstance(gc.lhs, Probability) \
                                    else (gc.lhs, gc.rhs)
                        gr = gp.args[0]
                        gset = Intersection(gr.as_set(), state_space)
                        gstate = list(gset)[0]
                        prob[gset] = p
                    else:
                        _, gstate = (gc.lhs.key, gc.rhs) if isinstance(gc.lhs, RandomIndexedSymbol) \
                                    else (gc.rhs.key, gc.lhs)

            if any((k not in self.index_set) for k in (rv.key, min_key_rv.key)):
                raise IndexError("The timestamps of the process are not in it's index set.")
            states = Intersection(states, state_space)
            for state in Union(states, FiniteSet(gstate)):
                if Ge(state, mat.shape[0]) == True:
                    raise IndexError("No information is available for (%s, %s) in "
                        "transition probabilities of shape, (%s, %s). "
                        "State space is zero indexed."
                        %(gstate, state, mat.shape[0], mat.shape[1]))
            if prob:
                gstates = Union(*prob.keys())
                if len(gstates) == 1:
                    gstate = list(gstates)[0]
                    gprob = list(prob.values())[0]
                    prob[gstates] = gprob
                elif len(gstates) == len(state_space) - 1:
                    gstate = list(state_space - gstates)[0]
                    gprob = S.One - sum(prob.values())
                    prob[state_space - gstates] = gprob
                else:
                    raise ValueError("Conflicting information.")
            else:
                gprob = S.One

            if min_key_rv == rv:
                return sum([prob[FiniteSet(state)] for state in states])
            if isinstance(self, ContinuousMarkovChain):
                return gprob * sum([trans_probs(rv.key - min_key_rv.key).__getitem__((gstate, state))
                                    for state in states])
            if isinstance(self, DiscreteMarkovChain):
                return gprob * sum([(trans_probs**(rv.key - min_key_rv.key)).__getitem__((gstate, state))
                                    for state in states])

        if isinstance(condition, Not):
            expr = condition.args[0]
            return S.One - self.probability(expr, given_condition, evaluate, **kwargs)

        if isinstance(condition, And):
            compute_later, state2cond, conds = [], dict(), condition.args
            for expr in conds:
                if isinstance(expr, Relational):
                    ris = list(expr.atoms(RandomIndexedSymbol))[0]
                    if state2cond.get(ris, None) is None:
                        state2cond[ris] = S.true
                    state2cond[ris] &= expr
                else:
                    compute_later.append(expr)
            ris = []
            for ri in state2cond:
                ris.append(ri)
                cset = Intersection(state2cond[ri].as_set(), state_space)
                if len(cset) == 0:
                    return S.Zero
                state2cond[ri] = cset.as_relational(ri)
            sorted_ris = sorted(ris, key=lambda ri: ri.key)
            prod = self.probability(state2cond[sorted_ris[0]], given_condition, evaluate, **kwargs)
            for i in range(1, len(sorted_ris)):
                ri, prev_ri = sorted_ris[i], sorted_ris[i-1]
                if not isinstance(state2cond[ri], Eq):
                    raise ValueError("The process is in multiple states at %s, unable to determine the probability."%(ri))
                mat_of = TransitionMatrixOf(self, mat) if isinstance(self, DiscreteMarkovChain) else GeneratorMatrixOf(self, mat)
                prod *= self.probability(state2cond[ri], state2cond[prev_ri]
                                 & mat_of
                                 & StochasticStateSpaceOf(self, state_space),
                                 evaluate, **kwargs)
            for expr in compute_later:
                prod *= self.probability(expr, given_condition, evaluate, **kwargs)
            return prod

        if isinstance(condition, Or):
            return sum([self.probability(expr, given_condition, evaluate, **kwargs)
                        for expr in condition.args])

        raise NotImplementedError("Mechanism for handling (%s, %s) queries hasn't been "
                                "implemented yet."%(condition, given_condition))

    def expectation(self, expr, condition=None, evaluate=True, **kwargs):
        """
        Handles expectation queries for markov process.

        Parameters
        ==========

        expr: RandomIndexedSymbol, Relational, Logic
            Condition for which expectation has to be computed. Must
            contain a RandomIndexedSymbol of the process.
        condition: Relational, Logic
            The given conditions under which computations should be done.

        Returns
        =======

        Expectation
            Unevaluated object if computations cannot be done due to
            insufficient information.
        Expr
            In all other cases when the computations are successful.

        Note
        ====

        Any information passed at the time of query overrides
        any information passed at the time of object creation like
        transition probabilities, state space.

        Pass the transition matrix using TransitionMatrixOf,
        generator matrix using GeneratorMatrixOf and state space
        using StochasticStateSpaceOf in given_condition using & or And.
        """

        check, mat, state_space, condition = \
            self._preprocess(condition, evaluate)

        if check:
            return Expectation(expr, condition)

        rvs = random_symbols(expr)
        if isinstance(expr, Expr) and isinstance(condition, Eq) \
            and len(rvs) == 1:
            # handle queries similar to E(f(X[i]), Eq(X[i-m], <some-state>))
            rv = list(rvs)[0]
            lhsg, rhsg = condition.lhs, condition.rhs
            if not isinstance(lhsg, RandomIndexedSymbol):
                lhsg, rhsg = (rhsg, lhsg)
            if rhsg not in self.state_space:
                raise ValueError("%s state is not in the state space."%(rhsg))
            if rv.key < lhsg.key:
                raise ValueError("Incorrect given condition is given, expectation "
                    "time %s < time %s"%(rv.key, rv.key))
            mat_of = TransitionMatrixOf(self, mat) if isinstance(self, DiscreteMarkovChain) else GeneratorMatrixOf(self, mat)
            cond = condition & mat_of & \
                    StochasticStateSpaceOf(self, state_space)
            func = lambda s: self.probability(Eq(rv, s), cond)*expr.subs(rv, s)
            return sum([func(s) for s in state_space])

        raise NotImplementedError("Mechanism for handling (%s, %s) queries hasn't been "
                                "implemented yet."%(expr, condition))


class CommunicationClass:
    def __init__(self, states: list, recurrent: bool = None, period: int = None):
        """Represents a single communication class for a Markov Chain.

        Parameters
        ==========

        states : list
            The list of integer states that are part of the equivalence class.
        recurrent : bool, optional
            Whether all the states in the class are recurrent or transient.
        period: int, optional
            The period of the states. None if it is transient.
        """
        self.states = states
        self.recurrent = recurrent
        self.period = period

    def __add__(self, other):
        """Finds the union of two communication classes"""
        states = list(set(self.states).union(other.states))
        if self.recurrent is not None:
            recurrent = self.recurrent
        elif other.recurrent is not None:
            recurrent = other.recurrent
        else:
            recurrent = None
        if self.period is not None:
            period = self.period
        elif other.period is not None:
            period = other.period
        else:
            period = None
        states.sort()
        return CommunicationClass(states, recurrent, period)

    def __contains__(self, item):
        """Checks if a state is in the communication class"""
        return item in self.states

    def __eq__(self, other):
        """Checks if two communication classes have the same states"""
        self.states.sort()
        other.states.sort()
        return self.states == other.states

    def __repr__(self):
        return str((self.states, self.recurrent, self.period))

class DiscreteMarkovChain(DiscreteTimeStochasticProcess, MarkovProcess):
    """
    Represents a finite discrete time-homogeneous Markov chain.

    This type of Markov Chain can be uniquely characterised by
    its (ordered) state space and its one-step transition probability
    matrix.

    Parameters
    ==========

    sym: Symbol/str
    state_space: Set, optional
        The default is ``S.Reals`` and becomes {0, 1, ..., n} if
        ``trans_probs`` is given. n is the number of states in ``trans_probs``.
    trans_probs: Matrix/ImmutableMatrix/MatrixSymbol/FunctionMatrix, optional
        The one-step transition probability matrix.

    Examples
    ========

    >>> from sympy.stats import DiscreteMarkovChain, TransitionMatrixOf
    >>> from sympy import Matrix, MatrixSymbol, Eq
    >>> from sympy.stats import P
    >>> T = Matrix([[0.5, 0.2, 0.3],[0.2, 0.5, 0.3],[0.2, 0.3, 0.5]])
    >>> Y = DiscreteMarkovChain("Y", [0, 1, 2], T)
    >>> YS = DiscreteMarkovChain("Y")
    >>> Y.state_space
    FiniteSet(0, 1, 2)
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

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Markov_chain#Discrete-time_Markov_chain
    .. [2] https://www.dartmouth.edu/~chance/teaching_aids/books_articles/probability_book/Chapter11.pdf
    """
    index_set = S.Naturals0

    def __new__(cls, sym, state_space=S.Reals, trans_probs=None):
        sym = _symbol_converter(sym)

        # maybe make trans_probs a MatrixSymbol if it is None
        if trans_probs is not None:
            trans_probs = _matrix_checks(trans_probs)

        if (state_space is S.Reals) and (trans_probs is not None):
            # handle symbolic-sized matrix
            n = trans_probs.shape[0]
            if isinstance(n, Symbol):
                state_space = Range(0, n)
            else:
                state_space = FiniteSet(*Range(0, n))
        state_space = _set_converter(state_space)

        return Basic.__new__(cls, sym, state_space, trans_probs)

    @property
    def transition_probabilities(self):
        """
        Transition probabilities of discrete Markov chain,
        either an instance of Matrix or MatrixSymbol.
        """
        return self.args[2]

    @property
    def num_states(self):
        """
        The number of states in the Markov Chain. Can be symbolic.
        """
        if self.transition_probabilities is None:
            return None
        else:
            return self.transition_probabilities.shape[0]

    @property
    def is_irreducible(self):
        """
        A Markov Chain is irreducible iff it consists of exactly
        one communication class.
        """
        trans_probs = self.transition_probabilities
        if isinstance(trans_probs, MatrixSymbol):
            return None  # cannot be known
        return len(self.communication_classes()) == 1


    def _transient2transient(self):
        """
        Computes the one step probabilities of transient
        states to transient states. Used in finding
        fundamental matrix, absorbing probabilities.

        This will be depreciated in the future.
        """
        trans_probs = self.transition_probabilities
        if not isinstance(trans_probs, ImmutableMatrix):
            return None

        m = trans_probs.shape[0]
        trans_states = [i for i in range(m) if trans_probs[i, i] != 1]
        t2t = [[trans_probs[si, sj] for sj in trans_states] for si in trans_states]

        return ImmutableMatrix(t2t)

    def _transient2absorbing(self):
        """
        Computes the one step probabilities of transient
        states to absorbing states. Used in finding
        fundamental matrix, absorbing probabilities.

        This will be depreciated in the future.
        """
        trans_probs = self.transition_probabilities
        if not isinstance(trans_probs, ImmutableMatrix):
            return None

        m, trans_states, absorb_states = \
            trans_probs.shape[0], [], []
        for i in range(m):
            if trans_probs[i, i] == 1:
                absorb_states.append(i)
            else:
                trans_states.append(i)

        if not absorb_states or not trans_states:
            return None

        t2a = [[trans_probs[si, sj] for sj in absorb_states]
               for si in trans_states]

        return ImmutableMatrix(t2a)

    def stationary_distribution(self):
        """
        The stationary distribution is a row vector, p, solves p = pP,
        is row stochastic and each element in p must be nonnegative.
        That means in matrix form: :math:`(P-I)^T p^T = 0` and
        :math:`(1, ..., 1) p = 1`
        where ``P`` is the one-step transition matrix.

        All time-homogeneous Markov Chains with a finite state space
        have at least one stationary distribution. In addition, if
        a finite time-homogeneous Markov Chain is irreducible, the
        stationary distribution is unique.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        An irreducible Markov Chain

        >>> T = Matrix([[S(1)/2, S(1)/2, 0],
        ...             [S(4)/5, S(1)/5, 0],
        ...             [1, 0, 0]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.stationary_distribution()
        Matrix([[8/13, 5/13, 0]])

        A reducible Markov Chain

        >>> T = Matrix([[S(1)/2, S(1)/2, 0],
        ...             [S(4)/5, S(1)/5, 0],
        ...             [0, 0, 1]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.stationary_distribution()
        Matrix([[8/13 - 8*tau0/13, 5/13 - 5*tau0/13, tau0]])

        References
        ==========

        .. [1] https://www.probabilitycourse.com/chapter11/11_2_6_stationary_and_limiting_distributions.php
        .. [2] https://galton.uchicago.edu/~yibi/teaching/stat317/2014/Lectures/Lecture4_6up.pdf

        See Also
        ========

        sympy.stats.stochastic_processes_types.DiscreteMarkovChain.limiting_distribution
        """
        trans_probs = self.transition_probabilities
        if trans_probs is None:
            return None
        n = self.num_states

        # matrix symbol version
        if isinstance(trans_probs, MatrixSymbol) or isinstance(n, Symbol):
            wm = MatrixSymbol('wm', 1, n)
            # the following throws an error when checking `w in set.subs(T_symbol, T_numeric)`
            # return ConditionSet(wm, Eq(wm*trans_probs, wm))  # and wm must be row stochastic
            if isinstance(trans_probs, FunctionMatrix):
                return Lambda(wm, Eq(wm * trans_probs, wm))
            return Lambda((wm, trans_probs), Eq(wm*trans_probs, wm))

        if n == 0:
            return Matrix([[]])

        # numeric matrix version
        a = Matrix(trans_probs - Identity(n)).T
        a[0, 0:n] = ones(rows=1, cols=n)

        b = zeros(rows=n, cols=1)
        b[0, 0] = 1

        try:
            pi_, params = a.gauss_jordan_solve(b)
            pi_ = pi_.T
        except ValueError:
            pi_ = None
        return pi_

    @property
    def limiting_distribution(self):
        """A wrapper property for the function ``limiting_dist()``."""
        return self.limiting_dist()

    def limiting_dist(self, p0=None):
        """
        The limiting distribution is the row vector
        equal to :math:`lim_{t \\rightarrow \\infty} p^{(n)}`.
        This property has been updated to a function.
        If you are exeriencing errors with old code,
        instead of ``X.limiting_distribution`` try
        ``X.limiting_distribution()`` instead. This was
        chosen to match ``ContinuousMarkovChain`` and
        to allow the argument ``p0``.

        where
        :math:`p^{(n)}` is the marginal state probability vector at time n.

        This is equal to the stationary distribution if
        the Markov Chain is irreducible, aperiodic and
        on a finite state space.

        Parameters
        ==========

        p0
            The inital state vector as a row Matrix. This
            is only needed for reducible or periodic
            Markov Chains. If None is given, it defaults
            to a ``MatrixSymbol`` of size 1 by n.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        An irreducible aperiodic finite Markov Chain

        >>> T = Matrix([[S(1)/2, S(1)/2, 0],
        ...             [S(2)/5, S(1)/5, S(2)/5],
        ...             [S(1)/5, S(3)/5, S(1)/5]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.limiting_distribution
        Matrix([[2/5, 2/5, 1/5]])

        The following has the same effect.
        This however lets you specify an argument ``p0`` if needed.

        >>> X.limiting_dist()
        Matrix([[2/5, 2/5, 1/5]])

        A reducible aperiodic Markov Chain
        >>> T = Matrix([[S(1)/2, S(1)/2, 0],
        ...             [S(4)/5, S(1)/5, 0],
        ...             [1, 0, 0]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.limiting_distribution
        Matrix([[8/13, 5/13, 0]])

        An irreducible prediodic Markov Chain if the limit cannot be found

        >>> T = Matrix([[0, 1],
        ...             [1, 0]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.limiting_dist(Matrix([[S(2)/3, S(1)/3]]))
        Matrix([[AccumBounds(1/3, 2/3), AccumBounds(1/3, 2/3)]])

        References
        ==========

        .. [1] https://www.probabilitycourse.com/chapter11/11_2_6_stationary_and_limiting_distributions.php
        .. [2] https://galton.uchicago.edu/~yibi/teaching/stat317/2014/Lectures/Lecture4_6up.pdf

        See Also
        ========

        sympy.stats.stochastic_processes_types.DiscreteMarkovChain.stationary_distribution
        """
        trans_probs = self.transition_probabilities
        n = self.num_states
        if (trans_probs is None) or (n is None):
            return None

        if isinstance(trans_probs, MatrixSymbol) or isinstance(n, Symbol):
            return self.stationary_distribution()

        comm_classes = self.communication_classes()

        # if irreducible and aperiodic
        if (len(comm_classes) == 1) and (comm_classes[0].period == 1):
            return self.stationary_distribution()

        # limits do not work well with matrices
        # I think this is still the stationary distribution
        # if all recurrent classes are aperiodic.
        # I have no sources to back this up though.
        if all([c.period == 1 for c in comm_classes if c.recurrent]):
            return self.stationary_distribution()

        # If the chain is reducible or periodic, the
        # limiting distribution might not be unique or exist.
        # It would depend on the time 0 marginal
        # state probability vector and we choose it
        # to simply be uniform.
        if trans_probs is None:
            return None

        if isinstance(trans_probs, MatrixSymbol):
            raise NotImplementedError("Limits on MatrixSymbol is not implemented.")

        if n == 0:  # a zero by zero matrix to any power is just itself
            return Matrix([[]])

        p0_is_given = p0 is not None  # whether or not the user gave their own p0.
        # using p0_is_given will give more refined AccumBounds in the example.
        # If we treat p0 the same way as whether it was given or not,
        # the lower and upper bounds will be wider than necessary for a given p0.
        # The last example in the docstring would give
        # Matrix([[AccumBounds(0, 1), AccumBounds(0, 1)]]).
        # so if p0 is given, we add it into the mix as soon as possible.

        if not p0_is_given:
            p0 = MatrixSymbol('p_0', 1, n)  # implicitly row-stochastic

        _n = Dummy("n", positive=True, integer=True)
        if p0_is_given:
            P_to_the_n = p0 * Matrix(trans_probs) ** _n
        else:
            P_to_the_n = Matrix(trans_probs) ** _n

        for row in range(P_to_the_n.shape[0]):
            for col in range(P_to_the_n.shape[1]):
                try:  # try to limit pn but it sometimes does not work
                    P_to_the_n[row, col] = limit_seq(P_to_the_n[row, col], _n)
                except NotImplementedError:
                    pass

        if p0_is_given:
            pn = P_to_the_n
        else:
            pn = p0 * P_to_the_n
        return pn

    def communication_classes(self):
        """
        Returns the list of communication classes that partition
        the states of the markov chain. It will return None if
        no transition matrix was given.

        Returns
        =======

        classes : List of CommunicationClass, optional
            The list of communication classes that make up the
            Markov Chain.

        Notes
        =====

        This method uses the following algorithm:

        To find the periods of each state, take the one-step transition
        matrix P, find P, P^2, ..., P^n where n is the number of states.
        Analyse the diagonals of these matrices. For state i, find the
        powers of P for which P^n[i, i] > 0. Find their greatest
        common divisor of this list for each state.

        To find the communication classes, one could take P^n and find
        which the states for which j is accessible from i and generate
        the classes. This, however, loses track of the states the process
        visited before the nth time step. Instead take X = (P + I)^n and
        one can build the classes this way. For example, if X[i, j] > 0
        then j is accessible from i.

        If n is large (about 30), the values in (P + I)^n may also become
        large (about 2^30). This may cause memory or overflow errors.
        Replacing I with 0.1*I in the formula could slow this process
        but it should not be necessary for most uses.
        """
        trans_probs = self.transition_probabilities
        n = self.num_states
        if (trans_probs is None) or (n is None):
            return None

        if isinstance(trans_probs, MatrixSymbol) or isinstance(n, Symbol):
            raise NotImplementedError("The transient and recurrent states cannot be determined.")

        temp_Pn = Identity(n)
        periods = [-1] * n
        for i in range(1, n + 1):
            temp_Pn = temp_Pn * trans_probs
            for diag in range(n):
                if (periods[diag] == -1) and (temp_Pn[diag, diag] != 0):
                    periods[diag] = i
                elif (periods[diag] != -1) and (temp_Pn[diag, diag] != 0) and (
                        i % periods[diag] != 0):
                    periods[diag] = 1

        # Normal power of P means that it will see periodic
        # states as being all in their own classes.
        # We must make a new P so that all states are aperiodic for this.
        # Note that P_prime does not have to be row stochastic
        P_prime = Matrix(trans_probs)
        for diag in range(n):
            P_prime[diag, diag] = 1

        Pn = P_prime ** n
        classes = [CommunicationClass([i], True, periods[i]) for i in range(n)]
        for row in range(0, n):
            for col in range(row + 1, n):
                # we use != instead of > in order to deal with symbols
                if (Pn[row, col] != 0) and (Pn[col, row] != 0):  # if row and col communicate
                    classes[row] = classes[row] + classes[col]
                    classes[col] = classes[row]
                elif (Pn[row, col] != 0) and (Pn[col, row] == 0):  # if row empties into column
                    classes[row].recurrent = False
                elif (Pn[row, col] == 0) and (Pn[col, row] != 0):  # if column empties into row
                    classes[col].recurrent = False

        new_classes = []
        to_remove = []
        for i in range(len(classes)):
            to_remove += classes[i].states[1:]
            if i not in to_remove:
                new_classes.append(classes[i])
        classes = new_classes
        return classes

    def canonical_form(self):
        """
        Reorders the one-step transition matrix
        so that recurrent states appear first and transient
        states appear last. Other notations include inserting
        transient states first and recurrent states last but
        that method creates n-step transition matrix with
        poor visual appeal.

        Returns
        =======

        (states, P_new)
            ``states`` is the list that describes the order of the
            new states in the matrix
            so that the ith element in ``states`` is the state of the
            ith row of A.
            ``P_new`` is the new transition matrix in canonical form.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        You can convert your chain into canonical form

        >>> T = Matrix([[S(1)/2, S(1)/2, 0, 0, 0],
        ...             [S(2)/5, S(1)/5, S(2)/5, 0, 0],
        ...             [0, 0, 1, 0, 0],
        ...             [0, 0, S(1)/2, S(1)/2, 0],
        ...             [S(1)/2, 0, 0, 0, S(1)/2]])
        >>> X = DiscreteMarkovChain('X', list(range(1, 6)), trans_probs=T)
        >>> X.canonical_form()
        ([3, 1, 2, 4, 5], Matrix([
        [  1,   0,   0,   0,   0],
        [  0, 1/2, 1/2,   0,   0],
        [2/5, 2/5, 1/5,   0,   0],
        [1/2,   0,   0, 1/2,   0],
        [  0, 1/2,   0,   0, 1/2]]))

        The new states are [3, 1, 2, 4, 5] and you can create a new chain with this

        >>> X = DiscreteMarkovChain('X', X.canonical_form()[0], X.canonical_form()[1])
        >>> X.canonical_form()
        ([1, 2, 3, 4, 5], Matrix([
        [  1,   0,   0,   0,   0],
        [  0, 1/2, 1/2,   0,   0],
        [2/5, 2/5, 1/5,   0,   0],
        [1/2,   0,   0, 1/2,   0],
        [  0, 1/2,   0,   0, 1/2]]))

        The output should be the same as the first call but setting a state space
        results in that state space being sorted instead of kept in place.

        See Also
        ========

        sympy.stats.stochastic_processes_types.DiscreteMarkovChain.decompose
        """
        trans_probs = self.transition_probabilities
        n = self.num_states
        if (trans_probs is None) or (n is None):
            return None

        if isinstance(trans_probs, MatrixSymbol) or isinstance(n, Symbol):
            raise NotImplementedError("The transient and recurrent states cannot be determined.")

        classes = self.communication_classes()
        r_states = []
        t_states = []
        for c in classes:
            if c.recurrent:
                r_states += c.states
            else:
                t_states += c.states

        states = r_states + t_states
        P_new = Matrix(n, n, lambda i, j: trans_probs[states[i], states[j]])

        # convert states to the user's states
        states = [self.state_space.args[state] for state in states]
        return states, P_new

    def decompose(self):
        """
        The transition matrix can be decomposed into 4 submatrices:

        - A - the submatrix from recurrent states to recurrent states
        - B - the submatrix from transient to recurrent states
        - C - the submatrix from transient to transient states
        - 0 - the submatrix of zeros for recurrent to transient states

        Returns
        =======

        (states, A, B, C)

            ``states`` - a list of state names with the first being
            the recurrent states and the last being
            the transient states in the order
            of the row names of A and then the row names of C.
            ``A`` - the submatrix from recurrent states to recurrent states.
            ``B`` - the submatrix from transient to recurrent states.
            ``C`` - the submatrix from transient to transient states.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        You can decompose this matrix for example

        >>> T = Matrix([[S(1)/2, S(1)/2, 0, 0, 0],
        ...             [S(2)/5, S(1)/5, S(2)/5, 0, 0],
        ...             [0, 0, 1, 0, 0],
        ...             [0, 0, S(1)/2, S(1)/2, 0],
        ...             [S(1)/2, 0, 0, 0, S(1)/2]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> states, A, B, C = X.decompose()
        >>> states
        [2, 0, 1, 3, 4]
        >>> A   # recurrent to recurrent
        Matrix([[1]])
        >>> B  # transient to recurrent
        Matrix([
        [  0],
        [2/5],
        [1/2],
        [  0]])
        >>> C  # transient to transient
        Matrix([
        [1/2, 1/2,   0,   0],
        [2/5, 1/5,   0,   0],
        [  0,   0, 1/2,   0],
        [1/2,   0,   0, 1/2]])

        This means that state 2 is the only absorbing state
        (since A is a 1x1 matrix). B is a 4x1 matrix since
        the 4 remaining transient states all merge into reccurent
        state 2. And C is the 4x4 matrix that shows how the
        transient states 0, 1, 3, 4 all interact.

        See Also
        ========

        sympy.stats.stochastic_processes_types.DiscreteMarkovChain.canonical_form

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Absorbing_Markov_chain
        .. [2] http://people.brandeis.edu/~igusa/Math56aS08/Math56a_S08_notes015.pdf

        """
        trans_probs = self.transition_probabilities
        if trans_probs is None:
            return None
        n = self.num_states

        if isinstance(trans_probs, MatrixSymbol) or isinstance(n, Symbol):
            raise NotImplementedError("The transient and recurrent states cannot be determined.")

        classes = self.communication_classes()
        r_states = []
        t_states = []
        for c in classes:
            if c.recurrent:
                r_states += c.states
            else:
                t_states += c.states

        states = r_states + t_states

        A = Matrix(len(r_states), len(r_states),
                   lambda i, j: trans_probs[states[i], states[j]])

        B = Matrix(len(t_states), len(r_states),
                   lambda i, j: trans_probs[states[len(r_states) + i], states[j]])

        C = Matrix(len(t_states), len(t_states),
                   lambda i, j: trans_probs[states[len(r_states) + i], states[len(r_states) + j]])

        return states, A, B, C

    def fundamental_matrix(self, C=None):
        """
        The fundamental matrix, :math:`M = (I - C)^{-1}` where C is
        the submatrix that takes transient states to transient states
        and I is the identity matrix. C can be obtained from the
        ``decompose`` method.

        Parameters
        ==========

        C
            The submatrix of the transition matrix that
            has probabilities of transitions from transient
            to transient states. This will be computed if
            not given.
        """
        if C is None:
            states, A, B, C = self.decompose()

        # expression make explicit since decompose needs an integer number of states
        I = Identity(C.shape[0]).as_explicit()
        if (I - C).det() == 0:  # .det() is not implemented for Identity
            raise ValueError("The the fundamental matrix does not exist.")
        M = (I - C)**-1
        return M

    def limiting_transient_matrix(self):
        """
        The limiting one-step sub-transition matrix
        for transitions from transient to recurrent states.
        It is given by :math:`MBQ` where M is the fundamental
        matrix, B is the submatrix for transitions from recurrent
        to transient states and Q is the limiting transition
        matrix for transitions from recurrent to recurrent
        states. Serves a similar purpose as absorbing probabilities.

        This is equal to the exit probability matrix
        if each state in the recurrent states is
        absorbing. That is, if each state,
        i, in the recurrent states has ``T[i, i] == 1`` where
        ``T`` is the one-step transition matrix.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        The recurrent states are not absorbing since they
        interact with one another.

        >>> T = Matrix([[S(1)/2, S(1)/2, 0, 0],
        ...             [S(4)/5, S(1)/5, 0, 0],
        ...             [S(1)/2, S(1)/3, S(1)/6, 0],
        ...             [0, S(1)/2, S(1)/2, 0]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.limiting_transient_matrix()
        Matrix([
        [8/13, 5/13],
        [8/13, 5/13]])

        Note that the 2 rows represent the two transient
        states and the two columns represent the 2 absorbing
        states.
        """
        states, A, B, C = self.decompose()

        M = Matrix(self.fundamental_matrix(C))

        _n = Dummy("n", positive=True, integer=True)
        Q_ = A**_n
        Q = Matrix(Q_.shape[0], Q_.shape[1],
                   lambda i, j: limit_seq(Q_[i, j], _n))
        return M*B*Q

    # Change to exit_probability_matrix in the future.
    # Typo in    `probabilities` for backwards compatibility.
    def absorbing_probabilites(self):
        """
        The exit probability matrix. The element :math:`e_{ij}` in E is the
        probability that, starting in state i, the process
        leaves the transient states and enters recurrent state j
        on the first step out of the transient states.
        It is given by :math:`MB` where M is the fundamental
        matrix and B is the submatrix for transitions from recurrent
        to transient states.
        Serves the same purpose as absorbing probabilities.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        The recurrent states are not absorbing since they
        interact with one another.

        >>> T = Matrix([[S(1)/2, S(1)/2, 0, 0],  # recurrent
        ...             [S(4)/5, S(1)/5, 0, 0],  # recurrent
        ...             [S(2)/3, S(1)/3, 0, 0],  # transient
        ...             [0, S(1)/2, S(1)/2, 0]])  # transient
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.absorbing_probabilites()
        Matrix([
        [2/3, 1/3],
        [1/3, 2/3]])

        Note that the 2 rows represent the two transient
        states and the two columns represent the 2 absorbing
        states. Notice how the first row of the output matches
        the third row of T since the process must exit the
        transient states.
        """
        states, A, B, C = self.decompose()
        new_from_old = states
        for m, n in enumerate(states):
            new_from_old[n] = m

        M = self.fundamental_matrix(C)
        E = M*B
        return E

    def expected_time_to_absorption(self):
        """
        Returns a column matrix where the element :math:`v_i` is
        the expected number of revisits (over
        the entire process) to transient states
        given that the process is currently in state i.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S

        The recurrent states are not absorbing since they
        interact with one another.

        >>> T = Matrix([[S(1)/2, S(1)/2, 0, 0],
        ...             [S(4)/5, S(1)/5, 0, 0],
        ...             [S(2)/3, S(1)/3, 0, 0],
        ...             [0, S(1)/2, S(1)/2, 0]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.expected_time_to_absorption()
        Matrix([
        [  0],
        [1/2]])

        Notice how the first element, the number of
        expected revisits to state 2 given that the
        process started in state 2 is 0 since it must
        leave the transient states immediatley. The last
        element is 1/2 since the process either stays
        0 or 1 time steps in the transient states.

        References
        ==========

        .. [1] https://www.stat.auckland.ac.nz/~fewster/325/notes/ch8.pdf
        """
        states, A, B, C = self.decompose()
        new_from_old = states
        for m, n in enumerate(states):
            new_from_old[n] = m

        M = self.fundamental_matrix(C)
        EV = M*ones(rows=M.shape[1], cols=1) - ones(rows=M.shape[1], cols=1)
        return EV

    def first_passage_matrix(self, t, i=None, j=None):
        """
        The first passage probability, :math:`f_{ij}^{(t)}` is the probability
        of transitioning from state i to
        state j for the first time in t number of steps. The first passage probability
        matrix, :math:`F^{(t)}` is a matrix with the (i,j)th element being
        :math:`f_{ij}^{(t)}`.
        This is a computationally expensive method especially if t is symbolic.

        Parameters
        ==========

        t : int, Symbol
            The number of time steps for which to calculate the
            first passage probability matrix. This can be an int or Symbol.

        i : int, optional
            The row index of the first passage probability matrix. This should
            be specified together with ``j`` if faster computation of a single
            value in the matrix is needed.

        j : int, optional
            The column index of the first passage probability matrix. This should
            be specified together with ``i`` if faster computation of a single
            value in the matrix is needed.

        Returns
        =======

        Ft : Matrix, Expr
            The first passage probability matrix for time step t. If i and j are
            specified, an expression is given intstead.

        Examples
        ========

        >>> from sympy.stats import DiscreteMarkovChain
        >>> from sympy import Matrix, S, symbols

        The following can be compared with the example from the References

        >>> T = Matrix([[S(2)/10, S(4)/10, S(4)/10],
        ...             [S(3)/10, S(3)/10, S(4)/10],
        ...             [S(5)/10, S(4)/10, S(1)/10]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.first_passage_matrix(2)
        Matrix([
        [  8/25, 6/25, 6/25],
        [29/100, 7/25, 6/25],
        [17/100, 6/25, 9/25]])

        The general version is the following. We restate ``T`` to floats
        so that an underlying diagonalization is done numerically. This
        is desired since this method is slow symbolically. The ``0**(t-1)``
        is 1 when t=1 and 0 when t>1.

        >>> T = Matrix([[0.6, 0.4],
        ...             [0.3, 0.7]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> t = symbols('t', integer=True, positive=True)
        >>> X.first_passage_matrix(t)
        Matrix([
        [0.428571428571428*0**(t - 1) + 0.171428571428571*0.7**(t - 1),                  0.4*0.6**(t - 1)],
        [                                             0.3*0.7**(t - 1), 0.5*0**(t - 1) + 0.2*0.6**(t - 1)]])

        The general 2 state model for t > 1

        >>> a, b = symbols('a b', positive=True)
        >>> T = Matrix([[1-a, a],
        ...             [b, 1-b]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> t = symbols('t', integer=True, positive=True)
        >>> X.first_passage_matrix(t)
        Matrix([
        [0**(t - 1)*(1 - a) + b*(0**(t - 1)*a/(b - 1) - a*(1 - b)**(t - 1)/(b - 1)),                                                         a*(1 - a)**(t - 1)],
        [                                                        b*(1 - b)**(t - 1), 0**(t - 1)*(1 - b) + a*(0**(t - 1)*b/(a - 1) - b*(1 - a)**(t - 1)/(a - 1))]])

        Specify i, j for quicker evaluation

        >>> T = Matrix([[S(1)/2, S(1)/2, 0, 0],
        ...             [S(4)/5, S(1)/5, 0, 0],
        ...             [S(2)/3, S(1)/3, 0, 0],
        ...             [0, S(1)/2, S(1)/2, 0]])
        >>> X = DiscreteMarkovChain('X', trans_probs=T)
        >>> X.first_passage_matrix(2, 1, 1)
        2/5

        References
        ==========

        .. [1] https://scholar.uwindsor.ca/cgi/viewcontent.cgi?article=1125&context=major-papers
        .. [2] http://maths.dur.ac.uk/stats/courses/ProbMC2H/_files/handouts/1516MarkovChains2H.pdf
        """
        P = self.transition_probabilities
        n = self.num_states
        if isinstance(n, Symbol):
            raise NotImplementedError("Cannot yet find Ft for symbolic sized transition matrix.")

        js_ = list(range(n))  # the columns to loop through
        calc_i_ne_j = True  # calculate the off-diagonals
        calc_i_eq_j = True  # calculate the diagonals
        if (i is not None) and (j is not None):
            js_ = [j]
            if i == j:
                calc_i_ne_j = False
            else:
                calc_i_eq_j = False

        Ft = zeros(rows=n, cols=n)  # empty matrix

        # if i != j
        if calc_i_ne_j:
            for j in js_:

                P0 = Matrix(P)
                P0[0:n, j] = zeros(rows=n, cols=1)
                F = P0**(t-1)*P
                Ft[0:n, j] = F[0:n, j]

        # if i == j
        if calc_i_eq_j:
            for j in js_:

                P_ = Matrix(P)
                P_[j, 0:n] = zeros(rows=1, cols=n)

                Pnew = zeros(rows=2*n, cols=2*n)
                Pnew[0:n, 0:n] = P
                Pnew[n:2*n, n:2*n] = P_
                Pnew[n+j, 0:n] = P[j, 0:n]

                P0 = Matrix(Pnew)
                P0[0:2*n, j] = zeros(rows=2*n, cols=1)

                F = P0**(t - 1)*Pnew

                Ft[j, j] = F[n+j, j]

        # there seem to be these popping up along the diagonals
        # if isinstance(t, Symbol):
            # Ft = Ft.replace(0**(t - 1), 0)

        if (i is not None) and (j is not None):
            return Ft[i, j]
        return Ft

    def is_regular(self):
        w = self.fixed_row_vector()
        if w is None or isinstance(w, (Lambda)):
            return None
        return all((wi > 0) == True for wi in w.row(0))

    def is_absorbing_state(self, state):
        """Checks whether a given state is absorbing (has ``P[state, state] == 1``)."""
        trans_probs = self.transition_probabilities
        if isinstance(trans_probs, ImmutableMatrix) and \
            state < trans_probs.shape[0]:
            return S(trans_probs[state, state]) is S.One

    def is_absorbing_chain(self):
        """Checks whether the Markov chain has at least one absorbing state."""
        trans_probs = self.transition_probabilities
        return any(self.is_absorbing_state(state) == True
                    for state in range(trans_probs.shape[0]))

    def fixed_row_vector(self):
        """
        The fixed row vector is the stationary
        distribution of a discrete Markov chain.
        """
        return self.stationary_distribution()

    def sample(self):
        """
        Returns
        =======

        sample: iterator object
            iterator object containing the sample

        """
        if not isinstance(self.transition_probabilities, (Matrix, ImmutableMatrix)):
            raise ValueError("Transition Matrix must be provided for sampling")
        Tlist = self.transition_probabilities.tolist()
        samps = [random.choice(list(self.state_space))]
        yield samps[0]
        time = 1
        densities = {}
        for state in self.state_space:
            states = list(self.state_space)
            densities[state] = {states[i]: Tlist[state][i]
                        for i in range(len(states))}
        while time < S.Infinity:
            samps.append((next(sample_iter(FiniteRV("_", densities[samps[time - 1]])))))
            yield samps[time]
            time += 1

class ContinuousMarkovChain(ContinuousTimeStochasticProcess, MarkovProcess):
    """
    Represents continuous time Markov chain.

    Parameters
    ==========

    sym: Symbol/str
    state_space: Set
        Optional, by default, S.Reals
    gen_mat: Matrix/ImmutableMatrix/MatrixSymbol
        Optional, by default, None

    Examples
    ========

    >>> from sympy.stats import ContinuousMarkovChain
    >>> from sympy import Matrix, S
    >>> G = Matrix([[-S(1), S(1)], [S(1), -S(1)]])
    >>> C = ContinuousMarkovChain('C', state_space=[0, 1], gen_mat=G)
    >>> C.limiting_distribution()
    Matrix([[1/2, 1/2]])

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Markov_chain#Continuous-time_Markov_chain
    .. [2] http://u.math.biu.ac.il/~amirgi/CTMCnotes.pdf
    """
    index_set = S.Reals

    def __new__(cls, sym, state_space=S.Reals, gen_mat=None):
        sym = _symbol_converter(sym)
        state_space = _set_converter(state_space)
        if gen_mat != None:
            gen_mat = _matrix_checks(gen_mat)
        return Basic.__new__(cls, sym, state_space, gen_mat)

    @property
    def generator_matrix(self):
        return self.args[2]

    @cacheit
    def transition_probabilities(self, gen_mat=None):
        t = Dummy('t')
        if isinstance(gen_mat, (Matrix, ImmutableMatrix)) and \
                gen_mat.is_diagonalizable():
            # for faster computation use diagonalized generator matrix
            Q, D = gen_mat.diagonalize()
            return Lambda(t, Q*exp(t*D)*Q.inv())
        if gen_mat != None:
            return Lambda(t, exp(t*gen_mat))

    def limiting_distribution(self):
        gen_mat = self.generator_matrix
        if gen_mat == None:
            return None
        if isinstance(gen_mat, MatrixSymbol):
            wm = MatrixSymbol('wm', 1, gen_mat.shape[0])
            return Lambda((wm, gen_mat), Eq(wm*gen_mat, wm))
        w = IndexedBase('w')
        wi = [w[i] for i in range(gen_mat.shape[0])]
        wm = Matrix([wi])
        eqs = (wm*gen_mat).tolist()[0]
        eqs.append(sum(wi) - 1)
        soln = list(linsolve(eqs, wi))[0]
        return ImmutableMatrix([[sol for sol in soln]])


class BernoulliProcess(DiscreteTimeStochasticProcess):
    """
    The Bernoulli process consists of repeated
    independent Bernoulli process trials with the same parameter `p`.
    It's assumed that the probability `p` applies to every
    trial and that the outcomes of each trial
    are independent of all the rest. Therefore Bernoulli Processs
    is Discrete State and Discrete Time Stochastic Process.

    Parameters
    ==========

    sym: Symbol/str
    success: Integer/str
            The event which is considered to be success, by default is 1.
    failure: Integer/str
            The event which is considered to be failure, by default is 0.
    p: Real Number between 0 and 1
            Represents the probability of getting success.

    Examples
    ========

    >>> from sympy.stats import BernoulliProcess, P, E
    >>> from sympy import Eq, Gt
    >>> B = BernoulliProcess("B", p=0.7, success=1, failure=0)
    >>> B.state_space
    FiniteSet(0, 1)
    >>> (B.p).round(2)
    0.70
    >>> B.success
    1
    >>> B.failure
    0
    >>> X = B[1] + B[2] + B[3]
    >>> P(Eq(X, 0)).round(2)
    0.03
    >>> P(Eq(X, 2)).round(2)
    0.44
    >>> P(Eq(X, 4)).round(2)
    0
    >>> P(Gt(X, 1)).round(2)
    0.78
    >>> P(Eq(B[1], 0) & Eq(B[2], 1) & Eq(B[3], 0) & Eq(B[4], 1)).round(2)
    0.04
    >>> B.joint_distribution(B[1], B[2])
    JointDistributionHandmade(Lambda((B[1], B[2]), Piecewise((0.7, Eq(B[1], 1)), (0.3, Eq(B[1], 0)), (0, True))*Piecewise((0.7, Eq(B[2], 1)), (0.3, Eq(B[2], 0)), (0, True))))
    >>> E(2*B[1] + B[2]).round(2)
    2.10
    >>> P(B[1] < 1).round(2)
    0.30

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Bernoulli_process
    .. [2] https://mathcs.clarku.edu/~djoyce/ma217/bernoulli.pdf

    """

    index_set = S.Naturals0

    def __new__(cls, sym, p, success=1, failure=0):
        _value_check(p >= 0 and p <= 1, 'Value of p must be between 0 and 1.')
        sym = _symbol_converter(sym)
        p = _sympify(p)
        success = _sym_sympify(success)
        failure = _sym_sympify(failure)
        return Basic.__new__(cls, sym, p, success, failure)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def p(self):
        return self.args[1]

    @property
    def success(self):
        return self.args[2]

    @property
    def failure(self):
        return self.args[3]

    @property
    def state_space(self):
        return _set_converter([self.success, self.failure])

    @property
    def distribution(self):
        return BernoulliDistribution(self.p)

    def simple_rv(self, rv):
        return Bernoulli(rv.name, p=self.p,
                succ=self.success, fail=self.failure)

    def expectation(self, expr, condition=None, evaluate=True, **kwargs):
        """
        Computes expectation.

        Parameters
        ==========

        expr: RandomIndexedSymbol, Relational, Logic
            Condition for which expectation has to be computed. Must
            contain a RandomIndexedSymbol of the process.
        condition: Relational, Logic
            The given conditions under which computations should be done.

        Returns
        =======

        Expectation of the RandomIndexedSymbol.

        """

        return _SubstituteRV._expectation(expr, condition, evaluate, **kwargs)

    def probability(self, condition, given_condition=None, evaluate=True, **kwargs):
        """
        Computes probability.

        Parameters
        ==========

        condition: Relational
                Condition for which probability has to be computed. Must
                contain a RandomIndexedSymbol of the process.
        given_condition: Relational/And
                The given conditions under which computations should be done.

        Returns
        =======

        Probability of the condition.

        """

        return _SubstituteRV._probability(condition, given_condition, evaluate, **kwargs)

    def density(self, x):
        return Piecewise((self.p, Eq(x, self.success)),
                         (1 - self.p, Eq(x, self.failure)),
                         (S.Zero, True))

class _SubstituteRV:
    """
    Internal class to handle the queries of expectation and probability
    by substitution.
    """

    @staticmethod
    def _rvindexed_subs(expr, condition=None):
        """
        Substitutes the RandomIndexedSymbol with the RandomSymbol with
        same name, distribution and probability as RandomIndexedSymbol.

        Parameters
        ==========

        expr: RandomIndexedSymbol, Relational, Logic
            Condition for which expectation has to be computed. Must
            contain a RandomIndexedSymbol of the process.
        condition: Relational, Logic
            The given conditions under which computations should be done.

        """

        rvs_expr = random_symbols(expr)
        if len(rvs_expr) != 0:
            swapdict_expr = {}
            for rv in rvs_expr:
                if isinstance(rv, RandomIndexedSymbol):
                    newrv = rv.pspace.process.simple_rv(rv) # substitute with equivalent simple rv
                    swapdict_expr[rv] = newrv
            expr = expr.subs(swapdict_expr)
        rvs_cond = random_symbols(condition)
        if len(rvs_cond)!=0:
            swapdict_cond = {}
            for rv in rvs_cond:
                if isinstance(rv, RandomIndexedSymbol):
                    newrv = rv.pspace.process.simple_rv(rv)
                    swapdict_cond[rv] = newrv
            condition = condition.subs(swapdict_cond)
        return expr, condition

    @classmethod
    def _expectation(self, expr, condition=None, evaluate=True, **kwargs):
        """
        Internal method for computing expectation of indexed RV.

        Parameters
        ==========

        expr: RandomIndexedSymbol, Relational, Logic
            Condition for which expectation has to be computed. Must
            contain a RandomIndexedSymbol of the process.
        condition: Relational, Logic
            The given conditions under which computations should be done.

        Returns
        =======

        Expectation of the RandomIndexedSymbol.

        """
        new_expr, new_condition = self._rvindexed_subs(expr, condition)

        if not is_random(new_expr):
            return new_expr
        new_pspace = pspace(new_expr)
        if new_condition is not None:
            new_expr = given(new_expr, new_condition)
        if new_expr.is_Add:  # As E is Linear
            return Add(*[new_pspace.compute_expectation(
                        expr=arg, evaluate=evaluate, **kwargs)
                        for arg in new_expr.args])
        return new_pspace.compute_expectation(
                new_expr, evaluate=evaluate, **kwargs)

    @classmethod
    def _probability(self, condition, given_condition=None, evaluate=True, **kwargs):
        """
        Internal method for computing probability of indexed RV

        Parameters
        ==========

        condition: Relational
                Condition for which probability has to be computed. Must
                contain a RandomIndexedSymbol of the process.
        given_condition: Relational/And
                The given conditions under which computations should be done.

        Returns
        =======

        Probability of the condition.

        """
        new_condition, new_givencondition = self._rvindexed_subs(condition, given_condition)

        if isinstance(new_givencondition, RandomSymbol):
            condrv = random_symbols(new_condition)
            if len(condrv) == 1 and condrv[0] == new_givencondition:
                return BernoulliDistribution(self._probability(new_condition), 0, 1)

            if any([dependent(rv, new_givencondition) for rv in condrv]):
                return Probability(new_condition, new_givencondition)
            else:
                return self._probability(new_condition)

        if new_givencondition is not None and \
                not isinstance(new_givencondition, (Relational, Boolean)):
            raise ValueError("%s is not a relational or combination of relationals"
                    % (new_givencondition))
        if new_givencondition == False or new_condition == False:
            return S.Zero
        if new_condition == True:
            return S.One
        if not isinstance(new_condition, (Relational, Boolean)):
            raise ValueError("%s is not a relational or combination of relationals"
                    % (new_condition))

        if new_givencondition is not None:  # If there is a condition
        # Recompute on new conditional expr
            return self._probability(given(new_condition, new_givencondition, **kwargs), **kwargs)
        result = pspace(new_condition).probability(new_condition, **kwargs)
        if evaluate and hasattr(result, 'doit'):
            return result.doit()
        else:
            return result

def get_timerv_swaps(expr, condition):
    """
    Finds the appropriate interval for each time stamp in expr by parsing
    the given condition and returns intervals for each timestamp and
    dictionary that maps variable time-stamped Random Indexed Symbol to its
    corresponding Random Indexed variable with fixed time stamp.

    Parameters
    ==========

    expr: Sympy Expression
        Expression containing Random Indexed Symbols with variable time stamps
    condition: Relational/Boolean Expression
        Expression containing time bounds of variable time stamps in expr

    Examples
    ========

    >>> from sympy.stats.stochastic_process_types import get_timerv_swaps, PoissonProcess
    >>> from sympy import symbols, Contains, Interval
    >>> x, t, d = symbols('x t d', positive=True)
    >>> X = PoissonProcess("X", 3)
    >>> get_timerv_swaps(x*X(t), Contains(t, Interval.Lopen(0, 1)))
    ([Interval.Lopen(0, 1)], {X(t): X(1)})
    >>> get_timerv_swaps((X(t)**2 + X(d)**2), Contains(t, Interval.Lopen(0, 1))
    ... & Contains(d, Interval.Ropen(1, 4))) # doctest: +SKIP
    ([Interval.Ropen(1, 4), Interval.Lopen(0, 1)], {X(d): X(3), X(t): X(1)})

    Returns
    =======

    intervals: list
        List of Intervals/FiniteSet on which each time stamp is defined
    rv_swap: dict
        Dictionary mapping variable time Random Indexed Symbol to constant time
        Random Indexed Variable

    """

    if not isinstance(condition, (Relational, Boolean)):
        raise ValueError("%s is not a relational or combination of relationals"
            % (condition))
    expr_syms = list(expr.atoms(RandomIndexedSymbol))
    if isinstance(condition, (And, Or)):
        given_cond_args = condition.args
    else: # single condition
        given_cond_args = (condition, )
    rv_swap = {}
    intervals = []
    for expr_sym in expr_syms:
        for arg in given_cond_args:
            if arg.has(expr_sym.key) and isinstance(expr_sym.key, Symbol):
                intv = _set_converter(arg.args[1])
                diff_key = intv._sup - intv._inf
                if diff_key == oo:
                    raise ValueError("%s should have finite bounds" % str(expr_sym.name))
                elif diff_key == S.Zero: # has singleton set
                    diff_key = intv._sup
                rv_swap[expr_sym] = expr_sym.subs({expr_sym.key: diff_key})
                intervals.append(intv)
    return intervals, rv_swap


class CountingProcess(ContinuousTimeStochasticProcess):
    """
    This class handles the common methods of the Counting Processes
    such as Poisson, Wiener and Gamma Processes
    """
    index_set = _set_converter(Interval(0, oo))

    @property
    def symbol(self):
        return self.args[0]

    def expectation(self, expr, condition=None, evaluate=True, **kwargs):
        """
        Computes expectation

        Parameters
        ==========

        expr: RandomIndexedSymbol, Relational, Logic
            Condition for which expectation has to be computed. Must
            contain a RandomIndexedSymbol of the process.
        condition: Relational, Boolean
            The given conditions under which computations should be done, i.e,
            the intervals on which each variable time stamp in expr is defined

        Returns
        =======

        Expectation of the given expr

        """
        if condition is not None:
            intervals, rv_swap = get_timerv_swaps(expr, condition)
             # they are independent when they have non-overlapping intervals
            if len(intervals) == 1 or all(Intersection(*intv_comb) == EmptySet
                for intv_comb in itertools.combinations(intervals, 2)):
                if expr.is_Add:
                    return Add.fromiter(self.expectation(arg, condition)
                            for arg in expr.args)
                expr = expr.subs(rv_swap)
            else:
                return Expectation(expr, condition)

        return _SubstituteRV._expectation(expr, evaluate=evaluate, **kwargs)

    def _solve_argwith_tworvs(self, arg):
        if arg.args[0].key >= arg.args[1].key or isinstance(arg, Eq):
            diff_key = abs(arg.args[0].key - arg.args[1].key)
            rv = arg.args[0]
            arg = arg.__class__(rv.pspace.process(diff_key), 0)
        else:
            diff_key = arg.args[1].key - arg.args[0].key
            rv = arg.args[1]
            arg = arg.__class__(rv.pspace.process(diff_key), 0)
        return arg

    def _solve_numerical(self, condition, given_condition=None):
        if isinstance(condition, And):
            args_list = list(condition.args)
        else:
            args_list = [condition]
        if given_condition is not None:
            if isinstance(given_condition, And):
                args_list.extend(list(given_condition.args))
            else:
                args_list.extend([given_condition])
        # sort the args based on timestamp to get the independent increments in
        # each segment using all the condition args as well as given_condition args
        args_list = sorted(args_list, key=lambda x: x.args[0].key)
        result = []
        cond_args = list(condition.args) if isinstance(condition, And) else [condition]
        if args_list[0] in cond_args and not (is_random(args_list[0].args[0])
                        and is_random(args_list[0].args[1])):
            result.append(_SubstituteRV._probability(args_list[0]))

        if is_random(args_list[0].args[0]) and is_random(args_list[0].args[1]):
            arg = self._solve_argwith_tworvs(args_list[0])
            result.append(_SubstituteRV._probability(arg))

        for i in range(len(args_list) - 1):
            curr, nex = args_list[i], args_list[i + 1]
            diff_key = nex.args[0].key - curr.args[0].key
            working_set = curr.args[0].pspace.process.state_space
            if curr.args[1] > nex.args[1]: #impossible condition so return 0
                result.append(0)
                break
            if isinstance(curr, Eq):
                working_set = Intersection(working_set, Interval.Lopen(curr.args[1], oo))
            else:
                working_set = Intersection(working_set, curr.as_set())
            if isinstance(nex, Eq):
                working_set = Intersection(working_set, Interval(-oo, nex.args[1]))
            else:
                working_set = Intersection(working_set, nex.as_set())
            if working_set == EmptySet:
                rv = Eq(curr.args[0].pspace.process(diff_key), 0)
                result.append(_SubstituteRV._probability(rv))
            else:
                if working_set.is_finite_set:
                    if isinstance(curr, Eq) and isinstance(nex, Eq):
                        rv = Eq(curr.args[0].pspace.process(diff_key), len(working_set))
                        result.append(_SubstituteRV._probability(rv))
                    elif isinstance(curr, Eq) ^ isinstance(nex, Eq):
                        result.append(Add.fromiter(_SubstituteRV._probability(Eq(
                        curr.args[0].pspace.process(diff_key), x))
                                for x in range(len(working_set))))
                    else:
                        n = len(working_set)
                        result.append(Add.fromiter((n - x)*_SubstituteRV._probability(Eq(
                        curr.args[0].pspace.process(diff_key), x)) for x in range(n)))
                else:
                    result.append(_SubstituteRV._probability(
                    curr.args[0].pspace.process(diff_key) <= working_set._sup - working_set._inf))
        return Mul.fromiter(result)


    def probability(self, condition, given_condition=None, evaluate=True, **kwargs):
        """
        Computes probability

        Parameters
        ==========

        condition: Relational
            Condition for which probability has to be computed. Must
            contain a RandomIndexedSymbol of the process.
        given_condition: Relational, Boolean
            The given conditions under which computations should be done, i.e,
            the intervals on which each variable time stamp in expr is defined

        Returns
        =======

        Probability of the condition

        """
        check_numeric = True
        if isinstance(condition, (And, Or)):
            cond_args = condition.args
        else:
            cond_args = (condition, )
        # check that condition args are numeric or not
        if not all(arg.args[0].key.is_number for arg in cond_args):
            check_numeric = False
        if given_condition is not None:
            check_given_numeric = True
            if isinstance(given_condition, (And, Or)):
                given_cond_args = given_condition.args
            else:
                given_cond_args = (given_condition, )
            # check that given condition args are numeric or not
            if given_condition.has(Contains):
                check_given_numeric = False
            # Handle numerical queries
            if check_numeric and check_given_numeric:
                res = []
                if isinstance(condition, Or):
                    res.append(Add.fromiter(self._solve_numerical(arg, given_condition)
                            for arg in condition.args))
                if isinstance(given_condition, Or):
                    res.append(Add.fromiter(self._solve_numerical(condition, arg)
                            for arg in given_condition.args))
                if res:
                    return Add.fromiter(res)
                return self._solve_numerical(condition, given_condition)

            # No numeric queries, go by Contains?... then check that all the
            # given condition are in form of `Contains`
            if not all(arg.has(Contains) for arg in given_cond_args):
                raise ValueError("If given condition is passed with `Contains`, then "
                "please pass the evaluated condition with its corresponding information "
                "in terms of intervals of each time stamp to be passed in given condition.")

            intervals, rv_swap = get_timerv_swaps(condition, given_condition)
            # they are independent when they have non-overlapping intervals
            if len(intervals) == 1 or all(Intersection(*intv_comb) == EmptySet
                for intv_comb in itertools.combinations(intervals, 2)):
                if isinstance(condition, And):
                    return Mul.fromiter(self.probability(arg, given_condition)
                            for arg in condition.args)
                elif isinstance(condition, Or):
                    return Add.fromiter(self.probability(arg, given_condition)
                            for arg in condition.args)
                condition = condition.subs(rv_swap)
            else:
                return Probability(condition, given_condition)
        if check_numeric:
            return self._solve_numerical(condition)
        return _SubstituteRV._probability(condition, evaluate=evaluate, **kwargs)

class PoissonProcess(CountingProcess):
    """
    The Poisson process is a counting process. It is usually used in scenarios
    where we are counting the occurrences of certain events that appear
    to happen at a certain rate, but completely at random.

    Parameters
    ==========

    sym: Symbol/str
    lamda: Positive number
        Rate of the process, ``lamda > 0``

    Examples
    ========

    >>> from sympy.stats import PoissonProcess, P, E
    >>> from sympy import symbols, Eq, Ne, Contains, Interval
    >>> X = PoissonProcess("X", lamda=3)
    >>> X.state_space
    Naturals0
    >>> X.lamda
    3
    >>> t1, t2 = symbols('t1 t2', positive=True)
    >>> P(X(t1) < 4)
    (9*t1**3/2 + 9*t1**2/2 + 3*t1 + 1)*exp(-3*t1)
    >>> P(Eq(X(t1), 2) | Ne(X(t1), 4), Contains(t1, Interval.Ropen(2, 4)))
    1 - 36*exp(-6)
    >>> P(Eq(X(t1), 2) & Eq(X(t2), 3), Contains(t1, Interval.Lopen(0, 2))
    ... & Contains(t2, Interval.Lopen(2, 4)))
    648*exp(-12)
    >>> E(X(t1))
    3*t1
    >>> E(X(t1)**2 + 2*X(t2),  Contains(t1, Interval.Lopen(0, 1))
    ... & Contains(t2, Interval.Lopen(1, 2)))
    18
    >>> P(X(3) < 1, Eq(X(1), 0))
    exp(-6)
    >>> P(Eq(X(4), 3), Eq(X(2), 3))
    exp(-6)
    >>> P(X(2) <= 3, X(1) > 1)
    5*exp(-3)

    Merging two Poisson Processes

    >>> Y = PoissonProcess("Y", lamda=4)
    >>> Z = X + Y
    >>> Z.lamda
    7

    Splitting a Poisson Process into two independent Poisson Processes

    >>> N, M = Z.split(l1=2, l2=5)
    >>> N.lamda, M.lamda
    (2, 5)

    References
    ==========

    .. [1] https://www.probabilitycourse.com
    .. [2] https://en.wikipedia.org/wiki/Poisson_point_process

    """

    def __new__(cls, sym, lamda):
        _value_check(lamda > 0, 'lamda should be a positive number.')
        sym = _symbol_converter(sym)
        lamda = _sympify(lamda)
        return Basic.__new__(cls, sym, lamda)

    @property
    def lamda(self):
        return self.args[1]

    @property
    def state_space(self):
        return S.Naturals0

    def distribution(self, rv):
        return PoissonDistribution(self.lamda*rv.key)

    def density(self, x):
        return (self.lamda*x.key)**x / factorial(x) * exp(-(self.lamda*x.key))

    def simple_rv(self, rv):
        return Poisson(rv.name, lamda=self.lamda*rv.key)

    def __add__(self, other):
        if not isinstance(other, PoissonProcess):
            raise ValueError("Only instances of Poisson Process can be merged")
        return PoissonProcess(Dummy(self.symbol.name + other.symbol.name),
                self.lamda + other.lamda)

    def split(self, l1, l2):
        if _sympify(l1 + l2) != self.lamda:
            raise ValueError("Sum of l1 and l2 should be %s" % str(self.lamda))
        return PoissonProcess(Dummy("l1"), l1), PoissonProcess(Dummy("l2"), l2)

class WienerProcess(CountingProcess):
    """
    The Wiener process is a real valued continuous-time stochastic process.
    In physics it is used to study Brownian motion and therefore also known as
    Brownian Motion.

    Parameters
    ==========

    sym: Symbol/str

    Examples
    ========

    >>> from sympy.stats import WienerProcess, P, E
    >>> from sympy import symbols, Contains, Interval
    >>> X = WienerProcess("X")
    >>> X.state_space
    Reals
    >>> t1, t2 = symbols('t1 t2', positive=True)
    >>> P(X(t1) < 7).simplify()
    erf(7*sqrt(2)/(2*sqrt(t1)))/2 + 1/2
    >>> P((X(t1) > 2) | (X(t1) < 4), Contains(t1, Interval.Ropen(2, 4))).simplify()
    -erf(1)/2 + erf(2)/2 + 1
    >>> E(X(t1))
    0
    >>> E(X(t1) + 2*X(t2),  Contains(t1, Interval.Lopen(0, 1))
    ... & Contains(t2, Interval.Lopen(1, 2)))
    0

    References
    ==========

    .. [1] https://www.probabilitycourse.com
    .. [2] https://en.wikipedia.org/wiki/Wiener_process

    """
    def __new__(cls, sym):
        sym = _symbol_converter(sym)
        return Basic.__new__(cls, sym)

    @property
    def state_space(self):
        return S.Reals

    def distribution(self, rv):
        return NormalDistribution(0, sqrt(rv.key))

    def density(self, x):
        return exp(-x**2/(2*x.key)) / (sqrt(2*pi)*sqrt(x.key))

    def simple_rv(self, rv):
        return Normal(rv.name, 0, sqrt(rv.key))


class GammaProcess(CountingProcess):
    """
    A Gamma process is a random process with independent gamma distributed
    increments.  It is a pure-jump increasing Levy process.

    Parameters
    ==========

    sym: Symbol/str
    lamda: Positive number
        Jump size of the process, ``lamda > 0``
    gamma: Positive number
        Rate of jump arrivals, ``gamma > 0``

    Examples
    ========

    >>> from sympy.stats import GammaProcess, E, P, variance
    >>> from sympy import symbols, Contains, Interval, Not
    >>> t, d, x, l, g = symbols('t d x l g', positive=True)
    >>> X = GammaProcess("X", l, g)
    >>> E(X(t))
    g*t/l
    >>> variance(X(t)).simplify()
    g*t/l**2
    >>> X = GammaProcess('X', 1, 2)
    >>> P(X(t) < 1).simplify()
    lowergamma(2*t, 1)/gamma(2*t)
    >>> P(Not((X(t) < 5) & (X(d) > 3)), Contains(t, Interval.Ropen(2, 4)) &
    ... Contains(d, Interval.Lopen(7, 8))).simplify()
    -4*exp(-3) + 472*exp(-8)/3 + 1
    >>> E(X(2) + x*E(X(5)))
    10*x + 4

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Gamma_process

    """
    def __new__(cls, sym, lamda, gamma):
        _value_check(lamda > 0, 'lamda should be a positive number')
        _value_check(gamma > 0, 'gamma should be a positive number')
        sym = _symbol_converter(sym)
        gamma = _sympify(gamma)
        lamda = _sympify(lamda)
        return Basic.__new__(cls, sym, lamda, gamma)

    @property
    def lamda(self):
        return self.args[1]

    @property
    def gamma(self):
        return self.args[2]

    @property
    def state_space(self):
        return _set_converter(Interval(0, oo))

    def distribution(self, rv):
        return GammaDistribution(self.gamma*rv.key, 1/self.lamda)

    def density(self, x):
        k = self.gamma*x.key
        theta = 1/self.lamda
        return x**(k - 1) * exp(-x/theta) / (gamma(k)*theta**k)

    def simple_rv(self, rv):
        return Gamma(rv.name, self.gamma*rv.key, 1/self.lamda)

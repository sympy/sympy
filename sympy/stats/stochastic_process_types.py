from __future__ import print_function, division

from sympy import (Symbol, Matrix, MatrixSymbol, S, Indexed, Basic,
                    Set, And, Tuple, Eq, FiniteSet, ImmutableMatrix,
                    nsimplify, Lambda, Mul, Sum, Dummy, Lt, IndexedBase,
                    linsolve, Piecewise, eye, Or, Ne, Not, Intersection,
                    Union, Expr, Function, sympify, Le, exp, cacheit)
from sympy.stats.rv import (RandomIndexedSymbol, random_symbols, RandomSymbol,
                            _symbol_converter)
from sympy.stats.joint_rv import JointDistributionHandmade, JointDistribution
from sympy.core.compatibility import string_types
from sympy.core.relational import Relational
from sympy.stats.symbolic_probability import Probability, Expectation
from sympy.stats.stochastic_process import StochasticPSpace
from sympy.core.logic import Logic
from sympy.logic.boolalg import Boolean

__all__ = [
    'StochasticProcess',
    'DiscreteTimeStochasticProcess',
    'DiscreteMarkovChain',
    'TransitionMatrixOf',
    'StochasticStateSpaceOf',
    'GeneratorMatrixOf',
    'ContinuousMarkovChain'
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
        # TODO: Add tests for the below part of the method, when implementation of Bernoulli Process
        # is completed
        pdf = Lambda(*[arg.name for arg in args],
                expr=Mul.fromiter(arg.pspace.distribution.pdf(arg) for arg in args))
        return JointDistributionHandmade(pdf)

    def expectation(self, condition, given_condition):
        raise NotImplementedError("Abstract method for expectation queries.")

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
        pspace_obj = StochasticPSpace(self.symbol, self)
        return RandomIndexedSymbol(func_obj, pspace_obj)

class StochasticProcessUtil(StochasticProcess):
    """
    Abstract class containining utility methods for
    processing information contained in queries for
    stochastic processes.
    """
    def _extract_information(self, given_condition):
        """
        Helper function to extract information, like,
        transition probabilities, state space, etc.
        """
        trans_probs, state_space = self.transition_probabilities, self.state_space
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
        if isinstance(given_condition, StochasticStateSpaceOf):
            state_space = given_condition.state_space
        return trans_probs, state_space, given_condition

    def _check_trans_probs(self, trans_probs):
        """
        Helper function for checking the validity of transition
        probabilities.
        """
        if not isinstance(trans_probs, MatrixSymbol):
            rows = trans_probs.tolist()
            for row in rows:
                if (sum(row) - 1) != 0:
                    raise ValueError("Probabilities in a row must sum to 1. "
                    "If you are using Float or floats then please use Rational.")

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
            raise ValueError("state space is not compatible with the transition probabilites.")
        state_space = FiniteSet(*[i for i in range(trans_probs.shape[0])])
        return state_space

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
            self._check_trans_probs(trans_probs)

            # working out state space
            state_space = self._work_out_state_space(state_space, given_condition, trans_probs)

        return is_insufficient, trans_probs, state_space, given_condition

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

class DiscreteMarkovChain(DiscreteTimeStochasticProcess, StochasticProcessUtil):
    """
    Represents discrete Markov chain.

    Parameters
    ==========

    sym: Symbol/string_types
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

    def _transient2transient(self):
        """
        Computes the one step probabilities of transient
        states to transient states. Used in finding
        fundamental matrix, absorbing probabilties.
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
        fundamental matrix, absorbing probabilties.
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

    def fundamental_matrix(self):
        Q = self._transient2transient()
        if Q == None:
            return None
        I = eye(Q.shape[0])
        if (I - Q).det() == 0:
            raise ValueError("Fundamental matrix doesn't exists.")
        return ImmutableMatrix((I - Q).inv().tolist())

    def absorbing_probabilites(self):
        """
        Computes the absorbing probabilities, i.e.,
        the ij-th entry of the matrix denotes the
        probability of Markov chain being absorbed
        in state j starting from state i.
        """
        R = self._transient2absorbing()
        N = self.fundamental_matrix()
        if R == None or N == None:
            return None
        return N*R

    def is_regular(self):
        w = self.fixed_row_vector()
        if w is None or isinstance(w, (Lambda)):
            return None
        return all((wi > 0) == True for wi in w.row(0))

    def is_absorbing_state(self, state):
        trans_probs = self.transition_probabilities
        if isinstance(trans_probs, ImmutableMatrix) and \
            state < trans_probs.shape[0]:
            return S(trans_probs[state, state]) == S.One

    def is_absorbing_chain(self):
        trans_probs = self.transition_probabilities
        return any(self.is_absorbing_state(state) == True
                    for state in range(trans_probs.shape[0]))

    def fixed_row_vector(self):
        trans_probs = self.transition_probabilities
        if trans_probs == None:
            return None
        if isinstance(trans_probs, MatrixSymbol):
            wm = MatrixSymbol('wm', 1, trans_probs.shape[0])
            return Lambda((wm, trans_probs), Eq(wm*trans_probs, wm))
        w = IndexedBase('w')
        wi = [w[i] for i in range(trans_probs.shape[0])]
        wm = Matrix([wi])
        eqs = (wm*trans_probs - wm).tolist()[0]
        eqs.append(sum(wi) - 1)
        soln = list(linsolve(eqs, wi))[0]
        return ImmutableMatrix([[sol for sol in soln]])

    @property
    def limiting_distribution(self):
        """
        The fixed row vector is the limiting
        distribution of a discrete Markov chain.
        """
        return self.fixed_row_vector()

    def probability(self, condition, given_condition=None, evaluate=True, **kwargs):
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

        check, trans_probs, state_space, given_condition = \
            self._preprocess(given_condition, evaluate)

        if check:
            return Probability(condition, given_condition)

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
            if Lt(stateg, trans_probs.shape[0]) == False or Lt(statec, trans_probs.shape[1]) == False:
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

        info = TransitionMatrixOf(self, trans_probs) & StochasticStateSpaceOf(self, state_space)
        new_gc = given_condition & info

        if isinstance(condition, And):
            # handle queries like,
            # P(Eq(X[i+k], s1) & Eq(X[i+m], s2) . . . & Eq(X[i], sn), Eq(P(Eq(X[i], si)), prob))
            conds = condition.args
            idx2state = dict()
            for cond in conds:
                idx, state = (cond.lhs, cond.rhs) if isinstance(cond.lhs, RandomIndexedSymbol) else \
                                (cond.rhs, cond.lhs)
                idx2state[idx] = cond if idx2state.get(idx, None) is None else \
                                           idx2state[idx] & cond
            if any(len(Intersection(idx2state[idx].as_set(), state_space)) != 1
                for idx in idx2state):
                return S.Zero # a RandomIndexedSymbol cannot go to different states simultaneously
            i, result = -1, 1
            conds = And.fromiter(Intersection(idx2state[idx].as_set(), state_space).as_relational(idx)
                                    for idx in idx2state)
            if not isinstance(conds, And):
                return self.probability(conds, new_gc)
            conds = conds.args
            while i > -len(conds):
                result *= self.probability(conds[i], conds[i-1] & info)
                i -= 1
            if isinstance(given_condition, (TransitionMatrixOf, StochasticStateSpaceOf)):
                return result * Probability(conds[i])
            if isinstance(given_condition, And):
                idx_sym = conds[i].atoms(RandomIndexedSymbol)
                prob, count = S(0), 0
                for gc in given_condition.args:
                    if gc.atoms(RandomIndexedSymbol) == idx_sym:
                        prob += gc.rhs if isinstance(gc.lhs, Probability) else gc.lhs
                        count += 1
                if isinstance(state_space, FiniteSet) and \
                    count == len(state_space) - 1:
                    given_condition = Eq(Probability(conds[i]), S(1) - prob)
            if isinstance(given_condition, Eq):
                if not isinstance(given_condition.lhs, Probability) or \
                    given_condition.lhs.args[0] != conds[i]:
                    raise ValueError("Probability for %s needed", conds[i])
                return result * given_condition.rhs

        if isinstance(condition, Or):
            conds, prob_sum = condition.args, S(0)
            idx2state = dict()
            for cond in conds:
                idx, state = (cond.lhs, cond.rhs) if isinstance(cond.lhs, RandomIndexedSymbol) else \
                                (cond.rhs, cond.lhs)
                idx2state[idx] = cond if idx2state.get(idx, None) is None else \
                                           idx2state[idx] | cond
            conds = Or.fromiter(Intersection(idx2state[idx].as_set(), state_space).as_relational(idx)
                        for idx in idx2state)
            if not isinstance(conds, Or):
                return self.probability(conds, new_gc)
            return sum([self.probability(cond, new_gc) for cond in conds.args])

        if isinstance(condition, Ne):
            prob = self.probability(Not(condition), new_gc)
            return S(1) - prob

        raise NotImplementedError("Mechanism for handling (%s, %s) queries hasn't been "
                                "implemented yet."%(condition, given_condition))

    def expectation(self, expr, condition=None, evaluate=True, **kwargs):
        """
        Handles expectation queries for discrete markov chains.

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
            In all other cases when the computations are successfull.

        Note
        ====

        Any information passed at the time of query overrides
        any information passed at the time of object creation like
        transition probabilities, state space.

        Pass the transition matrix using TransitionMatrixOf and state space
        using StochasticStateSpaceOf in given_condition using & or And.
        """

        check, trans_probs, state_space, condition = \
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
            cond = condition & TransitionMatrixOf(self, trans_probs) & \
                    StochasticStateSpaceOf(self, state_space)
            s = Dummy('s')
            func = Lambda(s, self.probability(Eq(rv, s), cond)*expr.subs(rv, s))
            return Sum(func(s), (s, state_space.inf, state_space.sup)).doit()

        raise NotImplementedError("Mechanism for handling (%s, %s) queries hasn't been "
                                "implemented yet."%(expr, condition))

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

class ContinuousMarkovChain(ContinuousTimeStochasticProcess, StochasticProcessUtil):
    """
    Represents continuous Markov chain.

    Parameters
    ==========

    sym: Symbol/string_types
    state_space: Set
        Optional, by default, S.Reals
    gen_mat: Matrix/ImmutableMatrix/MatrixSymbol
        Optional, by default, None

    Examples
    ========

    >>> from sympy.stats import ContinuousMarkovChain
    >>> from sympy import Matrix, S, MatrixSymbol
    >>> G = Matrix([[-S(1), S(1)], [S(1), -S(1)]])
    >>> C = ContinuousMarkovChain('C', state_space=[0, 1], gen_mat=G)
    >>> C.limiting_distribution()
    Matrix([[1/2, 1/2]])
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

    def probability(self, condition, given_condition=None, evaluate=True, **kwargs):
        check, gen_mat, state_space, given_condition = \
            self._preprocess(given_condition, evaluate)

        if check:
            return Probability(condition, given_condition)

        trans_probs = self.transition_probabilities(gen_mat)

        if isinstance(condition, tuple):
            # to be used by inner API for computations
            rv, states = condition
            gcs, grvs = given_condition.args, given_condition.atoms(RandomIndexedSymbol)
            for state in states:
                min_key_rv = None
                for grv in grvs:
                    if grv.key <= rv.key:
                        min_key_rv = grv
                if min_key_rv == None:
                    raise ValueError("Information insufficient for generating result.")

                prob, gstate = None, None
                for gc in gcs:
                    if gc.has(min_key_rv):
                        if gc.has(Probability):
                            prob = gc.lhs if isinstance(gc.lhs, Probability) \
                                        else gc.rhs
                        else:
                            key, gstate = gc.lhs.key, gc.rhs if isinstance(gc.lhs, RandomIndexedSymbol) \
                                        else gc.rhs.key, gc.lhs

                if state != None and key <= min_key_rv.key:
                    return trans_probs(min_key_rv.key - key).__getitem__(gstate, state)

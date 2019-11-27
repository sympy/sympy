from __future__ import print_function, division

from sympy import (Matrix, MatrixSymbol, S, Indexed, Basic,
                   Set, And, Eq, FiniteSet, ImmutableMatrix,
                   Lambda, Mul, Dummy, IndexedBase,
                   linsolve, eye, Or, Not, Intersection,
                   Union, Expr, Function, exp, cacheit,
                   Ge)
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean
from sympy.stats.joint_rv import JointDistributionHandmade, JointDistribution
from sympy.stats.rv import (RandomIndexedSymbol, random_symbols, RandomSymbol,
                            _symbol_converter)
from sympy.stats.stochastic_process import StochasticPSpace
from sympy.stats.symbolic_probability import Probability, Expectation

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
        raise TypeError("Transition probabilities either should "
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
            raise ValueError("state space is not compatible with the transition probabilites.")
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
                if grv.key <= rv.key:
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
                                "implemented yet."%(expr, condition))

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

class DiscreteMarkovChain(DiscreteTimeStochasticProcess, MarkovProcess):
    """
    Represents discrete time Markov chain.

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
            return S(trans_probs[state, state]) is S.One

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

class ContinuousMarkovChain(ContinuousTimeStochasticProcess, MarkovProcess):
    """
    Represents continuous time Markov chain.

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

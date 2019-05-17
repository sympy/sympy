from sympy import Symbol, Indexed, sympify, Basic, S
from sympy.stats.rv import RandomSymbol, _value_check,  NamedArgsMixin
from sympy.stats.joint_rv import JointPSpace
from sympy.core.compatibility import string_types

class StochasticRandomSymbol(RandomSymbol):

    def __getitem__(self, time):
        if not isinstance(self.pspace, StochasticProcess):
            raise ValueError("Currently %s process is not supported by SymPy."
                                %(str(self.pspace)))
        if time not in self.pspace.index_set:
            raise ValueError("The stochastic process doesn't have %s in its "
                                " index set"%(time))
        return Indexed(self, time)

    @property
    def process(self):
        return self.pspace.process

class StochasticProcess(JointPSpace):

    def __new__(cls, sym, process=None):
        if isinstance(sym, string_types):
            sym = Symbol(sym)
        if not isinstance(sym, Symbol):
            raise ValueError("The process name should be a String or Symbol.")
        if process == None:
            process = AbstractStochasticProcess
        return Basic.__new__(cls, sym, process)

    @property
    def index_set(self):
        return self.process.index_set

    @property
    def state_space(self):
        return self.process.state_space

    @property
    def symbol(self):
        return self.args[0]

    @property
    def process(self):
        return self.args[1]

    @property
    def value(self):
        return StochasticRandomSymbol(self.symbol, self)

    def joint_dist(self, n):
        return NotImplementedError()

    def first_passage_time(self, state):
        return NotImplementedError()

    def probability(self, initial_state):
        return NotImplementedError()

def stochastic_rv(sym, cls=None, *args):
    if cls is not None:
        process = cls(*args)
        process.check(*args)
    else:
        process = cls
    return StochasticProcess(sym, process).value

class AbstractStochasticProcess(Basic):

    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @staticmethod
    def check(*args):
        pass

    @property
    def index_set(self):
        raise NotImplementedError("Index set is not defined yet.")

    @property
    def state_space(self):
        raise NotImplementedError("State space is not defined yet.")


class MarkovChain(AbstractStochasticProcess, NamedArgsMixin):

    _argnames = ('trans_matrix', 'states')

    @staticmethod
    def check(trans_matrix, states):
        _value_check(trans_matrix.is_square, "The transition matrix must be "
                                              " square.")
        _value_check(trans_matrix.shape[0] == len(states), "The transition matrix "
                     "should have same number of rows as the number of states.")

    @property
    def index_set(self):
        return S.Naturals

    @property
    def state_space(self):
        return self.states

def MarkovChainRV(sym, trans_matrix, states):
    """
    Creates a stochastic random variable following markov chain.

    Parameters
    ==========
    trans_matrix : Matrix
                   Transition probabilites, must be sqaure
    states       : list/tuple
                   States that markov random variable can attain

    Returns
    =======

    A Stochastic Random Symbol

    Examples
    ========
    >>> from sympy.stats.stochastic_process import MarkovChainRV
    >>> from sympy import Matrix
    >>> tr = Matrix([[0.5, 0.2, 0.3], [0.2, 0.3, 0.5], [0.3, 0.5, 0.2]])
    >>> states = [1, 2, 3]
    >>> M = MarkovChainRV('M', tr, states)
    >>> M.process
    MarkovChain(Matrix([
    [0.5, 0.2, 0.3],
    [0.2, 0.3, 0.5],
    [0.3, 0.5, 0.2]]), [1, 2, 3])
    >>> M.process.index_set
    Naturals
    >>> M.process.state_space
    [1, 2, 3]
    >>> M[1]
    M[1]
    >>> M[1.5]
    Traceback (most recent call last):
        ...
    ValueError: The stochastic process doesn't have 1.5 in its  index set
    """
    return stochastic_rv(sym, MarkovChain, trans_matrix, states)

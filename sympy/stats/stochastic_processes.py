from sympy import Basic, S, Symbol, Indexed, Matrix, MatrixSymbol
from sympy.core.compatibility import string_types
from sympy.stats.rv import RandomSymbol
from sympy.stats.index_rv import IndexedRandomSymbol
from sympy.stats.joint_rv import JointPSpace
from sympy.core.mul import Mul

class StochasticProcess(JointPSpace):

    def __new__(cls, sym, *args, **kwargs):
        if isinstance(sym, string_types):
            sym = Symbol(sym)
        if not isinstance(sym, Symbol):
            raise ValueError("Name of StochasticProcess must be a string or Symbol.")
        obj = Basic.__new__(cls, sym, *args, **kwargs)
        obj.check(*obj.args[1:])
        return IndexedRandomSymbol(obj.symbol, obj)

    symbol = property(lambda self: self.args[0])
    name = property(lambda self: self.symbol.name)
    index_set = property(lambda self: self.index_set)
    state_space = property(lambda self: self.state_space)

    @staticmethod
    def check(*args):
        pass

    def joint_dist(self, n):
        raise NotImplementedError()

    def __getitem__(self, key):
        if key not in self.index_set:
            raise KeyError("The requested key is not in index set of %s"%(str(self)))
        return Indexed(self, key)

    def first_passage_time(self, state):
        raise NotImplementedError()

    def probability(self, *args, **kwargs):
        raise NotImplementedError()

class MarkovChain(StochasticProcess):

    trans_matrix = property(lambda self: self.args[1])
    index_set = S.Naturals0
    state_space = property(lambda self: self.args[2])

    def probability(self, condition, given):
        if not isinstance(given, (Indexed, Matrix)):
            raise ValueError("Expecting Indexed or Matrix, got %s"%(given))
        if isinstance(given, Matrix):
            return given * self.trans_matrix ** (kc - kg)
        kc, kg = condition.args[1], given.args[1]
        if kc < kg:
            raise ValueError("Cannot find probability of previous time step using"
                                " the future time step.")
        u = MatrixSymbol(given.name, 1, len(self.state_space))
        return u * (self.trans_matrix)**(kc - kg)

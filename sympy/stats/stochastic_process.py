from __future__ import print_function, division

from sympy import Basic, Symbol
from sympy.core.compatibility import string_types
from sympy.stats.rv import ProductDomain
from sympy.stats.joint_rv import ProductPSpace, JointRandomSymbol

class StochasticPSpace(ProductPSpace):
    """
    Represents probability space of stochastic processes
    and their random variables. Contains mechanics to do
    computations for queries of stochastic processes.

    Initialized by symbol and the specific process.
    """

    def __new__(cls, sym, process):
        if isinstance(sym, string_types):
            sym = Symbol(sym)
        if not isinstance(sym, Symbol):
            raise TypeError("Name of stochastic process should be either only "
                            "a string or Symbol.")
        from sympy.stats.stochastic_process_types import StochasticProcess
        if not isinstance(process, StochasticProcess):
            raise TypeError("`process` must be an instance of StochasticProcess.")
        return Basic.__new__(cls, sym, process)

    @property
    def process(self):
        """
        The associated stochastic process.
        """
        return self.args[1]

    @property
    def domain(self):
        return ProductDomain(self.process.index_set,
                             self.process.state_space)

    @property
    def symbol(self):
        return self.args[0]

    def probability(self, condition, given_condition=None, **kwargs):
        """
        Transfers the task of handling queries to the specific stochastic
        process because every process has their own logic of handling such
        queries.
        """
        return self.process.probability(condition, given_condition, **kwargs)

    def joint_distribution(self, *args):
        """
        Computes the joint distribution of the random indexed variables.

        Parameters
        ==========

        args: iterable
            The finite list of random indexed variables of a stochastic
            process whose joint distribution has to be computed.

        Returns
        =======

        JointDistribution
            The joint distribution of the list of random indexed variables.
            An unevaluated object is returned if it is not possible to
            compute the joint distribution.

        Raises
        ======

        ValueError: When the time/key of random indexed variables
                    is not in strictly increasing order.
        """
        NotImplementedError()

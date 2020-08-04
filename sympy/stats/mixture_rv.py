from sympy import Basic, ImmutableMatrix, S, Symbol, Lambda
from sympy.core.sympify import sympify
from sympy.stats.rv import _symbol_converter, RandomSymbol, NamedArgsMixin, PSpace
from sympy.stats.crv import SingleContinuousPSpace
from sympy.stats.frv import SingleFinitePSpace
from sympy.stats.drv import SingleDiscretePSpace
from sympy.stats.crv_types import ContinuousDistributionHandmade
from sympy.stats.drv_types import DiscreteDistributionHandmade
from sympy.stats.frv_types import FiniteDistributionHandmade


class MixturePSpace(PSpace):
    """
    A temporary Probability Space for the Mixture Distribution
    """

    def __new__(cls, s, distribution):
        s = _symbol_converter(s)
        if not isinstance(distribution, MixtureDistribution):
            raise ValueError("%s should be an isinstance of "
                        "MixtureDistribution"%(distribution))
        return Basic.__new__(cls, s, distribution)

    @property
    def value(self):
        return RandomSymbol(self.symbol, self)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def distribution(self):
        return self.args[1]

    @property
    def pdf(self):
        return self.distribution.pdf(self.symbol)

    @property
    def set(self):
        return self.distribution.set

    @property
    def domain(self):
        return self._get_newpspace().domain

    def _get_newpspace(self):
        new_pspace = self._transform_pspace()
        if new_pspace is not None:
            return new_pspace
        message = ("Mixture Distribution for %s is not implemeted yet" % str(self.distribution))
        raise NotImplementedError(message)

    def _transform_pspace(self):
        """
        This function returns the new pspace of the distribution using handmade
        Distributions and their corresponding pspace.
        """
        pdf = Lambda(self.symbol, self.pdf)
        _set = self.distribution.set
        if self.distribution.is_Continuous:
            return SingleContinuousPSpace(self.symbol, ContinuousDistributionHandmade(pdf, _set))
        elif self.distribution.is_Discrete:
            return SingleDiscretePSpace(self.symbol, DiscreteDistributionHandmade(pdf, _set))
        elif self.distribution.is_Finite:
            dens = dict((k, pdf(k)) for k in _set)
            return SingleFinitePSpace(self.symbol, FiniteDistributionHandmade(dens))

    def compute_density(self, expr, **kwargs):
        new_pspace = self._get_newpspace()
        expr = expr.subs({self.value: new_pspace.value})
        return new_pspace.compute_density(expr, **kwargs)

    def compute_cdf(self, expr, **kwargs):
        new_pspace = self._get_newpspace()
        expr = expr.subs({self.value: new_pspace.value})
        return new_pspace.compute_cdf(expr, **kwargs)

    def compute_expectation(self, expr, rvs=None, evaluate=False, **kwargs):
        new_pspace = self._get_newpspace()
        expr = expr.subs({self.value: new_pspace.value})
        if isinstance(new_pspace, SingleFinitePSpace):
            return new_pspace.compute_expectation(expr, rvs, **kwargs)
        return new_pspace.compute_expectation(expr, rvs, evaluate, **kwargs)

    def probability(self, condition, **kwargs):
        new_pspace = self._get_newpspace()
        condition = condition.subs({self.value: new_pspace.value})
        return new_pspace.probability(condition)

    def conditional_space(self, condition, **kwargs):
        new_pspace = self._get_newpspace()
        condition = condition.subs({self.value: new_pspace.value})
        return new_pspace.conditional_space(condition)

    def sample(self, size=(), library='scipy'):
        new_pspace = self._get_newpspace()
        samp = new_pspace.sample(size, library)
        return {self.value: samp[new_pspace.value]}

def rv(symbol, cls, args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    pspace = MixturePSpace(symbol, dist)
    return pspace.value


class MixtureDistribution(Basic, NamedArgsMixin):
    """Represents the Mixture distribution"""
    _argnames = ('wts', 'rvs')

    def __new__(cls, wts, rvs):
        wts = ImmutableMatrix(wts)
        rvs = ImmutableMatrix(rvs)
        return Basic.__new__(cls, wts, rvs)

    @property
    def set(self):
        return (self.rvs[0, 0]).pspace.distribution.set

    @property
    def is_Continuous(self):
        return (self.rvs[0, 0]).pspace.is_Continuous

    @property
    def is_Discrete(self):
        return (self.rvs[0, 0]).pspace.is_Discrete

    @property
    def is_Finite(self):
        return (self.rvs[0, 0]).pspace.is_Finite

    @staticmethod
    def check(wts, rvs):
        rvs, wts = list(rvs), list(wts)
        set_ = rvs[0].pspace.domain.set
        for rv in rvs:
            if not isinstance(rv, RandomSymbol):
                raise TypeError("Each of element should be a random variable")
            if rv.pspace.domain.set != set_:
                raise ValueError("Each random variable should be defined on same set")
        for wt in wts:
            if not wt.is_positive:
                raise ValueError("Weight of each random variable should be positive")
        if len(rvs) != len(wts):
            raise ValueError("Weights and RVs should be of same length")
        if sum(wts) != S.One:
            raise ValueError("Sum of the weights should be 1")

    def pdf(self, x):
        y = Symbol('y')
        rvs, wts = list(self.rvs), list(self.wts)
        pdf_ = S.Zero
        if self.is_Finite:
            for rv in range(len(rvs)):
                pdf_ = pdf_ + wts[rv]*rvs[rv].pspace.distribution.pmf(y)
        else:
            for rv in range(len(rvs)):
                pdf_ = pdf_ + wts[rv]*rvs[rv].pspace.distribution.pdf(y)
        return Lambda(y, pdf_)(x)

def Mixture(name, wts, rvs):
    """Creates a random variable with mixture distribution"""
    return rv(name, MixtureDistribution, (wts, rvs))

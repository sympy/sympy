from sympy import Basic, S, Symbol, Lambda, Tuple
from sympy.core.sympify import sympify
from sympy.stats.rv import _symbol_converter, RandomSymbol, NamedArgsMixin, PSpace
from sympy.stats.crv import SingleContinuousPSpace, ContinuousDistribution
from sympy.stats.frv import SingleFinitePSpace, SingleFiniteDistribution
from sympy.stats.drv import SingleDiscretePSpace, DiscreteDistribution
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
    _argnames = ('wts', 'dists')

    def __new__(cls, wts, dists):
        wts = Tuple(*list(wts))
        dists = Tuple(*list(dists))
        return Basic.__new__(cls, wts, dists)

    @property
    def set(self):
        return self.dists[0].set

    @property
    def is_Continuous(self):
        return isinstance(self.dists[0], ContinuousDistribution)

    @property
    def is_Discrete(self):
        return isinstance(self.dists[0], DiscreteDistribution)

    @property
    def is_Finite(self):
        return isinstance(self.dists[0], SingleFiniteDistribution)

    @staticmethod
    def check(wts, dists):
        dists, wts = list(dists), list(wts)
        set_ = dists[0].set
        if isinstance(dists[0], ContinuousDistribution):
            parent_dist = ContinuousDistribution
        elif isinstance(dists[0], DiscreteDistribution):
            parent_dist = DiscreteDistribution
        else:
            parent_dist = SingleFiniteDistribution
        for dist in dists:
            if not isinstance(dist, parent_dist):
                raise TypeError("Each of distribution should be an instance of %s"
                                % str(parent_dist))
            if dist.set != set_:
                raise ValueError("Each distribution should be defined on same set")
        for wt in wts:
            if not wt.is_positive:
                raise ValueError("Weight of each random variable should be positive")
        if len(dists) != len(wts):
            raise ValueError("Weights and RVs should be of same length")

    def pdf(self, x):
        y = Symbol('y')
        dists, wts = list(self.dists), list(self.wts)
        pdf_ = S.Zero
        tot_wt = sum(wts)
        wts = [wt/S(tot_wt) for wt in wts]
        if self.is_Finite:
            for dist in range(len(dists)):
                pdf_ = pdf_ + wts[dist]*dists[dist].pmf(y)
        else:
            for dist in range(len(dists)):
                pdf_ = pdf_ + wts[dist]*dists[dist].pdf(y)
        return Lambda(y, pdf_)(x)

def Mixture(name, wts, rvs):
    """Creates a random variable with mixture distribution"""
    return rv(name, MixtureDistribution, (wts, rvs))

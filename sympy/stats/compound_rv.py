from sympy import Basic, integrate, Sum, Dummy, Lambda
from sympy.stats.rv import (NamedArgsMixin, random_symbols, _symbol_converter,
                        PSpace, RandomSymbol, is_random)
from sympy.stats.crv import ContinuousDistribution, SingleContinuousPSpace
from sympy.stats.drv import DiscreteDistribution, SingleDiscretePSpace
from sympy.stats.frv import SingleFiniteDistribution, SingleFinitePSpace
from sympy.stats.crv_types import ContinuousDistributionHandmade
from sympy.stats.drv_types import DiscreteDistributionHandmade
from sympy.stats.frv_types import FiniteDistributionHandmade


class CompoundPSpace(PSpace):
    """
    A temporary Probability Space for the Compound Distribution. After
    Marginalization, this returns the corresponding Probability Space of the
    parent distribution.
    """

    is_Finite = None
    is_Continuous = None
    is_Discrete = None

    def __new__(cls, s, distribution):
        s = _symbol_converter(s)
        if isinstance(distribution, ContinuousDistribution):
            return SingleContinuousPSpace(s, distribution)
        if isinstance(distribution, DiscreteDistribution):
            return SingleDiscretePSpace(s, distribution)
        if isinstance(distribution, SingleFiniteDistribution):
            return SingleFinitePSpace(s, distribution)
        if not isinstance(distribution, CompoundDistribution):
            raise ValueError("%s should be an isinstance of "
                        "CompoundDistribution"%(distribution))
        cls.is_Finite = distribution.is_Finite
        cls.is_Continuous = distribution.is_Continuous
        cls.is_Discrete = distribution.is_Discrete
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
        x = Dummy('x')
        parent_dist = self.distribution.args[0]
        new_pspace = self._transform_pspace(self.symbol, parent_dist,
                            Lambda(x, self.distribution.pdf(x)))
        if new_pspace is not None:
            return new_pspace
        message = ("Compound Distribution for %s is not implemeted yet" % str(parent_dist))
        raise NotImplementedError(message)

    def _transform_pspace(self, sym, dist, pdf):
        """
        This function returns the new pspace of the distribution using handmade
        Distributions and their corresponding pspace.
        """
        pdf = Lambda(sym, pdf(sym))
        _set = dist.set
        if isinstance(dist, ContinuousDistribution):
            return SingleContinuousPSpace(sym, ContinuousDistributionHandmade(pdf, _set))
        elif isinstance(dist, DiscreteDistribution):
            return SingleDiscretePSpace(sym, DiscreteDistributionHandmade(pdf, _set))
        elif isinstance(dist, SingleFiniteDistribution):
            dens = dict((k, pdf(k)) for k in _set)
            return SingleFinitePSpace(sym, FiniteDistributionHandmade(dens))

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
        if rvs:
            rvs = rvs.subs({self.value: new_pspace.value})
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


class CompoundDistribution(Basic, NamedArgsMixin):
    """
    Class for Compound Distributions.

    Parameters
    ==========

    dist : Distribution
        Distribution must contain a random parameter

    Examples
    ========

    >>> from sympy.stats.compound_rv import CompoundDistribution
    >>> from sympy.stats.crv_types import NormalDistribution
    >>> from sympy.stats import Normal
    >>> from sympy.abc import x
    >>> X = Normal('X', 2, 4)
    >>> N = NormalDistribution(X, 4)
    >>> C = CompoundDistribution(N)
    >>> C.set
    Interval(-oo, oo)
    >>> C.pdf(x).simplify()
    exp(-x**2/64 + x/16 - 1/16)/(8*sqrt(pi))

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Compound_probability_distribution

    """

    is_Finite = None
    is_Continuous = None
    is_Discrete = None

    def __new__(cls, dist):
        if isinstance(dist, ContinuousDistribution):
            cls.is_Continuous = True
        elif isinstance(dist, DiscreteDistribution):
            cls.is_Discrete = True
        elif isinstance(dist, SingleFiniteDistribution):
            cls.is_Finite = True
        else:
            message = "Compound Distribution for %s is not implemeted yet" % str(dist)
            raise NotImplementedError(message)
        if not cls._compound_check(dist):
            return dist
        return Basic.__new__(cls, dist)

    @property
    def set(self):
        return self.args[0].set

    def pdf(self, x):
        dist = self.args[0]
        randoms = [rv for rv in dist.args if is_random(rv)]
        if isinstance(dist, SingleFiniteDistribution):
            y = Dummy('y', integer=True, negative=False)
            expr = dist.pmf(y)
        else:
            y = Dummy('y')
            expr = dist.pdf(y)
        for rv in randoms:
            expr = self._marginalise(expr, rv)
        return Lambda(y, expr)(x)

    def _marginalise(self, expr, rv):
        if isinstance(rv.pspace.distribution, SingleFiniteDistribution):
            rv_dens = rv.pspace.distribution.pmf(rv)
        else:
            rv_dens = rv.pspace.distribution.pdf(rv)
        rv_dom = rv.pspace.domain.set
        if rv.pspace.is_Discrete or rv.pspace.is_Finite:
            expr = Sum(expr*rv_dens, (rv, rv_dom._inf,
                    rv_dom._sup)).doit()
        else:
            expr = integrate(expr*rv_dens, (rv, rv_dom._inf,
                    rv_dom._sup))
        return expr

    @classmethod
    def _compound_check(self, dist):
        """
        Checks if the given distribution contains random parameters.
        """
        randoms = []
        for arg in dist.args:
            randoms.extend(random_symbols(arg))
        if len(randoms) == 0:
            return False
        return True

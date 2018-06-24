"""
Joint Random Variables Module

See Also
========
sympy.stats.rv
sympy.stats.frv
sympy.stats.crv
sympy.stats.drv
"""

from __future__ import print_function, division

# __all__ = ['marginal_distribution']

from sympy import Basic, Lambda, sympify
from sympy.concrete.summations import Sum
from sympy.integrals.integrals import Integral, integrate
from sympy.stats.rv import (ProductPSpace, NamedArgsMixin,
     ProductDomain)
from sympy.core.containers import Tuple

class JointPSpace(ProductPSpace):
    """
    Represents a joint probability space. Represented using symbols for
    each component and a distribution.
    """
    def __new__(cls, syms, dist, *args):
        syms = Tuple.fromiter(i for i in syms)
        return Basic.__new__(cls, syms, dist, *args)

    @property
    def set(self):
        return self.domain.set

    @property
    def symbols(self):
        return self.args[0]

    @property
    def distribution(self):
        return self.args[1]

    def pdf(self, *args):
        return self.distribution(*args)

    def marginal_distribution(self, indices, *sym):
        limits = list([i,] for i in self.symbols if i not in sym)
        for i in indices:
            limits[i].append(self.distribution.set.args[i])
            limits[i] = tuple(limits[i])
        limits = tuple(limits)
        if self.distribution.is_Continuous:
            return Lambda(sym, integrate(self.distribution(
                *self.symbols), limits))
        if self.distribution.is_Discrete:
            return Lambda(sym, integrate(self.distribution(
                *self.symbols), limits))

    @property
    def values(self):
        raise NotImplementedError

    def where(self, condition):
        raise NotImplementedError()

    def compute_density(self, expr):
        raise NotImplementedError()

    def sample(self):
        raise NotImplementedError()

    def probability(self, condition):
        raise NotImplementedError()

class JointDistribution(Basic, NamedArgsMixin):
    """
    Represented by the random variables part of the joint distribution.
    Contains methods for PDF, CDF, sampling, marginal densities, etc.
    """

    def __new__(cls, *args):
        args = list(map(sympify, args))
        for i in range(len(args)):
            if isinstance(args[i], list):
                args[i] = Tuple.fromiter(j for j in args[i])
        return Basic.__new__(cls, *args)

    @property
    def domain(self):
        return ProductDomain(self.symbols)

    @property
    def pdf(self, *args):
        return self.density.args[1]

    def cdf(self, other):
        assert isinstance(other, dict)
        rvs = other.keys()
        _set = self.domain.set
        expr = self.pdf(tuple(i.args[0] for i in self.symbols))
        for i in len(other):
            if rvs[i].is_Continuous:
                density = Integral(expr, (rvs[i], _set[i].inf,
                    other[rvs[i]]))
            elif rvs[i].is_Discrete:
                density = Sum(expr, (rvs[i], _set[i].inf,
                    other[rvs[i]]))
        return density

    def __call__(self, *args):
        return self.pdf(*args)

def marginal_distribution(prob_space, *sym):
    """
    Marginal density of a joint random variable.

    Parameters:
    ==========

    prob_space: A joint probability space.
    sym: list of symbols for whom the marginal density is to be
    calculated.

    Returns:
    =======
    A Lambda expression n `sym`.

    Example:
    =======
    >>> from sympy.stats.joint_rv_types import MultivariateNormal
    >>> from sympy.stats.joint_rv import marginal_distribution
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> m = MultivariateNormal(('x', 'y'), [1, 2], [[2, 1], [1, 2]])
    >>> marginal_distribution(m, x)(1)
    1/(2*sqrt(pi))

    """
    pspace_sym = prob_space.symbols
    if sym == ():
        raise ValueError(
            "At least one symbol for marginal density is needed.")
    assert all([i in pspace_sym for i in sym])
    indices = tuple((pspace_sym).index(s) for s in sym)
    if hasattr(prob_space.distribution, 'marginal_distribution'):
        return prob_space.distribution.marginal_distribution(indices, *sym)
    return prob_space.marginal_distribution(indices, *sym)

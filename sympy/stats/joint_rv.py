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

__all__ = ['marginal_density']

from sympy import Basic, Lambda, sympify, S
from sympy.concrete.summations import Sum
from sympy.integrals.integrals import Integral
from sympy.stats.rv import (ProductPSpace, RandomSymbol,
        NamedArgsMixin, ProductDomain)

class JointPSpace(ProductPSpace):
    """
    Represents a joint probability space. Represented using symbols for
    each component and a distribution.
    """
    def __new__(cls, syms, dist, **args):
        syms = tuple(map(sympify, syms))
        return Basic.__new__(cls, syms, dist, **args)

    @property
    def set(self):
        return self.domain.set

    @property
    def symbols(self):
        return tuple(self.args[0])

    @property
    def values(self):
        vls = []
        for sym in self.symbols:
            vls.extend(MarginalPSpace((sym,), self).values)
        return tuple(vls)

    @property
    def distribution(self):
        return self.args[1]

    def pdf(self, *args):
        return self.distribution(*args)

    def where(self, condition):
        raise NotImplementedError()

    def compute_density(self, expr):
        raise NotImplementedError()

    def sample(self):
        raise NotImplementedError()

    def probability(self, condition):
        raise NotImplementedError()


class MarginalPSpace(JointPSpace):
    """
    Space representing the marginal probaility space of a subset
    of random variables from a joint probability space.
    """
    def __new__(cls, sym, jpspace):
        return Basic.__new__(cls, sym, jpspace)

    @property
    def joint_space(self):
        return self.args[1]

    @property
    def symbols(self):
        return self.args[0]

    @property
    def density(self, *sym):
        space = self.joint_space
        sym_2 = tuple(i for i in space.symbols if i not in sym)
        if space.distributions.is_Continuous:
            limits = tuple((i, S.Reals) for i in self.symbols)
            return Lambda(sym_2, Integral(space.distribution(
                space.symbols), limits))
        if space.distributions.is_Discrete:
            limits = tuple((i, S.Integers) for i in self.symbols)
            return Lambda(sym_2, Integral(space.distribution(
                space.symbols), limits))

    @property
    def values(self):
        return tuple(RandomSymbol(sym, self) for sym in self.symbols)

class JointDistribution(Basic, NamedArgsMixin):
    """
    Represented by the random variables part of the joint distribution.
    Contains methods for PDF, CDF, sampling, marginal densities, etc.
    """

    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @property
    def domain(self):
        return ProductDomain(self.symbols)

    @property
    def symbols(self):
        return list(map(sympify, self.args[0]))

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

def marginal_density(prob_space, *sym):
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
    >>> from sympy.stats.joint_rv import marginal_density
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> m = MultivariateNormal(('x', 'y'), [1, 2], [[2, 1], [1, 2]])
    >>> marginal_density(m, x)(1)
    1/(2*sqrt(pi))

    """
    pspace_sym = prob_space.symbols
    if sym == ():
        raise ValueError(
            "At least one symbol for marginal density is needed.")
    indices = tuple((pspace_sym).index(s) for s in sym)
    if hasattr(prob_space.distribution, 'marginal_density'):
        return prob_space.distribution.marginal_density(indices, *sym)
    return MarginalPSpace(prob_space, *sym).density

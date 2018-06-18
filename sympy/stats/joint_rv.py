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

__all__ = ['Joint']

from sympy import Basic, Lambda, sympify
from sympy.core.symbol import Symbol
from sympy.core.compatibility import string_types
from sympy.concrete.summations import Sum
from sympy.integrals.integrals import Integral
from sympy.stats.rv import (ProductPSpace, RandomSymbol, pspace,
        NamedArgsMixin, ProductDomain)

class JointPSpace(ProductPSpace):
    """

    """
    def __new__(cls, name, syms, dist, **args):
        name = sympify(name)
        return Basic.__new__(cls, name, syms, dist, **args)

    @property
    def set(self):
        return self.domain.set

    @property
    def symbol(self):
        return self.args[0]

    @property
    def symbols(self):
        return tuple(self.args[1])

    @property
    def value(self):
        return RandomSymbol(self.symbol, self)

    @property
    def distribution(self):
        return self.args[2]

    # @property
    # def spaces(self):
    #     rvs = self.distribution.symbols
    #     return FiniteSet(*[rv.pspace for rv in rvs])

    def pdf(self, *args):
        return self.distribution(*args)


    # def marginal(self, rv):
    #     assert rv in self.distribution.symbols
    #     return self.distribution.marginal_density(rv)

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
        return Basic.__new__(cls, *args)

    @property
    def domain(self):
        return ProductDomain(self.symbols)

    @property
    def symbols(self):
        return list(self.args[0])

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

    def marginalise_out(self, expr, other, evaluate=False):
        if pspace(other).is_Continuous:
            result = Integral(expr, (other.args[0],
                other.pspace.set.inf, other.pspace.set.sup))
        elif pspace(other).is_Discrete:
            result = Sum(expr, (other.args[0],
                other.set.inf, other.set.sup))
        return result if not evaluate else result.doit()

    def marginal_density(self, rv, evaluate=False):
        if (rv not in self.symbols):
            raise ValueError("%s is not a part of the joint distribtion" %rv)
        if isinstance(rv, RandomSymbol):
            return Lambda(rv.symbol, pspace(rv).pdf)
        rvs = self.symbols
        rvs.remove(rv)
        expr = self.pdf
        for i in rvs:
            expr = self.marginalise_out(expr, i, evaluate)
            rvs.remove(i)
        return Lambda([i.symbol for i in rvs    ], expr)

    def __call__(self, *args):
        return self.pdf(*args)


def Joint(name, rvs, density=None):
    """
    Create a discrete random variable with a Joint probability distribution.

    The density of the distribution may/may not be provided by the user.

    Parameters
    ==========

    rvs: a tuple/list of Random symbols or symbols to be used as the constituents
    of the joint distribution.

    density: The probability density function of the joint distribution.
    This must be provided in case at least one of `rvs` is a `symbol`. The
    `_check` method of JointDistribution checks if the CDF over the entire
    domain evaluates to one, and raises an error if the condition evaluates
    to `False`.

    Returns
    =======

    A RandomSymbol.

    """
    if isinstance(name, string_types):
        name = Symbol(name)
    if not isinstance(name, Symbol):
        raise TypeError("Name must be a string or a symbol.")
    assert all([isinstance(rv, (RandomSymbol)) for rv in rvs])
    dist = JointDistribution(rvs, density)
    return JointPSpace(name, dist).value

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

from sympy import Basic, Eq, Mul, Lambda, S
from sympy.core.symbol import Dummy, Symbol
from sympy.core.compatibility import string_types
from sympy.concrete.summations import summation
from sympy.integrals.integrals import Integral, integrate
from sympy.stats.rv import (PSpace, RandomDomain, random_symbols,
        NamedArgsMixin, RandomSymbol, density)
from sympy.utilities.misc import filldedent
from sympy.sets.contains import Contains

class JointDomain(RandomDomain):
    """
    Represents the joint domain as a ProductSet.
    """
    @property
    def set(self):
        rv_sets = [rv.args[1].domain.set for rv in self.symbols]
        _set = rv_sets[0]
        for i in rv_sets[1:]:
            _set = _set*i
        return _set

    def as_boolean(self):
        return Contains(self.symbols, self.set)

    def __contains__(self, other):
        return other in self.set()

class JointPSpace(PSpace):
    """

    """
    @property
    def set(self):
        return self.domain.set

    @property
    def symbol(self):
        return self.args[0]

    @property
    def value(self):
        return RandomSymbol(self.symbol, self)

    @property
    def distribution(self):
        return self.args[1]

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

    def __new__(cls, _symbols, _density):
        if (_density!=None):
            jrv = Basic.__new__(cls, _symbols, _density)
            jrv._check()
        else:
            syms = [i.args[0] for i in _symbols]
            _density = Mul(*[density(_symbols[i])(syms[i]) for i in
                range(len(_symbols))])
            jrv = Basic.__new__(cls, _symbols, _density)
        return jrv

    def _check(self):
        x = Symbol('x', real=True)
        for sym in self.symbols:
            if isinstance(sym, RandomSymbol):
                assert Eq(density(sym)(x), self.marginal_density(sym)(x)) != False
        cdf_lim = {i:S.Infinity for i in self.symbols}
        if (self.cdf(cdf_lim) == 1) == False:
            raise ValueError("The density entered is invalid for the given variables")

    @property
    def domain(self):
        return JointDomain(self.symbols)

    @property
    def symbols(self):
        return list(self.args[0])

    @property
    def pdf(self):
        rvs = self.symbols
        z = [Dummy(str(i.args[0])) for i in rvs]
        _density = self.args[1]
        for i in range(len(rvs)):
            _density = _density.replace(rvs[i].symbol, z[i])
        return Lambda(z, _density)

    def cdf(self, other):
        assert isinstance(other, dict)
        rvs = other.keys()
        #all_syms = [Dummy(i.args[0]) for i in rvs]
        # _symbols = other.keys
        _set = self.domain.set
        expr = self.pdf(tuple(i.args[0] for i in self.symbols))
        for i in len(other):
            if rvs[i].is_Continuous:
                density = Integral(expr, (rvs[i], _set[i].inf,
                    other[rvs[i]]))
            elif rvs[i].is_Discrete:
                density = summation(expr, (rvs[i], _set[i].inf,
                    other[rvs[i]]))
        return density

    def marginalise_out(self, expr, other, evaluate):
        z = [Dummy(str(i.args[0])) for i in self.symbols]
        if evaluate:
            return Lambda(z, integrate(expr, (other.args[0],
             other.domain.set.inf, other.domain.set.sup)))
        else:
            return Lambda(z, Integral(expr, (other.args[0],
             other.set.inf, other.set.sup)))


    def marginal_density(self, rv, evaluate=True):
        if (rv not in self.symbols):
            raise ValueError("%s is not a part of the joint distribtion" %rv)
        rvs = self.symbols
        rvs.remove(rv)
        expr = self.pdf
        for i in rvs:
            expr = self.marginalise_out(expr, i, evaluate)
            rvs.remove(i)
        return expr

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
    assert all([isinstance(rv, (RandomSymbol, Symbol)) for rv in rvs])
    if any([isinstance(rv, Symbol) for rv in rvs]) and (density == None):
        raise ValueError(filldedent('''Density must be provided for user
            defined random distributions.'''))
    dist = JointDistribution(rvs, density)
    return JointPSpace(name, dist).value

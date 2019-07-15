"""
Compound Random Variables Module

See Also
========
sympy.stats.rv
sympy.stats.frv
sympy.stats.crv
sympy.stats.drv
sympy.stats.joint_rv
"""

from __future__ import print_function, division

from sympy import (Basic, Lambda, sympify, Indexed, Symbol, ProductSet, S,
                   Dummy)
from sympy.concrete.summations import Sum, summation
from sympy.concrete.products import Product
from sympy.core.compatibility import string_types, iterable
from sympy.core.containers import Tuple
from sympy.integrals.integrals import Integral, integrate
from sympy.matrices import ImmutableMatrix
from sympy.stats.crv import (ContinuousDistribution,
                             SingleContinuousDistribution, SingleContinuousPSpace)
from sympy.stats.drv import (DiscreteDistribution,
                             SingleDiscreteDistribution, SingleDiscretePSpace)
from sympy.stats.rv import (ProductPSpace, NamedArgsMixin,
                            ProductDomain, RandomSymbol, random_symbols, SingleDomain)
from sympy.utilities.misc import filldedent


class CompoundPSpace(ProductPSpace):
    """
    Represents a compound probability space. Represented using a distribution
    and parameter(s) of which one or more are distributed
    according to some distribution.
    """

    def __new__(cls, sym, dist):
        if isinstance(dist, SingleContinuousDistribution):
            return SingleContinuousPSpace(sym, dist)
        if isinstance(dist, SingleDiscreteDistribution):
            return SingleDiscretePSpace(sym, dist)
        if isinstance(sym, string_types):
            sym = Symbol(sym)
        if not isinstance(sym, Symbol):
            raise TypeError("s should have been string or Symbol")
        return Basic.__new__(cls, sym, dist)

    @property
    def set(self):
        return self.domain.set

    @property
    def symbol(self):
        return self.args[0]

    @property
    def distribution(self):
        return self.args[1]

    @property
    def component_count(self):
        _set = self.distribution.set
        if isinstance(_set, ProductSet):
            return S(len(_set.args))
        elif isinstance(_set, Product):
            return _set.limits[0][-1]
        return S(1)

    @property
    def pdf(self):
        sym = [Indexed(self.symbol, i) for i in range(self.component_count)]
        return self.distribution(*sym)

    @property
    def domain(self):
        rvs = random_symbols(self.distribution)
        if not rvs:
            return SingleDomain(self.symbol, self.distribution.set)
        return ProductDomain(*[rv.pspace.domain for rv in rvs])

    def component_domain(self, index):
        return self.set.args[index]

    def marginal_distribution(self, *indices):
        count = self.component_count
        if count.atoms(Symbol):
            raise ValueError("Marginal distributions cannot be computed "
                             "for symbolic dimensions. It is a work under progress.")
        orig = [Indexed(self.symbol, i) for i in range(count)]
        all_syms = [Symbol(str(i)) for i in orig]
        replace_dict = dict(zip(all_syms, orig))
        sym = [Symbol(str(Indexed(self.symbol, i))) for i in indices]
        limits = list([i, ] for i in all_syms if i not in sym)
        index = 0
        for i in range(count):
            if i not in indices:
                limits[index].append(self.distribution.set.args[i])
                limits[index] = tuple(limits[index])
                index += 1
        if self.distribution.is_Continuous:
            f = Lambda(sym, integrate(self.distribution(*all_syms), *limits))
        elif self.distribution.is_Discrete:
            f = Lambda(sym, summation(self.distribution(*all_syms), *limits))
        return f.xreplace(replace_dict)

    def compute_expectation(self, expr, rvs=None, evaluate=False, **kwargs):
        syms = tuple(self.value[i] for i in range(self.component_count))
        rvs = rvs or syms
        if not any([i in rvs for i in syms]):
            return expr
        expr = expr * self.pdf
        for rv in rvs:
            if isinstance(rv, Indexed):
                expr = expr.xreplace({rv: Indexed(str(rv.base), rv.args[1])})
            elif isinstance(rv, RandomSymbol):
                expr = expr.xreplace({rv: rv.symbol})
        if self.value in random_symbols(expr):
            raise NotImplementedError(filldedent('''
            Expectations of expression with unindexed joint random symbols
            cannot be calculated yet.'''))
        limits = tuple((Indexed(str(rv.base), rv.args[1]),
                        self.distribution.set.args[rv.args[1]]) for rv in syms)
        return Integral(expr, *limits)

    def where(self, condition):
        raise NotImplementedError()

    def compute_density(self, expr):
        raise NotImplementedError()

    def sample(self):
        raise NotImplementedError()

    def probability(self, condition):
        raise NotImplementedError()


class CompoundDistribution(Basic, NamedArgsMixin):
    """
    Represents a compound probability distribution.

    Constructed using a single probability distribution with a parameter
    distributed according to some given distribution.
    """

    def __new__(cls, dist):
        if not isinstance(dist, (ContinuousDistribution, DiscreteDistribution)):
            raise ValueError(filldedent('''CompoundDistribution can only be
             initialized from ContinuousDistribution or DiscreteDistribution
             '''))
        _args = dist.args
        if not any([isinstance(i, RandomSymbol) for i in _args]):
            return dist
        return Basic.__new__(cls, dist)

    @property
    def latent_distributions(self):
        return random_symbols(self.args[0])

    def pdf(self, *x):
        dist = self.args[0]
        z = Dummy('z')
        if isinstance(dist, ContinuousDistribution):
            rv = SingleContinuousPSpace(z, dist).value
        elif isinstance(dist, DiscreteDistribution):
            rv = SingleDiscretePSpace(z, dist).value
        return MarginalDistribution(self, (rv,)).pdf(*x)

    def set(self):
        return self.args[0].set

    def __call__(self, *args):
        return self.pdf(*args)

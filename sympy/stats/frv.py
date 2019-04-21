"""
Finite Discrete Random Variables Module

See Also
========
sympy.stats.frv_types
sympy.stats.rv
sympy.stats.crv
"""
from __future__ import print_function, division

from itertools import product

from sympy import (Basic, Symbol, cacheit, sympify, Mul,
        And, Or, Tuple, Piecewise, Eq, Lambda, exp, I, Dummy)
from sympy.sets.sets import FiniteSet
from sympy.stats.rv import (RandomDomain, ProductDomain, ConditionalDomain,
        PSpace, IndependentProductPSpace, SinglePSpace, random_symbols,
        sumsets, rv_subs, NamedArgsMixin)
from sympy.core.containers import Dict
import random

class FiniteDensity(dict):
    """
    A domain with Finite Density.

    """
    def __call__(self, item):
        """
        Make instance of a class callable.

        If item belongs to current instance of a class, return it.

        Otherwise, return 0.
        """
        item = sympify(item)
        if item in self:
            return self[item]
        else:
            return 0

    @property
    def dict(self):
        """Return item as dictionary"""
        return dict(self)

class FiniteDomain(RandomDomain):
    """
    A domain with discrete finite support

    Represented using a FiniteSet.
    """
    is_Finite = True

    @property
    def symbols(self):
        """Return set of symbols valid in the domain"""
        return FiniteSet(sym for sym, val in self.elements)

    @property
    def elements(self):
        """Return elements of the domain"""
        return self.args[0]

    @property
    def dict(self):
        """Return a dictionary of elements"""
        return FiniteSet(*[Dict(dict(el)) for el in self.elements])

    def __contains__(self, other):
        return other in self.elements

    def __iter__(self):
        """Return an iterator"""
        return self.elements.__iter__()

    def as_boolean(self):
        return Or(*[And(*[Eq(sym, val) for sym, val in item]) for item in self])


class SingleFiniteDomain(FiniteDomain):
    """
    A FiniteDomain over a single symbol/set

    Example
    =======

    The possibilities of a *single* die roll.
    """

    def __new__(cls, symbol, set):
        """
        Create  new Object of the SingleFiniteDomain Class

        Returns
        =======

        A new object
        """
        if not isinstance(set, FiniteSet):
            set = FiniteSet(*set)
        return Basic.__new__(cls, symbol, set)

    @property
    def symbol(self):
        """Return a symbol"""
        return self.args[0]
        return tuple(self.symbols)[0]

    @property
    def symbols(self):
        """Reurn a set of all symbols"""
        return FiniteSet(self.symbol)

    @property
    def set(self):
        """Return a set for the Domain"""
        return self.args[1]

    @property
    def elements(self):
        """Return a set of elements"""
        return FiniteSet(*[frozenset(((self.symbol, elem), )) for elem in self.set])

    def __iter__(self):
        """Return an iterator"""
        return (frozenset(((self.symbol, elem),)) for elem in self.set)

    def __contains__(self, other):
        sym, val = tuple(other)[0]
        return sym == self.symbol and val in self.set


class ProductFiniteDomain(ProductDomain, FiniteDomain):
    """
    A Finite domain consisting of several other FiniteDomains

    Example
    =======

    The possibilities of the rolls of three independent dice
    """

    def __iter__(self):
        """Create an iteraror"""
        proditer = product(*self.domains)
        return (sumsets(items) for items in proditer)

    @property
    def elements(self):
        """Return set of elements"""
        return FiniteSet(*self)


class ConditionalFiniteDomain(ConditionalDomain, ProductFiniteDomain):
    """
    A FiniteDomain that has been restricted by a condition

    Example
    =======

    The possibilities of a die roll under the condition that the
    roll is even.
    """

    def __new__(cls, domain, condition):
        """
        Create new instance of ConditionalFiniteDomain class

        Raises
        ======

        ValueError if condition contains foreign symbols

        Returns
        =======

        A new object
        """
        if condition is True:
            return domain
        cond = rv_subs(condition)
        # Check that we aren't passed a condition like die1 == z
        # where 'z' is a symbol that we don't know about
        # We will never be able to test this equality through iteration
        if not cond.free_symbols.issubset(domain.free_symbols):
            raise ValueError('Condition "%s" contains foreign symbols \n%s.\n' % (
                condition, tuple(cond.free_symbols - domain.free_symbols)) +
                "Will be unable to iterate using this condition")

        return Basic.__new__(cls, domain, cond)



    def _test(self, elem):
        """
        Test the value.

        If value is boolean, return it. If value is equality
        relational (two objects are equal), return it with left-hand side
        being equal to right-hand side. Otherwise, raise ValueError exception.

        Raises
        ======

        ValueError  if lhs == rhs
        """
        val = self.condition.xreplace(dict(elem))
        if val in [True, False]:
            return val
        elif val.is_Equality:
            return val.lhs == val.rhs
        raise ValueError("Undeciable if %s" % str(val))

    def __contains__(self, other):
        return other in self.fulldomain and self._test(other)

    def __iter__(self):
        """Create an iterator for the domain"""
        return (elem for elem in self.fulldomain if self._test(elem))

    @property
    def set(self):
        """
        Returns a set for the Conditional Finite Domain

        Raises
        ======

        NotImplementedError

        Returns
        =======

        A finite domain set
        """
        if isinstance(self.fulldomain, SingleFiniteDomain):
            return FiniteSet(*[elem for elem in self.fulldomain.set
                               if frozenset(((self.fulldomain.symbol, elem),)) in self])
        else:
            raise NotImplementedError(
                "Not implemented on multi-dimensional conditional domain")

    def as_boolean(self):
        return FiniteDomain.as_boolean(self)

class SingleFiniteDistribution(Basic, NamedArgsMixin):
    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @property
    @cacheit
    def dict(self):
        """Return a dictionary of self.pdf"""
        return dict((k, self.pdf(k)) for k in self.set)

    @property
    def pdf(self):
        """PDF for Single Finite Distribution"""
        x = Symbol('x')
        return Lambda(x, Piecewise(*(
            [(v, Eq(k, x)) for k, v in self.dict.items()] + [(0, True)])))

    @property
    def characteristic_function(self):
        """CDF for Single Finite Distribution"""
        t = Dummy('t', real=True)
        return Lambda(t, sum(exp(I*k*t)*v for k, v in self.dict.items()))

    @property
    def moment_generating_function(self):
        """Moment Generating Function for Single Finite Distribution"""
        t = Dummy('t', real=True)
        return Lambda(t, sum(exp(k * t) * v for k, v in self.dict.items()))

    @property
    def set(self):
        """Return a list of the set"""
        return list(self.dict.keys())

    values = property(lambda self: self.dict.values)
    items = property(lambda self: self.dict.items)
    __iter__ = property(lambda self: self.dict.__iter__)
    __getitem__ = property(lambda self: self.dict.__getitem__)

    __call__ = pdf

    def __contains__(self, other):
        return other in self.set


#=============================================
#=========  Probability Space  ===============
#=============================================


class FinitePSpace(PSpace):
    """
    A Finite Probability Space

    Represents the probabilities of a finite number of events.
    """
    is_Finite = True

    def __new__(cls, domain, density):
        """Create new object of the Finite Probability Space Class"""
        density = dict((sympify(key), sympify(val))
                for key, val in density.items())
        public_density = Dict(density)

        obj = PSpace.__new__(cls, domain, public_density)
        obj._density = density
        return obj

    def prob_of(self, elem):
        """PDF for Finite Probability Space Class"""
        elem = sympify(elem)
        return self._density.get(elem, 0)

    def where(self, condition):
        """Position inside the Finite Probability Space"""
        assert all(r.symbol in self.symbols for r in random_symbols(condition))
        return ConditionalFiniteDomain(self.domain, condition)

    def compute_density(self, expr):
        """Probability density for Finite Probability Space"""
        expr = expr.xreplace(dict(((rs, rs.symbol) for rs in self.values)))
        d = FiniteDensity()
        for elem in self.domain:
            val = expr.xreplace(dict(elem))
            prob = self.prob_of(elem)
            d[val] = d.get(val, 0) + prob
        return d

    @cacheit
    def compute_cdf(self, expr):
        """Cumulative Distribution Function for Finite Probability Space"""
        d = self.compute_density(expr)
        cum_prob = 0
        cdf = []
        for key in sorted(d):
            prob = d[key]
            cum_prob += prob
            cdf.append((key, cum_prob))

        return dict(cdf)

    @cacheit
    def sorted_cdf(self, expr, python_float=False):
        """Compute the sorted CDF for Finite Probability Space"""
        cdf = self.compute_cdf(expr)
        items = list(cdf.items())
        sorted_items = sorted(items, key=lambda val_cumprob: val_cumprob[1])
        if python_float:
            sorted_items = [(v, float(cum_prob))
                    for v, cum_prob in sorted_items]
        return sorted_items

    @cacheit
    def compute_characteristic_function(self, expr):
        """Characteristic Functionm for Finite Probability Space"""
        d = self.compute_density(expr)
        t = Dummy('t', real=True)

        return Lambda(t, sum(exp(I*k*t)*v for k,v in d.items()))

    @cacheit
    def compute_moment_generating_function(self, expr):
        """
        Compute Moment Generating Function
        for the given Finite Probability Space

        Returns
        =======
        A lambda with moment function
        """
        d = self.compute_density(expr)
        t = Dummy('t', real=True)

        return Lambda(t, sum(exp(k * t) * v for k, v in d.items()))

    def compute_expectation(self, expr, rvs=None, **kwargs):
        """
        Compute Expectation Valuw
        for the given Finite Probability Space

        Returns
        =======
        A lambda with expectation value
        """
        rvs = rvs or self.values
        expr = expr.xreplace(dict((rs, rs.symbol) for rs in rvs))
        return sum([expr.xreplace(dict(elem)) * self.prob_of(elem)
                for elem in self.domain])

    def probability(self, condition):
        """
        Compute the probability of following the provided condition

        Returns
        =======
        Probability of the condition being True
        """
        cond_symbols = frozenset(rs.symbol for rs in random_symbols(condition))
        assert cond_symbols.issubset(self.symbols)
        return sum(self.prob_of(elem) for elem in self.where(condition))

    def conditional_space(self, condition):
        """
        Create a conditional Space following given conditions

        Returns
        =======
        A FinitePSpace following given condition
        """
        domain = self.where(condition)
        prob = self.probability(condition)
        density = dict((key, val / prob)
                for key, val in self._density.items() if domain._test(key))
        return FinitePSpace(domain, density)

    def sample(self):
        """
        Internal sample method

        Returns
        =======
        A dictionary mapping RandomSymbol to realization value.
        """
        expr = Tuple(*self.values)
        cdf = self.sorted_cdf(expr, python_float=True)

        x = random.uniform(0, 1)
        # Find first occurrence with cumulative probability less than x
        # This should be replaced with binary search
        for value, cum_prob in cdf:
            if x < cum_prob:
                # return dictionary mapping RandomSymbols to values
                return dict(list(zip(expr, value)))

        assert False, "We should never have gotten to this point"


class SingleFinitePSpace(SinglePSpace, FinitePSpace):
    """
    A single finite probability space

    Represents the probabilities of a set of random events that can be
    attributed to a single variable/symbol.

    This class is implemented by many of the standard FiniteRV types such as
    Die, Bernoulli, Coin, etc....
    """
    @property
    def domain(self):
        return SingleFiniteDomain(self.symbol, self.distribution.set)

    @property
    @cacheit
    def _density(self):
        """
        Calculate the Probability Density
        for the given Single Finite PSpace


        Returns
        =======
        A dictionary with density expression
        """
        return dict((FiniteSet((self.symbol, val)), prob)
                    for val, prob in self.distribution.dict.items())


class ProductFinitePSpace(IndependentProductPSpace, FinitePSpace):
    """
    A collection of several independent finite probability spaces

    """
    @property
    def domain(self):
        """Return the domain of Single Finite PSpace"""
        return ProductFiniteDomain(*[space.domain for space in self.spaces])

    @property
    @cacheit
    def _density(self):
        """
        Calculate the Probability Density
        for the given Single Finite PSpace


        Returns
        =======
        A dictionary with density expression
        """
        proditer = product(*[iter(space._density.items())
            for space in self.spaces])
        d = {}
        for items in proditer:
            elems, probs = list(zip(*items))
            elem = sumsets(elems)
            prob = Mul(*probs)
            d[elem] = d.get(elem, 0) + prob
        return Dict(d)

    @property
    @cacheit
    def density(self):
        """Return a dictionary with density of PSpace"""
        return Dict(self._density)

    def probability(self, condition):
        """Return the probability of the condition being True inside PSpace"""
        return FinitePSpace.probability(self, condition)

    def compute_density(self, expr):
        """Return the density function for given expression"""
        return FinitePSpace.compute_density(self, expr)

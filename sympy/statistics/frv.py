from sympy import (And, Eq, Basic, S, Expr, Symbol, cacheit, sympify, Mul, Add,
        And, Or)
from sympy.core.sets import FiniteSet
from rv import (Domain,  ProductDomain, ConditionalDomain, PSpace,
        ProductPSpace, random_symbols, sumsets)
import itertools
from sympy.core.containers import Dict

class CountableDomain(Domain):
    def __iter__(self):
        raise NotImplementedError()
    def __contains__(self, other):
        other in self.__iter__()

class FiniteDomain(CountableDomain):
    is_finite = True
    def __new__(cls, elements):
        elements = FiniteSet(*elements)
        symbols = FiniteSet(sym for sym, val in elements)
        return Domain.__new__(cls, symbols, elements)

    @property
    def elements(self):
        return self.args[1]

    def __contains__(self, other):
        return other in self.elements

    def __iter__(self):
        return self.elements.__iter__()

    def as_boolean(self):
        return Or(*[And(*[Eq(sym, val) for sym, val in item]) for item in self])

class SingleFiniteDomain(FiniteDomain):
    def __new__(cls, symbol, set):
        return Domain.__new__(cls, (symbol, ), FiniteSet(*set))

    @property
    def symbol(self):
        return tuple(self.symbols)[0]
    @property
    def elements(self):
        return FiniteSet(frozenset(((self.symbol, elem), )) for elem in self.set)
    @property
    def set(self):
        return self.args[1]

    def __iter__(self):
        return (frozenset(((self.symbol, elem),)) for elem in self.set)

    def __contains__(self, other):
        sym, val = tuple(other)[0]
        return sym == self.symbol and val in self.set

class ProductFiniteDomain(ProductDomain, FiniteDomain):

    def __iter__(self):
        proditer = itertools.product(*self.domains)
        return (sumsets(items) for items in proditer)

    @property
    def elements(self):
        return FiniteSet(iter(self))

class ConditionalFiniteDomain(ConditionalDomain, ProductFiniteDomain):

    def _test(self, elem):
        val = self.condition.subs(dict(elem))
        assert val in [True, False]
        return val

    def __contains__(self, other):
        return other in self.fulldomain and self._test(other)

    def __iter__(self):
        return (elem for elem in self.fulldomain if self._test(elem))

    @property
    def set(self):
        return FiniteSet(elem for elem in self.fulldomain if elem in self)

#=============================================
#=========  Probability Space  ===============
#=============================================

class FiniteDensity(Dict):
    pass

class FinitePSpace(PSpace):

    is_finite = True
    def __new__(cls, domain, density):
        density = {sympify(key):sympify(val) for key, val in density.items()}
        public_density = Dict(density)

        obj = PSpace.__new__(cls, domain, public_density)
        obj._density = density
        return obj

    def prob_of(self, elem):
        return self._density.get(elem,0)

    def where(self, condition):
        assert all(r.symbol in self.symbols for r in random_symbols(condition))
        return ConditionalFiniteDomain(self.domain, condition)

    def compute_density(self, expr):
        expr = expr.subs({rs:rs.symbol for rs in self.values})
        d = {}
        for elem in self.domain:
            val = expr.subs(dict(elem))
            prob = self.prob_of(elem)
            d[val] = d.get(val, 0) + prob
        return d

    def integrate(self, expr, rvs=None):
        rvs = rvs or self.values
        expr = expr.subs({rs:rs.symbol for rs in rvs})
        return sum(expr.subs(dict(elem)) * self.prob_of(elem)
                for elem in self.domain)

    def P(self, condition):
        cond_symbols = frozenset(rs.symbol for rs in random_symbols(condition))
        assert cond_symbols.issubset(self.symbols)
        return sum(self.prob_of(elem) for elem in self.where(condition))

    def conditional_space(self, condition):
        domain = self.where(condition)
        prob = self.P(condition)
        density = {key: val / prob
                for key, val in self._density.items() if key in domain}
        return FinitePSpace(domain, density)

class SingleFinitePSpace(FinitePSpace):
    _count = 0
    _name = 'fx'

    @property
    def value(self):
        return tuple(self.values)[0]

def create_SingleFinitePSpace(density, symbol=None, cls = SingleFinitePSpace):
    symbol = symbol or cls.create_symbol()
    domain = SingleFiniteDomain(symbol, frozenset(density.keys()))
    density = dict((frozenset(((symbol, val),)) , prob)
            for val, prob in density.items())
    density = Dict(density)
    return FinitePSpace.__new__(cls, domain, density)

class ProductFinitePSpace(ProductPSpace, FinitePSpace):
    @property
    def domain(self):
        return ProductFiniteDomain(*[space.domain for space in self.spaces])
    @property
    @cacheit
    def _density(self):
        proditer = itertools.product(*[space._density.iteritems()
            for space in self.spaces])
        d = {}
        for items in proditer:
            elems, probs = zip(*items)
            elem = sumsets(elems)
            prob = Mul(*probs)
            d[elem] = d.get(elem, 0) + prob
        return Dict(d)

    @property
    @cacheit
    def density(self):
        return Dict(self._density)

def print_elem(elem):
    return And(*[Eq(sym, val) for sym, val in elem])

from rv import Domain, ProductDomain, PSpace, random_symbols, ProductPSpace
from sympy.functions.special.delta_functions import DiracDelta
from sympy import integrate, S, Interval, Dummy, FiniteSet, Mul

oo = S.Infinity

class ContinuousDomain(Domain):
    pass

class SingleContinuousDomain(ContinuousDomain):
    def __new__(cls, symbol, set):
        assert symbol.is_Symbol
        symbols = FiniteSet(symbol)
        return Domain.__new__(cls, symbols, set)

    @property
    def set(self):
        return self.args[1]
    @property
    def symbol(self):
        return tuple(self.symbols)[0]

    def __contains__(self, other):
        if len(other)!=1:
            return False
        sym, val = tuple(other)[0]
        return self.symbol == sym and val in self.set

    def integrate(self, expr, variables=None):
        assert not variables or frozenset(variables) == frozenset(self.symbols)
        return integrate(expr, (self.symbol, self.set)) # assumes only intervals

class ProductContinuousDomain(ProductDomain, ContinuousDomain):

    def integrate(self, expr, variables=None):
        variables = variables or self.symbols
        for domain in domains:
            expr = domain.integrate(expr, variables & domain.symbols)
        return expr

class ContinuousPSpace(PSpace):
    is_continuous = True

    def integrate(self, expr, rvs=None):
        rvs = rvs or self.values
        expr = expr.subs({rv:rv.symbol for rv in rvs})
        domain_symbols = frozenset(rv.symbol for rv in rvs)
        return self.domain.integrate(self.density * expr, domain_symbols)

    def computeDensity(self, expr):
        z = Dummy('z', real=True)
        return z, self.integrate(DiracDelta(expr - z))

class SingleContinuousPSpace(ContinuousPSpace):
    def __new__(cls, symbol, density, set=Interval(-oo, oo)):
        assert symbol.is_Symbol
        domain = SingleContinuousDomain(symbol, set)
        return ContinuousPSpace.__new__(cls, domain, density)

    @property
    def value(self):
        return tuple(self.values)[0]

class ProductContinuousPSpace(ProductPSpace, ContinuousPSpace):
    @property
    def density(self):
        return Mul(*[space.density for space in self.spaces])


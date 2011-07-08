from rv import (Domain, SingleDomain, ProductDomain, PSpace, random_symbols,
        ProductPSpace)
from sympy.functions.special.delta_functions import DiracDelta
from sympy import S, Interval, Dummy, FiniteSet, Mul, Integral
from sympy.solvers.inequalities import reduce_poly_inequalities
from sympy import integrate as sympy_integrate
oo = S.Infinity


def integrate(*args, **kwargs):
    lazy = kwargs.get('lazy', False)
    if not lazy:
        return sympy_integrate(*args)
    else:
        return Integral(*args)


class ContinuousDomain(Domain):
    is_continuous = True
    pass

class SingleContinuousDomain(ContinuousDomain, SingleDomain):
    def __new__(cls, symbol, set):
        assert symbol.is_Symbol
        symbols = FiniteSet(symbol)
        return Domain.__new__(cls, symbols, set)

    def integrate(self, expr, variables=None, **kwargs):
        if variables is None:
            variables = self.symbols
        if not variables:
            return expr
        assert frozenset(variables) == frozenset(self.symbols)
        # assumes only intervals
        return integrate(expr, (self.symbol, self.set), **kwargs)

class ProductContinuousDomain(ProductDomain, ContinuousDomain):

    def integrate(self, expr, variables=None, **kwargs):
        if variables is None:
            variables = self.symbols
        for domain in self.domains:
            domain_vars = frozenset(variables) & frozenset(domain.symbols)
            if domain_vars:
                expr = domain.integrate(expr, domain_vars, **kwargs)
        return expr

class ContinuousPSpace(PSpace):
    is_continuous = True

    def integrate(self, expr, rvs=None, **kwargs):
        if rvs == None:
            rvs = self.values
        else:
            rvs = frozenset(rvs)
        expr_rvs = frozenset(random_symbols(expr)) & rvs

        # Marginalize out every random variables that aren't in the expression
        # or aren't in the input list - marginalize every ignored variable
        inactive_rvs = self.values - expr_rvs # integrated rvs not in expr
        density = self.domain.integrate(self.density,
                {rv.symbol for rv in inactive_rvs}, **kwargs)

        expr = expr.subs({rv:rv.symbol for rv in rvs})

        # Now integrate over rvs found both in the expr and input rvs
        domain_symbols = frozenset(rv.symbol for rv in expr_rvs)

        return self.domain.integrate(density * expr, domain_symbols, **kwargs)

    def compute_density(self, expr, **kwargs):
        # Common case Density(X) where X in self.values
        if expr in self.values:
            # Marginalize all other random symbols out of the density
            density = self.domain.integrate(self.density, {rs.symbol
                for rs in self.values - frozenset((expr,))},  **kwargs)
            return expr.symbol, density

        z = Dummy('z', real=True)
        return z, self.integrate(DiracDelta(expr - z), **kwargs)

    def P(self, condition, **kwargs):
        # Univariate case can be handled by where
        try:
            domain = self.where(condition)
            rv = [rv for rv in self.values if rv.symbol == domain.symbol][0]
            # Integrate out all other random variables
            z, pdf = self.compute_density(rv, **kwargs)
            # Integrate out this last variable over the special domain
            return integrate(pdf, (z, domain.set), **kwargs)
        # Other cases can be turned into univariate case
        # by computing a density handled by density computation
        except NotImplementedError:
            expr = condition.lhs - condition.rhs
            val, density = self.compute_density(expr, **kwargs)
            # Turn problem into univariate case
            space = SingleContinuousPSpace(val, density)
            return space.P(condition.__class__(space.value, 0))


    def where(self, condition):
        rvs = frozenset(random_symbols(condition))
        if not (len(rvs)==1 and rvs.issubset(self.values)):
            raise NotImplementedError(
                    "Multiple continuous random variables not supported")
        rv = tuple(rvs)[0]
        interval = reduce_poly_inequalities([[condition]], rv, relational=False)
        return SingleContinuousDomain(rv.symbol, interval)

    def conditional_space(self, condition, normalize=True, **kwargs):
        condition = condition.subs({rv:rv.symbol for rv in self.values})
        if condition.is_Equality:
            domain = self.domain
            density = self.density * DiracDelta(condition.lhs-condition.rhs)
            if normalize:
                density = density / domain.integrate(density, **kwargs)
        else:
            NotImplementedError("Only Equality conditions supported")

        return ContinuousPSpace(domain, density)



class SingleContinuousPSpace(ContinuousPSpace):
    _count = 0
    _name = 'x'
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


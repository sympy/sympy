from rv import (Domain, SingleDomain, ConditionalDomain, ProductDomain, PSpace,
        random_symbols, ProductPSpace)
from sympy.functions.special.delta_functions import DiracDelta
from sympy import S, Interval, Dummy, FiniteSet, Mul, Integral, And
from sympy.solvers.inequalities import reduce_poly_inequalities
from sympy import integrate as sympy_integrate
oo = S.Infinity

def integrate(*args, **kwargs):
    """
    Wrap around sympy integrate function to include a lazy flag
    if lazy==True then just return the Integral object
    """
    lazy = kwargs.get('lazy', False)
    if not lazy:
        return sympy_integrate(*args)
    else:
        return Integral(*args)

class ContinuousDomain(Domain):
    """
    A domain with continuous support.
    Represented using symbols and Intervals
    """
    is_Continuous = True

    def as_boolean(self):
        return Or(*[And(*[Eq(sym, val) for sym, val in item]) for item in self])

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

    def as_boolean(self):
        return self.set.as_relational(self.symbol)


class ProductContinuousDomain(ProductDomain, ContinuousDomain):

    def integrate(self, expr, variables=None, **kwargs):
        if variables is None:
            variables = self.symbols
        for domain in self.domains:
            domain_vars = frozenset(variables) & frozenset(domain.symbols)
            if domain_vars:
                expr = domain.integrate(expr, domain_vars, **kwargs)
        return expr

class ConditionalContinuousDomain(ContinuousDomain, ConditionalDomain):


    def integrate(self, expr, variables=None, **kwargs):
        # Extract the full integral
        fullintegral = self.fulldomain.integrate(expr, variables, lazy=True)
        # separate into integrand and limits
        integrand, limits = fullintegral.function, fullintegral.limits

        conditions = [self.condition]
        while conditions:
            cond = conditions.pop()
            if cond.is_Boolean:
                if cond.is_And:
                    conditions.extend(cond.args)
                elif cond.is_Or:
                    raise NotImplementedError("Or not implemented here")
            elif cond.is_Relational:
                if cond.is_Equality:
                    # Add the appropriate Delta to the integrand
                    integrand *= DiracDelta(cond.lhs-cond.rhs)
                else:
                    raise NotImplementedError(
                            "Inequalities not yet implemented")
            else:
                raise ValueError(
                        "Condition %s is not a relational or Boolean"%cond)

        return integrate(integrand, *limits, **kwargs)



class ContinuousPSpace(PSpace):
    is_Continuous = True

    def integrate(self, expr, rvs=None, **kwargs):
        if rvs == None:
            rvs = self.values
        else:
            rvs = frozenset(rvs)

        expr = expr.subs(dict((rv, rv.symbol) for rv in rvs))

        domain_symbols = frozenset(rv.symbol for rv in rvs)

        return self.domain.integrate(self.density * expr,
                domain_symbols, **kwargs)

    def compute_density(self, expr, **kwargs):
        # Common case Density(X) where X in self.values
        if expr in self.values:
            # Marginalize all other random symbols out of the density
            density = self.domain.integrate(self.density, set(rs.symbol
                for rs in self.values - frozenset((expr,))),  **kwargs)
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
        interval = interval.intersect(self.domain.set)
        return SingleContinuousDomain(rv.symbol, interval)

    def conditional_space(self, condition, normalize=True, **kwargs):

        condition = condition.subs(dict((rv,rv.symbol) for rv in self.values))

        domain = ConditionalContinuousDomain(self.domain, condition)
        density = self.density
        if normalize:
            density = density / domain.integrate(density, **kwargs)

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


"""
Continuous Random Variables Module

See Also
========
sympy.stats.crv_types
sympy.stats.rv
sympy.stats.frv
"""

from __future__ import print_function, division

from sympy.stats.rv import (RandomDomain, SingleDomain, ConditionalDomain,
        ProductDomain, PSpace, SinglePSpace, random_symbols, ProductPSpace,
        NamedArgsMixin)
from sympy.functions.special.delta_functions import DiracDelta
from sympy import (Interval, Intersection, symbols, sympify, Dummy, Mul,
        Integral, And, Or, Piecewise, solve, cacheit, integrate, oo, Lambda,
        Basic, S)
from sympy.solvers.solveset import solveset
from sympy.solvers.inequalities import reduce_rational_inequalities
from sympy.polys.polyerrors import PolynomialError
import random


class ContinuousDomain(RandomDomain):
    """
    A domain with continuous support

    Represented using symbols and Intervals.
    """
    is_Continuous = True

    def as_boolean(self):
        raise NotImplementedError("Not Implemented for generic Domains")


class SingleContinuousDomain(ContinuousDomain, SingleDomain):
    """
    A univariate domain with continuous support

    Represented using a single symbol and interval.
    """
    def integrate(self, expr, variables=None, **kwargs):
        if variables is None:
            variables = self.symbols
        if not variables:
            return expr
        if frozenset(variables) != frozenset(self.symbols):
            raise ValueError("Values should be equal")
        # assumes only intervals
        return Integral(expr, (self.symbol, self.set), **kwargs)

    def as_boolean(self):
        return self.set.as_relational(self.symbol)


class ProductContinuousDomain(ProductDomain, ContinuousDomain):
    """
    A collection of independent domains with continuous support
    """

    def integrate(self, expr, variables=None, **kwargs):
        if variables is None:
            variables = self.symbols
        for domain in self.domains:
            domain_vars = frozenset(variables) & frozenset(domain.symbols)
            if domain_vars:
                expr = domain.integrate(expr, domain_vars, **kwargs)
        return expr

    def as_boolean(self):
        return And(*[domain.as_boolean() for domain in self.domains])


class ConditionalContinuousDomain(ContinuousDomain, ConditionalDomain):
    """
    A domain with continuous support that has been further restricted by a
    condition such as x > 3
    """

    def integrate(self, expr, variables=None, **kwargs):
        if variables is None:
            variables = self.symbols
        if not variables:
            return expr
        # Extract the full integral
        fullintgrl = self.fulldomain.integrate(expr, variables)
        # separate into integrand and limits
        integrand, limits = fullintgrl.function, list(fullintgrl.limits)

        conditions = [self.condition]
        while conditions:
            cond = conditions.pop()
            if cond.is_Boolean:
                if isinstance(cond, And):
                    conditions.extend(cond.args)
                elif isinstance(cond, Or):
                    raise NotImplementedError("Or not implemented here")
            elif cond.is_Relational:
                if cond.is_Equality:
                    # Add the appropriate Delta to the integrand
                    integrand *= DiracDelta(cond.lhs - cond.rhs)
                else:
                    symbols = cond.free_symbols & set(self.symbols)
                    if len(symbols) != 1:  # Can't handle x > y
                        raise NotImplementedError(
                            "Multivariate Inequalities not yet implemented")
                    # Can handle x > 0
                    symbol = symbols.pop()
                    # Find the limit with x, such as (x, -oo, oo)
                    for i, limit in enumerate(limits):
                        if limit[0] == symbol:
                            # Make condition into an Interval like [0, oo]
                            cintvl = reduce_rational_inequalities_wrap(
                                cond, symbol)
                            # Make limit into an Interval like [-oo, oo]
                            lintvl = Interval(limit[1], limit[2])
                            # Intersect them to get [0, oo]
                            intvl = cintvl.intersect(lintvl)
                            # Put back into limits list
                            limits[i] = (symbol, intvl.left, intvl.right)
            else:
                raise TypeError(
                    "Condition %s is not a relational or Boolean" % cond)

        return Integral(integrand, *limits, **kwargs)

    def as_boolean(self):
        return And(self.fulldomain.as_boolean(), self.condition)

    @property
    def set(self):
        if len(self.symbols) == 1:
            return (self.fulldomain.set & reduce_rational_inequalities_wrap(
                self.condition, tuple(self.symbols)[0]))
        else:
            raise NotImplementedError(
                "Set of Conditional Domain not Implemented")


class ContinuousDistribution(Basic):
    def __call__(self, *args):
        return self.pdf(*args)


class SingleContinuousDistribution(ContinuousDistribution, NamedArgsMixin):
    """ Continuous distribution of a single variable

    Serves as superclass for Normal/Exponential/UniformDistribution etc....

    Represented by parameters for each of the specific classes.  E.g
    NormalDistribution is represented by a mean and standard deviation.

    Provides methods for pdf, cdf, and sampling

    See Also:
        sympy.stats.crv_types.*
    """

    set = Interval(-oo, oo)

    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @staticmethod
    def check(*args):
        pass

    def sample(self):
        """ A random realization from the distribution """
        icdf = self._inverse_cdf_expression()
        return icdf(random.uniform(0, 1))

    @cacheit
    def _inverse_cdf_expression(self):
        """ Inverse of the CDF

        Used by sample
        """
        x, z = symbols('x, z', real=True, positive=True, cls=Dummy)
        # Invert CDF
        try:
            inverse_cdf = solveset(self.cdf(x) - z, x)
            if isinstance(inverse_cdf, Intersection) and inverse_cdf.args[0] == S.Reals:
                inverse_cdf = list(inverse_cdf.args[1])
        except NotImplementedError:
            inverse_cdf = None
        if not inverse_cdf or len(inverse_cdf) != 1:
            raise NotImplementedError("Could not invert CDF")

        return Lambda(z, inverse_cdf[0])

    @cacheit
    def compute_cdf(self, **kwargs):
        """ Compute the CDF from the PDF

        Returns a Lambda
        """
        x, z = symbols('x, z', real=True, finite=True, cls=Dummy)
        left_bound = self.set.start

        # CDF is integral of PDF from left bound to z
        pdf = self.pdf(x)
        cdf = integrate(pdf, (x, left_bound, z), **kwargs)
        # CDF Ensure that CDF left of left_bound is zero
        cdf = Piecewise((cdf, z >= left_bound), (0, True))
        return Lambda(z, cdf)

    def cdf(self, x, **kwargs):
        """ Cumulative density function """
        return self.compute_cdf(**kwargs)(x)

    def expectation(self, expr, var, evaluate=True, **kwargs):
        """ Expectation of expression over distribution """
        integral = Integral(expr * self.pdf(var), (var, self.set), **kwargs)
        return integral.doit() if evaluate else integral

class ContinuousDistributionHandmade(SingleContinuousDistribution):
    _argnames = ('pdf',)

    @property
    def set(self):
        return self.args[1]

    def __new__(cls, pdf, set=Interval(-oo, oo)):
        return Basic.__new__(cls, pdf, set)


class ContinuousPSpace(PSpace):
    """ Continuous Probability Space

    Represents the likelihood of an event space defined over a continuum.

    Represented with a ContinuousDomain and a PDF (Lambda-Like)
    """

    is_Continuous = True
    is_real = True

    @property
    def domain(self):
        return self.args[0]

    @property
    def density(self):
        return self.args[1]

    @property
    def pdf(self):
        return self.density(*self.domain.symbols)

    def integrate(self, expr, rvs=None, **kwargs):
        if rvs is None:
            rvs = self.values
        else:
            rvs = frozenset(rvs)

        expr = expr.xreplace(dict((rv, rv.symbol) for rv in rvs))

        domain_symbols = frozenset(rv.symbol for rv in rvs)

        return self.domain.integrate(self.pdf * expr,
                domain_symbols, **kwargs)

    def compute_density(self, expr, **kwargs):
        # Common case Density(X) where X in self.values
        if expr in self.values:
            # Marginalize all other random symbols out of the density
            randomsymbols = tuple(set(self.values) - frozenset([expr]))
            symbols = tuple(rs.symbol for rs in randomsymbols)
            pdf = self.domain.integrate(self.pdf, symbols, **kwargs)
            return Lambda(expr.symbol, pdf)

        z = Dummy('z', real=True, finite=True)
        return Lambda(z, self.integrate(DiracDelta(expr - z), **kwargs))

    @cacheit
    def compute_cdf(self, expr, **kwargs):
        if not self.domain.set.is_Interval:
            raise ValueError(
                "CDF not well defined on multivariate expressions")

        d = self.compute_density(expr, **kwargs)
        x, z = symbols('x, z', real=True, finite=True, cls=Dummy)
        left_bound = self.domain.set.start

        # CDF is integral of PDF from left bound to z
        cdf = integrate(d(x), (x, left_bound, z), **kwargs)
        # CDF Ensure that CDF left of left_bound is zero
        cdf = Piecewise((cdf, z >= left_bound), (0, True))
        return Lambda(z, cdf)

    def probability(self, condition, **kwargs):
        z = Dummy('z', real=True, finite=True)
        # Univariate case can be handled by where
        try:
            domain = self.where(condition)
            rv = [rv for rv in self.values if rv.symbol == domain.symbol][0]
            # Integrate out all other random variables
            pdf = self.compute_density(rv, **kwargs)
            # Integrate out the last variable over the special domain
            return Integral(pdf(z), (z, domain.set), **kwargs)

        # Other cases can be turned into univariate case
        # by computing a density handled by density computation
        except NotImplementedError:
            from sympy.stats.rv import density
            expr = condition.lhs - condition.rhs
            dens = density(expr, **kwargs)
            if not isinstance(dens, ContinuousDistribution):
                dens = ContinuousDistributionHandmade(dens)
            # Turn problem into univariate case
            space = SingleContinuousPSpace(z, dens)
            return space.probability(condition.__class__(space.value, 0))

    def where(self, condition):
        rvs = frozenset(random_symbols(condition))
        if not (len(rvs) == 1 and rvs.issubset(self.values)):
            raise NotImplementedError(
                "Multiple continuous random variables not supported")
        rv = tuple(rvs)[0]
        interval = reduce_rational_inequalities_wrap(condition, rv)
        interval = interval.intersect(self.domain.set)
        return SingleContinuousDomain(rv.symbol, interval)

    def conditional_space(self, condition, normalize=True, **kwargs):

        condition = condition.xreplace(dict((rv, rv.symbol) for rv in self.values))

        domain = ConditionalContinuousDomain(self.domain, condition)
        if normalize:
            pdf = self.pdf / domain.integrate(self.pdf, **kwargs)
            density = Lambda(domain.symbols, pdf)

        return ContinuousPSpace(domain, density)


class SingleContinuousPSpace(ContinuousPSpace, SinglePSpace):
    """
    A continuous probability space over a single univariate variable

    These consist of a Symbol and a SingleContinuousDistribution

    This class is normally accessed through the various random variable
    functions, Normal, Exponential, Uniform, etc....
    """

    @property
    def set(self):
        return self.distribution.set

    @property
    def domain(self):
        return SingleContinuousDomain(sympify(self.symbol), self.set)

    def sample(self):
        """
        Internal sample method

        Returns dictionary mapping RandomSymbol to realization value.
        """
        return {self.value: self.distribution.sample()}

    def integrate(self, expr, rvs=None, **kwargs):
        rvs = rvs or (self.value,)
        if self.value not in rvs:
            return expr

        expr = expr.xreplace(dict((rv, rv.symbol) for rv in rvs))

        x = self.value.symbol
        try:
            return self.distribution.expectation(expr, x, evaluate=False, **kwargs)
        except Exception:
            return Integral(expr * self.pdf, (x, self.set), **kwargs)

    def compute_cdf(self, expr, **kwargs):
        if expr == self.value:
            return self.distribution.compute_cdf(**kwargs)
        else:
            return ContinuousPSpace.compute_cdf(self, expr, **kwargs)

    def compute_density(self, expr, **kwargs):
        # http://en.wikipedia.org/wiki/Random_variable#Functions_of_random_variables
        if expr == self.value:
            return self.density
        y = Dummy('y')

        gs = solveset(expr - y, self.value)

        if isinstance(gs, Intersection) and gs.args[0] == S.Reals:
            gs = list(gs.args[1])

        if not gs:
            raise ValueError("Can not solve %s for %s"%(expr, self.value))
        fx = self.compute_density(self.value)
        fy = sum(fx(g) * abs(g.diff(y)) for g in gs)
        return Lambda(y, fy)


class ProductContinuousPSpace(ProductPSpace, ContinuousPSpace):
    """
    A collection of independent continuous probability spaces
    """
    @property
    def pdf(self):
        p = Mul(*[space.pdf for space in self.spaces])
        return p.subs(dict((rv, rv.symbol) for rv in self.values))

def _reduce_inequalities(conditions, var, **kwargs):
    try:
        return reduce_rational_inequalities(conditions, var, **kwargs)
    except PolynomialError:
        raise ValueError("Reduction of condition failed %s\n" % conditions[0])


def reduce_rational_inequalities_wrap(condition, var):
    if condition.is_Relational:
        return _reduce_inequalities([[condition]], var, relational=False)
    if condition.__class__ is Or:
        return _reduce_inequalities([list(condition.args)],
                var, relational=False)
    if condition.__class__ is And:
        intervals = [_reduce_inequalities([[arg]], var, relational=False)
            for arg in condition.args]
        I = intervals[0]
        for i in intervals:
            I = I.intersect(i)
        return I

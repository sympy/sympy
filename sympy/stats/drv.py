from __future__ import print_function, division

from sympy import (Basic, sympify, symbols, Dummy, Lambda, summation,
        Piecewise, S, cacheit, Sum, exp, I, oo)
from sympy.solvers.solveset import solveset
from sympy.solvers.inequalities import reduce_rational_inequalities
from sympy.stats.crv import (reduce_rational_inequalities_wrap,
        _reduce_inequalities)
from sympy.stats.rv import (NamedArgsMixin, SinglePSpace, SingleDomain,
        random_symbols, RandomDomain, PSpace)
from sympy.stats.symbolic_probability import Probability
from sympy.functions.elementary.integers import floor
from sympy.sets.fancysets import Range, FiniteSet
from sympy.sets.sets import Union
from sympy.utilities import filldedent
import random

class SingleDiscreteDistribution(Basic, NamedArgsMixin):
    """ Discrete distribution of a single variable

    Serves as superclass for PoissonDistribution etc....

    Provides methods for pdf, cdf, and sampling

    See Also:
        sympy.stats.crv_types.*
    """

    set = S.Integers

    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @staticmethod
    def check(*args):
        pass

    def sample(self):
        """ A random realization from the distribution """
        icdf = self._inverse_cdf_expression()
        while True:
            sample_ = floor(list(icdf(random.uniform(0, 1)))[0])
            if sample_ >= self.set.inf:
                return sample_

    @cacheit
    def _inverse_cdf_expression(self):
        """ Inverse of the CDF

        Used by sample
        """
        x = symbols('x', positive=True,
         integer=True, cls=Dummy)
        z = symbols('z', positive=True, cls=Dummy)
        cdf_temp = self.cdf(x)
        # Invert CDF
        try:
            inverse_cdf = solveset(cdf_temp - z, x, domain=S.Reals)
        except NotImplementedError:
            inverse_cdf = None
        if not inverse_cdf or len(inverse_cdf.free_symbols) != 1:
            raise NotImplementedError("Could not invert CDF")
        return Lambda(z, inverse_cdf)

    @cacheit
    def compute_cdf(self, **kwargs):
        """ Compute the CDF from the PDF

        Returns a Lambda
        """
        x, z = symbols('x, z', integer=True, finite=True, cls=Dummy)
        left_bound = self.set.inf

        # CDF is integral of PDF from left bound to z
        pdf = self.pdf(x)
        cdf = summation(pdf, (x, left_bound, z), **kwargs)
        # CDF Ensure that CDF left of left_bound is zero
        cdf = Piecewise((cdf, z >= left_bound), (0, True))
        return Lambda(z, cdf)

    def cdf(self, x, **kwargs):
        """ Cumulative density function """
        return self.compute_cdf(**kwargs)(x)

    @cacheit
    def compute_characteristic_function(self, **kwargs):
        """ Compute the characteristic function from the PDF

        Returns a Lambda
        """
        x, t = symbols('x, t', real=True, finite=True, cls=Dummy)
        pdf = self.pdf(x)
        cf = summation(exp(I*t*x)*pdf, (x, self.set.inf, self.set.sup))
        return Lambda(t, cf)

    def characteristic_function(self, t, **kwargs):
        """ Characteristic function """
        return self.compute_characteristic_function(**kwargs)(t)

    def expectation(self, expr, var, evaluate=True, **kwargs):
        """ Expectation of expression over distribution """
        # TODO: support discrete sets with non integer stepsizes
        if evaluate:
            return summation(expr * self.pdf(var),
                         (var, self.set.inf, self.set.sup), **kwargs)
        else:
            return Sum(expr * self.pdf(var),
                         (var, self.set.inf, self.set.sup), **kwargs)

    def __call__(self, *args):
        return self.pdf(*args)

class DiscreteDomain(RandomDomain):
    """
    A domain with discrete support

    Represented using an ImageSet or S.Integers, S.Naturals etc.
    """
    is_continuous = True

    def __contains__(self, other):
        return other is self.set

    def as_boolean(self):
        raise NotImplementedError

class SingleDiscreteDomain(DiscreteDomain, SingleDomain):
    """
    A domain with univariate discrete support

    Represented using an ImageSet or S.Integers, S.Naturals etc.
    """

    @property
    def symbol(self):
        return self.args[0]

    @property
    def set(self):
        return self.args[1]

    def as_boolean(self):
        raise NotImplementedError

# class ProductDiscretePSpace(ProductDomain, DiscreteDomain):
#     pass

class DiscretePSpace(PSpace):
    """
    Represents the likelihood of an event space defined on a countably
    infinite set.
    """

    is_real = True
    is_discrete = True

    @property
    def domain(self):
        return self.args[0]

    @property
    def density(self):
        self.args[1]

    def compute_density(self):
        pass

    def where(self, Condition):
        rvs = random_symbols(condition)
        assert all(r.symbol in self.symbols for r in rvs)
        if (len(rvs) > 1):
            raise NotImplementedError(filldedent('''Multivariate discrete
                random variables are not yet supported.'''))
        conditional_domain = reduce_rational_inequalities_wrap(condition,
            rvs[0])
        conditional_domain = conditional_domain.intersect(self.domain.set)
        return SingleDiscreteDomain(rvs, conditional_domain)

    def compute_cdf(self, expr):
        rvs = random_symbols(expr)
        if len(rvs) > 1:
            raise NotImplementedError(filldedent('''Multivariate discrete
                random variables are not yet supported.'''))
        density = self.pdf(expr)
        cdf = summation()

class SingleDiscretePSpace(SinglePSpace):
    """ Discrete probability space over a single univariate variable """
    is_real = True

    @property
    def set(self):
        return self.distribution.set

    @property
    def domain(self):
        return SingleDiscreteDomain(self.symbol, self.set)

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
            return self.distribution.expectation(expr, x, evaluate=False,
                    **kwargs)
        except NotImplementedError:
            return Sum(expr * self.pdf, (x, self.set.inf, self.set.sup),
                    **kwargs)

    def compute_cdf(self, expr, **kwargs):
        if expr == self.value:
            return self.distribution.compute_cdf(**kwargs)
        else:
            raise NotImplementedError()

    def compute_density(self, expr, **kwargs):
        if expr == self.value:
            return self.distribution
        raise NotImplementedError()

    def compute_characteristic_function(self, expr, **kwargs):
        if expr == self.value:
            return self.distribution.compute_characteristic_function(**kwargs)
        else:
            raise NotImplementedError()

    def restricted_domain(self, condition):
        rvs = random_symbols(condition)
        assert all(r.symbol in self.symbols for r in rvs)
        if (len(rvs) > 1):
            raise NotImplementedError(filldedent('''Multivariate discrete
            random variables are not yet supported.'''))
        conditional_domain = reduce_rational_inequalities_wrap(condition,
            rvs[0])
        conditional_domain = conditional_domain.intersect(self.domain.set)
        return conditional_domain

    def probability(self, condition):
        _domain = self.restricted_domain(condition)
        if condition == False or _domain is S.EmptySet:
            return S.Zero
        if condition == True or _domain == self.set:
            return S.One
        try:
            return self.eval_prob(_domain)
        except NotImplementedError:
            return Probability(condition)

    def eval_prob(self, _domain):
        if isinstance(_domain, Range):
            n = symbols('n')
            inf, sup, step = (r for r in _domain.args)
            summand = ((self.pdf).replace(
                self.symbol, inf + n*step))
            rv = summation(summand,
                (n, 0, floor((sup - inf)/step - 1))).doit()
            return rv
        elif isinstance(_domain, FiniteSet):
            pdf = Lambda(self.symbol, self.pdf)
            rv = sum(pdf(x) for x in _domain)
            return rv
        elif isinstance(_domain, Union):
            rv = sum(self.eval_prob(x) for x in _domain.args)
            return rv
        else:
            raise NotImplementedError(filldedent('''Probability for
                the domain %s cannot be calculated.'''%(_domain)))

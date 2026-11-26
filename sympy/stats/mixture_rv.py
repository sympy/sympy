from __future__ import annotations

from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.function import Lambda
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.sets.sets import Union
from sympy.core.symbol import Dummy
from sympy.stats.rv import (Distribution, NamedArgsMixin, RandomSymbol,
                            _symbol_converter, _value_check, PSpace)
from sympy.stats.crv import ContinuousDistribution, SingleContinuousPSpace
from sympy.stats.drv import DiscreteDistribution, SingleDiscretePSpace
from sympy.stats.frv import SingleFiniteDistribution, SingleFinitePSpace
from sympy.stats.crv_types import (ContinuousDistributionHandmade,
                                   NormalDistribution)
from sympy.stats.drv_types import DiscreteDistributionHandmade
from sympy.stats.frv_types import FiniteDistributionHandmade


class MixturePSpace(PSpace):
    """
    A probability space for a Mixture Distribution. Builds a handmade
    single-kind probability space from the weighted component densities and
    delegates the actual computations to it.
    """
    def __new__(cls, s, distribution):
        s = _symbol_converter(s)
        if not isinstance(distribution, MixtureDistribution):
            raise ValueError("%s should be an instance of MixtureDistribution" % distribution)
        return Basic.__new__(cls,s,distribution)

    @property
    def value(self):
        return RandomSymbol(self.symbol, self)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def is_Continuous(self):
        return self.distribution.is_Continuous

    @property
    def is_Finite(self):
        return self.distribution.is_Finite

    @property
    def is_Discrete(self):
        return self.distribution.is_Discrete

    @property
    def distribution(self):
        return self.args[1]

    @property
    def pdf(self):
        return self.distribution.pdf(self.symbol)

    @property
    def set(self):
        return self.distribution.set

    @property
    def domain(self):
        return self._get_newpspace().domain

    def _get_newpspace(self):
        x = Dummy('x')
        dist = self.distribution

        if dist.is_Continuous:
            pdf = Lambda(x, dist.pdf(x))
            return SingleContinuousPSpace(self.symbol,
                                          ContinuousDistributionHandmade(pdf, dist.set))

        if dist.is_Discrete:
            pdf = Lambda(x,dist.pdf(x))
            return SingleDiscretePSpace(self.symbol,
                                        DiscreteDistributionHandmade(pdf, dist.set))

        if dist.is_Finite:
            dens = {k: dist.pdf(k) for k in dist.set}
            return SingleFinitePSpace(self.symbol,
                                      FiniteDistributionHandmade(dens))

        raise NotImplementedError("Mixture pspace is not implemented for "
            "this distribution kind.")

    def _component_pspaces(self):
        pspaces = []
        for comp in self.distribution.components:
            if isinstance(comp, ContinuousDistribution):
                ps = SingleContinuousPSpace(self.symbol, comp)
            elif isinstance(comp, DiscreteDistribution):
                ps = SingleDiscretePSpace(self.symbol, comp)
            else:
                ps = SingleFinitePSpace(self.symbol, comp)
            pspaces.append(ps)
        return pspaces

    def compute_density(self, expr, **kwargs):
        new_pspace = self._get_newpspace()
        expr = expr.subs({self.value: new_pspace.value})
        return new_pspace.compute_density(expr, **kwargs)

    def compute_cdf(self, expr, **kwargs):
        dist = self.distribution
        # Continuous/discrete CDFs are the weighted sum of the component CDFs
        # (linearity). Integrating the combined PDF instead can introduce a
        # removable 0/0 singularity at a component's mean, so use the
        # distribution's own cdf here. Finite mixtures have no component cdf,
        # so they fall back to the handmade probability space.
        if expr == self.value and not dist.is_Finite:
            x = Dummy('x')
            return Lambda(x, dist.cdf(x))
        new_pspace = self._get_newpspace()
        expr = expr.subs({self.value: new_pspace.value})
        return new_pspace.compute_cdf(expr, **kwargs)

    def compute_expectation(self, expr, rvs=None, evaluate=False, **kwargs):
        # Linearity over components: E[g(M)] = sum_i w_i E[g(X_i)]. Integrating
        # the combined PDF over a union support (disjoint components) fails, so
        # delegate to each component's own probability space and weight.
        total = S.Zero
        for w, ps in zip(self.distribution.weights, self._component_pspaces()):
            e = expr.subs({self.value: ps.value})
            total += w*ps.compute_expectation(e, evaluate=evaluate, **kwargs)
        return total

    def probability(self, condition, **kwargs):
        # Same linearity argument: P(g(M)) = sum_i w_i P(g(X_i)).
        total = S.Zero
        for w, ps in zip(self.distribution.weights, self._component_pspaces()):
            cond = condition.subs({self.value: ps.value})
            total += w*ps.probability(cond)
        return total

    def compute_characteristic_function(self, expr, **kwargs):
        # Linearity: CF_M(t) = sum_i w_i * CF_{X_i}(t).  Each component pspace
        # returns a Lambda; evaluate them at a shared dummy and wrap the
        # weighted sum back into a Lambda so the result is callable.
        t = Dummy('t', real=True)
        total = S.Zero
        for w, ps in zip(self.distribution.weights, self._component_pspaces()):
            e = expr.subs({self.value: ps.value})
            total += w*ps.compute_characteristic_function(e, **kwargs)(t)
        return Lambda(t, total)

    def compute_moment_generating_function(self, expr, **kwargs):
        # Same linearity argument for the MGF.
        t = Dummy('t', real=True)
        total = S.Zero
        for w, ps in zip(self.distribution.weights, self._component_pspaces()):
            e = expr.subs({self.value: ps.value})
            total += w*ps.compute_moment_generating_function(e, **kwargs)(t)
        return Lambda(t, total)

    def sample(self, size=(), library='scipy', seed=None):

        import numpy
        if any(w.free_symbols for w in self.distribution.weights):
            raise ValueError("Cannot sample a mixture with symbolic weights.")
        if seed is None or isinstance(seed, int):
            rand_state = numpy.random.default_rng(seed=seed)
        else:
            rand_state = seed
        weights = [float(w) for w in self.distribution.weights]
        pspaces = self._component_pspaces()
        idx = rand_state.choice(len(pspaces), size=size, p=weights)
        # Draw every component at the requested size, threading the same
        # generator, then select element-wise by the chosen component index.
        draws = [list(ps.sample(size=size, library=library,
                                seed=rand_state).values())[0]
                 for ps in pspaces]
        stacked = numpy.stack(draws)
        result = numpy.take_along_axis(
            stacked, numpy.asarray(idx)[numpy.newaxis], axis=0)[0]
        return {self.value: result}

    def conditional_space(self, condition, **kwargs):
        new_pspace = self._get_newpspace()
        condition = condition.subs({self.value: new_pspace.value})
        return new_pspace.conditional_space(condition)


class MixtureDistribution(Distribution, NamedArgsMixin):
    """
    Class for finite Mixture Distributions.

    Parameters
    ==========

    weights : list of Expr
        The mixing weights, one per component. Must be non-negative; they are
        normalised to sum to one.
    components : list of Distribution
        The component distributions being mixed. They must all be of the same
        kind (all continuous, all discrete, or all finite).

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Mixture_distribution


    """

    _argnames = ('weights', 'components')

    def __new__(cls, weights, components):
        weights = [sympify(w) for w in weights]
        components = list(components)

        # one weight per component, and at least one component
        if len(weights) != len(components):
            raise ValueError("Got %d weights but %d components; they must "
                "match." % (len(weights), len(components)))
        if len(components) == 0:
            raise ValueError("A mixture needs at least one component.")

        # only the three univariate distribution kinds are supported
        for c in components:
            if not isinstance(c, (ContinuousDistribution, DiscreteDistribution,
                    SingleFiniteDistribution)):
                raise ValueError("Mixture component %s is not a supported "
                    "distribution." % c)

        # components must be homogeneous (all the same kind); supports are
        # combined as a union, so they need not coincide
        kinds = (ContinuousDistribution, DiscreteDistribution,
                 SingleFiniteDistribution)
        if not any(all(isinstance(c, k) for c in components) for k in kinds):
            raise ValueError("All mixture components must be of the same kind "
                "(all continuous, all discrete, or all finite).")

        # weights must be non-negative; then normalise them to sum to one
        for w in weights:
            _value_check(w >= 0, "Weights must be non-negative.")
        total = Add(*weights)
        if total == S.Zero:
            raise ValueError("Weights must not all be zero.")
        weights = Tuple(*[w/total for w in weights])
        components = Tuple(*components)

        return Basic.__new__(cls, weights, components)

    @property
    def is_Continuous(self):
        return all(isinstance(c, ContinuousDistribution)
                   for c in self.components)

    @property
    def is_Discrete(self):
        return all(isinstance(c, DiscreteDistribution)
                   for c in self.components)

    @property
    def is_Finite(self):
        return all(isinstance(c, SingleFiniteDistribution)
                   for c in self.components)

    @property
    def set(self):
        return Union(*[c.set for c in self.components])

    def pdf(self, x):
        if self.is_Finite:
            terms = [w*c.pmf(x) for w, c in zip(self.weights, self.components)]
        else:
            terms = [w*c.pdf(x) for w, c in zip(self.weights, self.components)]
        return Add(*terms)

    def cdf(self, x):
        if self.is_Finite:
            raise NotImplementedError("Finite mixture CDFs are computed via the "
                "probability space, not the distribution.")
        return Add(*[w*c.cdf(x) for w, c in zip(self.weights, self.components)])

    def _characteristic_function(self, t):
        # Linearity: CF_M(t) = sum_i w_i * CF_{X_i}(t).  Return None if any
        # component lacks a closed-form CF so the generic dispatcher in
        # ContinuousDistribution/DiscreteDistribution falls back to
        # integrating the pdf.
        parts = []
        for c in self.components:
            cf = getattr(c, '_characteristic_function', None)
            if cf is None:
                return None
            value = cf(t)
            if value is None:
                return None
            parts.append(value)
        return Add(*[w*p for w, p in zip(self.weights, parts)])

    def _moment_generating_function(self, t):
        # Same linearity argument for the MGF.
        parts = []
        for c in self.components:
            mgf = getattr(c, '_moment_generating_function', None)
            if mgf is None:
                return None
            value = mgf(t)
            if value is None:
                return None
            parts.append(value)
        return Add(*[w*p for w, p in zip(self.weights, parts)])


class _WeightedDistribution(Basic):
    """
    Internal helper used to support the ``w * Distribution + ...`` syntax for
    building two-component mixtures.  Pairs a sympified non-negative weight
    with a univariate component distribution.

    Summing two ``_WeightedDistribution``s, or adding a plain ``Distribution``
    to one, yields a :class:`MixtureDistribution`.  Three-or-more-component
    chaining is not supported by the operator syntax: use
    :func:`Mixture` directly for that.

    Not part of the public API.
    """

    def __new__(cls, weight, distribution):
        weight = sympify(weight)
        return Basic.__new__(cls, weight, distribution)

    @property
    def weight(self):
        return self.args[0]

    @property
    def distribution(self):
        return self.args[1]

    def __add__(self, other):
        if isinstance(other, _WeightedDistribution):
            return MixtureDistribution(
                [self.weight, other.weight],
                [self.distribution, other.distribution])
        if isinstance(other, (ContinuousDistribution, DiscreteDistribution,
                              SingleFiniteDistribution)):
            return MixtureDistribution(
                [self.weight, S.One],
                [self.distribution, other])
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, (ContinuousDistribution, DiscreteDistribution,
                              SingleFiniteDistribution)):
            return MixtureDistribution(
                [S.One, self.weight],
                [other, self.distribution])
        return NotImplemented


def Mixture(name, weights, components):
    """
    Create a random variable with a finite mixture distribution.

    A finite mixture draws its value from one of several component
    distributions, chosen according to fixed weights. Its density is the
    weighted sum of the component densities. All components must be of the same
    kind (all continuous, all discrete, or all finite). The weights need not
    sum to one; they are normalised automatically.

    Parameters
    ==========

    name : str
        Name of the random variable.
    weights : list of Expr
        Non-negative mixing weights, one per component.
    components : list of RandomSymbol or Distribution
        The component distributions being mixed. Random variables are accepted
        as a convenience; their distributions are used.

    Returns
    =======

    RandomSymbol

    Examples
    ========

    >>> from sympy.stats import Mixture, Die, E, P, density

    >>> M = Mixture('M', [1, 1], [Die('D1', 2), Die('D2', 4)])
    >>> E(M)
    2
    >>> P(M > 3)
    1/8
    >>> density(M).dict
    {1: 3/8, 2: 3/8, 3: 1/8, 4: 1/8}

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Mixture_distribution

    """
    dists = []
    for c in components:
        if isinstance(c, RandomSymbol):
            c = c.pspace.distribution
        dists.append(c)
    dist = MixtureDistribution(weights, dists)
    return MixturePSpace(name, dist).value


def GaussianMixture(name, weights, means, sigmas):
    """
    Create a random variable that is a mixture of univariate normals.

    Convenience constructor mirroring Mathematica's ``GaussianMixture``.
    Equivalent to ``Mixture(name, weights, [Normal(...), ...])`` but
    expressed in terms of the Gaussian parameters directly.

    Parameters
    ==========

    name : str
        Name of the random variable.
    weights : list of Expr
        Non-negative mixing weights, one per component; normalised to sum
        to one.
    means : list of Expr
        Means of the normal components.
    sigmas : list of Expr
        Standard deviations of the normal components (must be positive).

    Returns
    =======

    RandomSymbol

    Examples
    ========

    >>> from sympy.stats import GaussianMixture, E

    >>> G = GaussianMixture('G', [1, 1], [0, 5], [1, 1])
    >>> E(G)
    5/2

    References
    ==========

    .. [1] https://reference.wolfram.com/language/ref/method/GaussianMixture.html

    """
    if not (len(weights) == len(means) == len(sigmas)):
        raise ValueError("weights, means and sigmas must have the same length.")
    components = [NormalDistribution(m, s) for m, s in zip(means, sigmas)]
    return Mixture(name, weights, components)

from __future__ import annotations

from sympy.core.containers import Tuple
from sympy.core.function import expand_func
from sympy.core.numbers import oo
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.special.error_functions import erf
from sympy.sets.sets import Interval, Union
from sympy.simplify.simplify import simplify
from sympy.stats import (E, P, cdf, density, given, variance, sample, Die,
                         Mixture, GaussianMixture, Normal, Poisson,
                         characteristic_function, moment_generating_function)
from sympy.stats.crv_types import NormalDistribution, UniformDistribution
from sympy.stats.drv_types import PoissonDistribution
from sympy.stats.frv_types import BernoulliDistribution, DieDistribution
from sympy.stats.mixture_rv import (MixtureDistribution, MixturePSpace,
                                    _WeightedDistribution)
from sympy.external import import_module
from sympy.testing.pytest import raises, skip


def test_mixture_distribution_constructor():
    components = Tuple(NormalDistribution(0, 1), NormalDistribution(2, 3))
    dist = MixtureDistribution([3, 7], components)

    assert dist.args == (Tuple(S(3)/10, S(7)/10), components)
    assert dist.weights == Tuple(S(3)/10, S(7)/10)
    assert dist.components == components
    assert dist.is_Continuous
    assert not dist.is_Discrete
    assert not dist.is_Finite
    assert dist.set == Interval(-oo, oo)


def test_mixture_distribution_discrete_and_finite_types():
    discrete = MixtureDistribution([1, 1],
                                   [PoissonDistribution(1), PoissonDistribution(2)])
    finite = MixtureDistribution([1, 3],
                                 [DieDistribution(4),
                                  BernoulliDistribution(S.Half, 0, -1)])

    assert discrete.weights == Tuple(S.Half, S.Half)
    assert discrete.is_Discrete
    assert not discrete.is_Continuous
    assert not discrete.is_Finite
    assert discrete.set == S.Naturals0

    assert finite.weights == Tuple(S(1)/4, S(3)/4)
    assert finite.is_Finite
    assert not finite.is_Continuous
    assert not finite.is_Discrete
    assert finite.set == {-1, 0, 1, 2, 3, 4}


def test_mixture_distribution_set_union():
    dist = MixtureDistribution([1, 1],
                               [UniformDistribution(0, 1),
                                UniformDistribution(2, 3)])

    assert dist.set == Union(Interval(0, 1), Interval(2, 3))


def test_mixture_distribution_validation():
    normal = NormalDistribution(0, 1)
    other_normal = NormalDistribution(1, 2)

    raises(ValueError, lambda: MixtureDistribution([], []))
    raises(ValueError, lambda: MixtureDistribution([1], []))
    raises(ValueError, lambda: MixtureDistribution([1, 1], [normal]))
    raises(ValueError, lambda: MixtureDistribution([-1, 2],
                                                  [normal, other_normal]))
    raises(ValueError, lambda: MixtureDistribution([0, 0],
                                                  [normal, other_normal]))
    raises(ValueError, lambda: MixtureDistribution([1], [object()]))
    raises(ValueError, lambda: MixtureDistribution([1, 1],
                                                  [normal,
                                                   PoissonDistribution(1)]))


def test_mixture_pspace_continuous():
    x = Symbol('x')
    dist = MixtureDistribution([3, 7],
                               [NormalDistribution(0, 1),
                                NormalDistribution(5, 2)])
    M = MixturePSpace('M', dist).value
    N1, N2 = Normal('N1', 0, 1), Normal('N2', 5, 2)

    assert simplify(density(M)(x) -
                    (S(3)/10*density(N1)(x) + S(7)/10*density(N2)(x))) == 0
    # CDF is the weighted sum of the component CDFs, and stays finite at a
    # component mean (x = 5) where integrating the combined PDF would give 0/0.
    assert simplify(cdf(M)(x) -
                    (S(3)/10*cdf(N1)(x) + S(7)/10*cdf(N2)(x))) == 0
    assert cdf(M)(x).subs(x, 5) == S(3)*erf(5*sqrt(2)/2)/20 + S.Half
    assert simplify(E(M)) == S(7)/2
    # Var = sum w_i (Var_i + E_i**2) - E**2 = 3/10*1 + 7/10*29 - (7/2)**2
    assert simplify(expand_func(variance(M))) == S(167)/20


def test_mixture_pspace_discrete():
    x = Symbol('x')
    dist = MixtureDistribution([1, 1],
                               [PoissonDistribution(1), PoissonDistribution(2)])
    M = MixturePSpace('M', dist).value
    P1, P2 = Poisson('P1', 1), Poisson('P2', 2)

    # E = 1/2*1 + 1/2*2
    assert simplify(E(M)) == S(3)/2
    assert simplify(cdf(M)(x) - (S.Half*cdf(P1)(x) + S.Half*cdf(P2)(x))) == 0


def test_mixture_pspace_finite():
    dist = MixtureDistribution([1, 3],
                               [DieDistribution(4),
                                BernoulliDistribution(S.Half, 0, -1)])
    M = MixturePSpace('M', dist).value

    # E = 1/4*E(Die(4)) + 3/4*E(Bernoulli) = 1/4*5/2 + 3/4*(-1/2)
    assert simplify(E(M)) == S(1)/4
    # P(M > 0) = 1/4*P(Die(4) > 0) + 3/4*P(Bernoulli > 0) = 1/4*1 + 3/4*0
    assert simplify(P(M > 0)) == S(1)/4
    assert dict(cdf(M)) == {-1: S(3)/8, 0: S(3)/4, 1: S(13)/16,
                            2: S(7)/8, 3: S(15)/16, 4: S.One}


def test_mixture_disjoint_supports():
    # E and Var are computed by linearity over components, so each component
    # integrates over its own interval. Integrating the combined PDF over the
    # union support would raise NotImplementedError.
    dist = MixtureDistribution([1, 1],
                               [UniformDistribution(0, 1),
                                UniformDistribution(2, 3)])
    M = MixturePSpace('M', dist).value

    assert E(M) == S(3)/2
    assert simplify(variance(M)) == S(13)/12
    assert simplify(P(M > S(5)/2)) == S(1)/4


def test_mixture_conditional_space():
    # Conditioning needs MixturePSpace.conditional_space; finite and discrete
    # conditionals are exact (summation, not integration).
    MF = Mixture('MF', [1, 1], [Die('DF1', 2), Die('DF2', 4)])
    assert P(MF > 3, MF > 1) == S(1)/5
    assert E(MF, MF > 2) == S(7)/2

    MD = Mixture('MD', [1, 1], [Poisson('Q1', 1), Poisson('Q2', 2)])
    assert simplify(P(MD > 2, MD > 0) - P(MD > 2)/P(MD > 0)) == 0

    # Continuous conditioning is wired up (no AttributeError).
    Mc = Mixture('Mc', [S.Half, S.Half],
                 [Normal('Nc1', 0, 1), Normal('Nc2', 5, 2)])
    assert given(Mc, Mc > 1) is not None


def test_mixture_sample():
    numpy = import_module('numpy')
    scipy = import_module('scipy')
    if numpy is None or scipy is None:
        skip('numpy and scipy are required for sampling')
    M = Mixture('Ms', [1, 1], [Die('Da', 2), Die('Db', 4)])
    s = sample(M, size=(20,), seed=1)
    assert s.shape == (20,)
    # every draw lies in the union support {1, 2, 3, 4}
    assert set(map(int, s)).issubset({1, 2, 3, 4})
    # reproducible with the same seed
    assert (sample(M, size=(20,), seed=1) == s).all()


def test_mixture_characteristic_function():
    t = Symbol('t', real=True)

    # Continuous: CF of an equal-weight mixture of Normals matches the
    # weighted sum of the component CFs.
    M = Mixture('Mc', [1, 1], [Normal('a', 0, 1), Normal('b', 5, 1)])
    N1 = Normal('N1', 0, 1)
    N2 = Normal('N2', 5, 1)
    assert simplify(characteristic_function(M)(t) -
                    (S.Half*characteristic_function(N1)(t) +
                     S.Half*characteristic_function(N2)(t))) == 0

    # Discrete: same property for a Poisson mixture.
    MD = Mixture('Md', [1, 3], [Poisson('P1', 1), Poisson('P2', 2)])
    P1, P2 = Poisson('Q1', 1), Poisson('Q2', 2)
    assert simplify(characteristic_function(MD)(t) -
                    (S(1)/4*characteristic_function(P1)(t) +
                     S(3)/4*characteristic_function(P2)(t))) == 0


def test_mixture_moment_generating_function():
    t = Symbol('t', real=True)
    M = Mixture('Mg', [3, 7], [Normal('a', 0, 1), Normal('b', 2, 3)])
    N1 = Normal('N1', 0, 1)
    N2 = Normal('N2', 2, 3)
    assert simplify(moment_generating_function(M)(t) -
                    (S(3)/10*moment_generating_function(N1)(t) +
                     S(7)/10*moment_generating_function(N2)(t))) == 0


def test_gaussian_mixture():
    # GaussianMixture is equivalent to building the mixture by hand.
    G = GaussianMixture('G', [1, 1], [0, 5], [1, 1])
    M = Mixture('M', [1, 1], [Normal('a', 0, 1), Normal('b', 5, 1)])
    assert E(G) == E(M)
    assert simplify(variance(G) - variance(M)) == 0

    # Symbolic parameters are accepted (just like Normal).
    mu = Symbol('mu', real=True)
    GS = GaussianMixture('GS', [1, 1], [0, mu], [1, 1])
    assert E(GS) == mu/2

    # Length mismatch raises before any RV is built.
    raises(ValueError, lambda: GaussianMixture('X', [1], [0, 5], [1]))
    raises(ValueError, lambda: GaussianMixture('X', [1, 1], [0], [1, 1]))
    raises(ValueError, lambda: GaussianMixture('X', [1, 1], [0, 5], [1]))


def test_mixture_operator_overload():
    N1 = NormalDistribution(0, 1)
    N2 = NormalDistribution(5, 1)
    N3 = NormalDistribution(10, 1)

    # w * D returns a _WeightedDistribution.
    wd = S(3)/10 * N1
    assert isinstance(wd, _WeightedDistribution)
    assert wd.weight == S(3)/10
    assert wd.distribution == N1

    # The target case from issue #18730: w*D + w*D produces a mixture.
    m = S(3)/10*N1 + S(7)/10*N2
    assert isinstance(m, MixtureDistribution)
    assert m.weights == Tuple(S(3)/10, S(7)/10)
    assert tuple(m.components) == (N1, N2)

    # Plain D + D is an equal-weight two-component mixture.
    m2 = N1 + N2
    assert isinstance(m2, MixtureDistribution)
    assert m2.weights == Tuple(S.Half, S.Half)

    # Mixed forms: w*D + D and D + w*D treat the bare distribution as
    # weight one.
    m3 = S(1)/3*N1 + N2
    assert m3.weights == Tuple(S(1)/4, S(3)/4)
    m4 = N1 + S(2)/3*N2
    assert m4.weights == Tuple(S(3)/5, S(2)/5)

    # Three-or-more chaining is not supported by the operator syntax:
    # users should call Mixture(name, weights, components) directly.
    raises(TypeError, lambda: S(1)/3*N1 + S(1)/3*N2 + S(1)/3*N3)

from functools import singledispatch

from sympy.external import import_module
from sympy.stats.crv_types import BetaDistribution, ChiSquaredDistribution, ExponentialDistribution, GammaDistribution, \
    LogNormalDistribution, NormalDistribution, ParetoDistribution, UniformDistribution
from sympy.stats.drv_types import GeometricDistribution, PoissonDistribution, ZetaDistribution
from sympy.stats.frv_types import BinomialDistribution


numpy = import_module('numpy')


@singledispatch
def do_sample_numpy(dist, size):
    return None


# CRV:

@do_sample_numpy.register
def _(dist: BetaDistribution, size):
    return numpy.random.beta(a=float(dist.alpha), b=float(dist.beta), size=size)


@do_sample_numpy.register
def _(dist: ChiSquaredDistribution, size):
    return numpy.random.chisquare(df=float(dist.k), size=size)


@do_sample_numpy.register
def _(dist: ExponentialDistribution, size):
    return numpy.random.exponential(1 / float(dist.rate), size=size)


@do_sample_numpy.register
def _(dist: GammaDistribution, size):
    return numpy.random.gamma(float(dist.k), float(dist.theta), size=size)


@do_sample_numpy.register
def _(dist: LogNormalDistribution, size):
    return numpy.random.lognormal(float(dist.mean), float(dist.std), size=size)


@do_sample_numpy.register
def _(dist: NormalDistribution, size):
    return numpy.random.normal(float(dist.mean), float(dist.std), size=size)


@do_sample_numpy.register
def _(dist: ParetoDistribution, size):
    return (numpy.random.pareto(a=float(dist.alpha), size=size) + 1) * float(dist.xm)


@do_sample_numpy.register
def _(dist: UniformDistribution, size):
    return numpy.random.uniform(low=float(dist.left), high=float(dist.right), size=size)


# DRV:

@do_sample_numpy.register
def _(dist: GeometricDistribution, size):
    return numpy.random.geometric(p=float(dist.p), size=size)


@do_sample_numpy.register
def _(dist: PoissonDistribution, size):
    return numpy.random.poisson(lam=float(dist.lamda), size=size)


@do_sample_numpy.register
def _(dist: ZetaDistribution, size):
    return numpy.random.zipf(a=float(dist.s), size=size)


# FRV:

@do_sample_numpy.register
def _(dist: BinomialDistribution, size):
    return numpy.random.binomial(n=int(dist.n), p=float(dist.p), size=size)

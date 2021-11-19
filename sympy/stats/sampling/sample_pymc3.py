from functools import singledispatch
from sympy.external import import_module
from sympy.stats.crv_types import BetaDistribution, CauchyDistribution, ChiSquaredDistribution, ExponentialDistribution, \
    GammaDistribution, LogNormalDistribution, NormalDistribution, ParetoDistribution, UniformDistribution, \
    GaussianInverseDistribution
from sympy.stats.drv_types import PoissonDistribution, GeometricDistribution, NegativeBinomialDistribution
from sympy.stats.frv_types import BinomialDistribution, BernoulliDistribution


pymc3 = import_module('pymc3')


@singledispatch
def do_sample_pymc3(dist):
    return None


# CRV:

@do_sample_pymc3.register(BetaDistribution) # type: ignore
def _(dist: BetaDistribution):
    return pymc3.Beta('X', alpha=float(dist.alpha), beta=float(dist.beta))


@do_sample_pymc3.register(CauchyDistribution) # type: ignore
def _(dist: CauchyDistribution):
    return pymc3.Cauchy('X', alpha=float(dist.x0), beta=float(dist.gamma))


@do_sample_pymc3.register(ChiSquaredDistribution) # type: ignore
def _(dist: ChiSquaredDistribution):
    return pymc3.ChiSquared('X', nu=float(dist.k))


@do_sample_pymc3.register(ExponentialDistribution) # type: ignore
def _(dist: ExponentialDistribution):
    return pymc3.Exponential('X', lam=float(dist.rate))


@do_sample_pymc3.register(GammaDistribution) # type: ignore
def _(dist: GammaDistribution):
    return pymc3.Gamma('X', alpha=float(dist.k), beta=1 / float(dist.theta))


@do_sample_pymc3.register(LogNormalDistribution) # type: ignore
def _(dist: LogNormalDistribution):
    return pymc3.Lognormal('X', mu=float(dist.mean), sigma=float(dist.std))


@do_sample_pymc3.register(NormalDistribution) # type: ignore
def _(dist: NormalDistribution):
    return pymc3.Normal('X', float(dist.mean), float(dist.std))


@do_sample_pymc3.register(GaussianInverseDistribution) # type: ignore
def _(dist: GaussianInverseDistribution):
    return pymc3.Wald('X', mu=float(dist.mean), lam=float(dist.shape))


@do_sample_pymc3.register(ParetoDistribution) # type: ignore
def _(dist: ParetoDistribution):
    return pymc3.Pareto('X', alpha=float(dist.alpha), m=float(dist.xm))


@do_sample_pymc3.register(UniformDistribution) # type: ignore
def _(dist: UniformDistribution):
    return pymc3.Uniform('X', lower=float(dist.left), upper=float(dist.right))


# DRV:

@do_sample_pymc3.register(GeometricDistribution) # type: ignore
def _(dist: GeometricDistribution):
    return pymc3.Geometric('X', p=float(dist.p))


@do_sample_pymc3.register(NegativeBinomialDistribution) # type: ignore
def _(dist: NegativeBinomialDistribution):
    return pymc3.NegativeBinomial('X', mu=float((dist.p * dist.r) / (1 - dist.p)),
                                  alpha=float(dist.r))


@do_sample_pymc3.register(PoissonDistribution) # type: ignore
def _(dist: PoissonDistribution):
    return pymc3.Poisson('X', mu=float(dist.lamda))


# FRV:

@do_sample_pymc3.register(BernoulliDistribution) # type: ignore
def _(dist: BernoulliDistribution):
    return pymc3.Bernoulli('X', p=float(dist.p))


@do_sample_pymc3.register(BinomialDistribution) # type: ignore
def _(dist: BinomialDistribution):
    return pymc3.Binomial('X', n=int(dist.n), p=float(dist.p))

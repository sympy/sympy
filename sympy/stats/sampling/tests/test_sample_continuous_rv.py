from sympy import exp, Interval, oo, Symbol
from sympy.external import import_module
from sympy.stats import Beta, Chi, Normal, Gamma, Exponential, LogNormal, Pareto, ChiSquared, Uniform, sample, \
    BetaPrime, Cauchy, GammaInverse, GaussianInverse, StudentT, Weibull, density, ContinuousRV
from sympy.testing.pytest import skip, ignore_warnings, raises


def test_sample_numpy():
    distribs_numpy = [
        Beta("B", 1, 1),
        Normal("N", 0, 1),
        Gamma("G", 2, 7),
        Exponential("E", 2),
        LogNormal("LN", 0, 1),
        Pareto("P", 1, 1),
        ChiSquared("CS", 2),
        Uniform("U", 0, 1)
    ]
    size = 3
    numpy = import_module('numpy')
    if not numpy:
        skip('Numpy is not installed. Abort tests for _sample_numpy.')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            for X in distribs_numpy:
                samps = next(sample(X, size=size, library='numpy'))
                for sam in samps:
                    assert sam in X.pspace.domain.set
            raises(NotImplementedError,
                   lambda: next(sample(Chi("C", 1), library='numpy')))
    raises(NotImplementedError,
           lambda: Chi("C", 1).pspace.distribution.sample(library='tensorflow'))


def test_sample_scipy():
    distribs_scipy = [
        Beta("B", 1, 1),
        BetaPrime("BP", 1, 1),
        Cauchy("C", 1, 1),
        Chi("C", 1),
        Normal("N", 0, 1),
        Gamma("G", 2, 7),
        GammaInverse("GI", 1, 1),
        GaussianInverse("GUI", 1, 1),
        Exponential("E", 2),
        LogNormal("LN", 0, 1),
        Pareto("P", 1, 1),
        StudentT("S", 2),
        ChiSquared("CS", 2),
        Uniform("U", 0, 1)
    ]
    size = 3
    numsamples = 5
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests for _sample_scipy.')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            g_sample = list(sample(Gamma("G", 2, 7), size=size, numsamples=numsamples))
            assert len(g_sample) == numsamples
            for X in distribs_scipy:
                samps = next(sample(X, size=size, library='scipy'))
                samps2 = next(sample(X, size=(2, 2), library='scipy'))
                for sam in samps:
                    assert sam in X.pspace.domain.set
                for i in range(2):
                    for j in range(2):
                        assert samps2[i][j] in X.pspace.domain.set


def test_sample_pymc3():
    distribs_pymc3 = [
        Beta("B", 1, 1),
        Cauchy("C", 1, 1),
        Normal("N", 0, 1),
        Gamma("G", 2, 7),
        GaussianInverse("GI", 1, 1),
        Exponential("E", 2),
        LogNormal("LN", 0, 1),
        Pareto("P", 1, 1),
        ChiSquared("CS", 2),
        Uniform("U", 0, 1)
    ]
    size = 3
    pymc3 = import_module('pymc3')
    if not pymc3:
        skip('PyMC3 is not installed. Abort tests for _sample_pymc3.')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            for X in distribs_pymc3:
                samps = next(sample(X, size=size, library='pymc3'))
                for sam in samps:
                    assert sam in X.pspace.domain.set
            raises(NotImplementedError,
                   lambda: next(sample(Chi("C", 1), library='pymc3')))


def test_sampling_gamma_inverse():
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy not installed. Abort tests for sampling of gamma inverse.')
    X = GammaInverse("x", 1, 1)
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        assert next(sample(X)) in X.pspace.domain.set


def test_lognormal_sampling():
    # Right now, only density function and sampling works
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        for i in range(3):
            X = LogNormal('x', i, 1)
            assert next(sample(X)) in X.pspace.domain.set

    size = 5
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        samps = next(sample(X, size=size))
        for samp in samps:
            assert samp in X.pspace.domain.set


def test_sampling_gaussian_inverse():
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy not installed. Abort tests for sampling of Gaussian inverse.')
    X = GaussianInverse("x", 1, 1)
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        assert next(sample(X, library='scipy')) in X.pspace.domain.set


def test_prefab_sampling():
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    N = Normal('X', 0, 1)
    L = LogNormal('L', 0, 1)
    E = Exponential('Ex', 1)
    P = Pareto('P', 1, 3)
    W = Weibull('W', 1, 1)
    U = Uniform('U', 0, 1)
    B = Beta('B', 2, 5)
    G = Gamma('G', 1, 3)

    variables = [N, L, E, P, W, U, B, G]
    niter = 10
    size = 5
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        for var in variables:
            for _ in range(niter):
                assert next(sample(var)) in var.pspace.domain.set
                samps = next(sample(var, size=size))
                for samp in samps:
                    assert samp in var.pspace.domain.set


def test_sample_continuous():
    z = Symbol('z')
    Z = ContinuousRV(z, exp(-z), set=Interval(0, oo))
    assert density(Z)(-1) == 0

    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        assert next(sample(Z)) in Z.pspace.domain.set
    sym, val = list(Z.pspace.sample().items())[0]
    assert sym == Z and val in Interval(0, oo)

    libraries = ['scipy', 'numpy', 'pymc3']
    for lib in libraries:
        try:
            imported_lib = import_module(lib)
            if imported_lib:
                s0, s1, s2 = [], [], []
                s0 = list(sample(Z, numsamples=10, library=lib, seed=0))
                s1 = list(sample(Z, numsamples=10, library=lib, seed=0))
                s2 = list(sample(Z, numsamples=10, library=lib, seed=1))
                assert s0 == s1
                assert s1 != s2
        except NotImplementedError:
            continue

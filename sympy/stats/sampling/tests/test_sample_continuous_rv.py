from sympy.external import import_module
from sympy.stats import Beta, Chi, Normal, Gamma, Exponential, LogNormal, Pareto, ChiSquared, Uniform, sample, \
    BetaPrime, Cauchy, GammaInverse, GaussianInverse, StudentT
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

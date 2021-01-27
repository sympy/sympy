from sympy import S, Rational
from sympy.external import import_module
from sympy.stats import Binomial, sample, Die, FiniteRV, DiscreteUniform, Bernoulli, BetaBinomial, Hypergeometric, \
    Rademacher
from sympy.testing.pytest import skip, ignore_warnings, raises

def test_given_sample():
    X = Die('X', 6)
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        assert next(sample(X, X > 5)) == 6

def test_sample_numpy():
    distribs_numpy = [
        Binomial("B", 5, 0.4),
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
                   lambda: next(sample(Die("D"), library='numpy')))
    raises(NotImplementedError,
           lambda: Die("D").pspace.sample(library='tensorflow'))


def test_sample_scipy():
    distribs_scipy = [
        FiniteRV('F', {1: S.Half, 2: Rational(1, 4), 3: Rational(1, 4)}),
        DiscreteUniform("Y", list(range(5))),
        Die("D"),
        Bernoulli("Be", 0.3),
        Binomial("Bi", 5, 0.4),
        BetaBinomial("Bb", 2, 1, 1),
        Hypergeometric("H", 1, 1, 1),
        Rademacher("R")
    ]

    size = 3
    numsamples = 5
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy not installed. Abort tests for _sample_scipy.')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            h_sample = list(sample(Hypergeometric("H", 1, 1, 1), size=size, numsamples=numsamples))
            assert len(h_sample) == numsamples
            for X in distribs_scipy:
                samps = next(sample(X, size=size))
                samps2 = next(sample(X, size=(2, 2)))
                for sam in samps:
                    assert sam in X.pspace.domain.set
                for i in range(2):
                    for j in range(2):
                        assert samps2[i][j] in X.pspace.domain.set


def test_sample_pymc3():
    distribs_pymc3 = [
        Bernoulli('B', 0.2),
        Binomial('N', 5, 0.4)
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
                   lambda: next(sample(Die("D"), library='pymc3')))


def test_sample_seed():
    F = FiniteRV('F', {1: S.Half, 2: Rational(1, 4), 3: Rational(1, 4)})

    libraries = ['scipy', 'numpy', 'pymc3']
    for lib in libraries:
        try:
            imported_lib = import_module(lib)
            if imported_lib:
                s0 = list(sample(F, numsamples=10, library=lib, seed=0))
                s1 = list(sample(F, numsamples=10, library=lib, seed=0))
                s2 = list(sample(F, numsamples=10, library=lib, seed=1))
                assert s0 == s1
                assert s1 != s2
        except NotImplementedError:
            continue

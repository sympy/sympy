from __future__ import division
from sympy.stats import (P, E, where, density, variance, covariance, skewness,
                         given, pspace, cdf, ContinuousRV, sample,
                         Arcsin, Benini, Beta, BetaPrime, Cauchy,
                         Chi, ChiSquared,
                         ChiNoncentral, Dagum, Erlang, Exponential,
                         FDistribution, FisherZ, Frechet, Gamma, GammaInverse,
                         Gompertz, Gumbel, Kumaraswamy, Laplace, Logistic,
                         LogNormal, Maxwell, Nakagami, Normal, Pareto,
                         QuadraticU, RaisedCosine, Rayleigh, ShiftedGompertz,
                         StudentT, Triangular, Uniform, UniformSum,
                         VonMises, Weibull, WignerSemicircle, correlation,
                         moment, cmoment, smoment)

from sympy import (Symbol, Abs, exp, S, N, pi, simplify, Interval, erf, erfc,
                   Eq, log, lowergamma, Sum, symbols, sqrt, And, gamma, beta,
                   Piecewise, Integral, sin, cos, besseli, factorial, binomial,
                   floor, expand_func)


from sympy.stats.crv_types import NormalDistribution
from sympy.stats.rv import ProductPSpace

from sympy.utilities.pytest import raises, XFAIL, slow

from sympy.core.compatibility import range

oo = S.Infinity

x, y, z = map(Symbol, 'xyz')


def test_single_normal():
    mu = Symbol('mu', real=True, finite=True)
    sigma = Symbol('sigma', real=True, positive=True, finite=True)
    X = Normal('x', 0, 1)
    Y = X*sigma + mu

    assert simplify(E(Y)) == mu
    assert simplify(variance(Y)) == sigma**2
    pdf = density(Y)
    x = Symbol('x')
    assert (pdf(x) ==
            2**S.Half*exp(-(mu - x)**2/(2*sigma**2))/(2*pi**S.Half*sigma))

    assert P(X**2 < 1) == erf(2**S.Half/2)

    assert E(X, Eq(X, mu)) == mu


@XFAIL
def test_conditional_1d():
    X = Normal('x', 0, 1)
    Y = given(X, X >= 0)

    assert density(Y) == 2 * density(X)

    assert Y.pspace.domain.set == Interval(0, oo)
    assert E(Y) == sqrt(2) / sqrt(pi)

    assert E(X**2) == E(Y**2)


def test_ContinuousDomain():
    X = Normal('x', 0, 1)
    assert where(X**2 <= 1).set == Interval(-1, 1)
    assert where(X**2 <= 1).symbol == X.symbol
    where(And(X**2 <= 1, X >= 0)).set == Interval(0, 1)
    raises(ValueError, lambda: where(sin(X) > 1))

    Y = given(X, X >= 0)

    assert Y.pspace.domain.set == Interval(0, oo)


@slow
def test_multiple_normal():
    X, Y = Normal('x', 0, 1), Normal('y', 0, 1)

    assert E(X + Y) == 0
    assert variance(X + Y) == 2
    assert variance(X + X) == 4
    assert covariance(X, Y) == 0
    assert covariance(2*X + Y, -X) == -2*variance(X)
    assert skewness(X) == 0
    assert skewness(X + Y) == 0
    assert correlation(X, Y) == 0
    assert correlation(X, X + Y) == correlation(X, X - Y)
    assert moment(X, 2) == 1
    assert cmoment(X, 3) == 0
    assert moment(X + Y, 4) == 12
    assert cmoment(X, 2) == variance(X)
    assert smoment(X*X, 2) == 1
    assert smoment(X + Y, 3) == skewness(X + Y)
    assert E(X, Eq(X + Y, 0)) == 0
    assert variance(X, Eq(X + Y, 0)) == S.Half


@slow
def test_symbolic():
    mu1, mu2 = symbols('mu1 mu2', real=True, finite=True)
    s1, s2 = symbols('sigma1 sigma2', real=True, finite=True, positive=True)
    rate = Symbol('lambda', real=True, positive=True, finite=True)
    X = Normal('x', mu1, s1)
    Y = Normal('y', mu2, s2)
    Z = Exponential('z', rate)
    a, b, c = symbols('a b c', real=True, finite=True)

    assert E(X) == mu1
    assert E(X + Y) == mu1 + mu2
    assert E(a*X + b) == a*E(X) + b
    assert variance(X) == s1**2
    assert simplify(variance(X + a*Y + b)) == variance(X) + a**2*variance(Y)

    assert E(Z) == 1/rate
    assert E(a*Z + b) == a*E(Z) + b
    assert E(X + a*Z + b) == mu1 + a/rate + b


def test_cdf():
    X = Normal('x', 0, 1)

    d = cdf(X)
    assert P(X < 1) == d(1)
    assert d(0) == S.Half

    d = cdf(X, X > 0)  # given X>0
    assert d(0) == 0

    Y = Exponential('y', 10)
    d = cdf(Y)
    assert d(-5) == 0
    assert P(Y > 3) == 1 - d(3)

    raises(ValueError, lambda: cdf(X + Y))

    Z = Exponential('z', 1)
    f = cdf(Z)
    z = Symbol('z')
    assert f(z) == Piecewise((1 - exp(-z), z >= 0), (0, True))


def test_sample():
    z = Symbol('z')
    Z = ContinuousRV(z, exp(-z), set=Interval(0, oo))
    assert sample(Z) in Z.pspace.domain.set
    sym, val = list(Z.pspace.sample().items())[0]
    assert sym == Z and val in Interval(0, oo)


def test_ContinuousRV():
    x = Symbol('x')
    pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi))  # Normal distribution
    # X and Y should be equivalent
    X = ContinuousRV(x, pdf)
    Y = Normal('y', 0, 1)

    assert variance(X) == variance(Y)
    assert P(X > 0) == P(Y > 0)


def test_arcsin():
    a = Symbol("a", real=True)
    b = Symbol("b", real=True)

    X = Arcsin('x', a, b)
    assert density(X)(x) == 1/(pi*sqrt((-x + b)*(x - a)))


def test_benini():
    alpha = Symbol("alpha", positive=True)
    b = Symbol("beta", positive=True)
    sigma = Symbol("sigma", positive=True)

    X = Benini('x', alpha, b, sigma)
    assert density(X)(x) == ((alpha/x + 2*b*log(x/sigma)/x)
                          *exp(-alpha*log(x/sigma) - b*log(x/sigma)**2))


def test_beta():
    a, b = symbols('alpha beta', positive=True)

    B = Beta('x', a, b)

    assert pspace(B).domain.set == Interval(0, 1)

    dens = density(B)
    x = Symbol('x')
    assert dens(x) == x**(a - 1)*(1 - x)**(b - 1) / beta(a, b)

    # This is too slow
    # assert E(B) == a / (a + b)
    # assert variance(B) == (a*b) / ((a+b)**2 * (a+b+1))

    # Full symbolic solution is too much, test with numeric version
    a, b = 1, 2
    B = Beta('x', a, b)
    assert expand_func(E(B)) == a / S(a + b)
    assert expand_func(variance(B)) == (a*b) / S((a + b)**2 * (a + b + 1))


def test_betaprime():
    alpha = Symbol("alpha", positive=True)
    betap = Symbol("beta", positive=True)

    X = BetaPrime('x', alpha, betap)
    assert density(X)(x) == x**(alpha - 1)*(x + 1)**(-alpha - betap)/beta(alpha, betap)


def test_cauchy():
    x0 = Symbol("x0")
    gamma = Symbol("gamma", positive=True)

    X = Cauchy('x', x0, gamma)
    assert density(X)(x) == 1/(pi*gamma*(1 + (x - x0)**2/gamma**2))


def test_chi():
    k = Symbol("k", integer=True)

    X = Chi('x', k)
    assert density(X)(x) == 2**(-k/2 + 1)*x**(k - 1)*exp(-x**2/2)/gamma(k/2)

def test_chi_noncentral():
    k = Symbol("k", integer=True)
    l = Symbol("l")

    X = ChiNoncentral("x", k, l)
    assert density(X)(x) == (x**k*l*(x*l)**(-k/2)*
                          exp(-x**2/2 - l**2/2)*besseli(k/2 - 1, x*l))

def test_chi_squared():
    k = Symbol("k", integer=True)

    X = ChiSquared('x', k)
    assert density(X)(x) == 2**(-k/2)*x**(k/2 - 1)*exp(-x/2)/gamma(k/2)

def test_dagum():
    p = Symbol("p", positive=True)
    b = Symbol("b", positive=True)
    a = Symbol("a", positive=True)

    X = Dagum('x', p, a, b)
    assert density(X)(x) == a*p*(x/b)**(a*p)*((x/b)**a + 1)**(-p - 1)/x

def test_erlang():
    k = Symbol("k", integer=True, positive=True)
    l = Symbol("l", positive=True)

    X = Erlang("x", k, l)
    assert density(X)(x) == x**(k - 1)*l**k*exp(-x*l)/gamma(k)

def test_exponential():
    rate = Symbol('lambda', positive=True, real=True, finite=True)
    X = Exponential('x', rate)

    assert E(X) == 1/rate
    assert variance(X) == 1/rate**2
    assert skewness(X) == 2
    assert skewness(X) == smoment(X, 3)
    assert smoment(2*X, 4) == smoment(X, 4)
    assert moment(X, 3) == 3*2*1/rate**3
    assert P(X > 0) == S(1)
    assert P(X > 1) == exp(-rate)
    assert P(X > 10) == exp(-10*rate)

    assert where(X <= 1).set == Interval(0, 1)

def test_f_distribution():
    d1 = Symbol("d1", positive=True)
    d2 = Symbol("d2", positive=True)

    X = FDistribution("x", d1, d2)
    assert density(X)(x) == (d2**(d2/2)*sqrt((d1*x)**d1*(d1*x + d2)**(-d1 - d2))
                             /(x*beta(d1/2, d2/2)))

def test_fisher_z():
    d1 = Symbol("d1", positive=True)
    d2 = Symbol("d2", positive=True)

    X = FisherZ("x", d1, d2)
    assert density(X)(x) == (2*d1**(d1/2)*d2**(d2/2)*(d1*exp(2*x) + d2)
                             **(-d1/2 - d2/2)*exp(d1*x)/beta(d1/2, d2/2))

def test_frechet():
    a = Symbol("a", positive=True)
    s = Symbol("s", positive=True)
    m = Symbol("m", real=True)

    X = Frechet("x", a, s=s, m=m)
    assert density(X)(x) == a*((x - m)/s)**(-a - 1)*exp(-((x - m)/s)**(-a))/s

def test_gamma():
    k = Symbol("k", positive=True)
    theta = Symbol("theta", positive=True)

    X = Gamma('x', k, theta)
    assert density(X)(x) == x**(k - 1)*theta**(-k)*exp(-x/theta)/gamma(k)
    assert cdf(X, meijerg=True)(z) == Piecewise(
            (-k*lowergamma(k, 0)/gamma(k + 1) +
                k*lowergamma(k, z/theta)/gamma(k + 1), z >= 0),
            (0, True))
    # assert simplify(variance(X)) == k*theta**2  # handled numerically below
    assert E(X) == moment(X, 1)

    k, theta = symbols('k theta', real=True, finite=True, positive=True)
    X = Gamma('x', k, theta)
    assert simplify(E(X)) == k*theta
    # can't get things to simplify on this one so we use subs
    assert variance(X).subs(k, 5) == (k*theta**2).subs(k, 5)
    # The following is too slow
    # assert simplify(skewness(X)).subs(k, 5) == (2/sqrt(k)).subs(k, 5)

def test_gamma_inverse():
    a = Symbol("a", positive=True)
    b = Symbol("b", positive=True)

    X = GammaInverse("x", a, b)
    assert density(X)(x) == x**(-a - 1)*b**a*exp(-b/x)/gamma(a)

def test_gompertz():
    b = Symbol("b", positive=True)
    eta = Symbol("eta", positive=True)

    X = Gompertz("x", b, eta)
    assert density(X)(x) == b*eta*exp(eta)*exp(b*x)*exp(-eta*exp(b*x))

def test_gumbel():
    beta = Symbol("beta", positive=True)
    mu = Symbol("mu")
    x = Symbol("x")
    X = Gumbel("x", beta, mu)
    assert simplify(density(X)(x)) == exp((beta*exp((mu - x)/beta) + mu - x)/beta)/beta

def test_kumaraswamy():
    a = Symbol("a", positive=True)
    b = Symbol("b", positive=True)

    X = Kumaraswamy("x", a, b)
    assert density(X)(x) == x**(a - 1)*a*b*(-x**a + 1)**(b - 1)

def test_laplace():
    mu = Symbol("mu")
    b = Symbol("b", positive=True)

    X = Laplace('x', mu, b)
    assert density(X)(x) == exp(-Abs(x - mu)/b)/(2*b)

def test_logistic():
    mu = Symbol("mu", real=True)
    s = Symbol("s", positive=True)

    X = Logistic('x', mu, s)
    assert density(X)(x) == exp((-x + mu)/s)/(s*(exp((-x + mu)/s) + 1)**2)

def test_lognormal():
    mean = Symbol('mu', real=True, finite=True)
    std = Symbol('sigma', positive=True, real=True, finite=True)
    X = LogNormal('x', mean, std)
    # The sympy integrator can't do this too well
    #assert E(X) == exp(mean+std**2/2)
    #assert variance(X) == (exp(std**2)-1) * exp(2*mean + std**2)

    # Right now, only density function and sampling works
    # Test sampling: Only e^mean in sample std of 0
    for i in range(3):
        X = LogNormal('x', i, 0)
        assert S(sample(X)) == N(exp(i))
    # The sympy integrator can't do this too well
    #assert E(X) ==

    mu = Symbol("mu", real=True)
    sigma = Symbol("sigma", positive=True)

    X = LogNormal('x', mu, sigma)
    assert density(X)(x) == (sqrt(2)*exp(-(-mu + log(x))**2
                                    /(2*sigma**2))/(2*x*sqrt(pi)*sigma))

    X = LogNormal('x', 0, 1)  # Mean 0, standard deviation 1
    assert density(X)(x) == sqrt(2)*exp(-log(x)**2/2)/(2*x*sqrt(pi))

def test_maxwell():
    a = Symbol("a", positive=True)

    X = Maxwell('x', a)

    assert density(X)(x) == (sqrt(2)*x**2*exp(-x**2/(2*a**2))/
        (sqrt(pi)*a**3))
    assert E(X) == 2*sqrt(2)*a/sqrt(pi)
    assert simplify(variance(X)) == a**2*(-8 + 3*pi)/pi


def test_nakagami():
    mu = Symbol("mu", positive=True)
    omega = Symbol("omega", positive=True)

    X = Nakagami('x', mu, omega)
    assert density(X)(x) == (2*x**(2*mu - 1)*mu**mu*omega**(-mu)
                                *exp(-x**2*mu/omega)/gamma(mu))
    assert simplify(E(X, meijerg=True)) == (sqrt(mu)*sqrt(omega)
           *gamma(mu + S.Half)/gamma(mu + 1))
    assert simplify(variance(X, meijerg=True)) == (
    omega - omega*gamma(mu + S(1)/2)**2/(gamma(mu)*gamma(mu + 1)))


def test_pareto():
    xm, beta = symbols('xm beta', positive=True, finite=True)
    alpha = beta + 5
    X = Pareto('x', xm, alpha)

    dens = density(X)
    x = Symbol('x')
    assert dens(x) == x**(-(alpha + 1))*xm**(alpha)*(alpha)

    # These fail because SymPy can not deduce that 1/xm != 0
    # assert simplify(E(X)) == alpha*xm/(alpha-1)
    # assert simplify(variance(X)) == xm**2*alpha / ((alpha-1)**2*(alpha-2))


def test_pareto_numeric():
    xm, beta = 3, 2
    alpha = beta + 5
    X = Pareto('x', xm, alpha)

    assert E(X) == alpha*xm/S(alpha - 1)
    assert variance(X) == xm**2*alpha / S(((alpha - 1)**2*(alpha - 2)))
    # Skewness tests too slow. Try shortcutting function?


def test_raised_cosine():
    mu = Symbol("mu", real=True)
    s = Symbol("s", positive=True)

    X = RaisedCosine("x", mu, s)
    assert density(X)(x) == (Piecewise(((cos(pi*(x - mu)/s) + 1)/(2*s),
                          And(x <= mu + s, mu - s <= x)), (0, True)))


def test_rayleigh():
    sigma = Symbol("sigma", positive=True)

    X = Rayleigh('x', sigma)
    assert density(X)(x) ==  x*exp(-x**2/(2*sigma**2))/sigma**2
    assert E(X) == sqrt(2)*sqrt(pi)*sigma/2
    assert variance(X) == -pi*sigma**2/2 + 2*sigma**2

def test_shiftedgompertz():
    b = Symbol("b", positive=True)
    eta = Symbol("eta", positive=True)
    X = ShiftedGompertz("x", b, eta)
    assert density(X)(x) == b*(eta*(1 - exp(-b*x)) + 1)*exp(-b*x)*exp(-eta*exp(-b*x))

def test_studentt():
    nu = Symbol("nu", positive=True)

    X = StudentT('x', nu)
    assert density(X)(x) == (1 + x**2/nu)**(-nu/2 - 1/2)/(sqrt(nu)*beta(1/2, nu/2))


@XFAIL
def test_triangular():
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")

    X = Triangular('x', a, b, c)
    assert density(X)(x) == Piecewise(
                 ((2*x - 2*a)/((-a + b)*(-a + c)), And(a <= x, x < c)),
                 (2/(-a + b), x == c),
                 ((-2*x + 2*b)/((-a + b)*(b - c)), And(x <= b, c < x)),
                 (0, True))


def test_quadratic_u():
    a = Symbol("a", real=True)
    b = Symbol("b", real=True)

    X = QuadraticU("x", a, b)
    assert density(X)(x) == (Piecewise((12*(x - a/2 - b/2)**2/(-a + b)**3,
                          And(x <= b, a <= x)), (0, True)))

def test_uniform():
    l = Symbol('l', real=True, finite=True)
    w = Symbol('w', positive=True, finite=True)
    X = Uniform('x', l, l + w)

    assert simplify(E(X)) == l + w/2
    assert simplify(variance(X)) == w**2/12


    # With numbers all is well
    X = Uniform('x', 3, 5)
    assert P(X < 3) == 0 and P(X > 5) == 0
    assert P(X < 4) == P(X > 4) == S.Half


def test_uniform_P():
    """ This stopped working because SingleContinuousPSpace.compute_density no
    longer calls integrate on a DiracDelta but rather just solves directly.
    integrate used to call UniformDistribution.expectation which special-cased
    subsed out the Min and Max terms that Uniform produces

    I decided to regress on this class for general cleanliness (and I suspect
    speed) of the algorithm.
    """
    l = Symbol('l', real=True, finite=True)
    w = Symbol('w', positive=True, finite=True)
    X = Uniform('x', l, l + w)
    assert P(X < l) == 0 and P(X > l + w) == 0


@XFAIL
def test_uniformsum():
    n = Symbol("n", integer=True)
    _k = Symbol("k")

    X = UniformSum('x', n)
    assert density(X)(x) == (Sum((-1)**_k*(-_k + x)**(n - 1)
                        *binomial(n, _k), (_k, 0, floor(x)))/factorial(n - 1))


def test_von_mises():
    mu = Symbol("mu")
    k = Symbol("k", positive=True)

    X = VonMises("x", mu, k)
    assert density(X)(x) == exp(k*cos(x - mu))/(2*pi*besseli(0, k))


def test_weibull():
    a, b = symbols('a b', positive=True)
    X = Weibull('x', a, b)

    assert simplify(E(X)) == simplify(a * gamma(1 + 1/b))
    assert simplify(variance(X)) == simplify(a**2 * gamma(1 + 2/b) - E(X)**2)
    # Skewness tests too slow. Try shortcutting function?


def test_weibull_numeric():
    # Test for integers and rationals
    a = 1
    bvals = [S.Half, 1, S(3)/2, 5]
    for b in bvals:
        X = Weibull('x', a, b)
        assert simplify(E(X)) == simplify(a * gamma(1 + 1/S(b)))
        assert simplify(variance(X)) == simplify(
            a**2 * gamma(1 + 2/S(b)) - E(X)**2)
        # Not testing Skew... it's slow with int/frac values > 3/2


def test_wignersemicircle():
    R = Symbol("R", positive=True)

    X = WignerSemicircle('x', R)
    assert density(X)(x) == 2*sqrt(-x**2 + R**2)/(pi*R**2)
    assert E(X) == 0


def test_prefab_sampling():
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
    for var in variables:
        for i in range(niter):
            assert sample(var) in var.pspace.domain.set


def test_input_value_assertions():
    a, b = symbols('a b')
    p, q = symbols('p q', positive=True)
    m, n = symbols('m n', positive=False, real=True)

    raises(ValueError, lambda: Normal('x', 3, 0))
    raises(ValueError, lambda: Normal('x', m, n))
    Normal('X', a, p)  # No error raised
    raises(ValueError, lambda: Exponential('x', m))
    Exponential('Ex', p)  # No error raised
    for fn in [Pareto, Weibull, Beta, Gamma]:
        raises(ValueError, lambda: fn('x', m, p))
        raises(ValueError, lambda: fn('x', p, n))
        fn('x', p, q)  # No error raised


@XFAIL
def test_unevaluated():
    X = Normal('x', 0, 1)
    assert E(X, evaluate=False) == (
        Integral(sqrt(2)*x*exp(-x**2/2)/(2*sqrt(pi)), (x, -oo, oo)))

    assert E(X + 1, evaluate=False) == (
        Integral(sqrt(2)*x*exp(-x**2/2)/(2*sqrt(pi)), (x, -oo, oo)) + 1)

    assert P(X > 0, evaluate=False) == (
        Integral(sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)), (x, 0, oo)))

    assert P(X > 0, X**2 < 1, evaluate=False) == (
        Integral(sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)*
            Integral(sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)),
                (x, -1, 1))), (x, 0, 1)))


def test_probability_unevaluated():
    T = Normal('T', 30, 3)
    assert type(P(T > 33, evaluate=False)) == Integral

def test_density_unevaluated():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 2)
    assert isinstance(density(X+Y, evaluate=False)(z), Integral)


def test_NormalDistribution():
    nd = NormalDistribution(0, 1)
    x = Symbol('x')
    assert nd.cdf(x) == (1 - erfc(sqrt(2)*x/2))/2 + S.One/2
    assert isinstance(nd.sample(), float) or nd.sample().is_Number
    assert nd.expectation(1, x) == 1
    assert nd.expectation(x, x) == 0
    assert nd.expectation(x**2, x) == 1

def test_random_parameters():
    mu = Normal('mu', 2, 3)
    meas = Normal('T', mu, 1)
    assert density(meas, evaluate=False)(z)
    assert isinstance(pspace(meas), ProductPSpace)
    #assert density(meas, evaluate=False)(z) == Integral(mu.pspace.pdf *
    #        meas.pspace.pdf, (mu.symbol, -oo, oo)).subs(meas.symbol, z)

def test_random_parameters_given():
    mu = Normal('mu', 2, 3)
    meas = Normal('T', mu, 1)
    assert given(meas, Eq(mu, 5)) == Normal('T', 5, 1)

def test_conjugate_priors():
    mu = Normal('mu', 2, 3)
    x = Normal('x', mu, 1)
    assert isinstance(simplify(density(mu, Eq(x, y), evaluate=False)(z)),
            Integral)

def test_difficult_univariate():
    """ Since using solve in place of deltaintegrate we're able to perform
    substantially more complex density computations on single continuous random
    variables """
    x = Normal('x', 0, 1)
    assert density(x**3)
    assert density(exp(x**2))
    assert density(log(x))


def test_issue_10003():
    X = Exponential('x', 3)
    G = Gamma('g', 1, 2)
    assert P(X < -1) == S.Zero
    assert P(G < -1) == S.Zero

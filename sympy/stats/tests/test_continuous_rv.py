from sympy.stats import (P, E, where, density, variance, covariance, skewness,
                         given, pspace, cdf, ContinuousRV, sample)
from sympy.stats import (Arcsin, Benini, Beta, BetaPrime, Cauchy, Chi, Dagum,
                         Exponential, Gamma, Laplace, Logistic, LogNormal,
                         Maxwell, Nakagami, Normal, Pareto, Rayleigh, StudentT,
                         Triangular, Uniform, UniformSum, Weibull,
                         WignerSemicircle)
from sympy import (Symbol, Dummy, Abs, exp, S, N, pi, simplify, Interval, erf,
                   Eq, log, lowergamma, Sum, symbols, sqrt, And, gamma, beta,
                   Piecewise, Integral, sin, Lambda, factorial, binomial, floor)
from sympy.utilities.pytest import raises, XFAIL

oo = S.Infinity

_x = Dummy("x")
_z = Dummy("z")

def test_single_normal():
    mu = Symbol('mu', real=True, bounded=True)
    sigma = Symbol('sigma', real=True, positive=True, bounded=True)
    X = Normal('x', 0,1)
    Y = X*sigma + mu

    assert simplify(E(Y)) == mu
    assert simplify(variance(Y)) == sigma**2
    pdf = density(Y)
    x = Symbol('x')
    assert (pdf(x) ==
            2**S.Half*exp(-(x - mu)**2/(2*sigma**2))/(2*pi**S.Half*sigma))

    assert P(X**2 < 1) == erf(2**S.Half/2)

    assert E(X, Eq(X, mu)) == mu

@XFAIL
def test_conditional_1d():
    X = Normal('x', 0,1)
    Y = given(X, X>=0)

    assert density(Y) == 2 * density(X)

    assert Y.pspace.domain.set == Interval(0, oo)
    assert E(Y) == sqrt(2) / sqrt(pi)

    assert E(X**2) == E(Y**2)

def test_ContinuousDomain():
    X = Normal('x', 0,1)
    assert where(X**2<=1).set == Interval(-1,1)
    assert where(X**2<=1).symbol == X.symbol
    where(And(X**2<=1, X>=0)).set == Interval(0,1)
    raises(ValueError, lambda: where(sin(X)>1))

    Y = given(X, X>=0)

    assert Y.pspace.domain.set == Interval(0, oo)

def test_multiple_normal():
    X, Y = Normal('x', 0,1), Normal('y', 0,1)

    assert E(X+Y) == 0
    assert variance(X+Y) == 2
    assert variance(X+X) == 4
    assert covariance(X, Y) == 0
    assert covariance(2*X + Y, -X) == -2*variance(X)

    assert E(X, Eq(X+Y, 0)) == 0
    assert variance(X, Eq(X+Y, 0)) == S.Half

def test_symbolic():
    mu1, mu2 = symbols('mu1 mu2', real=True, bounded=True)
    s1, s2 = symbols('sigma1 sigma2', real=True, bounded=True, positive=True)
    rate = Symbol('lambda', real=True, positive=True, bounded=True)
    X = Normal('x', mu1, s1)
    Y = Normal('y', mu2, s2)
    Z = Exponential('z', rate)
    a, b, c = symbols('a b c', real=True, bounded=True)

    assert E(X) == mu1
    assert E(X+Y) == mu1+mu2
    assert E(a*X+b) == a*E(X)+b
    assert variance(X) == s1**2
    assert simplify(variance(X+a*Y+b)) == variance(X) + a**2*variance(Y)

    assert E(Z) == 1/rate
    assert E(a*Z+b) == a*E(Z)+b
    assert E(X+a*Z+b) == mu1 + a/rate + b

def test_cdf():
    X = Normal('x', 0,1)

    d = cdf(X)
    assert P(X<1) == d(1)
    assert d(0) == S.Half

    d = cdf(X, X>0) # given X>0
    assert d(0) == 0

    Y = Exponential('y', 10)
    d = cdf(Y)
    assert d(-5) == 0
    assert P(Y > 3) == 1 - d(3)

    raises(ValueError, lambda: cdf(X+Y))

    Z = Exponential('z', 1)
    f = cdf(Z)
    z = Symbol('z')
    assert f(z) == Piecewise((1 - exp(-z), z >= 0), (0, True))

def test_sample():
    z = Symbol('z')
    Z = ContinuousRV(z, exp(-z), set=Interval(0,oo))
    assert sample(Z) in Z.pspace.domain.set
    sym, val = Z.pspace.sample().items()[0]
    assert sym == Z and val in Interval(0, oo)

def test_ContinuousRV():
    x = Symbol('x')
    pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)) # Normal distribution
    # X and Y should be equivalent
    X = ContinuousRV(x, pdf)
    Y = Normal('y', 0, 1)

    assert variance(X) == variance(Y)
    assert P(X>0) == P(Y>0)


def test_arcsin():
    a = Symbol("a", real=True)
    b = Symbol("b", real=True)

    X = Arcsin('x', a, b)
    assert density(X) == Lambda(_x, 1/(pi*sqrt((-_x + b)*(_x - a))))


def test_benini():
    alpha = Symbol("alpha", positive=True)
    b = Symbol("beta", positive=True)
    sigma = Symbol("sigma", positive=True)

    X = Benini('x', alpha, b, sigma)
    assert density(X) == (Lambda(_x, (alpha/_x + 2*b*log(_x/sigma)/_x)
                          *exp(-alpha*log(_x/sigma) - b*log(_x/sigma)**2)))


def test_beta():
    a, b = symbols('alpha beta', positive=True)

    B = Beta('x', a, b)

    assert pspace(B).domain.set == Interval(0, 1)

    dens = density(B)
    x = Symbol('x')
    assert dens(x) == x**(a-1)*(1-x)**(b-1) / beta(a,b)

    # This is too slow
    # assert E(B) == a / (a + b)
    # assert variance(B) == (a*b) / ((a+b)**2 * (a+b+1))

    # Full symbolic solution is too much, test with numeric version
    a, b = 1, 2
    B = Beta('x', a, b)
    assert E(B) == a / S(a + b)
    assert variance(B) == (a*b) / S((a+b)**2 * (a+b+1))


def test_betaprime():
    alpha = Symbol("alpha", positive=True)
    beta = Symbol("beta", positive=True)

    X = BetaPrime('x', alpha, beta)
    assert density(X) == (Lambda(_x, _x**(alpha - 1)*(_x + 1)**(-alpha - beta)
                          *gamma(alpha + beta)/(gamma(alpha)*gamma(beta))))


def test_cauchy():
    x0 = Symbol("x0")
    gamma = Symbol("gamma", positive=True)

    X = Cauchy('x', x0, gamma)
    assert density(X) == Lambda(_x, 1/(pi*gamma*(1 + (_x - x0)**2/gamma**2)))


def test_chi():
    k = Symbol("k", integer=True)

    X = Chi('x', k)
    assert density(X) == (Lambda(_x, 2**(-k/2 + 1)*_x**(k - 1)
                          *exp(-_x**2/2)/gamma(k/2)))


def test_dagum():
    p = Symbol("p", positive=True)
    b = Symbol("b", positive=True)
    a = Symbol("a", positive=True)

    X = Dagum('x', p, a, b)
    assert density(X) == Lambda(_x,
                                a*p*(_x/b)**(a*p)*((_x/b)**a + 1)**(-p - 1)/_x)


def test_exponential():
    rate = Symbol('lambda', positive=True, real=True, bounded=True)
    X = Exponential('x', rate)

    assert E(X) == 1/rate
    assert variance(X) == 1/rate**2
    assert skewness(X) == 2
    assert P(X>0) == S(1)
    assert P(X>1) == exp(-rate)
    assert P(X>10) == exp(-10*rate)

    assert where(X<=1).set == Interval(0,1)


def test_gamma():
    k = Symbol("k", positive=True)
    theta = Symbol("theta", positive=True)

    X = Gamma('x', k, theta)
    assert density(X) == Lambda(_x,
                                _x**(k - 1)*theta**(-k)*exp(-_x/theta)/gamma(k))
    assert cdf(X, meijerg=True) == Lambda(_z, Piecewise(
    (-k*lowergamma(k, 0)/gamma(k + 1) + k*lowergamma(k, _z/theta)/gamma(k + 1), _z >= 0), (0, True)))
    assert variance(X) == (-theta**2*gamma(k + 1)**2/gamma(k)**2 +
           theta*theta**(-k)*theta**(k + 1)*gamma(k + 2)/gamma(k))

    k, theta = symbols('k theta', real=True, bounded=True, positive=True)
    X = Gamma('x', k, theta)
    assert simplify(E(X)) == k*theta
    # can't get things to simplify on this one so we use subs
    assert variance(X).subs(k,5) == (k*theta**2).subs(k, 5)
    # The following is too slow
    # assert simplify(skewness(X)).subs(k, 5) == (2/sqrt(k)).subs(k, 5)


def test_laplace():
    mu = Symbol("mu")
    b = Symbol("b", positive=True)

    X = Laplace('x', mu, b)
    assert density(X) == Lambda(_x, exp(-Abs(_x - mu)/b)/(2*b))


def test_logistic():
    mu = Symbol("mu", real=True)
    s = Symbol("s", positive=True)

    X = Logistic('x', mu, s)
    assert density(X) == Lambda(_x,
                                exp((-_x + mu)/s)/(s*(exp((-_x + mu)/s) + 1)**2))


def test_lognormal():
    mean = Symbol('mu', real=True, bounded=True)
    std = Symbol('sigma', positive=True, real=True, bounded=True)
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
    assert density(X) == (Lambda(_x, sqrt(2)*exp(-(-mu + log(_x))**2
                                    /(2*sigma**2))/(2*_x*sqrt(pi)*sigma)))

    X = LogNormal('x', 0, 1) # Mean 0, standard deviation 1
    assert density(X) == Lambda(_x, sqrt(2)*exp(-log(_x)**2/2)/(2*_x*sqrt(pi)))


def test_maxwell():
    a = Symbol("a", positive=True)

    X = Maxwell('x', a)

    assert density(X) == (Lambda(_x, sqrt(2)*_x**2*exp(-_x**2/(2*a**2))/
        (sqrt(pi)*a**3)))
    assert E(X) == 2*sqrt(2)*a/sqrt(pi)
    assert simplify(variance(X)) == a**2*(-8 + 3*pi)/pi

def test_nakagami():
    mu = Symbol("mu", positive=True)
    omega = Symbol("omega", positive=True)

    X = Nakagami('x', mu, omega)
    assert density(X) == (Lambda(_x, 2*_x**(2*mu - 1)*mu**mu*omega**(-mu)
                                *exp(-_x**2*mu/omega)/gamma(mu)))
    assert simplify(E(X, meijerg=True)) == (sqrt(mu)*sqrt(omega)
           *gamma(mu + S.Half)/gamma(mu + 1))
    assert (simplify(variance(X, meijerg=True)) ==
                            (omega*(gamma(mu)*gamma(mu + 1)
                          - gamma(mu + S.Half)**2)/(gamma(mu)*gamma(mu + 1))))

def test_pareto():
    xm, beta = symbols('xm beta', positive=True, bounded=True)
    alpha = beta + 5
    X = Pareto('x', xm, alpha)

    dens = density(X)
    x = Symbol('x')
    assert dens(x) == x**(-(alpha+1))*xm**(alpha)*(alpha)

    # These fail because SymPy can not deduce that 1/xm != 0
    # assert simplify(E(X)) == alpha*xm/(alpha-1)
    # assert simplify(variance(X)) == xm**2*alpha / ((alpha-1)**2*(alpha-2))

def test_pareto_numeric():
    xm, beta = 3, 2
    alpha = beta + 5
    X = Pareto('x', xm, alpha)

    assert E(X) == alpha*xm/S(alpha-1)
    assert variance(X) == xm**2*alpha / S(((alpha-1)**2*(alpha-2)))
    # Skewness tests too slow. Try shortcutting function?

def test_rayleigh():
    sigma = Symbol("sigma", positive=True)

    X = Rayleigh('x', sigma)
    assert density(X) == Lambda(_x, _x*exp(-_x**2/(2*sigma**2))/sigma**2)
    assert E(X) == sqrt(2)*sqrt(pi)*sigma/2
    assert variance(X) == -pi*sigma**2/2 + 2*sigma**2

def test_studentt():
    nu = Symbol("nu", positive=True)

    X = StudentT('x', nu)
    assert density(X) == (Lambda(_x, (_x**2/nu + 1)**(-nu/2 - S.Half)
                        *gamma(nu/2 + S.Half)/(sqrt(pi)*sqrt(nu)*gamma(nu/2))))

@XFAIL
def test_triangular():
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")

    X = Triangular('x', a,b,c)
    assert Density(X) == Lambda(_x,
             Piecewise(((2*_x - 2*a)/((-a + b)*(-a + c)), And(a <= _x, _x < c)),
                       (2/(-a + b), _x == c),
                       ((-2*_x + 2*b)/((-a + b)*(b - c)), And(_x <= b, c < _x)),
                       (0, True)))

def test_uniform():
    l = Symbol('l', real=True, bounded=True)
    w = Symbol('w', positive=True, bounded=True)
    X = Uniform('x', l, l+w)

    assert simplify(E(X)) == l + w/2
    assert simplify(variance(X)) == w**2/12

    assert P(X<l) == 0 and P(X>l+w) == 0

    # With numbers all is well
    X = Uniform('x', 3, 5)
    assert P(X<3) == 0 and P(X>5) == 0
    assert P(X<4) == P(X>4) == S.Half


@XFAIL
def test_uniformsum():
    n = Symbol("n", integer=True)
    _k = Symbol("k")

    X = UniformSum('x', n)
    assert density(X) == (Lambda(_x, Sum((-1)**_k*(-_k + _x)**(n - 1)
                        *binomial(n, _k), (_k, 0, floor(_x)))/factorial(n - 1)))

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
    assert density(X) == Lambda(_x, 2*sqrt(-_x**2 + R**2)/(pi*R**2))
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

    variables = [N,L,E,P,W,U,B,G]
    niter = 10
    for var in variables:
        for i in xrange(niter):
            assert sample(var) in var.pspace.domain.set

def test_input_value_assertions():
    a, b = symbols('a b')
    p, q = symbols('p q', positive=True)

    raises(ValueError, lambda: Normal('x', 3, 0))
    raises(ValueError, lambda: Normal('x', a, b))
    Normal('X', a, p) # No error raised
    raises(ValueError, lambda: Exponential('x', a))
    Exponential('Ex', p) # No error raised
    for fn in [Pareto, Weibull, Beta, Gamma]:
        raises(ValueError, lambda: fn('x', a, p))
        raises(ValueError, lambda: fn('x', p, a))
        fn('x', p, q) # No error raised

@XFAIL
def test_unevaluated():
    X = Normal('x', 0,1)
    assert E(X, evaluate=False) == (
            Integral(sqrt(2)*x*exp(-x**2/2)/(2*sqrt(pi)), (x, -oo, oo)))

    assert E(X+1, evaluate=False) == (
            Integral(sqrt(2)*x*exp(-x**2/2)/(2*sqrt(pi)), (x, -oo, oo)) + 1)

    assert P(X>0, evaluate=False) == (
            Integral(sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)), (x, 0, oo)))

    assert P(X>0, X**2<1, evaluate=False) == (
            Integral(sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)*
            Integral(sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)),
                (x, -1, 1))), (x, 0, 1)))

def test_probability_unevaluated():
     T = Normal('T', 30, 3)
     assert type(P(T>33, evaluate=False)) == Integral

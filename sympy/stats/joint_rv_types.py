from sympy import (sympify, S, pi, sqrt, exp, Lambda, Indexed, Gt, IndexedBase,
                    besselk, gamma, Interval, Range, factorial, Mul, Integer,
                    Add, rf, Eq, Piecewise, Symbol, imageset, Intersection)
from sympy.matrices import ImmutableMatrix
from sympy.matrices.expressions.determinant import det
from sympy.stats.joint_rv import (JointDistribution, JointPSpace,
    JointDistributionHandmade, MarginalDistribution)
from sympy.stats.rv import _value_check, random_symbols

# __all__ = ['MultivariateNormal',
# 'MultivariateLaplace',
# 'MultivariateT',
# 'NormalGamma'
# ]

def multivariate_rv(cls, sym, *args):
    args = list(map(sympify, args))
    dist = cls(*args)
    args = dist.args
    dist.check(*args)
    return JointPSpace(sym, dist).value

def JointRV(symbol, pdf, _set=None):
    """
    Create a Joint Random Variable where each of its component is conitinuous,
    given the following:

    -- a symbol
    -- a PDF in terms of indexed symbols of the symbol given
     as the first argument

    NOTE: As of now, the set for each component for a `JointRV` is
    equal to the set of all integers, which can not be changed.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy import symbols, exp, pi, Indexed, S
    >>> from sympy.stats import density
    >>> from sympy.stats.joint_rv_types import JointRV

    >>> x1, x2 = (Indexed('x', i) for i in (1, 2))
    >>> pdf = exp(-x1**2/2 + x1 - x2**2/2 - S(1)/2)/(2*pi)

    >>> N1 = JointRV('x', pdf) #Multivariate Normal distribution
    >>> density(N1)(1, 2)
    exp(-2)/(2*pi)
    """
    #TODO: Add support for sets provided by the user
    symbol = sympify(symbol)
    syms = list(i for i in pdf.free_symbols if isinstance(i, Indexed)
        and i.base == IndexedBase(symbol))
    syms.sort(key = lambda index: index.args[1])
    _set = S.Reals**len(syms)
    pdf = Lambda(syms, pdf)
    dist = JointDistributionHandmade(pdf, _set)
    jrv = JointPSpace(symbol, dist).value
    rvs = random_symbols(pdf)
    if len(rvs) != 0:
        dist = MarginalDistribution(dist, (jrv,))
        return JointPSpace(symbol, dist).value
    return jrv

#-------------------------------------------------------------------------------
# Multivariate Normal distribution ---------------------------------------------------------

class MultivariateNormalDistribution(JointDistribution):
    _argnames = ['mu', 'sigma']

    is_Continuous=True

    @property
    def set(self):
        k = len(self.mu)
        return S.Reals**k

    def check(self, mu, sigma):
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the mean vector and covariance matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([Gt(i, 0) != False for i in sigma.eigenvals().keys()]),
            "The covariance matrix must be positive definite. ")

    def pdf(self, *args):
        mu, sigma = self.mu, self.sigma
        k = len(mu)
        args = ImmutableMatrix(args)
        x = args - mu
        return  S(1)/sqrt((2*pi)**(k)*det(sigma))*exp(
            -S(1)/2*x.transpose()*(sigma.inv()*\
                x))[0]

    def marginal_distribution(self, indices, sym):
        sym = ImmutableMatrix([Indexed(sym, i) for i in indices])
        _mu, _sigma = self.mu, self.sigma
        k = len(self.mu)
        for i in range(k):
            if i not in indices:
                _mu = _mu.row_del(i)
                _sigma = _sigma.col_del(i)
                _sigma = _sigma.row_del(i)
        return Lambda(sym, S(1)/sqrt((2*pi)**(len(_mu))*det(_sigma))*exp(
            -S(1)/2*(_mu - sym).transpose()*(_sigma.inv()*\
                (_mu - sym)))[0])

#-------------------------------------------------------------------------------
# Multivariate Laplace distribution ---------------------------------------------------------

class MultivariateLaplaceDistribution(JointDistribution):
    _argnames = ['mu', 'sigma']
    is_Continuous=True

    @property
    def set(self):
        k = len(self.mu)
        return S.Reals**k

    def check(self, mu, sigma):
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the mean vector and covariance matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([Gt(i, 0) != False for i in sigma.eigenvals().keys()]),
            "The covariance matrix must be positive definite. ")

    def pdf(self, *args):
        mu, sigma = self.mu, self.sigma
        mu_T = mu.transpose()
        k = S(len(mu))
        sigma_inv = sigma.inv()
        args = ImmutableMatrix(args)
        args_T = args.transpose()
        x = (mu_T*sigma_inv*mu)[0]
        y = (args_T*sigma_inv*args)[0]
        v = 1 - k/2
        return S(2)/((2*pi)**(S(k)/2)*sqrt(det(sigma)))\
        *(y/(2 + x))**(S(v)/2)*besselk(v, sqrt((2 + x)*(y)))\
        *exp((args_T*sigma_inv*mu)[0])


#-------------------------------------------------------------------------------
# Multivariate StudentT distribution ---------------------------------------------------------

class MultivariateTDistribution(JointDistribution):
    _argnames = ['mu', 'shape_mat', 'dof']
    is_Continuous=True

    @property
    def set(self):
        k = len(self.mu)
        return S.Reals**k

    def check(self, mu, sigma, v):
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the location vector and shape matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([Gt(i, 0) != False for i in sigma.eigenvals().keys()]),
            "The shape matrix must be positive definite. ")

    def pdf(self, *args):
        mu, sigma = self.mu, self.shape_mat
        v = S(self.dof)
        k = S(len(mu))
        sigma_inv = sigma.inv()
        args = ImmutableMatrix(args)
        x = args - mu
        return gamma((k + v)/2)/(gamma(v/2)*(v*pi)**(k/2)*sqrt(det(sigma)))\
        *(1 + 1/v*(x.transpose()*sigma_inv*x)[0])**((-v - k)/2)

def MultivariateT(syms, mu, sigma, v):
    """
    Creates a joint random variable with multivariate T-distribution.

    Parameters
    ==========

    syms: list/tuple/set of symbols for identifying each component
    mu: A list/tuple/set consisting of k means,represents a k
        dimensional location vector
    sigma: The shape matrix for the distribution

    Returns
    =======

    A random symbol
    """
    return multivariate_rv(MultivariateTDistribution, syms, mu, sigma, v)


#-------------------------------------------------------------------------------
# Multivariate Normal Gamma distribution ---------------------------------------------------------

class NormalGammaDistribution(JointDistribution):

    _argnames = ['mu', 'lamda', 'alpha', 'beta']
    is_Continuous=True

    def check(self, mu, lamda, alpha, beta):
        _value_check(mu.is_real, "Location must be real.")
        _value_check(lamda > 0, "Lambda must be positive")
        _value_check(alpha > 0, "alpha must be positive")
        _value_check(beta > 0, "beta must be positive")

    @property
    def set(self):
        return S.Reals*Interval(0, S.Infinity)

    def pdf(self, x, tau):
        beta, alpha, lamda = self.beta, self.alpha, self.lamda
        mu = self.mu

        return beta**alpha*sqrt(lamda)/(gamma(alpha)*sqrt(2*pi))*\
        tau**(alpha - S(1)/2)*exp(-1*beta*tau)*\
        exp(-1*(lamda*tau*(x - mu)**2)/S(2))

    def marginal_distribution(self, indices, *sym):
        if len(indices) == 2:
            return self.pdf(*sym)
        if indices[0] == 0:
            #For marginal over `x`, return non-standardized Student-T's
            #distribution
            x = sym[0]
            v, mu, sigma = self.alpha - S(1)/2, self.mu, \
                S(self.beta)/(self.lamda * self.alpha)
            return Lambda(sym, gamma((v + 1)/2)/(gamma(v/2)*sqrt(pi*v)*sigma)*\
                (1 + 1/v*((x - mu)/sigma)**2)**((-v -1)/2))
        #For marginal over `tau`, return Gamma distribution as per construction
        from sympy.stats.crv_types import GammaDistribution
        return Lambda(sym, GammaDistribution(self.alpha, self.beta)(sym[0]))

def NormalGamma(syms, mu, lamda, alpha, beta):
    """
    Creates a bivariate joint random variable with multivariate Normal gamma
    distribution.

    Parameters
    ==========

    syms: list/tuple/set of two symbols for identifying each component
    mu: A real number, as the mean of the normal distribution
    alpha: a positive integer
    beta: a positive integer
    lamda: a positive integer

    Returns
    =======

    A random symbol
    """
    return multivariate_rv(NormalGammaDistribution, syms, mu, lamda, alpha, beta)

#-------------------------------------------------------------------------------
# Multivariate Beta/Dirichlet distribution ---------------------------------------------------------

class MultivariateBetaDistribution(JointDistribution):

    _argnames = ['alpha']
    is_Continuous = True

    def check(self, alpha):
        _value_check(len(alpha) >= 2, "At least two categories should be passed.")
        for a_k in alpha:
            _value_check((a_k > 0) != False, "Each concentration parameter"
                                            " should be positive.")

    @property
    def set(self):
        k = len(self.alpha)
        return Interval(0, 1)**k

    def pdf(self, *syms):
        alpha = self.alpha
        B = Mul.fromiter(map(gamma, alpha))/gamma(Add(*alpha))
        return Mul.fromiter([sym**(a_k - 1) for a_k, sym in zip(alpha, syms)])/B

def MultivariateBeta(syms, *alpha):
    """
    Creates a continuous random variable with Dirichlet/Multivariate Beta
    Distribution.

    The density of the dirichlet distribution can be found at [1].

    Parameters
    ==========

    alpha: positive real numbers signifying concentration numbers.

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import density
    >>> from sympy.stats.joint_rv import marginal_distribution
    >>> from sympy.stats.joint_rv_types import MultivariateBeta
    >>> from sympy import Symbol
    >>> a1 = Symbol('a1', positive=True)
    >>> a2 = Symbol('a2', positive=True)
    >>> B = MultivariateBeta('B', [a1, a2])
    >>> C = MultivariateBeta('C', a1, a2)
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> density(B)(x, y)
    x**(a1 - 1)*y**(a2 - 1)*gamma(a1 + a2)/(gamma(a1)*gamma(a2))
    >>> marginal_distribution(C, C[0])(x)
    x**(a1 - 1)*gamma(a1 + a2)/(a2*gamma(a1)*gamma(a2))

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Dirichlet_distribution
    .. [2] http://mathworld.wolfram.com/DirichletDistribution.html

    """
    if not isinstance(alpha[0], list):
        alpha = (list(alpha),)
    return multivariate_rv(MultivariateBetaDistribution, syms, alpha[0])

Dirichlet = MultivariateBeta

#-------------------------------------------------------------------------------
# Multivariate Ewens distribution ---------------------------------------------------------

class MultivariateEwensDistribution(JointDistribution):

    _argnames = ['n', 'theta']
    is_Discrete = True
    is_Continuous = False

    def check(self, n, theta):
        _value_check(isinstance(n, Integer) and (n > 0) == True,
                        "sample size should be positive integer.")
        _value_check(theta.is_positive, "mutation rate should be positive.")

    @property
    def set(self):
        prod_set = Range(0, self.n//1 + 1)
        for i in range(2, self.n + 1):
            prod_set *= Range(0, self.n//i + 1)
        return prod_set

    def pdf(self, *syms):
        n, theta = self.n, self.theta
        term_1 = factorial(n)/rf(theta, n)
        term_2 = Mul.fromiter([theta**syms[j]/((j+1)**syms[j]*factorial(syms[j]))
                            for j in range(n)])
        cond = Eq(sum([(k+1)*syms[k] for k in range(n)]), n)
        return Piecewise((term_1 * term_2, cond), (0, True))

def MultivariateEwens(syms, n, theta):
    """
    Creates a discrete random variable with Multivariate Ewens
    Distribution.

    The density of the said distribution can be found at [1].

    Parameters
    ==========

    n: postive integer of class Integer,
            size of the sample or the integer whose partitions are considered
    theta: mutation rate, must be positive real number.

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import density
    >>> from sympy.stats.joint_rv import marginal_distribution
    >>> from sympy.stats.joint_rv_types import MultivariateEwens
    >>> from sympy import Symbol
    >>> a1 = Symbol('a1', positive=True)
    >>> a2 = Symbol('a2', positive=True)
    >>> ed = MultivariateEwens('E', 2, 1)
    >>> density(ed)(a1, a2)
    Piecewise((2**(-a2)/(factorial(a1)*factorial(a2)), Eq(a1 + 2*a2, 2)), (0, True))
    >>> marginal_distribution(ed, ed[0])(a1)
    Piecewise((1/factorial(a1), Eq(a1, 2)), (0, True))

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Ewens%27s_sampling_formula
    .. [2] http://www.stat.rutgers.edu/home/hcrane/Papers/STS529.pdf

    """
    return multivariate_rv(MultivariateEwensDistribution, syms, n, theta)

#-------------------------------------------------------------------------------
# Multinomial distribution ---------------------------------------------------------

class MultinomialDistribution(JointDistribution):

    _argnames = ['n', 'p']
    is_Continuous=False
    is_Discrete = True

    def check(self, n, p):
        _value_check(n > 0,
                        "number of trials must be a positve integer")
        for p_k in p:
            _value_check((p_k >= 0, p_k <= 1),
                        "probability must be at least a positive symbol.")
        _value_check(Eq(sum(p), 1),
                        "probabilities must sum to 1")

    @property
    def set(self):
        return Intersection(S.Naturals0, Interval(0, self.n))**len(self.p)

    def pdf(self, *x):
        n, p = self.n, self.p
        term_1 = factorial(n)/Mul.fromiter([factorial(x_k) for x_k in x])
        term_2 = Mul.fromiter([p_k**x_k for p_k, x_k in zip(p, x)])
        return Piecewise((term_1 * term_2, Eq(sum(x), n)), (0, True))

def Multinomial(syms, n, *p):
    """
    Creates a discrete random variable with Multinomial Distribution.

    The density of the said distribution can be found at [1].

    Parameters
    ==========
    n: postive integer of class Integer,
       number of trials
    p: event probabilites, >= 0 and <= 1

    Returns
    =======
    A RandomSymbol.

    Examples
    ========
    >>> from sympy.stats import density
    >>> from sympy.stats.joint_rv import marginal_distribution
    >>> from sympy.stats.joint_rv_types import Multinomial
    >>> from sympy import symbols
    >>> x1, x2, x3 = symbols('x1, x2, x3', nonnegative=True, integer=True)
    >>> p1, p2, p3 = symbols('p1, p2, p3', positive=True)
    >>> M = Multinomial('M', 3, p1, p2, p3)
    >>> density(M)(x1, x2, x3)
    Piecewise((6*p1**x1*p2**x2*p3**x3/(factorial(x1)*factorial(x2)*factorial(x3)),
    Eq(x1 + x2 + x3, 3)), (0, True))
    >>> marginal_distribution(M, M[0])(x1).subs(x1, 1)
    3*p1*p2**2 + 6*p1*p2*p3 + 3*p1*p3**2

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Multinomial_distribution
    .. [2] http://mathworld.wolfram.com/MultinomialDistribution.html
    """
    if not isinstance(p[0], list):
        p = (list(p), )
    return multivariate_rv(MultinomialDistribution, syms, n, p[0])

#-------------------------------------------------------------------------------
# Negative Multinomial Distribution ---------------------------------------------------------

class NegativeMultinomialDistribution(JointDistribution):

    _argnames = ['k0', 'p']
    is_Continuous=False
    is_Discrete = True

    def check(self, k0, p):
        _value_check(((k0 > 0) != False),
                        "number of failures must be a positve integer")
        for p_k in p:
            _value_check((p_k >= 0) != False,
                        "probability must be at least a positive symbol.")
        _value_check((sum(p) <= 1) != False,
                        "success probabilities must not be greater than 1.")

    @property
    def set(self):
        return Range(0, S.Infinity)**len(self.p)

    def pdf(self, *k):
        k0, p = self.k0, self.p
        term_1 = (gamma(k0 + sum(k))*(1 - sum(p))**k0)/gamma(k0)
        term_2 = Mul.fromiter([pi**ki/factorial(ki) for pi, ki in zip(p, k)])
        return term_1 * term_2

def NegativeMultinomial(syms, k0, *p):
    """
    Creates a discrete random variable with Negative Multinomial Distribution.

    The density of the said distribution can be found at [1].

    Parameters
    ==========
    k0: postive integer of class Integer,
        number of failures before the experiment is stopped
    p: event probabilites, >= 0 and <= 1

    Returns
    =======
    A RandomSymbol.

    Examples
    ========
    >>> from sympy.stats import density
    >>> from sympy.stats.joint_rv import marginal_distribution
    >>> from sympy.stats.joint_rv_types import NegativeMultinomial
    >>> from sympy import symbols
    >>> x1, x2, x3 = symbols('x1, x2, x3', nonnegative=True, integer=True)
    >>> p1, p2, p3 = symbols('p1, p2, p3', positive=True)
    >>> N = NegativeMultinomial('M', 3, p1, p2, p3)
    >>> N_c = NegativeMultinomial('M', 3, 0.1, 0.1, 0.1)
    >>> density(N)(x1, x2, x3)
    p1**x1*p2**x2*p3**x3*(-p1 - p2 - p3 + 1)**3*gamma(x1 + x2 +
    x3 + 3)/(2*factorial(x1)*factorial(x2)*factorial(x3))
    >>> marginal_distribution(N_c, N_c[0])(1).evalf().round(2)
    0.25


    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Negative_multinomial_distribution
    .. [2] http://mathworld.wolfram.com/NegativeBinomialDistribution.html
    """
    if not isinstance(p[0], list):
        p = (list(p), )
    return multivariate_rv(NegativeMultinomialDistribution, syms, k0, p[0])

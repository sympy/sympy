from sympy import sympify, S, pi, sqrt, exp, Lambda
from sympy.stats.rv import _value_check
from sympy.stats.joint_rv import JointDistribution, JointPSpace
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions.determinant import det
from sympy.core.containers import Tuple

__all__ = ['MultivariateNormal',
'MultivariateLaplace',
'MultivariateT',
'NormalGamma'
]

def rv(cls, syms, *args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return JointPSpace(syms, dist)

#-------------------------------------------------------------------------------
# Multivariate Normal distribution ---------------------------------------------------------

class MultivariateNormalDistribution(JointDistribution):
    _argnames = ['mu', 'sigma']

    is_Continuous=True

    def set(mu):
        k = len(mu)
        return S.Reals**k

    def check(self, mu, sigma):
        mu, sigma = Matrix([mu]), Matrix(sigma)
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the mean vector and covariance matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([i > 0 for i in sigma.eigenvals().keys()]),
            "The covariance matrix must be positive definite. ")

    def pdf(self, *args):
        mu, sigma = Matrix(self.mu), Matrix(self.sigma)
        k = len(mu)
        args = Matrix(args)
        return  S(1)/sqrt((2*pi)**(k)*det(sigma))*exp(
            -S(1)/2*(mu - args).transpose()*(sigma**(-1)*\
                (mu - args)))[0]

    def marginal_density(self, indices, *sym):
        _mu, _sigma = Matrix(self.mu), Matrix(self.sigma)
        sym = Matrix(list(sym))
        k = len(self.mu)
        for i in range(k):
            if i not in indices:
                _mu.row_del(i)
                _sigma.col_del(i)
                _sigma.row_del(i)
        return Lambda(sym, S(1)/sqrt((2*pi)**(len(_mu))*det(_sigma))*exp(
            -S(1)/2*(_mu - sym).transpose()*(_sigma**(-1)*\
                (_mu - sym)))[0])


def MultivariateNormal(syms, mu, sigma):
    """
    Creates a joint random variable with multivariate noramal distribution.

    Parameters:
    ==========

    syms: list/tuple/set of symbols for identifying each component
    mu: A list/tuple/set consisting of k means,represents a k
        dimensional location vector
    covariance_mat: The covariance matrix for the distribution

    Returns:
    =======

    A random symbol
    """
    return rv(MultivariateNormalDistribution, syms, mu, sigma)

#-------------------------------------------------------------------------------
# Multivariate laplace distribution ---------------------------------------------------------

class MultivariateLaplaceDistribution(JointDistribution):

    _argnames = ['mu', 'sigma']
    is_Continuous=True

    def set(mu):
        k = len(mu)
        return S.Reals**k

    def check(self, mu, sigma):
        mu, sigma = Matrix([mu]), Matrix(sigma)
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the mean vector and covariance matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([i > 0 for i in sigma.eigenvals().keys()]),
            "The covariance matrix must be positive definite. ")

    def pdf(self, *args):
        from sympy.functions.special.bessel import besselk
        mu, sigma = Matrix(self.mu), Matrix(self.sigma)
        k = len(mu)
        sigma_inv = sigma**(-1)
        args = Matrix(args)
        v = 1 - S(k)/2
        return S(2)/((2*pi)**(S(k)/2)*sqrt(det(sigma))) \
        *((args.transpose()*sigma_inv*args)[0]/(2 + (mu.transpose()*sigma_inv*mu)[0])) \
        **(S(v)/2)*besselk(v, sqrt((2 + (mu.transpose()*sigma_inv*mu)[0])
            *((args.transpose()*sigma_inv*args)[0])))\
        *exp((args.transpose()*sigma_inv*mu)[0])

def MultivariateLaplace(syms, mu, sigma):
    """
    Creates a joint random variable with multivariate laplace distribution.

    Parameters:
    ==========

    syms: list/tuple/set of symbols for identifying each component
    mean_vector: A list/tuple/set consisting of k means,represents a k
        dimensional mean vector
    covariance_mat: The covariance matrix for the distribution

    Returns:
    =======

    A random symbol
    """
    return rv(MultivariateLaplaceDistribution, syms, mu, sigma)

#-------------------------------------------------------------------------------
# Multivariate T-distribution ---------------------------------------------------------

class MultivariateTDistribution(JointDistribution):

    _argnames = ['mu', 'shape_mat', 'dof']
    is_Continuous=True

    def set(mu):
        k = len(mu)
        return S.Reals**k

    def check(self, mu, sigma, v):
        mu, sigma = Matrix([mu]), Matrix(sigma)
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the location vector and shape matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([i > 0 for i in sigma.eigenvals().keys()]),
            "The shape matrix must be positive definite. ")

    def pdf(self, *args):
        from sympy.functions.special.gamma_functions import gamma
        mu, sigma = Matrix(self.mu), Matrix(self.shape_mat)
        v = S(self.dof)
        k = S(len(mu))
        sigma_inv = sigma**(-1)
        args = Matrix(args)
        return gamma((k + v)/2)/(gamma(v/2)*(v*pi)**(k/2)*sqrt(det(sigma)))\
        *(1 + 1/v*((args - mu).transpose()*sigma_inv*(args - mu))[0])**((-v - k)/2)

def MultivariateT(syms, mu, sigma, v):
    """
    Creates a joint random variable with multivariate T-distribution.

    Parameters:
    ==========

    syms: list/tuple/set of symbols for identifying each component
    mu: A list/tuple/set consisting of k means,represents a k
        dimensional location vector
    sigma: The shape matrix for the distribution

    Returns:
    =======

    A random symbol
    """
    return rv(MultivariateTDistribution, syms, mu, sigma, v)

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

    def set(self):
        return S.Reals**2

    def pdf(self, x, tau):
        from sympy.functions.special.gamma_functions import gamma
        beta, alpha, lamda = self.beta, self.alpha, self.lamda
        mu = self.mu

        return beta**alpha*sqrt(lamda)/(gamma(alpha)*sqrt(2*pi))*\
        tau**(alpha - S(1)/2)*exp(-1*beta*tau)*\
        exp(-1*(lamda*tau*(x - mu)**2)/S(2))

def NormalGamma(syms, mu, lamda, alpha, beta):
    """
    Creates a bivariate joint random variable with multivariate Normal gamma
    distribution.

    Parameters:
    ==========

    syms: list/tuple/set of two symbols for identifying each component
    mu: A real number, as the mean of the normal distribution
    alpha: a positive integer
    beta: a positive integer
    lamda: a positive integer

    Returns:
    =======

    A random symbol
    """
    return rv(NormalGammaDistribution, syms, mu, lamda, alpha, beta)

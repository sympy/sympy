from sympy import sympify, S, pi, sqrt, exp, Lambda
from sympy.stats.rv import _value_check
from sympy.stats.joint_rv import JointDistribution, JointPSpace
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions.determinant import det

def rv(cls, syms, *args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return JointPSpace(syms, dist)


def MultivariateNormal(syms, mu, sigma):
    """
    Creates a joint random variable with multivariate noramal idstribution.

    Parameters:
    ==========

    syms: list/tuple/set of symbols for identifying each component
    k: An integer greater than one, number of components
    mean_vector: A list/tuple/set consisting of k means,represents a k
        dimensional mean vector
    covariance_mat: The covariance matrix for the distribution

    Returns:
    =======

    A random symbol
    """
    return rv(MultivariateNormalDistribution, syms, mu, sigma)

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
        mu, sigma = Matrix([self.mu]), Matrix(self.sigma)
        k = len(mu)
        args = Matrix([args])
        return  S(1)/sqrt((2*pi)**(k)*det(sigma))*exp(
            -S(1)/2*(mu- args)*(sigma**(-1)*\
                (mu - args).transpose()))[0]

    def marginal_density(self, indices, *sym):
        _mu, _sigma = Matrix([self.mu]), Matrix(self.sigma)
        sym = Matrix([list(sym)])
        k = len(self.mu)
        for i in range(k):
            if i not in indices:
                _mu.col_del(i)
                _sigma.col_del(i)
                _sigma.row_del(i)
        return Lambda(sym, S(1)/sqrt((2*pi)**(len(_mu))*det(_sigma))*exp(
            -S(1)/2*(_mu - sym)*(_sigma**(-1)*\
                (_mu - sym).transpose()))[0])

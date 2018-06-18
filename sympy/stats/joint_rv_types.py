from sympy import sympify, S, pi, sqrt, exp
from sympy.stats.rv import _value_check
from sympy.stats.joint_rv import JointDistribution, JointPSpace
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions.determinant import det

def rv(name, cls, syms, *args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return JointPSpace(name, syms, dist).value


def MultivariateNormal(name, syms, mu, sigma):
    """
    Creates a joint random variable with multivariate noramal idstribution.

    Parameters:
    ==========

    name: A string or a symbol
    k: An integer greater than one, number of components
    mean_vector: A list/tuple consisting of k means, represents a k
        dimensional mean vector
    covariance_mat: The covariance matrix for the distribution

    Returns:
    =======

    A random symbol
    """
    mu, sigma = Matrix([mu]), Matrix(sigma)
    return rv(name, MultivariateNormalDistribution, syms, mu, sigma)

class MultivariateNormalDistribution(JointDistribution):
    _argnames = ['mu', 'sigma']

    def set(mu):
        k = len(mu)
        return S.Reals**k

    def check(self, mu, sigma):
        _value_check(len(mu) == len(sigma.col(0)),
            "Size of the mean vector and covariance matrix are incorrect.")
        #check if covariance matrix is positive definite or not.
        _value_check(all([i > 0 for i in sigma.eigenvals().keys()]),
            "The covariance matrix must be positive definite. ")

    def pdf(self, *args):
        k = len(self.mu)
        args = Matrix([args])
        return  1/sqrt((2*pi)**(k)*det(self.sigma))*exp(
            -1/2*(self.mu- args)*(self.sigma**(-1)*\
                (self.mu - args).transpose()))[0]

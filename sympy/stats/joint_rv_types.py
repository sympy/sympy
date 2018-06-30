from sympy import sympify, S, pi, sqrt, exp, Lambda, Indexed, Symbol
from sympy.stats.rv import _value_check
from sympy.stats.joint_rv import JointDistribution, JointPSpace
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions.determinant import det

# __all__ = ['MultivariateNormal',
# 'MultivariateLaplace',
# 'MultivariateT',
# 'NormalGamma'
# ]

def multivariate_rv(cls, sym, *args):
    sym = sympify(sym)
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return JointPSpace(sym, dist).value

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

    def marginal_distribution(self, indices, sym):
        sym = Matrix([Symbol(str(Indexed(sym, i))) for i in indices])
        _mu, _sigma = Matrix(self.mu), Matrix(self.sigma)
        k = len(self.mu)
        for i in range(k):
            if i not in indices:
                _mu.row_del(i)
                _sigma.col_del(i)
                _sigma.row_del(i)
        return Lambda(sym, S(1)/sqrt((2*pi)**(len(_mu))*det(_sigma))*exp(
            -S(1)/2*(_mu - sym).transpose()*(_sigma**(-1)*\
                (_mu - sym)))[0])

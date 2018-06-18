from sympy import S, symbols, sqrt, exp, pi
from sympy.stats import Poisson, Geometric, density
from sympy.stats.rv import pspace
from sympy.stats.joint_rv import Joint
from sympy.stats.joint_rv_types import MultivariateNormal
from sympy.stats.crv_types import Normal, Exponential
from sympy.utilities.pytest import raises
x, y, n = symbols('x y n')

def test_multivariate_normal():
    m = MultivariateNormal('m', ('m1', 'm2'), [1, 2], [[1, 0], [0, 1]])
    assert density(m)(1, 2) == 1/(2*pi)
    raises (ValueError,\
        lambda: MultivariateNormal('m', ('m1', 'm2'),[1, 2], [[0, 0], [0, 1]]))
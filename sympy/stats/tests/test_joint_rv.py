from sympy import symbols, pi
from sympy.stats import density
from sympy.stats.joint_rv_types import MultivariateNormal
from sympy.stats.joint_rv import marginal_density
from sympy.utilities.pytest import raises
x, y, n = symbols('x y z')

def test_multivariate_normal():
    m = MultivariateNormal(('m1', 'm2'), [1, 2], [[1, 0], [0, 1]])
    assert density(m)(1, 2) == 1/(2*pi)
    raises (ValueError,\
        lambda: MultivariateNormal(('m1', 'm2'),[1, 2], [[0, 0], [0, 1]]))
    n = MultivariateNormal(('x', 'y', 'z'), [1, 2, 3], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert marginal_density(n, x, y)(1, 2) == 1/(2*pi)

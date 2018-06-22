from sympy import symbols, pi, oo, sqrt, exp
from sympy.stats import density
from sympy.stats.joint_rv import marginal_density
from sympy.utilities.pytest import raises
from sympy.integrals.integrals import integrate
x, y, z, a, b = symbols('x y z a b')

def test_MultivariateNormal():
    from sympy.stats.joint_rv_types import MultivariateNormal
    m = MultivariateNormal((x, y), [1, 2], [[1, 0], [0, 1]])
    assert density(m)(1, 2) == 1/(2*pi)
    raises (ValueError,\
        lambda: MultivariateNormal(('m1', 'm2'),[1, 2], [[0, 0], [0, 1]]))
    n = MultivariateNormal(('x', 'y', 'z'), [1, 2, 3], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert marginal_density(n, x, y)(1, 2) == 1/(2*pi)
    assert integrate(density(m)(x, y), (x, -oo, oo), (y, -oo, oo)).evalf() == 1
    raises (ValueError, lambda: MultivariateNormal(
        ('m1', 'm2'), [1, 2], [[1, 1], [1, -1]]))

def test_MultivariateTDist():
    from sympy.stats.joint_rv_types import MultivariateT
    t1 = MultivariateT((x, y), [0, 0], [[1, 0], [0, 1]], 2)
    assert(density(t1))(1, 1) == 1/(8*pi)
    assert integrate(density(t1)(x, y), (x, -oo, oo), \
        (y, -oo, oo)).evalf() == 1
    raises(ValueError, lambda: MultivariateT(
        ('t1', 't2'), [1, 2], [[1, 1], [1, -1]], 1))

def test_NormalGamma():
    from sympy.stats.joint_rv_types import NormalGamma
    ng = NormalGamma(('x', 'y'), 1, 2, 3, 4)
    assert density(ng)(1, 1) == 32*exp(-4)/sqrt(pi)
    raises(ValueError, lambda:NormalGamma(('x', 'y'), 1, 2, 3, -1))

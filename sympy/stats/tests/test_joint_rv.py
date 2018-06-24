from sympy import S, symbols, pi, oo, sqrt, exp
from sympy.stats import density
from sympy.stats.joint_rv import marginal_distribution
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
    assert marginal_distribution(n, x, y)(1, 2) == 1/(2*pi)
    assert integrate(density(m)(x, y), (x, -oo, oo), (y, -oo, oo)).evalf() == 1
    raises (ValueError, lambda: MultivariateNormal(
        ('m1', 'm2'), [1, 2], [[1, 1], [1, -1]]))

def test_MultivariateTDist():
    from sympy.stats.joint_rv_types import MultivariateT
    t1 = MultivariateT(('x', 'y'), [0, 0], [[1, 0], [0, 1]], 2)
    assert(density(t1))(1, 1) == 1/(8*pi)
    assert integrate(density(t1)(x, y), (x, -oo, oo), \
        (y, -oo, oo)).evalf() == 1
    raises(ValueError, lambda: MultivariateT(
        ('t1', 't2'), [1, 2], [[1, 1], [1, -1]], 1))

def test_NormalGamma():
    from sympy.stats.joint_rv_types import NormalGamma
    from sympy import gamma
    ng = NormalGamma(('x', 'y'), 1, 2, 3, 4)
    assert density(ng)(1, 1) == 32*exp(-4)/sqrt(pi)
    raises(ValueError, lambda:NormalGamma(('x', 'y'), 1, 2, 3, -1))
    assert marginal_distribution(ng, x)(1) == \
        3*sqrt(10)*gamma(S(7)/4)/(10*sqrt(pi)*gamma(S(5)/4))
    assert marginal_distribution(ng, y)(1) == exp(-S(1)/4)/128

def test_JointPSpace_margial_distribution():
    from sympy.stats.joint_rv_types import MultivariateT
    from sympy import polar_lift
    T = MultivariateT(('x', 'y'), [0, 0], [[1, 0], [0, 1]], 2)
    assert marginal_distribution(T, x)(x) == sqrt(2)*(x**2 + 2)/(
        8*polar_lift(x**2/2 + 1)**(5/2))
    assert integrate(marginal_distribution(T, x)(x), (x, -oo, oo)) == 1
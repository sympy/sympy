from sympy import S, symbols, pi, oo, sqrt, exp
from sympy.stats import density
from sympy.stats.joint_rv import marginal_distribution
from sympy.utilities.pytest import raises
from sympy.integrals.integrals import integrate
x, y, z, a, b = symbols('x y z a b')

def test_MultivariateNormal():
    from sympy.stats.joint_rv_types import MultivariateNormal
    m = MultivariateNormal('A', [1, 2], [[1, 0], [0, 1]])
    assert density(m)(1, 2) == 1/(2*pi)
    raises (ValueError,\
        lambda: MultivariateNormal('M',[1, 2], [[0, 0], [0, 1]]))
    n = MultivariateNormal('B', [1, 2, 3], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert marginal_distribution(n, 0, 1)(1, 2) == 1/(2*pi)
    assert integrate(density(m)(x, y), (x, -oo, oo), (y, -oo, oo)).evalf() == 1
    raises (ValueError, lambda: MultivariateNormal('M', [1, 2], [[1, 1], [1, -1]]))

def test_MultivariateTDist():
    from sympy.stats.joint_rv_types import MultivariateT
    t1 = MultivariateT('T', [0, 0], [[1, 0], [0, 1]], 2)
    assert(density(t1))(1, 1) == 1/(8*pi)
    assert integrate(density(t1)(x, y), (x, -oo, oo), \
        (y, -oo, oo)).evalf() == 1
    raises(ValueError, lambda: MultivariateT('T', [1, 2], [[1, 1], [1, -1]], 1))

def test_NormalGamma():
    from sympy.stats.joint_rv_types import NormalGamma
    from sympy import gamma
    ng = NormalGamma('G', 1, 2, 3, 4)
    assert density(ng)(1, 1) == 32*exp(-4)/sqrt(pi)
    raises(ValueError, lambda:NormalGamma('G', 1, 2, 3, -1))
    assert marginal_distribution(ng, 0)(1) == \
        3*sqrt(10)*gamma(S(7)/4)/(10*sqrt(pi)*gamma(S(5)/4))
    assert marginal_distribution(ng, y)(1) == exp(-S(1)/4)/128

def test_JointPSpace_margial_distribution():
    from sympy.stats.joint_rv_types import MultivariateT
    from sympy import polar_lift
    T = MultivariateT('T', [0, 0], [[1, 0], [0, 1]], 2)
    assert marginal_distribution(T, T[1])(x) == sqrt(2)*(x**2 + 2)/(
        8*polar_lift(x**2/2 + 1)**(5/2))
    assert integrate(marginal_distribution(T, 1)(x), (x, -oo, oo)) == 1

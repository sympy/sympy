from sympy import (symbols, pi, oo, S, exp, sqrt, besselk, Indexed, Sum, simplify,
                    Mul, Rational, Integral)
from sympy.stats import density
from sympy.stats.joint_rv import marginal_distribution
from sympy.stats.joint_rv_types import JointRV
from sympy.stats.crv_types import Normal
from sympy.utilities.pytest import raises, XFAIL
from sympy.integrals.integrals import integrate
from sympy.matrices import Matrix
x, y, z, a, b = symbols('x y z a b')

def test_Normal():
    m = Normal('A', [1, 2], [[1, 0], [0, 1]])
    assert density(m)(1, 2) == 1/(2*pi)
    raises (ValueError, lambda:m[2])
    raises (ValueError,\
        lambda: Normal('M',[1, 2], [[0, 0], [0, 1]]))
    n = Normal('B', [1, 2, 3], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    p = Normal('C',  Matrix([1, 2]), Matrix([[1, 0], [0, 1]]))
    assert density(m)(x, y) == density(p)(x, y)
    assert marginal_distribution(n, 0, 1)(1, 2) == 1/(2*pi)
    assert integrate(density(m)(x, y), (x, -oo, oo), (y, -oo, oo)).evalf() == 1
    N = Normal('N', [1, 2], [[x, 0], [0, y]])
    assert density(N)(0, 0) == exp(-2/y - 1/(2*x))/(2*pi*sqrt(x*y))

    raises (ValueError, lambda: Normal('M', [1, 2], [[1, 1], [1, -1]]))

def test_MultivariateTDist():
    from sympy.stats.joint_rv_types import MultivariateT
    t1 = MultivariateT('T', [0, 0], [[1, 0], [0, 1]], 2)
    assert(density(t1))(1, 1) == 1/(8*pi)
    assert integrate(density(t1)(x, y), (x, -oo, oo), \
        (y, -oo, oo)).evalf() == 1
    raises(ValueError, lambda: MultivariateT('T', [1, 2], [[1, 1], [1, -1]], 1))
    t2 = MultivariateT('t2', [1, 2], [[x, 0], [0, y]], 1)
    assert density(t2)(1, 2) == 1/(2*pi*sqrt(x*y))

def test_multivariate_laplace():
    from sympy.stats.crv_types import Laplace
    raises(ValueError, lambda: Laplace('T', [1, 2], [[1, 2], [2, 1]]))
    L = Laplace('L', [1, 0], [[1, 2], [0, 1]])
    assert density(L)(2, 3) == exp(2)*besselk(0, sqrt(3))/pi
    L1 = Laplace('L1', [1, 2], [[x, 0], [0, y]])
    assert density(L1)(0, 1) == \
        exp(2/y)*besselk(0, sqrt((2 + 4/y + 1/x)/y))/(pi*sqrt(x*y))

def test_NormalGamma():
    from sympy.stats.joint_rv_types import NormalGamma
    from sympy import gamma
    ng = NormalGamma('G', 1, 2, 3, 4)
    assert density(ng)(1, 1) == 32*exp(-4)/sqrt(pi)
    raises(ValueError, lambda:NormalGamma('G', 1, 2, 3, -1))
    assert marginal_distribution(ng, 0)(1) == \
        3*sqrt(10)*gamma(S(7)/4)/(10*sqrt(pi)*gamma(S(5)/4))
    assert marginal_distribution(ng, y)(1) == exp(-S(1)/4)/128

def test_GeneralizedMultivariateLogGammaDistribution():
    from sympy.stats.joint_rv_types import GMVLG
    from sympy import gamma
    omega = Matrix([[1, 0.5, 0.5, 0.5],
                     [0.5, 1, 0.5, 0.5],
                     [0.5, 0.5, 1, 0.5],
                     [0.5, 0.5, 0.5, 1]])
    v, l, mu = (4, [1, 2, 3, 4], [1, 2, 3, 4])
    y_1, y_2, y_3, y_4 = symbols('y_1:5', real=True)
    n = symbols('n', negative=False, integer=True)
    G = GMVLG('G', omega, v, l, mu)
    den = "5*2**(2/3)*5**(1/3)*Sum(4*24**(-n - 4)*(-2**(2/3)*5**(1/3)/4 + 1)**n*"+\
         "exp((n + 4)*(y_1 + 2*y_2 + 3*y_3 + 4*y_4) - exp(y_1) - exp(2*y_2)/2 - "+\
         "exp(3*y_3)/3 - exp(4*y_4)/4)/(gamma(n + 1)*gamma(n + 4)**3), (n, 0, oo))/64"
    assert str(density(G)(y_1, y_2, y_3, y_4)) == den
    marg = "5*2**(2/3)*5**(1/3)*exp(4*y_1)*exp(-exp(y_1))*Integral(exp(-exp(4*G[3])"+\
           "/4)*exp(16*G[3])*Integral(exp(-exp(3*G[2])/3)*exp(12*G[2])*Integral(exp("+\
           "-exp(2*G[1])/2)*exp(8*G[1])*Sum((-1/4)**n*24**(-n)*(-4 + 2**(2/3)*5**(1/3"+\
           "))**n*exp(n*y_1)*exp(2*n*G[1])*exp(3*n*G[2])*exp(4*n*G[3])/(gamma(n + 1)"+\
           "*gamma(n + 4)**3), (n, 0, oo)), (G[1], -oo, oo)), (G[2], -oo, oo)), (G[3]"+\
           ", -oo, oo))/5308416"
    assert str(marginal_distribution(G, G[0])(y_1)) == marg
    omega_f1 = Matrix([[1, 0.5, 0.5]])
    omega_f2 = Matrix([[1, 0.5, 0.5, 0.5],
                     [0.5, 1, 2, 0.5],
                     [0.5, 0.5, 1, 0.5],
                     [0.5, 0.5, 0.5, 1]])
    omega_f3 = Matrix([[6, 0.5, 0.5, 0.5],
                     [0.5, 1, 2, 0.5],
                     [0.5, 0.5, 1, 0.5],
                     [0.5, 0.5, 0.5, 1]])
    v_f = symbols("v_f", positive=False)
    l_f = [1, 2, v_f, 4]
    m_f = [v_f, 2, 3, 4]
    omega_f4 = Matrix([[1, 0.5, 0.5, 0.5, 0.5],
                     [0.5, 1, 0.5, 0.5, 0.5],
                     [0.5, 0.5, 1, 0.5, 0.5],
                     [0.5, 0.5, 0.5, 1, 0.5],
                     [0.5, 0.5, 0.5, 0.5, 1]])
    l_f1 = [1, 2, 3, 4, 5]
    raises(ValueError, lambda: GMVLG('G', omega_f1, v, l, mu))
    raises(ValueError, lambda: GMVLG('G', omega_f2, v, l, mu))
    raises(ValueError, lambda: GMVLG('G', omega_f3, v, l, mu))
    raises(ValueError, lambda: GMVLG('G', omega, v_f, l, mu))
    raises(ValueError, lambda: GMVLG('G', omega, v, l_f, mu))
    raises(ValueError, lambda: GMVLG('G', omega, v, l, m_f))
    raises(ValueError, lambda: GMVLG('G', omega_f4, v, l, mu))
    raises(ValueError, lambda: GMVLG('G', omega, v, l_f1, mu))

def test_JointPSpace_margial_distribution():
    from sympy.stats.joint_rv_types import MultivariateT
    from sympy import polar_lift
    T = MultivariateT('T', [0, 0], [[1, 0], [0, 1]], 2)
    assert marginal_distribution(T, T[1])(x) == sqrt(2)*(x**2 + 2)/(
        8*polar_lift(x**2/2 + 1)**(S(5)/2))
    assert integrate(marginal_distribution(T, 1)(x), (x, -oo, oo)) == 1
    t = MultivariateT('T', [0, 0, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 3)
    assert marginal_distribution(t, 0)(1).evalf().round(1) == 0.2


def test_JointRV():
    from sympy.stats.joint_rv import JointDistributionHandmade
    x1, x2 = (Indexed('x', i) for i in (1, 2))
    pdf = exp(-x1**2/2 + x1 - x2**2/2 - S(1)/2)/(2*pi)
    X = JointRV('x', pdf)
    assert density(X)(1, 2) == exp(-2)/(2*pi)
    assert isinstance(X.pspace.distribution, JointDistributionHandmade)
    assert marginal_distribution(X, 0)(2) == sqrt(2)*exp(-S(1)/2)/(2*sqrt(pi))

def test_expectation():
    from sympy import simplify
    from sympy.stats import E
    m = Normal('A', [x, y], [[1, 0], [0, 1]])
    assert simplify(E(m[1])) == y

@XFAIL
def test_joint_vector_expectation():
    from sympy.stats import E
    m = Normal('A', [x, y], [[1, 0], [0, 1]])
    assert E(m) == (x, y)

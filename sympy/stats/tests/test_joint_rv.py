from sympy import (symbols, pi, oo, S, exp, sqrt, besselk, Indexed, Rational,
                    simplify, Piecewise, factorial, Eq, gamma, Sum)
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

def test_MultivariateBeta():
    from sympy.stats.joint_rv_types import MultivariateBeta
    from sympy import gamma
    a1, a2 = symbols('a1, a2', positive=True)
    a1_f, a2_f = symbols('a1, a2', positive=False)
    mb = MultivariateBeta('B', [a1, a2])
    mb_c = MultivariateBeta('C', a1, a2)
    assert density(mb)(1, 2) == S(2)**(a2 - 1)*gamma(a1 + a2)/\
                                (gamma(a1)*gamma(a2))
    assert marginal_distribution(mb_c, 0)(3) == S(3)**(a1 - 1)*gamma(a1 + a2)/\
                                                (a2*gamma(a1)*gamma(a2))
    raises(ValueError, lambda: MultivariateBeta('b1', [a1_f, a2]))
    raises(ValueError, lambda: MultivariateBeta('b2', [a1, a2_f]))
    raises(ValueError, lambda: MultivariateBeta('b3', [0, 0]))
    raises(ValueError, lambda: MultivariateBeta('b4', [a1_f, a2_f]))

def test_MultivariateEwens():
    from sympy.stats.joint_rv_types import MultivariateEwens
    n, theta = symbols('n theta', positive=True)
    theta_f = symbols('t_f', negative=True)
    a = symbols('a_1:4', positive = True, integer = True)
    ed = MultivariateEwens('E', 3, theta)
    assert density(ed)(a[0], a[1], a[2]) == Piecewise((6*2**(-a[1])*3**(-a[2])*
                                            theta**a[0]*theta**a[1]*theta**a[2]/
                                            (theta*(theta + 1)*(theta + 2)*
                                            factorial(a[0])*factorial(a[1])*
                                            factorial(a[2])), Eq(a[0] + 2*a[1] +
                                            3*a[2], 3)), (0, True))
    assert marginal_distribution(ed, ed[1])(a[1]) == Piecewise((6*2**(-a[1])*
                                                    theta**a[1]/((theta + 1)*
                                                    (theta + 2)*factorial(a[1])),
                                                    Eq(2*a[1] + 1, 3)), (0, True))
    raises(ValueError, lambda: MultivariateEwens('e1', 5, theta_f))
    raises(ValueError, lambda: MultivariateEwens('e1', n, theta))

def test_Multinomial():
    from sympy.stats.joint_rv_types import Multinomial
    n, x1, x2, x3, x4 = symbols('n, x1, x2, x3, x4', nonnegative=True, integer=True)
    p1, p2, p3, p4 = symbols('p1, p2, p3, p4', positive=True)
    p1_f, n_f = symbols('p1_f, n_f', negative=True)
    M = Multinomial('M', n, [p1, p2, p3, p4])
    C = Multinomial('C', 3, p1, p2, p3)
    f = factorial
    assert density(M)(x1, x2, x3, x4) == Piecewise((p1**x1*p2**x2*p3**x3*p4**x4*
                                            f(n)/(f(x1)*f(x2)*f(x3)*f(x4)),
                                            Eq(n, x1 + x2 + x3 + x4)), (0, True))
    assert marginal_distribution(C, C[0])(x1).subs(x1, 1) ==\
                                                            3*p1*p2**2 +\
                                                            6*p1*p2*p3 +\
                                                            3*p1*p3**2
    raises(ValueError, lambda: Multinomial('b1', 5, [p1, p2, p3, p1_f]))
    raises(ValueError, lambda: Multinomial('b2', n_f, [p1, p2, p3, p4]))
    raises(ValueError, lambda: Multinomial('b3', n, 0.5, 0.4, 0.3, 0.1))

def test_NegativeMultinomial():
    from sympy.stats.joint_rv_types import NegativeMultinomial
    k0, x1, x2, x3, x4 = symbols('k0, x1, x2, x3, x4', nonnegative=True, integer=True)
    p1, p2, p3, p4 = symbols('p1, p2, p3, p4', positive=True)
    p1_f = symbols('p1_f', negative=True)
    N = NegativeMultinomial('N', 4, [p1, p2, p3, p4])
    C = NegativeMultinomial('C', 4, 0.1, 0.2, 0.3)
    g = gamma
    f = factorial
    assert simplify(density(N)(x1, x2, x3, x4) -
            p1**x1*p2**x2*p3**x3*p4**x4*(-p1 - p2 - p3 - p4 + 1)**4*g(x1 + x2 +
            x3 + x4 + 4)/(6*f(x1)*f(x2)*f(x3)*f(x4))) == S(0)
    assert marginal_distribution(C, C[0])(1).evalf().round(2) == 0.33
    raises(ValueError, lambda: NegativeMultinomial('b1', 5, [p1, p2, p3, p1_f]))
    raises(ValueError, lambda: NegativeMultinomial('b2', k0, 0.5, 0.4, 0.3, 0.4))

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

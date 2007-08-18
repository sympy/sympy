from sympy import *
from sympy.specfun.factorials import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

fs = factorial_simplify
fac = factorial

def test_factorial1():
    assert [fac(t) for t in [0,1,2,3,4]] == [1,1,2,6,24]
    assert (fac(-1) == oo) == True
    assert fac(Rational(1,2)) == Rational(1,2)*sqrt(pi)
    assert fac(Rational(3,2)) == Rational(3,4)*sqrt(pi)
    assert fac(Rational(5,2)) == Rational(15,8)*sqrt(pi)
    assert fac(Rational(-1,2)) == sqrt(pi)
    assert fac(Rational(-3,2)) == -2*sqrt(pi)
    assert fac(Rational(-5,2)) == Rational(4,3)*sqrt(pi)
    assert fac(Rational(-17,2)) == Rational(256,2027025)*sqrt(pi)
    #assert latex(fac(x, evaluate=False)) == "$x!$"
    #assert latex(fac(-4, evaluate=False)) == "$(-4)!$"
    #assert latex(fac(-x, evaluate=False)) == "$(- x)!$"

def test_factorial_evalf():
    def relcmp(a,b):
        return abs(a-b)/abs(a)
    for i in range(20):
        assert relcmp(unfac(i).evalf(), fac(i)) < Real(1e-10)
    for i in range(-7, 7, 2):
        x = Rational(i)/2
        assert relcmp(unfac(x).evalf(), fac(x).evalf()) < Real(1e-10)
    assert relcmp(unfac(165).evalf(), fac(165)) < Real(1e-10)
    assert relcmp(fac(-0.75).evalf(), 3.6256099082219083) < Real(1e-10)
    assert relcmp((1-fac(1e-8).evalf())*10**8, 0.5772157) < Real(1e-5)
    assert relcmp((fac(-1e-8).evalf()-1)*10**8, 0.5772157) < Real(1e-5)
    #re, im = fac(7+8*I).evalf().as_real_imag()
    #assert relcmp(re, 1.84428481562553) < Real(1e-10)
    #assert relcmp(im, -125.96060801751909) < Real(1e-10)

def test_multifactorials():
    # OEIS: A006882 
    assert [factorial(n, 2) for n in range(10)] == [1, 1, 2, 3, 8, 15, 48, 105, 384, 945]

    assert (factorial(-2, 2) == oo) == True
    assert (factorial(-4, 2) == oo) == True
    assert factorial(-1, 2) == 1
    assert factorial(-3, 2) == -1
    assert factorial(-7, 2) == Rational(-1,15)
    assert factorial(-9, 2) == Rational(1,105)
    #assert latex(factorial2(x, evaluate=False)) == "$x!!$"
    #assert latex(factorial2(-4, evaluate=False)) == "$(-4)!!$"
    #assert latex(factorial2(-x, evaluate=False)) == "$(- x)!!$"

    # OEIS: A007661
    assert [factorial(n, 3) for n in range(10)] == [1, 1, 2, 3, 4, 10, 18, 28, 80, 162]
    # OEIS: A007662
    assert [factorial(n, 4) for n in range(10)] == [1, 1, 2, 3, 4, 5, 12, 21, 32, 45]

def test_factorial_simplify():
    assert fs(fac(x+5)/fac(x+5)) == 1
    assert fs(fac(x+1)/fac(x)) == 1+x
    assert fs(fac(x+2)/fac(x)) == (1+x)*(2+x)
    assert fs(fac(x+3)/fac(x)) == (1+x)*(2+x)*(3+x)
    assert fs(fac(x-1)/fac(x)) == (1/x)
    assert fs(fac(x-2)/fac(x)) == 1/(x*(-1+x))
    assert fs(fac(x-3)/fac(x)) == 1/(x*(-1+x)*(-2+x))
    assert fs(fac(x)*(x+1)*(x+2)) == fac(x+2)
    assert fs(fac(x)/x/(x-1)) == fac(x-2)
    assert fs(x*(x-1)/fac(x)) == 1/fac(x-2)
    assert fs((1/(x+1))/fac(x)) == 1/fac(x+1)
    assert fs(fac(x)*fac(y-2)*fac(z+2)/fac(z)/fac(y+1)) == fac(x)*(z+1)*(z+2)/(y-1)/y/(y+1)
    assert fs(fac(x)*fac(y+1)*fac(z+2)/fac(z)/fac(y-2)) == fac(x)*(z+1)*(z+2)*(y-1)*y*(y+1)

def test_rising_falling():
    assert rising_factorial(x, 0) == 1
    assert rising_factorial(x, 1) == x
    assert rising_factorial(x, 2) == x*(x+1)
    assert rising_factorial(x, 3) == x*(x+1)*(x+2)
    assert falling_factorial(x, 0) == 1
    assert falling_factorial(x, 1) == x
    assert falling_factorial(x, 2) == x*(x-1)
    assert falling_factorial(x, 3) == x*(x-1)*(x-2)
    assert rising_factorial(1, x) == fac(x)
    assert falling_factorial(1, x) == 1/fac(1-x)
    assert falling_factorial(15, 8) == 259459200
    assert falling_factorial(-3,4) == 360
    assert falling_factorial(3,-4) == Rational(1,840)
    assert rising_factorial(15, 8) == 12893126400
    assert rising_factorial(-3,4) == 0
    #n = Symbol('n')
    #assert latex(rising_factorial(x, n, evaluate=False)) == "${(x)}^{(n)}$"
    #assert latex(falling_factorial(x, n, evaluate=False)) == "${(x)}_{(n)}$"

def test_binomial2():
    assert binomial2(x, 0) == 1
    assert binomial2(x, x) == 1
    assert binomial2(0, 0) == 1
    assert binomial2(0, 1) == 0
    assert binomial2(0, x) == sin(pi*x)/(pi*x)
    assert binomial2(x, 1) == x
    assert binomial2(x, 2) == Rational(1,2)*x*(x-1)
    assert binomial2(x, 3) == Rational(1,6)*x*(x-1)*(x-2)
    assert [binomial2(4,k) for k in range(5)] == [1,4,6,4,1]
    assert [binomial2(5,k) for k in range(6)] == [1,5,10,10,5,1]
    assert sum(binomial2(20, k) for k in range(21)) == 2**20
    assert binomial2(10**20, 10**20 - 2) == \
        4999999999999999999950000000000000000000
    #assert latex(binomial(8,3,evaluate=False)) == r"${{8}\choose{3}}$"

#def test_gamma():
#    assert gamma(0) == oo
#    assert gamma(1) == 1
#    assert gamma(2) == 1
#    assert gamma(3) == 2
#    assert gamma(Rational(1,2)) == sqrt(pi)
#    #assert latex(gamma(3+x)) == "$\Gamma(3+x)$"
#    assert upper_gamma(4,0) == 6
#    assert upper_gamma(3,oo) == 0
#    assert lower_gamma(1,x) + upper_gamma(1,x) == gamma(1)
#    assert lower_gamma(5,x) + upper_gamma(5,x) == gamma(5)

#def test_derivatives():
#    x = Symbol('x')
#    from sympy.specfun.zeta_functions import polygamma
#    assert gamma(4*x).diff(x) == 4*gamma(4*x)*polygamma(0, 4*x)
#    assert fac(4*x).diff(x) == 4*gamma(1+4*x)*polygamma(0, 1+4*x)
#    assert log(gamma(x)).diff(x, 5) == polygamma(4, x)
#    assert fac(x).diff(x).subs(x, 0) == -EulerGamma
#    assert fac(x).diff(x).subs(x, 4) == 24*(Rational(25, 12)-EulerGamma)

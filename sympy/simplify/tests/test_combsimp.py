from sympy import (
    Rational, combsimp, factorial, gamma, binomial, Symbol, pi, S,
    sin, exp, powsimp, sqrt, sympify, FallingFactorial, RisingFactorial,
    simplify, symbols, cos)

from sympy.abc import x, y, z, t, a, b, c, d, e, f, g, h, i, k


def test_combsimp():
    from sympy.abc import n, k

    assert combsimp(factorial(n)) == factorial(n)
    assert combsimp(binomial(n, k)) == binomial(n, k)

    assert combsimp(factorial(n)/factorial(n - 3)) == n*(-1 + n)*(-2 + n)
    assert combsimp(binomial(n + 1, k + 1)/binomial(n, k)) == (1 + n)/(1 + k)

    assert combsimp(binomial(3*n + 4, n + 1)/binomial(3*n + 1, n)) == \
        S(3)/2*((3*n + 2)*(3*n + 4)/((n + 1)*(2*n + 3)))

    assert combsimp(factorial(n)**2/factorial(n - 3)) == \
        factorial(n)*n*(-1 + n)*(-2 + n)
    assert combsimp(factorial(n)*binomial(n + 1, k + 1)/binomial(n, k)) == \
        factorial(n)*(1 + n)/(1 + k)

    assert combsimp(binomial(n - 1, k)) == -((-n + k)*binomial(n, k))/n

    assert combsimp(binomial(n + 2, k + S(1)/2)) == 4*((n + 1)*(n + 2) *
        binomial(n, k + S(1)/2))/((2*k - 2*n - 1)*(2*k - 2*n - 3))
    assert combsimp(binomial(n + 2, k + 2.0)) == \
        -((1.0*n + 2.0)*binomial(n + 1.0, k + 2.0))/(k - n)

    # coverage tests
    assert combsimp(factorial(n*(1 + n) - n**2 - n)) == 1
    assert combsimp(binomial(n + k - 2, n)) == \
        k*(k - 1)*binomial(n + k, n)/((n + k)*(n + k - 1))
    i = Symbol('i', integer=True)
    e = gamma(i + 3)
    assert combsimp(e) == e
    e = gamma(exp(i))
    assert combsimp(e) == e
    e = gamma(n + S(1)/3)*gamma(n + S(2)/3)
    assert combsimp(e) == e
    assert combsimp(gamma(4*n + S(1)/2)/gamma(2*n - S(3)/4)) == \
        2**(4*n - S(5)/2)*(8*n - 3)*gamma(2*n + S(3)/4)/sqrt(pi)

    assert combsimp(6*FallingFactorial(-4, n)/factorial(n)) == \
        (-1)**n*(n + 1)*(n + 2)*(n + 3)
    assert combsimp(6*FallingFactorial(-4, n - 1)/factorial(n - 1)) == \
        (-1)**(n - 1)*n*(n + 1)*(n + 2)
    assert combsimp(6*FallingFactorial(-4, n - 3)/factorial(n - 3)) == \
        (-1)**(n - 3)*n*(n - 1)*(n - 2)
    assert combsimp(6*FallingFactorial(-4, -n - 1)/factorial(-n - 1)) == \
        -(-1)**(-n - 1)*n*(n - 1)*(n - 2)

    assert combsimp(6*RisingFactorial(4, n)/factorial(n)) == \
        (n + 1)*(n + 2)*(n + 3)
    assert combsimp(6*RisingFactorial(4, n - 1)/factorial(n - 1)) == \
        n*(n + 1)*(n + 2)
    assert combsimp(6*RisingFactorial(4, n - 3)/factorial(n - 3)) == \
        n*(n - 1)*(n - 2)
    assert combsimp(6*RisingFactorial(4, -n - 1)/factorial(-n - 1)) == \
        -n*(n - 1)*(n - 2)


def test_combsimp_gamma():
    from sympy.abc import x, y
    R = Rational

    assert combsimp(gamma(x)) == gamma(x)
    assert combsimp(gamma(x + 1)/x) == gamma(x)
    assert combsimp(gamma(x)/(x - 1)) == gamma(x - 1)
    assert combsimp(x*gamma(x)) == gamma(x + 1)
    assert combsimp((x + 1)*gamma(x + 1)) == gamma(x + 2)
    assert combsimp(gamma(x + y)*(x + y)) == gamma(x + y + 1)
    assert combsimp(x/gamma(x + 1)) == 1/gamma(x)
    assert combsimp((x + 1)**2/gamma(x + 2)) == (x + 1)/gamma(x + 1)
    assert combsimp(x*gamma(x) + gamma(x + 3)/(x + 2)) == \
        (x + 2)*gamma(x + 1)

    assert combsimp(gamma(2*x)*x) == gamma(2*x + 1)/2
    assert combsimp(gamma(2*x)/(x - S(1)/2)) == 2*gamma(2*x - 1)

    assert combsimp(gamma(x)*gamma(1 - x)) == pi/sin(pi*x)
    assert combsimp(gamma(x)*gamma(-x)) == -pi/(x*sin(pi*x))
    assert combsimp(1/gamma(x + 3)/gamma(1 - x)) == \
        sin(pi*x)/(pi*x*(x + 1)*(x + 2))

    assert powsimp(combsimp(
        gamma(x)*gamma(x + S(1)/2)*gamma(y)/gamma(x + y))) == \
        2**(-2*x + 1)*sqrt(pi)*gamma(2*x)*gamma(y)/gamma(x + y)
    assert combsimp(1/gamma(x)/gamma(x - S(1)/3)/gamma(x + S(1)/3)) == \
        3**(3*x - S(3)/2)/(2*pi*gamma(3*x - 1))
    assert simplify(
        gamma(S(1)/2 + x/2)*gamma(1 + x/2)/gamma(1 + x)/sqrt(pi)*2**x) == 1
    assert combsimp(gamma(S(-1)/4)*gamma(S(-3)/4)) == 16*sqrt(2)*pi/3

    assert powsimp(combsimp(gamma(2*x)/gamma(x))) == \
        2**(2*x - 1)*gamma(x + S(1)/2)/sqrt(pi)

    # issue 6792
    e = (-gamma(k)*gamma(k + 2) + gamma(k + 1)**2)/gamma(k)**2
    assert combsimp(e) == -k
    assert combsimp(1/e) == -1/k
    e = (gamma(x) + gamma(x + 1))/gamma(x)
    assert combsimp(e) == x + 1
    assert combsimp(1/e) == 1/(x + 1)
    e = (gamma(x) + gamma(x + 2))*(gamma(x - 1) + gamma(x))/gamma(x)
    assert combsimp(e) == (x**2 + x + 1)*gamma(x + 1)/(x - 1)
    e = (-gamma(k)*gamma(k + 2) + gamma(k + 1)**2)/gamma(k)**2
    assert combsimp(e**2) == k**2
    assert combsimp(e**2/gamma(k + 1)) == k/gamma(k)
    a = R(1, 2) + R(1, 3)
    b = a + R(1, 3)
    assert combsimp(gamma(2*k)/gamma(k)*gamma(k + a)*gamma(k + b))
    3*2**(2*k + 1)*3**(-3*k - 2)*sqrt(pi)*gamma(3*k + R(3, 2))/2

    A, B = symbols('A B', commutative=False)
    assert combsimp(e*B*A) == combsimp(e)*B*A

    # check iteration
    assert combsimp(gamma(2*k)/gamma(k)*gamma(-k - R(1, 2))) == (
        -2**(2*k + 1)*sqrt(pi)/(2*((2*k + 1)*cos(pi*k))))
    assert combsimp(
        gamma(k)*gamma(k + R(1, 3))*gamma(k + R(2, 3))/gamma(3*k/2)) == (
        3*2**(3*k + 1)*3**(-3*k - S.Half)*sqrt(pi)*gamma(3*k/2 + S.Half)/2)

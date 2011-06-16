from sympy.integrals.transforms import (mellin_transform,
                                        inverse_mellin_transform)
from sympy import (gamma, exp, oo, Heaviside, symbols, re, factorial, pi,
                   cos, S, And)
from sympy.abc import x, s, a
nu, beta, rho = symbols('nu beta rho')

def test_undefined_function():
    from sympy import Function, MellinTransform
    f = Function('f')
    assert mellin_transform(f(x), x, s) == MellinTransform(f(x), x, s)
    assert mellin_transform(f(x) + exp(-x), x, s) == \
           (MellinTransform(f(x), x, s) + gamma(s), (0, oo), True)

def test_free_symbols():
    from sympy import Function
    f = Function('f')
    assert mellin_transform(f(x), x, s).free_symbols == set([s])
    assert mellin_transform(f(x)*a, x, s).free_symbols == set([s, a])

def test_as_integral():
    from sympy import Function, Integral
    f = Function('f')
    assert mellin_transform(f(x), x, s).rewrite('Integral') == \
           Integral(x**(s - 1)*f(x), (x, 0, oo))

def test_mellin_transform():
    MT = mellin_transform

    # 8.4.2
    assert MT(x**nu*Heaviside(x - 1), x, s) \
           == (1/(-nu - s), (-oo, -re(nu)), True)
    assert MT(x**nu*Heaviside(1 - x), x, s) \
           == (1/(nu + s), (-re(nu), oo), True)

    assert MT((1-x)**(beta - 1)*Heaviside(1-x), x, s) \
           == (gamma(beta)*gamma(s)/gamma(beta + s),
               (0, oo), re(-beta) < 0)
    assert MT((x-1)**(beta - 1)*Heaviside(x-1), x, s) \
           == (gamma(beta)*gamma(1 - beta - s)/gamma(1 - s),
               (-oo, -re(beta) + 1), re(-beta) < 0)

    assert MT((1+x)**(-rho), x, s) == (gamma(s)*gamma(rho-s)/gamma(rho),
                                       (0, re(rho)), True)

    # TODO also the conditions should be simplified
    assert MT(abs(1-x)**(-rho), x, s) == \
        (cos(pi*rho/2 - pi*s)*gamma(s)*gamma(rho-s)/(cos(pi*rho/2)*gamma(rho)),\
         (0, re(rho)), And(re(rho) - 1 < 0, re(rho) < 1))

    mt = MT((1-x)**(beta-1)*Heaviside(1-x)
            + a*(x-1)**(beta-1)*Heaviside(x-1), x, s)
    assert mt[1], mt[2] == ((0, -re(beta) + 1), True)
    # TODO ...

    # 8.4.2
    assert MT(exp(-x), x, s) == (gamma(s), (0, oo), True)
    assert MT(exp(-1/x), x, s) == (gamma(-s), (-oo, 0), True)

def test_inverse_mellin_transform():
    from sympy import sin, simplify, expand_func, powsimp
    IMT = inverse_mellin_transform

    assert IMT(gamma(s), s, x, (0, oo)) == exp(-x)
    assert IMT(gamma(-s), s, x, (-oo, 0)) == exp(-1/x)
    assert IMT(s/(2*s**2 - 2), s, x, (2, oo)) \
           == (x**2 + 1)*Heaviside(1 - x)/(4*x)

    # test passing "None"
    assert IMT(1/(s**2 - 1), s, x, (-1, None)) \
           == -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)
    assert IMT(1/(s**2 - 1), s, x, (None, 1)) \
           == -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)

    # test expansion of sums
    assert IMT(gamma(s) + gamma(s-1), s, x, (1, oo)) == (x + 1)*exp(-x)/x

    # test factorisation of polys
    assert simplify(expand_func(IMT(1/(s**2 + 1), s, exp(-x),
                                    (None, oo))).rewrite(sin)) \
           == sin(x)*Heaviside(1 - exp(-x))

    # test multiplicative substitution
    a, b = symbols('a b', positive=True)
    assert IMT(b**(-s/a)*factorial(s/a)/s, s, x, (0, oo)) == exp(-b*x**a)
    assert IMT(factorial(a/b + s/b)/(a+ s), s, x, (-a, oo)) == x**a*exp(-x**b)

    def simp_pows(expr): return simplify(powsimp(expr, force=True))

    # Now test the inverses of all direct transforms tested above
    assert IMT(-1/(nu + s), s, x, (-oo, None)) == x**nu*Heaviside(x - 1)
    assert IMT(1/(nu + s), s, x, (None, oo)) == x**nu*Heaviside(1 - x)
    assert IMT(gamma(beta)*gamma(s)/gamma(s + beta), s, x, (0, oo)) \
           == (1 - x)**(beta - 1)*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(1-beta-s)/gamma(1-s),
                         s, x, (-oo, None))) \
           == (x - 1)**(beta - 1)*Heaviside(x - 1)
    assert simp_pows(IMT(gamma(s)*gamma(rho-s)/gamma(rho), s, x, (0, None))) \
           == (1/(x + 1))**rho
    # TODO when better combsimp is in place, test abs(1-x)**(-rho)

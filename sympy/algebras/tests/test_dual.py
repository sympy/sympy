from sympy import symbols, simplify, exp, ln
from sympy.algebras.dual import Dual

a, b, c, d = symbols('a:d')

def test_dual_construction():
    x = Dual(a, b)
    assert x + x == Dual(2*a, 2*b)

def test_dual_power():
    x = Dual(a, b)
    assert x**2 == Dual(a**2, 2*a*b)
    assert x**4 == Dual(a**4, 4*b*a**3)

def test_dual_add():
    d1 = Dual(a, b)
    d2 = Dual(c, d)
    assert d1 + d2 == Dual(a+c, b+d)
    assert d1 + c == Dual(a+c, b)

def test_dual_mul():
    d1 = Dual(a, b)
    d2 = Dual(c, d)
    assert d1 * d2 == Dual(a*c, a*d+b*c)
    assert d1 * c == Dual(a*c, b*c)

def test_dual_div():
    d1 = Dual(a, b)
    d2 = Dual(c, d)

    # Inversion
    assert d1.inverse() == Dual(1/a, -b/a**2)
    assert d1*d1.inverse() == Dual(1, 0)
    # Regular division
    assert simplify(d1/d2 - Dual(a/c, (b*c-a*d)/c**2)) == Dual(0, 0)
    # Normalization
    assert d1.normalize() == Dual(1, b/a)

def test_dual_exp():
    x = Dual(a, b)
    assert x.exp() == Dual(exp(a), b*exp(a))

def test_dual_ln():
    x = Dual(a, b)
    assert x._ln().exp() == x

def test_nested_dual():
    d1 = Dual(a, b)
    d2 = Dual(c, d)
    z = Dual(d1, d2)
    assert z**2 == Dual(Dual(a**2, 2*a*b), Dual(2*a*c, 2*a*d + 2*b*c))

def test_autodiff():
    """
    The main reason for introducing dual numbers is automatic differentiation.
    Using nested dual numbers we can automatically compute derivatives of
    arbitrary order with a single evaluation. In this test we compute the first
    and second derivative of a function.
    """
    # Function contains both log and exp to really test it. If these work all
    # trig should also work.
    f = lambda x: x**5 + 2*x**2.3 + (x.exp() if isinstance(x, Dual) else exp(x)) +\
                  (x._ln() if isinstance(x, Dual) else ln(x))
    y = f(Dual(Dual(a, 1), Dual(1, 0)))
    assert simplify(y.a.a - f(a)) == 0
    assert simplify(y.b.a - f(a).diff()) == 0
    assert simplify(y.b.b - f(a).diff().diff()) == 0

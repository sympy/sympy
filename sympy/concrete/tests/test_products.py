from sympy import (symbols, product, factorial, rf, sqrt, cos,
                   Function, Product, Rational, Sum, oo)
from sympy.utilities.pytest import raises

a, k, n, m, x = symbols('a,k,n,m,x', integer=True)
f = Function('f')

def test_simple_products():
    assert product(2, (k, a, n)) == 2**(n-a+1)
    assert product(k, (k, 1, n)) == factorial(n)
    assert product(k**3, (k, 1, n)) == factorial(n)**3

    assert product(k+1, (k, 0, n-1)) == factorial(n)
    assert product(k+1, (k, a, n-1)) == rf(1+a, n-a)

    assert product(cos(k), (k, 0, 5)) == cos(1)*cos(2)*cos(3)*cos(4)*cos(5)
    assert product(cos(k), (k, 3, 5)) == cos(3)*cos(4)*cos(5)
    assert product(cos(k), (k, 1, Rational(5, 2))) != cos(1)*cos(2)

    assert isinstance(product(k**k, (k, 1, n)), Product)

    raises(ValueError, lambda: Product(n))
    raises(ValueError, lambda: Product(n, k))
    raises(ValueError, lambda: Product(n, k, 1))
    raises(ValueError, lambda: Product(n, k, 1, 10))
    raises(ValueError, lambda: Product(n, (k, 1)))

def test_multiple_products():
    assert product(x, (n, 1, k), (k, 1, m)) == x**(m**2/2 + m/2)
    assert product(f(n), (n, 1, m), (m, 1, k)) == Product(f(n), (n, 1, m), (m, 1, k)).doit()
    assert Product(f(n), (m, 1, k), (n, 1, k)).doit() == \
        Product(Product(f(n), (m, 1, k)), (n, 1, k)).doit() == \
        product(f(n), (m, 1, k), (n, 1, k)) == \
        product(product(f(n), (m, 1, k)), (n, 1, k)) == \
        Product(f(n)**k, (n, 1, k))
    assert Product(x, (x, 1, k), (k, 1, n)).doit() == Product(factorial(k), (k, 1, n))

def test_rational_products():
    assert product(1+1/k, (k, 1, n)) == rf(2, n)/factorial(n)

def test_special_products():
    # Wallis product
    assert product((4*k)**2 / (4*k**2-1), (k, 1, n)) == \
        4**n*factorial(n)**2/rf(Rational(1, 2), n)/rf(Rational(3, 2), n)

    # Euler's product formula for sin
    assert product(1 + a/k**2, (k, 1, n)) == \
        rf(1 - sqrt(-a), n)*rf(1 + sqrt(-a), n)/factorial(n)**2

def test__eval_product():
    from sympy.abc import i, n
    # 1710
    a = Function('a')
    assert product(2*a(i), (i, 1, n)) == 2**n * Product(a(i), (i, 1, n))
    # 1711
    assert product(2**i, (i, 1, n)) == 2**(n/2 + n**2/2)

def test_product_pow():
    # Issue 1718
    assert product(2**f(k), (k, 1, n)) == 2**Sum(f(k), (k, 1, n))
    assert product(2**(2*f(k)), (k, 1, n)) == 2**Sum(2*f(k), (k, 1, n))

def test_infinite_product():
    # Issue 2638
    assert isinstance(Product(2**(1/factorial(n)), (n, 0, oo)), Product)

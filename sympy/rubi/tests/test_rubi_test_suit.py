from sympy.core.symbol import symbols
from sympy.rubi.rubi import rubi_integrate
from sympy.functions import log, sqrt

a, b, c, d, e, f, m, x, u = symbols('a b c d e f m x u')

def test_rubi_algebriac_1_2():

    # Integrands of the form x**m
    assert rubi_integrate(x**m, x) == x**(1 + m)/(1 + m)
    assert rubi_integrate(x**100, x) == 1/101*x**101
    assert rubi_integrate(x**3, x) == 1/4*x**4
    assert rubi_integrate(x**2, x) == 1/3*x**3
    assert rubi_integrate(x, x) == 1/2*x**2

    #assert rubi_integrate(1, x) == x
    assert rubi_integrate(1/x, x) == log(x)
    assert rubi_integrate(1/x**2, x) == (-1)/x
    assert rubi_integrate(1/x**3, x) == (-1/2)/x**2
    assert rubi_integrate(1/x**4, x) == (-1/3)/x**3
    assert rubi_integrate(1/x**100, x) == (-1/99)/x**99

    # Integrands of the form x**(m/2)
    assert rubi_integrate(x**(5/2), x) == 2/7*x**(7/2)
    assert rubi_integrate(x**(3/2), x) == 2/5*x**(5/2)
    assert rubi_integrate(x**(1/2), x) == 2/3*x**(3/2)
    assert rubi_integrate(1/x**(1/2), x) == 2*sqrt(x)
    assert rubi_integrate(1/x**(3/2), x) == (-2)/sqrt(x)
    assert rubi_integrate(1/x**(5/2), x) == (-2/3)/x**(3/2)

    '''
    # Integrands of the form x**(m/3
    assert rubi_integrate(x**(5/3), x) == 3/8*x**(8/3)
    assert rubi_integrate(x**(4/3), x) == 3/7*x**(7/3)
    assert rubi_integrate(x**(2/3), x) == 3/5*x**(5/3)
    assert rubi_integrate(x**(1/3), x) == 3/4*x**(4/3)
    assert rubi_integrate(1/x**(1/3), x) == 3/2*x**(2/3)
    assert rubi_integrate(1/x**(2/3), x) == 3*x**(1/3)
    assert rubi_integrate(1/x**(4/3), x) == (-3)/x**(1/3)
    assert rubi_integrate(1/x**(5/3), x) == (-3/2)/x**(2/3)
    '''

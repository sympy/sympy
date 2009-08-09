from sympy import symbols, log
x, y, z = symbols('x y z')

def test_expand_no_log():
    assert ((1+log(x**4))**2).expand(log=False) == 1 + 2*log(x**4) + log(x**4)**2
    assert ((1+log(x**4))*(1+log(x**3))).expand(log=False) == 1 + log(x**4) + log(x**3) + log(x**4)*log(x**3)

def test_expand_no_multinomial():
    assert ((1+x)*(1+(1+x)**4)).expand(multinomial=False) == 1 + x + (1+x)**4 + x*(1+x)**4

def test_expand_negative_integer_powers():
    expr = (x+y)**(-2)
    assert expr.expand() == 1 / (2*x*y + x**2 + y**2)
    assert expr.expand(multinomial=False) == (x+y)**(-2)
    expr = (x+y)**(-3)
    assert expr.expand() == 1 / (3*x*x*y + 3*x*y*y + x**3 + y**3)
    assert expr.expand(multinomial=False) == (x+y)**(-3)
    expr = (x+y)**(2) * (x+y)**(-4)
    assert expr.expand() == 1 / (2*x*y + x**2 + y**2)
    assert expr.expand(multinomial=False) == (x+y)**(-2)

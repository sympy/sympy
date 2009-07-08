from sympy import Symbol, log
x = Symbol('x')

def test_expand_no_log():
    assert ((1+log(x**4))**2).expand(log=False) == 1 + 2*log(x**4) + log(x**4)**2
    assert ((1+log(x**4))*(1+log(x**3))).expand(log=False) == 1 + log(x**4) + log(x**3) + log(x**4)*log(x**3)

def test_expand_no_multinomial():
    assert ((1+x)*(1+(1+x)**4)).expand(multinomial=False) == 1 + x + (1+x)**4 + x*(1+x)**4

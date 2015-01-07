from sympy.parsing.mathematica import mathematica
from sympy import sympify


def test_mathematica():
    d = {
        '- 6x': '-6*x',
        'Sin[x]^2': 'sin(x)**2',
        '2(x-1)': '2*(x-1)',
        '3y+8': '3*y+8',
        'Arcsin[2x+9(4-x)^2]/x': 'asin(2*x+9*(4-x)**2)/x',
        'x+y': 'x+y',
        '355/113': '355/113',
        '2.718281828': '2.718281828',
        'Sin[12]': 'sin(12)',
        'Exp[Log[4]]': 'exp(log(4))',
        '(x+1)(x+3)': '(x+1)*(x+3)',
        'Cos[Arccos[3.6]]': 'cos(acos(3.6))',
        'Cos[x]==Sin[y]': 'cos(x)==sin(y)',
        '2*Sin[x+y]': '2*sin(x+y)',
        'Sin[x]+Cos[y]': 'sin(x)+cos(y)',
        'Sin[Cos[x]]': 'sin(cos(x))',
        '2*Sqrt[x+y]': '2*sqrt(x+y)',   # Test case from the issue 4259
        'x y': 'x*y',
        'x Sin[x]': 'x*sin(x)',
        'x Sin[1/x]': 'x*sin(1/x)',
        'Sin[1/x] x': 'x*sin(1/x)',
        'x*Sin[1/x]': 'x*sin(1/x)',   # Test case from the issue 8501
        'Sin[1/x]*x': 'x*sin(1/x)'}
    for e in d:
        assert mathematica(e) == sympify(d[e])

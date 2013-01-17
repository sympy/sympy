from sympy.parsing.mathematica import mathematica
from sympy import sympify


def test_mathematica():
    d = {'Sin[x]^2': 'sin(x)**2',
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
        'Cos[x]==Sin[y]': 'cos(x)==sin(y)'}
    for e in d:
        assert mathematica(e) == sympify(d[e])

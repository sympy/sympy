from sympy import Integral, pretty, latex, Function
from sympy.vector import CoordSysCartesian, Vector, Dyadic, express
from sympy.abc import a, b, c
from sympy.core.compatibility import u


#Initialize the basic and tedious vector/dyadic expressions
#needed for testing
N = CoordSysCartesian('N')
C = N.orient_new_axis('C', a, N.k)
v = []
d = []
v.append(Vector.zero)
v.append(N.i)
v.append(-N.i)
v.append(N.i + N.j)
v.append(a*N.i)
v.append(a*N.i - b*N.j)
v.append((a**2 + N.x)*N.i + N.k)
v.append((a**2 + b)*N.i + 3*(C.y - c)*N.k)
f = Function('f')
v.append(N.j - (Integral(f(b)) - C.x**2)*N.k)
v.append(N.i + C.k)
v.append(express(N.i, C))
v.append((a**2 + b)*N.i + (Integral(f(b)))*N.k)
for x in v:
    d.append(x | N.k)
s = 3*N.x**2*C.y


def test_str_printing():
    assert str(v[0]) == '0'
    assert str(v[1]) == 'N.i'
    assert str(v[2]) == '(-1)*N.i'
    assert str(v[3]) == 'N.i + N.j'
    assert str(v[8]) == 'N.j + (C.x**2 - Integral(f(b), b))*N.k'
    assert str(v[9]) == 'C.k + N.i'
    assert str(s) == '3*C.y*N.x**2'
    assert str(d[0]) == '0'
    assert str(d[1]) == '(N.i|N.k)'
    assert str(d[4]) == 'a*(N.i|N.k)'
    assert str(d[5]) == 'a*(N.i|N.k) + (-b)*(N.j|N.k)'
    assert str(d[8]) == ('(N.j|N.k) + (C.x**2 - ' +
                         'Integral(f(b), b))*(N.k|N.k)')


def test_pretty_printing():
    assert pretty(v[0]) == u('0')
    assert pretty(v[1]) == u('N_i')
    assert pretty(v[5]) == u('(a) N_i + (-b) N_j')
    assert pretty(v[8]) == u('N_j + \u239b   2   \u2320        \u239e '+
                             'N_k\n      \u239cC_x  - \u23ae f(b) db\u2' +
                             '39f    \n      \u239d       \u2321        ' +
                             '\u23a0    ')
    assert pretty(v[2]) == u('(-1) N_i')
    assert pretty(v[11]) == u('\u239b 2    \u239e N_i + \u239b\u2320     '+
                              '   \u239e N_k\n\u239da  + b\u23a0       \u2' +
                              '39c\u23ae f(b) db\u239f    \n         ' +
                              '      \u239d\u2321        \u23a0    ')
    assert pretty(s) == u('         2\n3\u22c5C_y\u22c5N_x ')
    assert pretty(d[0]) == u('(0|0)')
    assert pretty(d[5]) == u('(a) (N_i|N_k) + (-b) (N_j|N_k)')
    assert pretty(d[7]) == u('\u239b 2    \u239e (N_i|N_k) + (3\u22c5'+
                             'C_y - 3\u22c5c) (N_k|N_k)\n\u239da  + '+
                             'b\u23a0          ')
    assert pretty(d[10]) == u('(cos(a)) (C_i|N_k) + (-sin(a)) (C_j|N_k)')


def test_latex_printing():
    assert latex(v[0]) == '\\mathbf{\\hat{0}}'
    assert latex(v[1]) == '\\mathbf{\\hat{i}_{N}}'
    assert latex(v[2]) == '- \\mathbf{\\hat{i}_{N}}'
    assert latex(v[5]) == ('(a)\\mathbf{\\hat{i}_{N}} + ' +
                           '(- b)\\mathbf{\\hat{j}_{N}}')
    assert latex(v[6]) == ('(\\mathbf{{x}_{N}} + a^{2})\\mathbf{\\' +
                           'hat{i}_{N}} + \\mathbf{\\hat{k}_{N}}')
    assert latex(v[8]) == ('\\mathbf{\\hat{j}_{N}} + (\\mathbf{{x}_' +
                           '{C}}^{2} - \\int f{\\left (b \\right )}\\,' +
                           ' db)\\mathbf{\\hat{k}_{N}}')
    assert latex(s) == '3 \\mathbf{{y}_{C}} \\mathbf{{x}_{N}}^{2}'
    assert latex(d[0]) == '(\\mathbf{\\hat{0}}|\\mathbf{\\hat{0}})'
    assert latex(d[4]) == ('(a)(\\mathbf{\\hat{i}_{N}}{|}\\mathbf' +
                           '{\\hat{k}_{N}})')
    assert latex(d[9]) == ('(\\mathbf{\\hat{k}_{C}}{|}\\mathbf{\\' +
                           'hat{k}_{N}}) + (\\mathbf{\\hat{i}_{N}}{|' +
                           '}\\mathbf{\\hat{k}_{N}})')
    assert latex(d[11]) == ('(a^{2} + b)(\\mathbf{\\hat{i}_{N}}{|}\\' +
                            'mathbf{\\hat{k}_{N}}) + (\\int f{\\left (' +
                            'b \\right )}\\, db)(\\mathbf{\\hat{k}_{N}' +
                            '}{|}\\mathbf{\\hat{k}_{N}})')

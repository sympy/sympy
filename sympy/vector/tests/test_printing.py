# -*- coding: utf-8 -*-
from sympy import Integral, latex, Function
from sympy import pretty as xpretty
from sympy.vector import CoordSys3D, Vector, express
from sympy.abc import a, b, c
from sympy.core.compatibility import u_decode as u
from sympy.testing.pytest import XFAIL


def pretty(expr):
    """ASCII pretty-printing"""
    return xpretty(expr, use_unicode=False, wrap_line=False)


def upretty(expr):
    """Unicode pretty-printing"""
    return xpretty(expr, use_unicode=True, wrap_line=False)


# Initialize the basic and tedious vector/dyadic expressions
# needed for testing.
# Some of the pretty forms shown denote how the expressions just
# above them should look with pretty printing.
N = CoordSys3D('N')
C = N.orient_new_axis('C', a, N.k)  # type: ignore
v = []
d = []
v.append(Vector.zero)
v.append(N.i)  # type: ignore
v.append(-N.i)  # type: ignore
v.append(N.i + N.j)  # type: ignore
v.append(a*N.i)  # type: ignore
v.append(a*N.i - b*N.j)  # type: ignore
v.append((a**2 + N.x)*N.i + N.k)  # type: ignore
v.append((a**2 + b)*N.i + 3*(C.y - c)*N.k)  # type: ignore
f = Function('f')
v.append(N.j - (Integral(f(b)) - C.x**2)*N.k)  # type: ignore
upretty_v_8 = u(
"""\
      ⎛   2   ⌠        ⎞    \n\
j_N + ⎜x_C  - ⎮ f(b) db⎟ k_N\n\
      ⎝       ⌡        ⎠    \
""")
pretty_v_8 = u(
    """\
j_N + /         /       \\\n\
      |   2    |        |\n\
      |x_C  -  | f(b) db|\n\
      |        |        |\n\
      \\       /         / \
""")

v.append(N.i + C.k)  # type: ignore
v.append(express(N.i, C))  # type: ignore
v.append((a**2 + b)*N.i + (Integral(f(b)))*N.k)  # type: ignore
upretty_v_11 = u(
"""\
⎛ 2    ⎞        ⎛⌠        ⎞    \n\
⎝a  + b⎠ i_N  + ⎜⎮ f(b) db⎟ k_N\n\
                ⎝⌡        ⎠    \
""")
pretty_v_11 = u(
"""\
/ 2    \\ + /  /       \\\n\
\\a  + b/ i_N| |        |\n\
           | | f(b) db|\n\
           | |        |\n\
           \\/         / \
""")

for x in v:
    d.append(x | N.k)  # type: ignore
s = 3*N.x**2*C.y  # type: ignore
upretty_s = u(
"""\
         2\n\
3⋅y_C⋅x_N \
""")
pretty_s = u(
"""\
         2\n\
3*y_C*x_N \
""")

# This is the pretty form for ((a**2 + b)*N.i + 3*(C.y - c)*N.k) | N.k
upretty_d_7 = u(
"""\
⎛ 2    ⎞                                     \n\
⎝a  + b⎠ (i_N|k_N)  + (3⋅y_C - 3⋅c) (k_N|k_N)\
""")
pretty_d_7 = u(
"""\
/ 2    \\ (i_N|k_N) + (3*y_C - 3*c) (k_N|k_N)\n\
\\a  + b/                                    \
""")


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


@XFAIL
def test_pretty_printing_ascii():
    assert pretty(v[0]) == '0'
    assert pretty(v[1]) == 'i_N'
    assert pretty(v[5]) == '(a) i_N + (-b) j_N'
    assert pretty(v[8]) == pretty_v_8
    assert pretty(v[2]) == '(-1) i_N'
    assert pretty(v[11]) == pretty_v_11
    assert pretty(s) == pretty_s
    assert pretty(d[0]) == '(0|0)'
    assert pretty(d[5]) == '(a) (i_N|k_N) + (-b) (j_N|k_N)'
    assert pretty(d[7]) == pretty_d_7
    assert pretty(d[10]) == '(cos(a)) (i_C|k_N) + (-sin(a)) (j_C|k_N)'


def test_pretty_print_unicode_v():
    assert upretty(v[0]) == '0'
    assert upretty(v[1]) == 'i_N'
    assert upretty(v[5]) == '(a) i_N + (-b) j_N'
    # Make sure the printing works in other objects
    assert upretty(v[5].args) == '((a) i_N, (-b) j_N)'
    assert upretty(v[8]) == upretty_v_8
    assert upretty(v[2]) == '(-1) i_N'
    assert upretty(v[11]) == upretty_v_11
    assert upretty(s) == upretty_s
    assert upretty(d[0]) == '(0|0)'
    assert upretty(d[5]) == '(a) (i_N|k_N) + (-b) (j_N|k_N)'
    assert upretty(d[7]) == upretty_d_7
    assert upretty(d[10]) == '(cos(a)) (i_C|k_N) + (-sin(a)) (j_C|k_N)'


def test_latex_printing():
    assert latex(v[0]) == '\\mathbf{\\hat{0}}'
    assert latex(v[1]) == '\\mathbf{\\hat{i}_{N}}'
    assert latex(v[2]) == '- \\mathbf{\\hat{i}_{N}}'
    assert latex(v[5]) == ('(a)\\mathbf{\\hat{i}_{N}} + ' +
                           '(- b)\\mathbf{\\hat{j}_{N}}')
    assert latex(v[6]) == ('(\\mathbf{{x}_{N}} + a^{2})\\mathbf{\\hat{i}_' +
                          '{N}} + \\mathbf{\\hat{k}_{N}}')
    assert latex(v[8]) == ('\\mathbf{\\hat{j}_{N}} + (\\mathbf{{x}_' +
                           '{C}}^{2} - \\int f{\\left(b \\right)}\\,' +
                           ' db)\\mathbf{\\hat{k}_{N}}')
    assert latex(s) == '3 \\mathbf{{y}_{C}} \\mathbf{{x}_{N}}^{2}'
    assert latex(d[0]) == '(\\mathbf{\\hat{0}}|\\mathbf{\\hat{0}})'
    assert latex(d[4]) == ('(a)(\\mathbf{\\hat{i}_{N}}{|}\\mathbf' +
                           '{\\hat{k}_{N}})')
    assert latex(d[9]) == ('(\\mathbf{\\hat{k}_{C}}{|}\\mathbf{\\' +
                           'hat{k}_{N}}) + (\\mathbf{\\hat{i}_{N}}{|' +
                           '}\\mathbf{\\hat{k}_{N}})')
    assert latex(d[11]) == ('(a^{2} + b)(\\mathbf{\\hat{i}_{N}}{|}\\' +
                            'mathbf{\\hat{k}_{N}}) + (\\int f{\\left(' +
                            'b \\right)}\\, db)(\\mathbf{\\hat{k}_{N}' +
                            '}{|}\\mathbf{\\hat{k}_{N}})')


def test_custom_names():
    A = CoordSys3D('A', vector_names=['x', 'y', 'z'],
                   variable_names=['i', 'j', 'k'])
    assert A.i.__str__() == 'A.i'
    assert A.x.__str__() == 'A.x'
    assert A.i._pretty_form == 'i_A'
    assert A.x._pretty_form == 'x_A'
    assert A.i._latex_form == r'\mathbf{{i}_{A}}'
    assert A.x._latex_form == r"\mathbf{\hat{x}_{A}}"

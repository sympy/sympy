from sympy import symbols, sin, cos, sqrt
from sympy.physics.vector import dynamicsymbols, ReferenceFrame
from sympy.physics.vector.printers import VectorPrettyPrinter, VectorLatexPrinter

a, b, c = symbols('a, b, c')
alpha, omega, beta = dynamicsymbols('alpha, omega, beta')

N = ReferenceFrame('N')

v = a ** 2 * N.x + b * N.y + c * sin(alpha) * N.z
w = alpha * N.x + sin(omega) * N.y + alpha * beta * N.z

# TODO : The division does not print correctly here:
# w = alpha * N.x + sin(omega) * N.y + alpha / beta * N.z

u = a ** 2 * (N.x | N.y) + b * (N.y | N.y) + c * sin(alpha) * (N.z | N.y)
x = alpha * (N.x | N.x) + sin(omega) * (N.y | N.z) + alpha * beta * (N.z | N.x)


def test_vector_pretty_print():

    pp = VectorPrettyPrinter()

    expected = (u' 2\na *\x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m + '
                u'b*\x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + '
                u'c\u22c5sin(\u03b1)*\x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m')

    assert expected == pp.doprint(v)

    expected = (u'\u03b1*\x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m + '
                u'sin(\u03c9)*\x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + '
                u'\u03b1\u22c5\u03b2*\x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m')

    assert expected == pp.doprint(w)


def test_vector_latex():

    N = ReferenceFrame('N')

    a, b, c, d, omega = symbols('a, b, c, d, omega')

    v = (a ** 2 + b / c) * N.x + sqrt(d) * N.y + cos(omega) * N.z

    expected = ('(a^{2} + \\frac{b}{c})\\mathbf{\\hat{n}_x} + '
                '\\sqrt{d}\\mathbf{\\hat{n}_y} + '
                '\\operatorname{cos}\\left(\\omega\\right)\\mathbf{\\hat{n}_z}')

    assert v._latex() == expected
    lp = VectorLatexPrinter()
    assert lp.doprint(v) == expected

    N = ReferenceFrame('N', latexs=(r'\hat{i}', r'\hat{j}', r'\hat{k}'))

    v = (a ** 2 + b / c) * N.x + sqrt(d) * N.y + cos(omega) * N.z

    expected = ('(a^{2} + \\frac{b}{c})\\hat{i} + '
                '\\sqrt{d}\\hat{j} + '
                '\\operatorname{cos}\\left(\\omega\\right)\\hat{k}')
    assert v._latex() == expected


def test_dyadic_pretty_print():

    expected = (u' 2\na  \x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m\u2a02 \x1b[9'
                u'4m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + b \x1b[94m\x1b[1mn_y'
                u'\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b'
                u'[0;0m + c\u22c5sin(\u03b1) \x1b[94m\x1b[1mn_z\x1b[0;0m'
                u'\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m')

    result = u._pretty().render()

    assert expected == result

    expected = (u'\u03b1 \x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m\u2a02 \x1b[9'
                u'4m\x1b[1mn_x\x1b[0;0m\x1b[0;0m + sin(\u03c9) \x1b[94m'
                u'\x1b[1mn_y\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_z'
                u'\x1b[0;0m\x1b[0;0m + \u03b1\u22c5\u03b2 \x1b[94m\x1b[1mn'
                u'_z\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_x\x1b[0;0m'
                u'\x1b[0;0m')

    result = x._pretty().render()

    assert expected == result

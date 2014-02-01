from sympy import symbols, sin
from sympy.physics.vector import dynamicsymbols, ReferenceFrame
from sympy.physics.vector.printers import VectorPrettyPrinter

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

    expected = (u' 2\na *\x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m + ' +
                u'b*\x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + ' +
                u'c\u22c5sin(\u03b1)*\x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m')

    assert expected == pp.doprint(v)

    expected = (u'\u03b1*\x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m + ' +
                u'sin(\u03c9)*\x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + ' +
                u'\u03b1\u22c5\u03b2*\x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m')

    assert expected == pp.doprint(w)


def test_dyadic_pretty_print():

    expected = (u' 2\na  \x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + ' +
                u'b \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m + ' +
                u'c\u22c5sin(\u03b1) \x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m')

    result = u._pretty().render()

    assert expected == result

    expected = (u'\u03b1 \x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m + ' +
                u'sin(\u03c9) \x1b[94m\x1b[1mn_y\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m + ' +
                u'\u03b1\u22c5\u03b2 \x1b[94m\x1b[1mn_z\x1b[0;0m\x1b[0;0m\u2a02 \x1b[94m\x1b[1mn_x\x1b[0;0m\x1b[0;0m')

    result = x._pretty().render()

    assert expected == result

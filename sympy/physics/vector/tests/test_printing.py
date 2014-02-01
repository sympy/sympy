from sympy import symbols, sin, cos, sqrt, Function
from sympy.physics.vector import ReferenceFrame, dynamicsymbols
from sympy.physics.vector.printing import VectorLatexPrinter

A = ReferenceFrame('A')


def test_latex_printer():
    r = Function('r')('t')
    assert VectorLatexPrinter().doprint(r**2) == "r^{2}"


def test_vector_latex():

    a, b, c, d, omega = symbols('a, b, c, d, omega')

    v = (a ** 2 + b / c) * A.x + sqrt(d) * A.y + cos(omega) * A.z

    assert v._latex() == ('(a^{2} + \\frac{b}{c})\\mathbf{\\hat{a}_x} + ' +
        '\\sqrt{d}\\mathbf{\\hat{a}_y} + ' +
        '\\operatorname{cos}\\left(\\omega\\right)\\mathbf{\\hat{a}_z}')

    theta, omega, alpha, q = dynamicsymbols('theta, omega, alpha, q')

    v = theta * A.x + omega * omega * A.y + (q * alpha) * A.z

    # NOTE : "gamma" is not an accepted latex symbol it is always overidden
    # by SymPy's printing for Gamma functions. It would be best to not use
    # "gamma" as a symbol or dynamicsymbol name when using the mechanics
    # package. I noticed this when writing this test which originally
    # included a "gamma" dynamicsymbol.

    assert v._latex() == ('\\theta\\mathbf{\\hat{a}_x} + ' +
                          '\\omega^{2}\\mathbf{\\hat{a}_y} + ' +
                          '\\alpha q\\mathbf{\\hat{a}_z}')

    phi1, phi2, phi3 = dynamicsymbols('phi1, phi2, phi3')
    theta1, theta2, theta3 = symbols('theta1, theta2, theta3')

    v = (sin(theta1) * A.x +
         cos(phi1) * cos(phi2) * A.y +
         cos(theta1 + phi3) * A.z)

    assert v._latex() == ('\\operatorname{sin}\\left(\\theta_{1}\\right)\\mathbf{\\hat{a}_x} + ' +
                          '\\operatorname{cos}\\left(\\phi_{1}\\right) \\operatorname{cos}\\left(\\phi_{2}\\right)\\mathbf{\\hat{a}_y} + ' +
                          '\\operatorname{cos}\\left(\\theta_{1} + \\phi_{3}\\right)\\mathbf{\\hat{a}_z}')

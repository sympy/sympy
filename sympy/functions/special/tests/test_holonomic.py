from sympy.functions.special.holonomic import DifferentialOperator, HoloFunc, DiffOperatorAlgebra
from sympy import symbols
from sympy import ZZ, QQ, RR


def test_DifferentialOperator():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(QQ.old_poly_ring(x), 'Dx')
    assert Dx == R.derivative_operator
    assert Dx == DifferentialOperator([R.base.zero, R.base.one], R)
    assert x * Dx + x**2 * Dx**2 == DifferentialOperator([0, x, x**2], R)
    assert (x**2 + 1) + Dx + x * \
        Dx**5 == DifferentialOperator([x**2 + 1, 1, 0, 0, 0, x], R)
    assert (x * Dx + x**2 + 1 - Dx * (x**3 + x))**3 == (-48 * x**6) + \
        (-57 * x**7) * Dx + (-15 * x**8) * Dx**2 + (-x**9) * Dx**3
    p = (x * Dx**2 + (x**2 + 3) * Dx**5) * (Dx + x**2)
    q = (2 * x) + (4 * x**2) * Dx + (x**3) * Dx**2 + \
        (20 * x**2 + x + 60) * Dx**3 + (10 * x**3 + 30 * x) * Dx**4 + \
        (x**4 + 3 * x**2) * Dx**5 + (x**2 + 3) * Dx**6
    assert p == q


def test_HoloFunc_addition():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    p = HoloFunc(Dx**2 * x, x)
    q = HoloFunc((2) * Dx + (x) * Dx**2, x)
    assert p == q
    p = HoloFunc(x * Dx + 1, x)
    q = HoloFunc(Dx + 1, x)
    r = HoloFunc((x - 2) + (x**2 - 2) * Dx + (x**2 - x) * Dx**2, x)
    assert p + q == r
    p = HoloFunc(x * Dx + Dx**2 * (x**2 + 2), x)
    q = HoloFunc(Dx - 3, x)
    r = HoloFunc((-54 * x**2 - 126 * x - 150) + (-135 * x**3 - 252 * x**2 - 270 * x + 140) * Dx +\
                 (-27 * x**4 - 24 * x**2 + 14 * x - 150) * Dx**2 + \
                 (9 * x**4 + 15 * x**3 + 38 * x**2 + 30 * x +40) * Dx**3, x)
    assert p + q == r
    p = HoloFunc(Dx**5 - 1, x)
    q = HoloFunc(x**3 + Dx, x)
    r = HoloFunc((-x**18 + 45*x**14 - 525*x**10 + 1575*x**6 - x**3 - 630*x**2) + \
        (-x**15 + 30*x**11 - 195*x**7 + 210*x**3 - 1)*Dx + (x**18 - 45*x**14 + 525*x**10 - \
        1575*x**6 + x**3 + 630*x**2)*Dx**5 + (x**15 - 30*x**11 + 195*x**7 - 210*x**3 + \
        1)*Dx**6, x)
    assert p+q == r

def test_HoloFunc_multiplication():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    p = HoloFunc(Dx+x+x*Dx**2,x)
    q = HoloFunc(x*Dx+Dx*x+Dx**2,x)
    r = HoloFunc((8*x**6 + 4*x**4 + 6*x**2 + 3) + (24*x**5 - 4*x**3 + 24*x)*Dx + \
        (8*x**6 + 20*x**4 + 12*x**2 + 2)*Dx**2 + (8*x**5 + 4*x**3 + 4*x)*Dx**3 + \
        (2*x**4 + x**2)*Dx**4, x)
    assert p*q == r
    p = HoloFunc(Dx**2+1,x)
    q = HoloFunc(Dx-1,x)
    r = HoloFunc((2) + (-2)*Dx + (1)*Dx**2, x)
    assert p*q == r
    p = HoloFunc(Dx**2+1+x+Dx,x)
    q = HoloFunc((Dx*x-1)**2,x)
    r = HoloFunc((4*x**7 + 11*x**6 + 16*x**5 + 4*x**4 - 6*x**3 - 7*x**2 - 8*x - 2) + \
        (8*x**6 + 26*x**5 + 24*x**4 - 3*x**3 - 11*x**2 - 6*x - 2)*Dx + \
        (8*x**6 + 18*x**5 + 15*x**4 - 3*x**3 - 6*x**2 - 6*x - 2)*Dx**2 + (8*x**5 + \
            10*x**4 + 6*x**3 - 2*x**2 - 4*x)*Dx**3 + (4*x**5 + 3*x**4 - x**2)*Dx**4, x)
    assert p*q == r
    p = HoloFunc(x*Dx**2-1,x)
    q = HoloFunc(Dx*x-x,x)
    r = HoloFunc((x - 3) + (-2*x + 2)*Dx + (x)*Dx**2, x)
    assert p*q == r

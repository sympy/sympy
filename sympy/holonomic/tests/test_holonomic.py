from sympy.holonomic.holonomic import DifferentialOperator, HolonomicFunction, DiffOperatorAlgebra
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


def test_HolonomicFunction_addition():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx**2 * x, x)
    q = HolonomicFunction((2) * Dx + (x) * Dx**2, x)
    assert p == q
    p = HolonomicFunction(x * Dx + 1, x)
    q = HolonomicFunction(Dx + 1, x)
    r = HolonomicFunction((x - 2) + (x**2 - 2) * Dx + (x**2 - x) * Dx**2, x)
    assert p + q == r
    p = HolonomicFunction(x * Dx + Dx**2 * (x**2 + 2), x)
    q = HolonomicFunction(Dx - 3, x)
    r = HolonomicFunction((-54 * x**2 - 126 * x - 150) + (-135 * x**3 - 252 * x**2 - 270 * x + 140) * Dx +\
                 (-27 * x**4 - 24 * x**2 + 14 * x - 150) * Dx**2 + \
                 (9 * x**4 + 15 * x**3 + 38 * x**2 + 30 * x +40) * Dx**3, x)
    assert p + q == r
    p = HolonomicFunction(Dx**5 - 1, x)
    q = HolonomicFunction(x**3 + Dx, x)
    r = HolonomicFunction((-x**18 + 45*x**14 - 525*x**10 + 1575*x**6 - x**3 - 630*x**2) + \
        (-x**15 + 30*x**11 - 195*x**7 + 210*x**3 - 1)*Dx + (x**18 - 45*x**14 + 525*x**10 - \
        1575*x**6 + x**3 + 630*x**2)*Dx**5 + (x**15 - 30*x**11 + 195*x**7 - 210*x**3 + \
        1)*Dx**6, x)
    assert p+q == r

def test_HolonomicFunction_multiplication():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx+x+x*Dx**2, x)
    q = HolonomicFunction(x*Dx+Dx*x+Dx**2, x)
    r = HolonomicFunction((8*x**6 + 4*x**4 + 6*x**2 + 3) + (24*x**5 - 4*x**3 + 24*x)*Dx + \
        (8*x**6 + 20*x**4 + 12*x**2 + 2)*Dx**2 + (8*x**5 + 4*x**3 + 4*x)*Dx**3 + \
        (2*x**4 + x**2)*Dx**4, x)
    assert p*q == r
    p = HolonomicFunction(Dx**2+1, x)
    q = HolonomicFunction(Dx-1, x)
    r = HolonomicFunction((2) + (-2)*Dx + (1)*Dx**2, x)
    assert p*q == r
    p = HolonomicFunction(Dx**2+1+x+Dx, x)
    q = HolonomicFunction((Dx*x-1)**2, x)
    r = HolonomicFunction((4*x**7 + 11*x**6 + 16*x**5 + 4*x**4 - 6*x**3 - 7*x**2 - 8*x - 2) + \
        (8*x**6 + 26*x**5 + 24*x**4 - 3*x**3 - 11*x**2 - 6*x - 2)*Dx + \
        (8*x**6 + 18*x**5 + 15*x**4 - 3*x**3 - 6*x**2 - 6*x - 2)*Dx**2 + (8*x**5 + \
            10*x**4 + 6*x**3 - 2*x**2 - 4*x)*Dx**3 + (4*x**5 + 3*x**4 - x**2)*Dx**4, x)
    assert p*q == r
    p = HolonomicFunction(x*Dx**2-1, x)
    q = HolonomicFunction(Dx*x-x, x)
    r = HolonomicFunction((x - 3) + (-2*x + 2)*Dx + (x)*Dx**2, x)
    assert p*q == r

def test_addition_initial_condition():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(QQ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx-1, x, 0, 3)
    q = HolonomicFunction(Dx**2+1, x, 0, [1, 0])
    r = HolonomicFunction(-1 + Dx - Dx**2 + Dx**3, x, 0, [4, 3, 2])
    assert p + q == r
    p = HolonomicFunction(Dx - x + Dx**2, x, 0, [1, 2])
    q = HolonomicFunction(Dx**2 + x, x, 0, [1, 0])
    r = HolonomicFunction((-4*x**4 - x**3 - 4*x**2 + 1) + (4*x**3 + x**2 + 3*x + 4)*Dx + \
        (-6*x + 7)*Dx**2 + (4*x**2 - 7*x + 1)*Dx**3 + (4*x**2 + x + 2)*Dx**4, x, 0, [2, 2, -2, 2])
    assert p + q == r
    p = HolonomicFunction(Dx**2 + 4*x*Dx + x**2, x, 0, [3, 4])
    q = HolonomicFunction(Dx**2 + 1, x, 0, [1, 1])
    r = HolonomicFunction((x**6 + 2*x**4 - 5*x**2 - 6) + (4*x**5 + 36*x**3 - 32*x)*Dx + \
         (x**6 + 3*x**4 + 5*x**2 - 9)*Dx**2 + (4*x**5 + 36*x**3 - 32*x)*Dx**3 + (x**4 + \
            10*x**2 - 3)*Dx**4, x, 0, [4, 5, -1, -17])
    q = HolonomicFunction(Dx**3 + x, x, 2, [3, 0, 1])
    p = HolonomicFunction(Dx - 1, x, 2, [1])
    r = HolonomicFunction((-x**2 - x + 1) + (x**2 + x)*Dx + (-x - 2)*Dx**3 + \
        (x + 1)*Dx**4, x, 2, [4, 1, 2, -5 ])
    assert p + q == r

def test_multiplication_initial_condition():
    x = symbols('x')
    R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx**2 + x*Dx - 1, x, 0, [3, 1])
    q = HolonomicFunction(Dx**2 + 1, x, 0, [1, 1])
    r = HolonomicFunction((x**4 + 14*x**2 + 60) + 4*x*Dx + (x**4 + 9*x**2 + 20)*Dx**2 + \
        (2*x**3 + 18*x)*Dx**3 + (x**2 + 10)*Dx**4, x, 0, [3, 4, 2, 3])
    assert p * q == r
    p = HolonomicFunction(Dx**2 + x, x, 0, [1, 0])
    q = HolonomicFunction(Dx**3 - x**2, x, 0, [3, 3, 3])
    r = HolonomicFunction((27*x**8 - 37*x**7 - 10*x**6 - 492*x**5 - 552*x**4 + 160*x**3 + \
        1212*x**2 + 216*x + 360) + (162*x**7 - 384*x**6 - 294*x**5 - 84*x**4 + 24*x**3 + \
        756*x**2 + 120*x - 1080)*Dx + (81*x**6 - 246*x**5 + 228*x**4 + 36*x**3 + \
        660*x**2 - 720*x)*Dx**2 + (-54*x**6 + 128*x**5 - 18*x**4 - 240*x**2 + 600)*Dx**3 + \
        (81*x**5 - 192*x**4 - 84*x**3 + 162*x**2 - 60*x - 180)*Dx**4 + (-108*x**3 + \
        192*x**2 + 72*x)*Dx**5 + (27*x**4 - 64*x**3 - 36*x**2 + 60)*Dx**6, x, 0, [3, 3, 3, -3, -12, -24])
    assert p * q == r
    p = HolonomicFunction(Dx - 1, x, 0, [2])
    q = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])
    r = HolonomicFunction(2 -2*Dx + Dx**2, x, 0, [0, 2])
    assert p * q == r
    q = HolonomicFunction(x*Dx**2+1+2*Dx,x,0,[0,1])
    r = HolonomicFunction((x - 1) + (-2*x + 2)*Dx + x*Dx**2, x, 0, [0, 2])
    assert p * q == r
    p = HolonomicFunction(Dx**2 - 1, x, 0, [1, 3])
    q = HolonomicFunction(Dx**3 + 1, x, 0, [1, 2, 1])
    r = HolonomicFunction(6*Dx + 3*Dx**2 + 2*Dx**3 - 3*Dx**4 + Dx**6, x, 0, [1, 5, 14, 17, 17, 2])
    assert p * q == r

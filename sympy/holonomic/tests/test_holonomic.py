from sympy.holonomic import (DifferentialOperator, HolonomicFunction,
    DifferentialOperators, from_hyper, from_meijerg, from_sympy)
from sympy.holonomic.recurrence import RecurrenceOperators, HolonomicSequence
from sympy import (symbols, hyper, S, sqrt, pi, exp, erf, erfc, sstr,
    O, I, meijerg, sin, cos, log, cosh, besselj, hyperexpand, Ci, EulerGamma, Si)
from sympy import ZZ, QQ


def test_DifferentialOperator():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
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
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
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

    p = x**2 + 3*x + 8
    q = x**3 - 7*x + 5
    p = p*Dx - p.diff()
    q = q*Dx - q.diff()
    r = HolonomicFunction(p, x) + HolonomicFunction(q, x)
    s = HolonomicFunction((6*x**2 + 18*x + 14) + (-4*x**3 - 18*x**2 - 62*x + 10)*Dx +\
        (x**4 + 6*x**3 + 31*x**2 - 10*x - 71)*Dx**2, x)
    assert r == s


def test_HolonomicFunction_multiplication():
    x = symbols('x')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
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
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx-1, x, 0, 3)
    q = HolonomicFunction(Dx**2+1, x, 0, [1, 0])
    r = HolonomicFunction(-1 + Dx - Dx**2 + Dx**3, x, 0, [4, 3, 2])
    assert p + q == r
    p = HolonomicFunction(Dx - x + Dx**2, x, 0, [1, 2])
    q = HolonomicFunction(Dx**2 + x, x, 0, [1, 0])
    r = HolonomicFunction((-x**4 - x**3/4 - x**2 + 1/4) + (x**3 + x**2/4 + 3*x/4 + 1)*Dx + \
        (-3*x/2 + 7/4)*Dx**2 + (x**2 - 7*x/4 + 1/4)*Dx**3 + (x**2 + x/4 + 1/2)*Dx**4, x, 0, [2, 2, -2, 2])
    assert p + q == r
    p = HolonomicFunction(Dx**2 + 4*x*Dx + x**2, x, 0, [3, 4])
    q = HolonomicFunction(Dx**2 + 1, x, 0, [1, 1])
    r = HolonomicFunction((x**6 + 2*x**4 - 5*x**2 - 6) + (4*x**5 + 36*x**3 - 32*x)*Dx + \
         (x**6 + 3*x**4 + 5*x**2 - 9)*Dx**2 + (4*x**5 + 36*x**3 - 32*x)*Dx**3 + (x**4 + \
            10*x**2 - 3)*Dx**4, x, 0, [4, 5, -1, -17])
    assert p + q == r
    q = HolonomicFunction(Dx**3 + x, x, 2, [3, 0, 1])
    p = HolonomicFunction(Dx - 1, x, 2, [1])
    r = HolonomicFunction((-x**2 - x + 1) + (x**2 + x)*Dx + (-x - 2)*Dx**3 + \
        (x + 1)*Dx**4, x, 2, [4, 1, 2, -5 ])
    assert p + q == r
    p = from_sympy(sin(x))
    q = from_sympy(1/x)
    r = HolonomicFunction((x**2 + 6) + (x**3 + 2*x)*Dx + (x**2 + 6)*Dx**2 + (x**3 + 2*x)*Dx**3, \
        x, 1, [sin(1) + 1, -1 + cos(1), -sin(1) + 2])
    assert p + q == r

def test_multiplication_initial_condition():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx**2 + x*Dx - 1, x, 0, [3, 1])
    q = HolonomicFunction(Dx**2 + 1, x, 0, [1, 1])
    r = HolonomicFunction((x**4 + 14*x**2 + 60) + 4*x*Dx + (x**4 + 9*x**2 + 20)*Dx**2 + \
        (2*x**3 + 18*x)*Dx**3 + (x**2 + 10)*Dx**4, x, 0, [3, 4, 2, 3])
    assert p * q == r
    p = HolonomicFunction(Dx**2 + x, x, 0, [1, 0])
    q = HolonomicFunction(Dx**3 - x**2, x, 0, [3, 3, 3])
    r = HolonomicFunction((x**8 - 37*x**7/27 - 10*x**6/27 - 164*x**5/9 - 184*x**4/9 + \
        160*x**3/27 + 404*x**2/9 + 8*x + 40/3) + (6*x**7 - 128*x**6/9 - 98*x**5/9 - 28*x**4/9 + \
        8*x**3/9 + 28*x**2 + 40*x/9 - 40)*Dx + (3*x**6 - 82*x**5/9 + 76*x**4/9 + 4*x**3/3 + \
        220*x**2/9 - 80*x/3)*Dx**2 + (-2*x**6 + 128*x**5/27 - 2*x**4/3 -80*x**2/9 + 200/9)*Dx**3 + \
        (3*x**5 - 64*x**4/9 - 28*x**3/9 + 6*x**2 - 20*x/9 - 20/3)*Dx**4 + (-4*x**3 + 64*x**2/9 + \
            8*x/3)*Dx**5 + (x**4 - 64*x**3/27 - 4*x**2/3 + 20/9)*Dx**6, x, 0, [3, 3, 3, -3, -12, -24])
    assert p * q == r
    p = HolonomicFunction(Dx - 1, x, 0, [2])
    q = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])
    r = HolonomicFunction(2 -2*Dx + Dx**2, x, 0, [0, 2])
    assert p * q == r
    q = HolonomicFunction(x*Dx**2 + 1 + 2*Dx, x, 0,[0, 1])
    r = HolonomicFunction((x - 1) + (-2*x + 2)*Dx + x*Dx**2, x, 0, [0, 2])
    assert p * q == r
    p = HolonomicFunction(Dx**2 - 1, x, 0, [1, 3])
    q = HolonomicFunction(Dx**3 + 1, x, 0, [1, 2, 1])
    r = HolonomicFunction(6*Dx + 3*Dx**2 + 2*Dx**3 - 3*Dx**4 + Dx**6, x, 0, [1, 5, 14, 17, 17, 2])
    assert p * q == r
    p = from_sympy(sin(x))
    q = from_sympy(1/x)
    r = HolonomicFunction(x + 2*Dx + x*Dx**2, x, 1, [sin(1), -sin(1) + cos(1)])
    assert p * q == r

def test_HolonomicFunction_composition():
    x = symbols('x')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx-1, x).composition(x**2+x)
    r = HolonomicFunction((-2*x - 1) + Dx, x)
    assert p == r
    p = HolonomicFunction(Dx**2+1, x).composition(x**5+x**2+1)
    r = HolonomicFunction((125*x**12 + 150*x**9 + 60*x**6 + 8*x**3) + (-20*x**3 - 2)*Dx + \
        (5*x**4 + 2*x)*Dx**2, x)
    assert p == r
    p = HolonomicFunction(Dx**2*x+x, x).composition(2*x**3+x**2+1)
    r = HolonomicFunction((216*x**9 + 324*x**8 + 180*x**7 + 152*x**6 + 112*x**5 + \
        36*x**4 + 4*x**3) + (24*x**4 + 16*x**3 + 3*x**2 - 6*x - 1)*Dx + (6*x**5 + 5*x**4 + \
        x**3 + 3*x**2 + x)*Dx**2, x)
    assert p == r
    p = HolonomicFunction(Dx**2+1, x).composition(1-x**2)
    r = HolonomicFunction((4*x**3) - Dx + x*Dx**2, x)
    assert p == r
    p = HolonomicFunction(Dx**2+1, x).composition(x - 2/(x**2 + 1))
    r = HolonomicFunction((x**12 + 6*x**10 + 12*x**9 + 15*x**8 + 48*x**7 + 68*x**6 + \
        72*x**5 + 111*x**4 + 112*x**3 + 54*x**2 + 12*x + 1) + (12*x**8 + 32*x**6 + \
        24*x**4 - 4)*Dx + (x**12 + 6*x**10 + 4*x**9 + 15*x**8 + 16*x**7 + 20*x**6 + 24*x**5+ \
        15*x**4 + 16*x**3 + 6*x**2 + 4*x + 1)*Dx**2, x)
    assert p == r

def test_from_hyper():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    p = hyper([1, 1], [S(3)/2], x**2/4)
    q = HolonomicFunction((4*x) + (5*x**2 - 8)*Dx + (x**3 - 4*x)*Dx**2, x, 1, [2*sqrt(3)*pi/9, -4*sqrt(3)*pi/27 + 4/3])
    r = from_hyper(p)
    assert r == q
    p = from_hyper(hyper([1], [S(3)/2], x**2/4))
    q = HolonomicFunction(-x + (-x**2/2 + 2)*Dx + x*Dx**2, x)
    x0 = 1
    y0 = '[sqrt(pi)*exp(1/4)*erf(1/2), -sqrt(pi)*exp(1/4)*erf(1/2)/2 + 1]'
    assert sstr(p.y0) == y0
    assert q.annihilator == p.annihilator

def test_from_meijerg():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    p = from_meijerg(meijerg(([], [S(3)/2]), ([S(1)/2], [S(1)/2, 1]), x))
    q = HolonomicFunction(x/2 - 1/4 + (-x**2 + x/4)*Dx + x**2*Dx**2 + x**3*Dx**3, x, 1, \
        [1/sqrt(pi), 1/(2*sqrt(pi)), -1/(4*sqrt(pi))])
    assert p == q
    p = from_meijerg(meijerg(([], []), ([0], []), x))
    q = HolonomicFunction(1 + Dx, x, 0, 1)
    assert p == q
    p = from_meijerg(meijerg(([1], []), ([S(1)/2], [0]), x))
    q = HolonomicFunction((x + 1/2)*Dx + x*Dx**2, x, 1, [sqrt(pi)*erf(1), exp(-1)])
    assert p == q
    p = from_meijerg(meijerg(([0], [1]), ([0], []), 2*x**2))
    q = HolonomicFunction((3*x**2 - 1)*Dx + x**3*Dx**2, x, 1, [-exp(-S(1)/2) + 1, -exp(-S(1)/2)])
    assert p == q

def test_to_Sequence():
    x = symbols('x')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
    n = symbols('n', integer=True)
    _, Sn = RecurrenceOperators(ZZ.old_poly_ring(n), 'Sn')
    p = HolonomicFunction(x**2*Dx**4 + x + Dx, x).to_sequence()
    q = (HolonomicSequence(1 + (n + 2)*Sn**2 + (n**4 + 6*n**3 + 11*n**2 + 6*n)*Sn**3), 1)
    assert p == q
    p = HolonomicFunction(x**2*Dx**4 + x**3 + Dx**2, x).to_sequence()
    q = (HolonomicSequence(1 + (n**4 + 14*n**3 + 72*n**2 + 163*n + 140)*Sn**5), 0)
    assert p == q
    p = HolonomicFunction(x**3*Dx**4 + 1 + Dx**2, x).to_sequence()
    q = (HolonomicSequence(1 + (n**4 - 2*n**3 - n**2 + 2*n)*Sn + (n**2 + 3*n + 2)*Sn**2), 3)
    assert p == q
    p = HolonomicFunction(3*x**3*Dx**4 + 2*x*Dx + x*Dx**3, x).to_sequence()
    q = (HolonomicSequence(2*n + (3*n**4 - 6*n**3 - 3*n**2 + 6*n)*Sn + (n**3 + 3*n**2 + 2*n)*Sn**2), 3)
    assert p == q

def test_to_Sequence_Initial_Coniditons():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    n = symbols('n', integer=True)
    _, Sn = RecurrenceOperators(QQ.old_poly_ring(n), 'Sn')
    p = HolonomicFunction(Dx - 1, x, 0, 1).to_sequence()
    q = (HolonomicSequence(-1 + (n + 1)*Sn, 1), 0)
    assert p == q
    p = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1]).to_sequence()
    q = (HolonomicSequence(1 + (n**2 + 3*n + 2)*Sn**2, [0, 1]), 0)
    assert p == q
    p = HolonomicFunction(Dx**2 + 1 + x**3*Dx, x, 0, [2, 3]).to_sequence()
    q = (HolonomicSequence(n + Sn**2 + (n**2 + 7*n + 12)*Sn**4, [2, 3, -1, -1/2, 1/12]), 1)
    assert p == q
    p = HolonomicFunction(x**3*Dx**5 + 1 + Dx, x).to_sequence()
    q = (HolonomicSequence(1 + (n + 1)*Sn + (n**5 - 5*n**3 + 4*n)*Sn**2), 3)
    assert p == q
    C_1, C_2, C_3 = symbols('C_1, C_2, C_3')
    p = from_sympy(log(1+x**2))
    q = (HolonomicSequence(n**2 + (n**2 + 2*n)*Sn**2, [0, 0, C_2, 0]), 2)
    assert p.to_sequence() == q
    p = p.diff()
    q = (HolonomicSequence((n + 1) + (n + 1)*Sn**2, [0, C_1, 0]), 1)
    assert p.to_sequence() == q
    p = from_sympy(erf(x) + x).to_sequence()
    q = (HolonomicSequence((2*n**2 - 2*n) + (n**3 + 2*n**2 - n - 2)*Sn**2, [0, 1 + 2/sqrt(pi), 0, C_3]), 2)
    assert p == q

def test_series():
    x = symbols('x')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx**2 + 2*x*Dx, x, 0, [0, 1]).series(n=10)
    q = x - x**3/3 + x**5/10 - x**7/42 + x**9/216 + O(x**10)
    assert p == q
    p = HolonomicFunction(Dx - 1, x).composition(x**2, 0, 1)  # e^(x**2)
    q = HolonomicFunction(Dx**2 + 1, x, 0, [1, 0])  # cos(x)
    r = (p * q).series(n=10)  # expansion of cos(x) * exp(x**2)
    s = 1 + x**2/2 + x**4/24 - 31*x**6/720 - 179*x**8/8064 + O(x**10)
    assert r == s
    t = HolonomicFunction((1 + x)*Dx**2 + Dx, x, 0, [0, 1])  # log(1 + x)
    r = (p * t + q).series(n=10)
    s = 1 + x - x**2 + 4*x**3/3 - 17*x**4/24 + 31*x**5/30 - 481*x**6/720 +\
     71*x**7/105 - 20159*x**8/40320 + 379*x**9/840 + O(x**10)
    assert r == s
    p = HolonomicFunction((6+6*x-3*x**2) - (10*x-3*x**2-3*x**3)*Dx + \
        (4-6*x**3+2*x**4)*Dx**2, x, 0, [0, 1]).series(n=7)
    q = x + x**3/6 - 3*x**4/16 + x**5/20 - 23*x**6/960 + O(x**7)
    assert p == q
    p = HolonomicFunction((6+6*x-3*x**2) - (10*x-3*x**2-3*x**3)*Dx + \
        (4-6*x**3+2*x**4)*Dx**2, x, 0, [1, 0]).series(n=7)
    q = 1 - 3*x**2/4 - x**3/4 - 5*x**4/32 - 3*x**5/40 - 17*x**6/384 + O(x**7)
    assert p == q
    p = from_sympy(erf(x) + x).series(n=10)
    C_3 = symbols('C_3')
    q = (erf(x) + x).series(n=10)
    assert p.subs(C_3, -2/(3*sqrt(pi))) == q

def test_evalf_euler():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')

    # log(1+x)
    p = HolonomicFunction((1 + x)*Dx**2 + Dx, x, 0, [0, 1])

    # path taken is a straight line from 0 to 1, on the real axis
    r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    s = '0.699525841805253'  # approx. equal to log(2) i.e. 0.693147180559945
    assert sstr(p.evalf(r, method='Euler')[-1]) == s

    # path taken is a traingle 0-->1+i-->2
    r = [0.1 + 0.1*I]
    for i in range(9):
        r.append(r[-1]+0.1+0.1*I)
    for i in range(10):
        r.append(r[-1]+0.1-0.1*I)

    # close to the exact solution 1.09861228866811
    # imaginary part also close to zero
    s = '1.07530466271334 - 0.0251200594793912*I'
    assert sstr(p.evalf(r, method='Euler')[-1]) == s

    # sin(x)
    p = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])
    s = '0.905546532085401 - 6.93889390390723e-18*I'
    assert sstr(p.evalf(r, method='Euler')[-1]) == s

    # computing sin(pi/2) using this method
    # using a linear path from 0 to pi/2
    r = [0.1]
    for i in range(14):
        r.append(r[-1] + 0.1)
    r.append(pi/2)
    s = '1.08016557252834' # close to 1.0 (exact solution)
    assert sstr(p.evalf(r, method='Euler')[-1]) == s

    # trying different path, a rectangle (0-->i-->pi/2 + i-->pi/2)
    # computing the same value sin(pi/2) using different path
    r = [0.1*I]
    for i in range(9):
        r.append(r[-1]+0.1*I)
    for i in range(15):
        r.append(r[-1]+0.1)
    r.append(pi/2+I)
    for i in range(10):
        r.append(r[-1]-0.1*I)

    # close to 1.0
    s = '0.976882381836257 - 1.65557671738537e-16*I'
    assert sstr(p.evalf(r, method='Euler')[-1]) == s

    # cos(x)
    p = HolonomicFunction(Dx**2 + 1, x, 0, [1, 0])
    # compute cos(pi) along 0-->pi
    r = [0.05]
    for i in range(61):
        r.append(r[-1]+0.05)
    r.append(pi)
    # close to -1 (exact answer)
    s = '-1.08140824719196'
    assert sstr(p.evalf(r, method='Euler')[-1]) == s

    # a rectangular path (0 -> i -> 2+i -> 2)
    r = [0.1*I]
    for i in range(9):
        r.append(r[-1]+0.1*I)
    for i in range(20):
        r.append(r[-1]+0.1)
    for i in range(10):
        r.append(r[-1]-0.1*I)

    p = HolonomicFunction(Dx**2 + 1, x, 0, [1,1]).evalf(r, method='Euler')
    s = '0.501421652861245 - 3.88578058618805e-16*I'
    assert sstr(p[-1]) == s

def test_evalf_rk4():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')

    # log(1+x)
    p = HolonomicFunction((1 + x)*Dx**2 + Dx, x, 0, [0, 1])

    # path taken is a straight line from 0 to 1, on the real axis
    r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    s = '0.693146363174626'  # approx. equal to log(2) i.e. 0.693147180559945
    assert sstr(p.evalf(r)[-1]) == s

    # path taken is a traingle 0-->1+i-->2
    r = [0.1 + 0.1*I]
    for i in range(9):
        r.append(r[-1]+0.1+0.1*I)
    for i in range(10):
        r.append(r[-1]+0.1-0.1*I)

    # close to the exact solution 1.09861228866811
    # imaginary part also close to zero
    s = '1.09861574485151 + 1.36082967699958e-7*I'
    assert sstr(p.evalf(r)[-1]) == s

    # sin(x)
    p = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])
    s = '0.90929463522785 + 1.52655665885959e-16*I'
    assert sstr(p.evalf(r)[-1]) == s

    # computing sin(pi/2) using this method
    # using a linear path from 0 to pi/2
    r = [0.1]
    for i in range(14):
        r.append(r[-1] + 0.1)
    r.append(pi/2)
    s = '0.999999895088917' # close to 1.0 (exact solution)
    assert sstr(p.evalf(r)[-1]) == s

    # trying different path, a rectangle (0-->i-->pi/2 + i-->pi/2)
    # computing the same value sin(pi/2) using different path
    r = [0.1*I]
    for i in range(9):
        r.append(r[-1]+0.1*I)
    for i in range(15):
        r.append(r[-1]+0.1)
    r.append(pi/2+I)
    for i in range(10):
        r.append(r[-1]-0.1*I)

    # close to 1.0
    s = '1.00000003415141 + 6.11940487991086e-16*I'
    assert sstr(p.evalf(r)[-1]) == s

    # cos(x)
    p = HolonomicFunction(Dx**2 + 1, x, 0, [1, 0])
    # compute cos(pi) along 0-->pi
    r = [0.05]
    for i in range(61):
        r.append(r[-1]+0.05)
    r.append(pi)
    # close to -1 (exact answer)
    s = '-0.999999993238714'
    assert sstr(p.evalf(r)[-1]) == s

    # a rectangular path (0 -> i -> 2+i -> 2)
    r = [0.1*I]
    for i in range(9):
        r.append(r[-1]+0.1*I)
    for i in range(20):
        r.append(r[-1]+0.1)
    for i in range(10):
        r.append(r[-1]-0.1*I)

    p = HolonomicFunction(Dx**2 + 1, x, 0, [1,1]).evalf(r)
    s = '0.493152791638442 - 1.41553435639707e-15*I'
    assert sstr(p[-1]) == s

def test_from_sympy():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    p = from_sympy((sin(x)/x)**2)
    q = HolonomicFunction(8*x + (4*x**2 + 6)*Dx + 6*x*Dx**2 + x**2*Dx**3, x, 0, \
        [1, 0, -2/3])
    assert p == q
    p = from_sympy(1/(1+x**2)**2)
    q = HolonomicFunction(4*x + (x**2 + 1)*Dx, x, 0, 1)
    assert p == q
    p = from_sympy(exp(x)*sin(x)+x*log(1+x))
    q = HolonomicFunction((2*x**3 + 10*x**2 + 20*x + 18) + (-2*x**4 - 10*x**3 - 20*x**2 \
        - 18*x)*Dx + (2*x**5 + 6*x**4 + 7*x**3 + 8*x**2 + 10*x - 4)*Dx**2 + \
        (-2*x**5 - 5*x**4 - 2*x**3 + 2*x**2 - x + 4)*Dx**3 + (x**5 + 2*x**4 - x**3 - \
        7*x**2/2 + x + 5/2)*Dx**4, x, 0, [0, 1, 4, -1])
    assert p == q
    p = from_sympy(x*exp(x)+cos(x)+1)
    q = HolonomicFunction((-x - 3)*Dx + (x + 2)*Dx**2 + (-x - 3)*Dx**3 + (x + 2)*Dx**4, x, \
        0, [2, 1, 1, 3])
    assert p == q
    assert (x*exp(x)+cos(x)+1).series(n=10) == p.series(n=10)
    p = from_sympy(log(1 + x)**2 + 1)
    q = HolonomicFunction(Dx + (3*x + 3)*Dx**2 + (x**2 + 2*x + 1)*Dx**3, x, 0, [1, 0, 2])
    assert p == q
    p = from_sympy(erf(x)**2 + x)
    q = HolonomicFunction((8*x**4 - 2*x**2 + 2)*Dx**2 + (6*x**3 - x/2)*Dx**3 + \
        (x**2+ 1/4)*Dx**4, x, 0, [0, 1, 8/pi, 0])
    assert p == q
    p = from_sympy(cosh(x)*x)
    q = HolonomicFunction((-x**2 + 2) -2*x*Dx + x**2*Dx**2, x, 0, [0, 1])
    assert p == q
    p = from_sympy(besselj(2, x))
    q = HolonomicFunction((x**2 - 4) + x*Dx + x**2*Dx**2, x, 0, [0, 0])
    assert p == q
    p = from_sympy(besselj(0, x) + exp(x))
    q = HolonomicFunction((-x**2 - x/2 + 1/2) + (x**2 - x/2 - 3/2)*Dx + (-x**2 + x/2 + 1)*Dx**2 +\
        (x**2 + x/2)*Dx**3, x, 0, [2, 1, 1/2])
    assert p == q
    p = from_sympy(sin(x)**2/x)
    q = HolonomicFunction(4 + 4*x*Dx + 3*Dx**2 + x*Dx**3, x, 0, [0, 1, 0])
    assert p == q
    p = from_sympy(sin(x)**2/x, x0=2)
    q = HolonomicFunction((4) + (4*x)*Dx + (3)*Dx**2 + (x)*Dx**3, x, 2, [sin(2)**2/2,
        sin(2)*cos(2) - sin(2)**2/4, -3*sin(2)**2/4 + cos(2)**2 - sin(2)*cos(2)])
    assert p == q
    p = from_sympy(log(x)/2 - Ci(2*x)/2 + Ci(2)/2)
    q = HolonomicFunction(4*Dx + 4*x*Dx**2 + 3*Dx**3 + x*Dx**4, x, 0, \
        [-log(2)/2 - EulerGamma/2 + Ci(2)/2, 0, 1, 0])
    assert p == q
    p = p.to_sympy()
    q = log(x)/2 - Ci(2*x)/2 + Ci(2)/2
    assert p == q


def test_to_hyper():
    x = symbols('x')
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx - 2, x, 0, 3).to_hyper()
    q = 3 * hyper([], [], 2*x)
    assert p == q
    p = hyperexpand(HolonomicFunction((1 + x) * Dx - 3, x, 0, 2).to_hyper()).expand()
    q = 2*x**3 + 6*x**2 + 6*x + 2
    assert p == q
    p = HolonomicFunction((1 + x)*Dx**2 + Dx, x, 0, [0, 1]).to_hyper()
    q = -x**2*hyper((2, 2, 1), (2, 3), -x)/2 + x
    assert p == q
    p = HolonomicFunction(2*x*Dx + Dx**2, x, 0, [0, 2/sqrt(pi)]).to_hyper()
    q = 2*x*hyper((1/2,), (3/2,), -x**2)/sqrt(pi)
    assert p == q
    p = hyperexpand(HolonomicFunction(2*x*Dx + Dx**2, x, 0, [1, -2/sqrt(pi)]).to_hyper())
    q = erfc(x)
    assert p.rewrite(erfc) == q
    p =  hyperexpand(HolonomicFunction((x**2 - 1) + x*Dx + x**2*Dx**2,
        x, 0, [0, S(1)/2]).to_hyper())
    q = besselj(1, x)
    assert p == q
    p = hyperexpand(HolonomicFunction(x*Dx**2 + Dx + x, x, 0, [1, 0]).to_hyper())
    q = besselj(0, x)
    assert p == q

def test_to_sympy():
    x = symbols('x')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(Dx - 1, x, 0, 1).to_sympy()
    q = exp(x)
    assert p == q
    p = HolonomicFunction(Dx**2 + 1, x, 0, [1, 0]).to_sympy()
    q = cos(x)
    assert p == q
    p = HolonomicFunction(Dx**2 - 1, x, 0, [1, 0]).to_sympy()
    q = cosh(x)
    assert p == q
    p = HolonomicFunction(2 + (4*x - 1)*Dx + \
        (x**2 - x)*Dx**2, x, 0, [1, 2]).to_sympy()
    q = 1/(x**2 - 2*x + 1)
    assert p == q
    p = from_sympy(sin(x)**2/x).integrate((x, 0, x)).to_sympy()
    q = (sin(x)**2/x).integrate((x, 0, x))
    assert p == q
    C_1, C_2, C_3 = symbols('C_1, C_2, C_3')
    p = from_sympy(log(1+x**2)).to_sympy()
    q = C_2*log(x**2 + 1)
    assert p == q
    p = from_sympy(log(1+x**2)).diff().to_sympy()
    q = C_1*x/(x**2 + 1)
    assert p == q
    p = from_sympy(erf(x) + x).to_sympy()
    q = 3*C_3*x - 3*sqrt(pi)*C_3*erf(x)/2 + x + 2*x/sqrt(pi)
    assert p == q

def test_integrate():
    x = symbols('x')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
    p = from_sympy(sin(x)**2/x, x0=1).integrate((x, 2, 3))
    q = '0.166270406994788'
    assert sstr(p) == q
    p = from_sympy(sin(x)).integrate((x, 0, x)).to_sympy()
    q = 1 - cos(x)
    assert p == q
    p = from_sympy(sin(x)).integrate((x, 0, 3))
    q = '1.98999246812687'
    assert sstr(p) == q
    p = from_sympy(sin(x)/x, x0=1).integrate((x, 1, 2))
    q = '0.659329913368450'
    assert sstr(p) == q
    p = from_sympy(sin(x)**2/x, x0=1).integrate((x, 1, 0))
    q = '-0.423690480850035'
    assert sstr(p) == q

def test_diff():
    x, y = symbols('x, y')
    R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
    p = HolonomicFunction(x*Dx**2 + 1, x, 0, [0, 1])
    assert p.diff().to_sympy() == p.to_sympy().diff().simplify()
    p = HolonomicFunction(Dx**2 - 1, x, 0, [1, 0])
    assert p.diff(x, 2).to_sympy() == p.to_sympy()
    p = from_sympy(Si(x))
    assert p.diff().to_sympy() == sin(x)/x
    assert p.diff(y) == 0
    C_1, C_2, C_3 = symbols('C_1, C_2, C_3')
    q = Si(x)
    assert p.diff(x).to_sympy() == q.diff()
    assert p.diff(x, 2).to_sympy().subs(C_1, -S(1)/3) == q.diff(x, 2).simplify()
    assert p.diff(x, 3).series().subs(C_2, S(1)/10) == q.diff(x, 3).series()

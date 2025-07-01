
from sympy.core.symbol import symbols
from sympy.physics.control.routh_table import (RouthHurwitz,
                                            negative_real_part_conditions)
from sympy.matrices.dense import Matrix
from sympy.polys import Poly, PurePoly
from sympy.logic.boolalg import true, false
from sympy import Symbol, Add, Mul, Integer, Pow

s = symbols('s')

def test_table():
    b0, b1, b2, b3, b4 = symbols('b_0 b_1 b_2 b_3 b_4')
    p1 = b4 * s**4 + b3 * s**3 + b2 * s**2 + b1 * s + b0

    # generic polynomial tests
    t1 = RouthHurwitz(p1, s)
    expected1 = Matrix([
        [b4, b2, b0], [b3, b1, 0], [-b1*b4/b3 + b2, b0, 0],
        [(b0*b3**2 + b1*(b1*b4 - b2*b3))/(b1*b4 - b2*b3), 0, 0],
        [b0, 0, 0]])

    assert t1.equals(expected1)

    expected1_1 = Matrix([
        [b4], [b3], [-b1*b4/b3 + b2],
        [(b0*b3**2 + b1*(b1*b4 - b2*b3))/(b1*b4 - b2*b3)],
        [b0]])

    assert t1[:, 0].equals(expected1_1)
    assert t1.zero_row_case is False
    assert t1.auxiliary_polynomials is None

    # zero in the first column test case
    p2 = s**4 + s**3 + 3*s**2 + 3*s + 3
    t2 = RouthHurwitz(p2, s)

    expected2 = Matrix([
        [1, 3, 3], [1, 3, 0],[-3, 3, 0],
        [4, 0, 0], [3, 0, 0]])

    assert t2.equals(expected2)
    assert t2.zero_row_case is False

    # zero row test case
    p3 = s**6 + 2*s**5 + 8*s**4 + 12*s**3 + 20*s**2 + 16*s + 16
    t3 = RouthHurwitz(p3, s)

    expected3 = Matrix([
        [1, 8, 20, 16], [2, 12, 16, 0], [2, 12, 16, 0],
        [8, 24, 0, 0], [6, 16, 0, 0], [8/3, 0, 0, 0], [16, 0, 0, 0]])

    assert t3.equals(expected3)
    assert t3.zero_row_case is True
    assert t3.auxiliary_polynomials == [Poly(2*s**4 + 12*s**2 + 16, s)]

def test_negative_real_part_conditions():
    b0, b1, b2, b3, b4 = symbols('b_0 b_1 b_2 b_3 b_4')
    p1 = b4 * s**4 + b3 * s**3 + b2 * s**2 + b1 * s + b0

    conds = negative_real_part_conditions(p1, s)
    assert conds == [
        b3*b4 > 0, b3**2*(-b1*b4 + b2*b3) > 0,
        (-b0*b3**3 + b1*b3*(-b1*b4 + b2*b3))*(-b1*b4 + b2*b3)**2 > 0,
        b0*b3*(-b0*b3**3 + b1*b3*(-b1*b4 + b2*b3))**3*(-b1*b4 + b2*b3) > 0]

    p2 = -3*s**2 - 2*s - b0
    assert negative_real_part_conditions(p2, s) == [true, 8 * b0 > 0]

    a = symbols('a', nonpositive = True)

    p4 = b0*s**2 + a*s + 3
    assert negative_real_part_conditions(p4, s) == [a * b0 > 0, false]

    p5 = b0*s**2 + a*s - 3
    assert negative_real_part_conditions(p5, s) == [a * b0 > 0, -3 * a**3 > 0]

    p6 = b0 + b1*s**2 + b1*s + b3*s**4 + b3*s**3
    expected6 = [b3**2 > 0, false]

    assert negative_real_part_conditions(p6, s) == expected6

    # test for issue https://github.com/sympy/sympy/issues/28010
    # In that test we want to be sure that negative_real_conditions works fast
    I_L = Symbol('I_L', nonnegative=True, real=True)
    I_T = Symbol('I_T', nonnegative=True, real=True)
    d_L = Symbol('d_L', nonnegative=True, real=True)
    d_T = Symbol('d_T', nonnegative=True, real=True)
    l_L = Symbol('l_L', nonnegative=True, real=True)
    m_L = Symbol('m_L', nonnegative=True, real=True)
    m_T = Symbol('m_T', nonnegative=True, real=True)
    g = Symbol('g', nonnegative=True, real=True)
    k_00, k_01, k_02, k_03 = symbols('k_00 k_01 k_02 k_03')
    k_10, k_11, k_12, k_13 = symbols('k_10 k_11 k_12 k_13')
    s_00, s_01, s_02, s_03 = symbols('s_00 s_01 s_02 s_03')
    s_10, s_11, s_12, s_13 = symbols('s_10 s_11 s_12 s_13')

    common_denom = (I_L*I_T + I_L*d_T**2*m_T + I_T*d_L**2*m_L +
                    I_T*l_L**2*m_T + d_L**2*d_T**2*m_L*m_T)

    p7 = PurePoly(
        s**4 + s**3/common_denom*(
            I_L*k_13*s_13 + I_T*k_02*s_02- I_T*k_03*s_03 - I_T*k_12*s_12 +
            I_T*k_13*s_13 + d_L**2*k_13*m_L*s_13 + d_T**2*k_02*m_T*s_02 -
            d_T**2*k_03*m_T*s_03 - d_T**2*k_12*m_T*s_12 + d_T**2*k_13*m_T*s_13 -
            d_T*k_03*l_L*m_T*s_03 - d_T*k_12*l_L*m_T*s_12 +
            2*d_T*k_13*l_L*m_T*s_13 + k_13*l_L**2*m_T*s_13) +
        s**2/common_denom*(
            I_L*d_T*g*m_T - I_L*k_11*s_11 + I_T*d_L*g*m_L + I_T*g*l_L*m_T -
            I_T*k_00*s_00 + I_T*k_01*s_01 + I_T*k_10*s_10 - I_T*k_11*s_11 +
            d_L**2*d_T*g*m_L*m_T - d_L**2*k_11*m_L*s_11 + d_L*d_T**2*g*m_L*m_T +
            d_T**2*g*l_L*m_T**2 - d_T**2*k_00*m_T*s_00 + d_T**2*k_01*m_T*s_01 +
            d_T**2*k_10*m_T*s_10 - d_T**2*k_11*m_T*s_11 + d_T*g*l_L**2*m_T**2 +
            d_T*k_01*l_L*m_T*s_01 + d_T*k_10*l_L*m_T*s_10 -
            2*d_T*k_11*l_L*m_T*s_11 + k_02*k_13*s_02*s_13 - k_03*k_12*s_03*s_12
            - k_11*l_L**2*m_T*s_11) +
        s/common_denom*(
            d_L*g*k_13*m_L*s_13 +d_T*g*k_02*m_T*s_02 - d_T*g*k_03*m_T*s_03 -
            d_T*g*k_12*m_T*s_12 + d_T*g*k_13*m_T*s_13 + g*k_13*l_L*m_T*s_13 -
            k_00*k_13*s_00*s_13 + k_01*k_12*s_01*s_12 - k_02*k_11*s_02*s_11 +
            k_03*k_10*s_03*s_10) +
            1/common_denom*(d_L*d_T*g**2*m_L*m_T - d_L*g*k_11*m_L*s_11 +
            d_T*g**2*l_L*m_T**2 - d_T*g*k_00*m_T*s_00 + d_T*g*k_01*m_T*s_01 +
            d_T*g*k_10*m_T*s_10 - d_T*g*k_11*m_T*s_11 - g*k_11*l_L*m_T*s_11 +
            k_00*k_11*s_00*s_11 - k_01*k_10*s_01*s_10), s)

    negative_real_part_conditions(p7, s)

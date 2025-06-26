
from sympy.core.symbol import symbols
from sympy.physics.control.routh_table import (RouthHurwitz,
                                            negative_real_root_conditions)
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

def test_negative_real_root_conditions():
    b0, b1, b2, b3, b4 = symbols('b_0 b_1 b_2 b_3 b_4')
    p1 = b4 * s**4 + b3 * s**3 + b2 * s**2 + b1 * s + b0

    conds = negative_real_root_conditions(p1, s)
    assert conds == [
        b3*b4 > 0, b3**2*(-b1*b4 + b2*b3) > 0,
        (-b0*b3**3 + b1*b3*(-b1*b4 + b2*b3))*(-b1*b4 + b2*b3)**2 > 0,
        b0*b3*(-b0*b3**3 + b1*b3*(-b1*b4 + b2*b3))**3*(-b1*b4 + b2*b3) > 0]

    p2 = -3*s**2 - 2*s - b0
    assert negative_real_root_conditions(p2, s) == [true, 8 * b0 > 0]

    a = symbols('a', nonpositive = True)

    p4 = b0*s**2 + a*s + 3
    assert negative_real_root_conditions(p4, s) == [a * b0 > 0, false]

    p5 = b0*s**2 + a*s - 3
    assert negative_real_root_conditions(p5, s) == [a * b0 > 0, -3 * a**3 > 0]

    p6 = b0 + b1*s**2 + b1*s + b3*s**4 + b3*s**3
    expected6 = [b3**2 > 0, false]

    assert negative_real_root_conditions(p6, s) == expected6

    # test for issue https://github.com/sympy/sympy/issues/28010
    # In that test we want to be sure that negative_real_conditions works fast
    p7 = PurePoly(Add(Pow(Symbol('s'), Integer(4)), Mul(Pow(Symbol('s'),
                Integer(3)),
                Pow(Add(Mul(Symbol('I_L', nonnegative=True, real=True), Symbol('I_T',
                nonnegative=True, real=True)), Mul(Symbol('I_L', nonnegative=True, real=True),
                Pow(Symbol('d_T', nonnegative=True, real=True), Integer(2)), Symbol('m_T',
                nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True, real=True),
                Pow(Symbol('d_L', nonnegative=True, real=True), Integer(2)), Symbol('m_L',
                nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True, real=True),
                Pow(Symbol('l_L', nonnegative=True, real=True), Integer(2)), Symbol('m_T',
                nonnegative=True, real=True)), Mul(Pow(Symbol('d_L', nonnegative=True,
                real=True), Integer(2)), Pow(Symbol('d_T', nonnegative=True, real=True),
                Integer(2)), Symbol('m_L', nonnegative=True, real=True), Symbol('m_T',
                nonnegative=True, real=True))), Integer(-1)), Add(Mul(Symbol('I_L',
                nonnegative=True, real=True), Symbol('k_13'), Symbol('s_13')), Mul(Symbol('I_T',
                nonnegative=True, real=True), Symbol('k_02'), Symbol('s_02')), Mul(Integer(-1),
                Symbol('I_T', nonnegative=True, real=True), Symbol('k_03'), Symbol('s_03')),
                Mul(Integer(-1), Symbol('I_T', nonnegative=True, real=True), Symbol('k_12'),
                Symbol('s_12')), Mul(Symbol('I_T', nonnegative=True, real=True), Symbol('k_13'),
                Symbol('s_13')), Mul(Pow(Symbol('d_L', nonnegative=True, real=True),
                Integer(2)), Symbol('k_13'), Symbol('m_L', nonnegative=True, real=True),
                Symbol('s_13')), Mul(Pow(Symbol('d_T', nonnegative=True, real=True),
                Integer(2)), Symbol('k_02'), Symbol('m_T', nonnegative=True, real=True),
                Symbol('s_02')), Mul(Integer(-1), Pow(Symbol('d_T', nonnegative=True,
                real=True), Integer(2)), Symbol('k_03'), Symbol('m_T', nonnegative=True,
                real=True), Symbol('s_03')), Mul(Integer(-1), Pow(Symbol('d_T',
                nonnegative=True, real=True), Integer(2)), Symbol('k_12'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_12')), Mul(Pow(Symbol('d_T',
                nonnegative=True, real=True), Integer(2)), Symbol('k_13'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_13')), Mul(Integer(-1), Symbol('d_T',
                nonnegative=True, real=True), Symbol('k_03'), Symbol('l_L', nonnegative=True,
                real=True), Symbol('m_T', nonnegative=True, real=True), Symbol('s_03')),
                Mul(Integer(-1), Symbol('d_T', nonnegative=True, real=True), Symbol('k_12'),
                Symbol('l_L', nonnegative=True, real=True), Symbol('m_T', nonnegative=True,
                real=True), Symbol('s_12')), Mul(Integer(2), Symbol('d_T', nonnegative=True,
                real=True), Symbol('k_13'), Symbol('l_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True), Symbol('s_13')), Mul(Symbol('k_13'),
                Pow(Symbol('l_L', nonnegative=True, real=True), Integer(2)), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_13')))), Mul(Pow(Symbol('s'),
                Integer(2)), Pow(Add(Mul(Symbol('I_L', nonnegative=True, real=True),
                Symbol('I_T', nonnegative=True, real=True)), Mul(Symbol('I_L', nonnegative=True,
                real=True), Pow(Symbol('d_T', nonnegative=True, real=True), Integer(2)),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True,
                real=True), Pow(Symbol('d_L', nonnegative=True, real=True), Integer(2)),
                Symbol('m_L', nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True,
                real=True), Pow(Symbol('l_L', nonnegative=True, real=True), Integer(2)),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Pow(Symbol('d_L',
                nonnegative=True, real=True), Integer(2)), Pow(Symbol('d_T', nonnegative=True,
                real=True), Integer(2)), Symbol('m_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True))), Integer(-1)),
                Add(Mul(Symbol('I_L', nonnegative=True, real=True), Symbol('d_T',
                nonnegative=True, real=True), Symbol('g', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Integer(-1), Symbol('I_L',
                nonnegative=True, real=True), Symbol('k_11'), Symbol('s_11')), Mul(Symbol('I_T',
                nonnegative=True, real=True), Symbol('d_L', nonnegative=True, real=True),
                Symbol('g', nonnegative=True, real=True), Symbol('m_L', nonnegative=True,
                real=True)), Mul(Symbol('I_T', nonnegative=True, real=True), Symbol('g',
                nonnegative=True, real=True), Symbol('l_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Integer(-1), Symbol('I_T',
                nonnegative=True, real=True), Symbol('k_00'), Symbol('s_00')), Mul(Symbol('I_T',
                nonnegative=True, real=True), Symbol('k_01'), Symbol('s_01')), Mul(Symbol('I_T',
                nonnegative=True, real=True), Symbol('k_10'), Symbol('s_10')), Mul(Integer(-1),
                Symbol('I_T', nonnegative=True, real=True), Symbol('k_11'), Symbol('s_11')),
                Mul(Pow(Symbol('d_L', nonnegative=True, real=True), Integer(2)), Symbol('d_T',
                nonnegative=True, real=True), Symbol('g', nonnegative=True, real=True),
                Symbol('m_L', nonnegative=True, real=True), Symbol('m_T', nonnegative=True,
                real=True)), Mul(Integer(-1), Pow(Symbol('d_L', nonnegative=True, real=True),
                Integer(2)), Symbol('k_11'), Symbol('m_L', nonnegative=True, real=True),
                Symbol('s_11')), Mul(Symbol('d_L', nonnegative=True, real=True),
                Pow(Symbol('d_T', nonnegative=True, real=True), Integer(2)), Symbol('g',
                nonnegative=True, real=True), Symbol('m_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Pow(Symbol('d_T',
                nonnegative=True, real=True), Integer(2)), Symbol('g', nonnegative=True,
                real=True), Symbol('l_L', nonnegative=True, real=True), Pow(Symbol('m_T',
                nonnegative=True, real=True), Integer(2))), Mul(Integer(-1), Pow(Symbol('d_T',
                nonnegative=True, real=True), Integer(2)), Symbol('k_00'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_00')), Mul(Pow(Symbol('d_T',
                nonnegative=True, real=True), Integer(2)), Symbol('k_01'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_01')), Mul(Pow(Symbol('d_T',
                nonnegative=True, real=True), Integer(2)), Symbol('k_10'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_10')), Mul(Integer(-1),
                Pow(Symbol('d_T', nonnegative=True, real=True), Integer(2)), Symbol('k_11'),
                Symbol('m_T', nonnegative=True, real=True), Symbol('s_11')), Mul(Symbol('d_T',
                nonnegative=True, real=True), Symbol('g', nonnegative=True, real=True),
                Pow(Symbol('l_L', nonnegative=True, real=True), Integer(2)), Pow(Symbol('m_T',
                nonnegative=True, real=True), Integer(2))), Mul(Symbol('d_T', nonnegative=True,
                real=True), Symbol('k_01'), Symbol('l_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True), Symbol('s_01')), Mul(Symbol('d_T',
                nonnegative=True, real=True), Symbol('k_10'), Symbol('l_L', nonnegative=True,
                real=True), Symbol('m_T', nonnegative=True, real=True), Symbol('s_10')),
                Mul(Integer(-1), Integer(2), Symbol('d_T', nonnegative=True, real=True),
                Symbol('k_11'), Symbol('l_L', nonnegative=True, real=True), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_11')), Mul(Symbol('k_02'),
                Symbol('k_13'), Symbol('s_02'), Symbol('s_13')), Mul(Integer(-1),
                Symbol('k_03'), Symbol('k_12'), Symbol('s_03'), Symbol('s_12')),
                Mul(Integer(-1), Symbol('k_11'), Pow(Symbol('l_L', nonnegative=True, real=True),
                Integer(2)), Symbol('m_T', nonnegative=True, real=True), Symbol('s_11')))),
                Mul(Symbol('s'), Pow(Add(Mul(Symbol('I_L', nonnegative=True, real=True),
                Symbol('I_T', nonnegative=True, real=True)), Mul(Symbol('I_L', nonnegative=True,
                real=True), Pow(Symbol('d_T', nonnegative=True, real=True), Integer(2)),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True,
                real=True), Pow(Symbol('d_L', nonnegative=True, real=True), Integer(2)),
                Symbol('m_L', nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True,
                real=True), Pow(Symbol('l_L', nonnegative=True, real=True), Integer(2)),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Pow(Symbol('d_L',
                nonnegative=True, real=True), Integer(2)), Pow(Symbol('d_T', nonnegative=True,
                real=True), Integer(2)), Symbol('m_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True))), Integer(-1)),
                Add(Mul(Symbol('d_L', nonnegative=True, real=True), Symbol('g',
                nonnegative=True, real=True), Symbol('k_13'), Symbol('m_L', nonnegative=True,
                real=True), Symbol('s_13')), Mul(Symbol('d_T', nonnegative=True, real=True),
                Symbol('g', nonnegative=True, real=True), Symbol('k_02'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_02')), Mul(Integer(-1), Symbol('d_T',
                nonnegative=True, real=True), Symbol('g', nonnegative=True, real=True),
                Symbol('k_03'), Symbol('m_T', nonnegative=True, real=True), Symbol('s_03')),
                Mul(Integer(-1), Symbol('d_T', nonnegative=True, real=True), Symbol('g',
                nonnegative=True, real=True), Symbol('k_12'), Symbol('m_T', nonnegative=True,
                real=True), Symbol('s_12')), Mul(Symbol('d_T', nonnegative=True, real=True),
                Symbol('g', nonnegative=True, real=True), Symbol('k_13'), Symbol('m_T',
                nonnegative=True, real=True), Symbol('s_13')), Mul(Symbol('g', nonnegative=True,
                real=True), Symbol('k_13'), Symbol('l_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True), Symbol('s_13')), Mul(Integer(-1),
                Symbol('k_00'), Symbol('k_13'), Symbol('s_00'), Symbol('s_13')),
                Mul(Symbol('k_01'), Symbol('k_12'), Symbol('s_01'), Symbol('s_12')),
                Mul(Integer(-1), Symbol('k_02'), Symbol('k_11'), Symbol('s_02'),
                Symbol('s_11')), Mul(Symbol('k_03'), Symbol('k_10'), Symbol('s_03'),
                Symbol('s_10')))), Mul(Pow(Add(Mul(Symbol('I_L', nonnegative=True, real=True),
                Symbol('I_T', nonnegative=True, real=True)), Mul(Symbol('I_L', nonnegative=True,
                real=True), Pow(Symbol('d_T', nonnegative=True, real=True), Integer(2)),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True,
                real=True), Pow(Symbol('d_L', nonnegative=True, real=True), Integer(2)),
                Symbol('m_L', nonnegative=True, real=True)), Mul(Symbol('I_T', nonnegative=True,
                real=True), Pow(Symbol('l_L', nonnegative=True, real=True), Integer(2)),
                Symbol('m_T', nonnegative=True, real=True)), Mul(Pow(Symbol('d_L',
                nonnegative=True, real=True), Integer(2)), Pow(Symbol('d_T', nonnegative=True,
                real=True), Integer(2)), Symbol('m_L', nonnegative=True, real=True),
                Symbol('m_T', nonnegative=True, real=True))), Integer(-1)),
                Add(Mul(Symbol('d_L', nonnegative=True, real=True), Symbol('d_T',
                nonnegative=True, real=True), Pow(Symbol('g', nonnegative=True, real=True),
                Integer(2)), Symbol('m_L', nonnegative=True, real=True), Symbol('m_T',
                nonnegative=True, real=True)), Mul(Integer(-1), Symbol('d_L', nonnegative=True,
                real=True), Symbol('g', nonnegative=True, real=True), Symbol('k_11'),
                Symbol('m_L', nonnegative=True, real=True), Symbol('s_11')), Mul(Symbol('d_T',
                nonnegative=True, real=True), Pow(Symbol('g', nonnegative=True, real=True),
                Integer(2)), Symbol('l_L', nonnegative=True, real=True), Pow(Symbol('m_T',
                nonnegative=True, real=True), Integer(2))), Mul(Integer(-1), Symbol('d_T',
                nonnegative=True, real=True), Symbol('g', nonnegative=True, real=True),
                Symbol('k_00'), Symbol('m_T', nonnegative=True, real=True), Symbol('s_00')),
                Mul(Symbol('d_T', nonnegative=True, real=True), Symbol('g', nonnegative=True,
                real=True), Symbol('k_01'), Symbol('m_T', nonnegative=True, real=True),
                Symbol('s_01')), Mul(Symbol('d_T', nonnegative=True, real=True), Symbol('g',
                nonnegative=True, real=True), Symbol('k_10'), Symbol('m_T', nonnegative=True,
                real=True), Symbol('s_10')), Mul(Integer(-1), Symbol('d_T', nonnegative=True,
                real=True), Symbol('g', nonnegative=True, real=True), Symbol('k_11'),
                Symbol('m_T', nonnegative=True, real=True), Symbol('s_11')), Mul(Integer(-1),
                Symbol('g', nonnegative=True, real=True), Symbol('k_11'), Symbol('l_L',
                nonnegative=True, real=True), Symbol('m_T', nonnegative=True, real=True),
                Symbol('s_11')), Mul(Symbol('k_00'), Symbol('k_11'), Symbol('s_00'),
                Symbol('s_11')), Mul(Integer(-1), Symbol('k_01'), Symbol('k_10'),
                Symbol('s_01'), Symbol('s_10'))))), Symbol('s'))

    negative_real_root_conditions(p7, s)
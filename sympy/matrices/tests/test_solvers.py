import random

from sympy import (
    Abs, Add, E, Float, I, Integer, Max, Min, N, Poly, Pow, PurePoly, Rational,
    S, Symbol, cos, exp, log, expand_mul, oo, pi, signsimp, simplify, sin,
    sqrt, symbols, sympify, trigsimp, tan, sstr, diff, Function, expand)
from sympy.matrices.matrices import (ShapeError, MatrixError,
    NonSquareMatrixError, DeferredVector, _find_reasonable_pivot_naive,
    _simplify)
from sympy.matrices import (
    GramSchmidt, ImmutableMatrix, ImmutableSparseMatrix, Matrix,
    SparseMatrix, casoratian, diag, eye, hessian,
    matrix_multiply_elementwise, ones, randMatrix, rot_axis1, rot_axis2,
    rot_axis3, wronskian, zeros, MutableDenseMatrix, ImmutableDenseMatrix, MatrixSymbol)
from sympy.core.compatibility import iterable, Hashable
from sympy.core import Tuple, Wild
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.utilities.iterables import flatten, capture
from sympy.testing.pytest import raises, XFAIL, skip, warns_deprecated_sympy
from sympy.solvers import solve
from sympy.assumptions import Q
from sympy.tensor.array import Array
from sympy.matrices.expressions import MatPow

from sympy.abc import a, b, c, d, x, y, z, t

# don't re-order this list
classes = (Matrix, SparseMatrix, ImmutableMatrix, ImmutableSparseMatrix)

def test_issue_17247_expression_blowup_29():
    M = Matrix(S('''[
        [             -3/4,       45/32 - 37*I/16,                   0,                     0],
        [-149/64 + 49*I/32, -177/128 - 1369*I/128,                   0, -2063/256 + 541*I/128],
        [                0,         9/4 + 55*I/16, 2473/256 + 137*I/64,                     0],
        [                0,                     0,                   0, -177/128 - 1369*I/128]]'''))
    assert M.gauss_jordan_solve(ones(4, 1)) == (Matrix(S('''[
        [                          -32549314808672/3306971225785 - 17397006745216*I/3306971225785],
        [                               67439348256/3306971225785 - 9167503335872*I/3306971225785],
        [-15091965363354518272/21217636514687010905 + 16890163109293858304*I/21217636514687010905],
        [                                                          -11328/952745 + 87616*I/952745]]''')), Matrix(0, 1, []))

@XFAIL # dotprodsimp is not on by default in this function
def test_issue_17247_expression_blowup_30():
    M = Matrix(S('''[
        [             -3/4,       45/32 - 37*I/16,                   0,                     0],
        [-149/64 + 49*I/32, -177/128 - 1369*I/128,                   0, -2063/256 + 541*I/128],
        [                0,         9/4 + 55*I/16, 2473/256 + 137*I/64,                     0],
        [                0,                     0,                   0, -177/128 - 1369*I/128]]'''))
    assert M.cholesky_solve(ones(4, 1)) == Matrix(S('''[
        [                          -32549314808672/3306971225785 - 17397006745216*I/3306971225785],
        [                               67439348256/3306971225785 - 9167503335872*I/3306971225785],
        [-15091965363354518272/21217636514687010905 + 16890163109293858304*I/21217636514687010905],
        [                                                          -11328/952745 + 87616*I/952745]]'''))

# This test is commented out because without dotprodsimp this calculation hangs.
# @XFAIL # dotprodsimp is not on by default in this function
# def test_issue_17247_expression_blowup_31():
#     M = Matrix([
#         [x + 1, 1 - x,     0,     0],
#         [1 - x, x + 1,     0, x + 1],
#         [    0, 1 - x, x + 1,     0],
#         [    0,     0,     0, x + 1]])
#     assert M.LDLsolve(ones(4, 1)) == Matrix([
#         [(x + 1)/(4*x)],
#         [(x - 1)/(4*x)],
#         [(x + 1)/(4*x)],
#         [    1/(x + 1)]])

@XFAIL # dotprodsimp is not on by default in this function
def test_issue_17247_expression_blowup_32():
    M = Matrix([
        [x + 1, 1 - x,     0,     0],
        [1 - x, x + 1,     0, x + 1],
        [    0, 1 - x, x + 1,     0],
        [    0,     0,     0, x + 1]])
    assert M.LUsolve(ones(4, 1)) == Matrix([
        [(x + 1)/(4*x)],
        [(x - 1)/(4*x)],
        [(x + 1)/(4*x)],
        [    1/(x + 1)]])

def test_issue_18531():
    # solve_linear_system still needs fixing but the rref works.
    M = Matrix([
        [1, 1, 1, 1, 1, 0, 1, 0, 0],
        [1 + sqrt(2), -1 + sqrt(2), 1 - sqrt(2), -sqrt(2) - 1, 1, 1, -1, 1, 1],
        [-5 + 2*sqrt(2), -5 - 2*sqrt(2), -5 - 2*sqrt(2), -5 + 2*sqrt(2), -7, 2, -7, -2, 0],
        [-3*sqrt(2) - 1, 1 - 3*sqrt(2), -1 + 3*sqrt(2), 1 + 3*sqrt(2), -7, -5, 7, -5, 3],
        [7 - 4*sqrt(2), 4*sqrt(2) + 7, 4*sqrt(2) + 7, 7 - 4*sqrt(2), 7, -12, 7, 12, 0],
        [-1 + 3*sqrt(2), 1 + 3*sqrt(2), -3*sqrt(2) - 1, 1 - 3*sqrt(2), 7, -5, -7, -5, 3],
        [-3 + 2*sqrt(2), -3 - 2*sqrt(2), -3 - 2*sqrt(2), -3 + 2*sqrt(2), -1, 2, -1, -2, 0],
        [1 - sqrt(2), -sqrt(2) - 1, 1 + sqrt(2), -1 + sqrt(2), -1, 1, 1, 1, 1]
        ])
    assert M.rref() == (Matrix([
        [1, 0, 0, 0, 0, 0, 0, 0,  1/2],
        [0, 1, 0, 0, 0, 0, 0, 0, -1/2],
        [0, 0, 1, 0, 0, 0, 0, 0,  1/2],
        [0, 0, 0, 1, 0, 0, 0, 0, -1/2],
        [0, 0, 0, 0, 1, 0, 0, 0,    0],
        [0, 0, 0, 0, 0, 1, 0, 0, -1/2],
        [0, 0, 0, 0, 0, 0, 1, 0,    0],
        [0, 0, 0, 0, 0, 0, 0, 1, -1/2]]), (0, 1, 2, 3, 4, 5, 6, 7))

def test_LUsolve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = Matrix([[2, 1], [1, 0], [1, 0]])   # issue 14548
    b = Matrix([3, 1, 1])
    assert A.LUsolve(b) == Matrix([1, 1])
    b = Matrix([3, 1, 2])                  # inconsistent
    raises(ValueError, lambda: A.LUsolve(b))
    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4],
                [2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix([2, 1, -4])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = Matrix([[0, -1, 2], [5, 10, 7]])  # underdetermined
    x = Matrix([-1, 2, 0])
    b = A*x
    raises(NotImplementedError, lambda: A.LUsolve(b))

    A = Matrix(4, 4, lambda i, j: 1/(i+j+1) if i != 3 else 0)
    b = Matrix.zeros(4, 1)
    raises(NotImplementedError, lambda: A.LUsolve(b))

def test_QRsolve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x
    x = Matrix([[1, 2], [3, 4], [5, 6]])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x

    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x
    x = Matrix([[7, 8], [9, 10], [11, 12]])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x

def test_cholesky_solve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x
    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x
    A = Matrix(((1, 5), (5, 1)))
    x = Matrix((4, -3))
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x
    A = Matrix(((9, 3*I), (-3*I, 5)))
    x = Matrix((-2, 1))
    b = A*x
    soln = A.cholesky_solve(b)
    assert expand_mul(soln) == x
    A = Matrix(((9*I, 3), (-3 + I, 5)))
    x = Matrix((2 + 3*I, -1))
    b = A*x
    soln = A.cholesky_solve(b)
    assert expand_mul(soln) == x
    a00, a01, a11, b0, b1 = symbols('a00, a01, a11, b0, b1')
    A = Matrix(((a00, a01), (a01, a11)))
    b = Matrix((b0, b1))
    x = A.cholesky_solve(b)
    assert simplify(A*x) == b

def test_LDLsolve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.LDLsolve(b)
    assert soln == x

    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.LDLsolve(b)
    assert soln == x

    A = Matrix(((9, 3*I), (-3*I, 5)))
    x = Matrix((-2, 1))
    b = A*x
    soln = A.LDLsolve(b)
    assert expand_mul(soln) == x

    A = Matrix(((9*I, 3), (-3 + I, 5)))
    x = Matrix((2 + 3*I, -1))
    b = A*x
    soln = A.LDLsolve(b)
    assert expand_mul(soln) == x

    A = Matrix(((9, 3), (3, 9)))
    x = Matrix((1, 1))
    b = A * x
    soln = A.LDLsolve(b)
    assert expand_mul(soln) == x

    A = Matrix([[-5, -3, -4], [-3, -7, 7]])
    x = Matrix([[8], [7], [-2]])
    b = A * x
    raises(NotImplementedError, lambda: A.LDLsolve(b))

def test_lower_triangular_solve():

    raises(NonSquareMatrixError,
        lambda: Matrix([1, 0]).lower_triangular_solve(Matrix([0, 1])))
    raises(ShapeError,
        lambda: Matrix([[1, 0], [0, 1]]).lower_triangular_solve(Matrix([1])))
    raises(ValueError,
        lambda: Matrix([[2, 1], [1, 2]]).lower_triangular_solve(
            Matrix([[1, 0], [0, 1]])))

    A = Matrix([[1, 0], [0, 1]])
    B = Matrix([[x, y], [y, x]])
    C = Matrix([[4, 8], [2, 9]])

    assert A.lower_triangular_solve(B) == B
    assert A.lower_triangular_solve(C) == C

def test_upper_triangular_solve():

    raises(NonSquareMatrixError,
        lambda: Matrix([1, 0]).upper_triangular_solve(Matrix([0, 1])))
    raises(ShapeError,
        lambda: Matrix([[1, 0], [0, 1]]).upper_triangular_solve(Matrix([1])))
    raises(TypeError,
        lambda: Matrix([[2, 1], [1, 2]]).upper_triangular_solve(
            Matrix([[1, 0], [0, 1]])))

    A = Matrix([[1, 0], [0, 1]])
    B = Matrix([[x, y], [y, x]])
    C = Matrix([[2, 4], [3, 8]])

    assert A.upper_triangular_solve(B) == B
    assert A.upper_triangular_solve(C) == C


def test_diagonal_solve():
    raises(TypeError, lambda: Matrix([1, 1]).diagonal_solve(Matrix([1])))
    A = Matrix([[1, 0], [0, 1]])*2
    B = Matrix([[x, y], [y, x]])
    assert A.diagonal_solve(B) == B/2

    A = Matrix([[1, 0], [1, 2]])
    raises(TypeError, lambda: A.diagonal_solve(B))

def test_pinv_solve():
    # Fully determined system (unique result, identical to other solvers).
    A = Matrix([[1, 5], [7, 9]])
    B = Matrix([12, 13])
    assert A.pinv_solve(B) == A.cholesky_solve(B)
    assert A.pinv_solve(B) == A.LDLsolve(B)
    assert A.pinv_solve(B) == Matrix([sympify('-43/26'), sympify('71/26')])
    assert A * A.pinv() * B == B
    # Fully determined, with two-dimensional B matrix.
    B = Matrix([[12, 13, 14], [15, 16, 17]])
    assert A.pinv_solve(B) == A.cholesky_solve(B)
    assert A.pinv_solve(B) == A.LDLsolve(B)
    assert A.pinv_solve(B) == Matrix([[-33, -37, -41], [69, 75, 81]]) / 26
    assert A * A.pinv() * B == B
    # Underdetermined system (infinite results).
    A = Matrix([[1, 0, 1], [0, 1, 1]])
    B = Matrix([5, 7])
    solution = A.pinv_solve(B)
    w = {}
    for s in solution.atoms(Symbol):
        # Extract dummy symbols used in the solution.
        w[s.name] = s
    assert solution == Matrix([[w['w0_0']/3 + w['w1_0']/3 - w['w2_0']/3 + 1],
                               [w['w0_0']/3 + w['w1_0']/3 - w['w2_0']/3 + 3],
                               [-w['w0_0']/3 - w['w1_0']/3 + w['w2_0']/3 + 4]])
    assert A * A.pinv() * B == B
    # Overdetermined system (least squares results).
    A = Matrix([[1, 0], [0, 0], [0, 1]])
    B = Matrix([3, 2, 1])
    assert A.pinv_solve(B) == Matrix([3, 1])
    # Proof the solution is not exact.
    assert A * A.pinv() * B != B

def test_gauss_jordan_solve():

    # Square, full rank, unique solution
    A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
    b = Matrix([3, 6, 9])
    sol, params = A.gauss_jordan_solve(b)
    assert sol == Matrix([[-1], [2], [0]])
    assert params == Matrix(0, 1, [])

    # Square, full rank, unique solution, B has more columns than rows
    A = eye(3)
    B = Matrix([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    sol, params = A.gauss_jordan_solve(B)
    assert sol == B
    assert params == Matrix(0, 4, [])

    # Square, reduced rank, parametrized solution
    A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    b = Matrix([3, 6, 9])
    sol, params, freevar = A.gauss_jordan_solve(b, freevar=True)
    w = {}
    for s in sol.atoms(Symbol):
        # Extract dummy symbols used in the solution.
        w[s.name] = s
    assert sol == Matrix([[w['tau0'] - 1], [-2*w['tau0'] + 2], [w['tau0']]])
    assert params == Matrix([[w['tau0']]])
    assert freevar == [2]

    # Square, reduced rank, parametrized solution, B has two columns
    A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    B = Matrix([[3, 4], [6, 8], [9, 12]])
    sol, params, freevar = A.gauss_jordan_solve(B, freevar=True)
    w = {}
    for s in sol.atoms(Symbol):
        # Extract dummy symbols used in the solution.
        w[s.name] = s
    assert sol == Matrix([[w['tau0'] - 1, w['tau1'] - Rational(4, 3)],
                          [-2*w['tau0'] + 2, -2*w['tau1'] + Rational(8, 3)],
                          [w['tau0'], w['tau1']],])
    assert params == Matrix([[w['tau0'], w['tau1']]])
    assert freevar == [2]

    # Square, reduced rank, parametrized solution
    A = Matrix([[1, 2, 3], [2, 4, 6], [3, 6, 9]])
    b = Matrix([0, 0, 0])
    sol, params = A.gauss_jordan_solve(b)
    w = {}
    for s in sol.atoms(Symbol):
        w[s.name] = s
    assert sol == Matrix([[-2*w['tau0'] - 3*w['tau1']],
                         [w['tau0']], [w['tau1']]])
    assert params == Matrix([[w['tau0']], [w['tau1']]])

    # Square, reduced rank, parametrized solution
    A = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    b = Matrix([0, 0, 0])
    sol, params = A.gauss_jordan_solve(b)
    w = {}
    for s in sol.atoms(Symbol):
        w[s.name] = s
    assert sol == Matrix([[w['tau0']], [w['tau1']], [w['tau2']]])
    assert params == Matrix([[w['tau0']], [w['tau1']], [w['tau2']]])

    # Square, reduced rank, no solution
    A = Matrix([[1, 2, 3], [2, 4, 6], [3, 6, 9]])
    b = Matrix([0, 0, 1])
    raises(ValueError, lambda: A.gauss_jordan_solve(b))

    # Rectangular, tall, full rank, unique solution
    A = Matrix([[1, 5, 3], [2, 1, 6], [1, 7, 9], [1, 4, 3]])
    b = Matrix([0, 0, 1, 0])
    sol, params = A.gauss_jordan_solve(b)
    assert sol == Matrix([[Rational(-1, 2)], [0], [Rational(1, 6)]])
    assert params == Matrix(0, 1, [])

    # Rectangular, tall, full rank, unique solution, B has less columns than rows
    A = Matrix([[1, 5, 3], [2, 1, 6], [1, 7, 9], [1, 4, 3]])
    B = Matrix([[0,0], [0, 0], [1, 2], [0, 0]])
    sol, params = A.gauss_jordan_solve(B)
    assert sol == Matrix([[Rational(-1, 2), Rational(-2, 2)], [0, 0], [Rational(1, 6), Rational(2, 6)]])
    assert params == Matrix(0, 2, [])

    # Rectangular, tall, full rank, no solution
    A = Matrix([[1, 5, 3], [2, 1, 6], [1, 7, 9], [1, 4, 3]])
    b = Matrix([0, 0, 0, 1])
    raises(ValueError, lambda: A.gauss_jordan_solve(b))

    # Rectangular, tall, full rank, no solution, B has two columns (2nd has no solution)
    A = Matrix([[1, 5, 3], [2, 1, 6], [1, 7, 9], [1, 4, 3]])
    B = Matrix([[0,0], [0, 0], [1, 0], [0, 1]])
    raises(ValueError, lambda: A.gauss_jordan_solve(B))

    # Rectangular, tall, full rank, no solution, B has two columns (1st has no solution)
    A = Matrix([[1, 5, 3], [2, 1, 6], [1, 7, 9], [1, 4, 3]])
    B = Matrix([[0,0], [0, 0], [0, 1], [1, 0]])
    raises(ValueError, lambda: A.gauss_jordan_solve(B))

    # Rectangular, tall, reduced rank, parametrized solution
    A = Matrix([[1, 5, 3], [2, 10, 6], [3, 15, 9], [1, 4, 3]])
    b = Matrix([0, 0, 0, 1])
    sol, params = A.gauss_jordan_solve(b)
    w = {}
    for s in sol.atoms(Symbol):
        w[s.name] = s
    assert sol == Matrix([[-3*w['tau0'] + 5], [-1], [w['tau0']]])
    assert params == Matrix([[w['tau0']]])

    # Rectangular, tall, reduced rank, no solution
    A = Matrix([[1, 5, 3], [2, 10, 6], [3, 15, 9], [1, 4, 3]])
    b = Matrix([0, 0, 1, 1])
    raises(ValueError, lambda: A.gauss_jordan_solve(b))

    # Rectangular, wide, full rank, parametrized solution
    A = Matrix([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 1, 12]])
    b = Matrix([1, 1, 1])
    sol, params = A.gauss_jordan_solve(b)
    w = {}
    for s in sol.atoms(Symbol):
        w[s.name] = s
    assert sol == Matrix([[2*w['tau0'] - 1], [-3*w['tau0'] + 1], [0],
                         [w['tau0']]])
    assert params == Matrix([[w['tau0']]])

    # Rectangular, wide, reduced rank, parametrized solution
    A = Matrix([[1, 2, 3, 4], [5, 6, 7, 8], [2, 4, 6, 8]])
    b = Matrix([0, 1, 0])
    sol, params = A.gauss_jordan_solve(b)
    w = {}
    for s in sol.atoms(Symbol):
        w[s.name] = s
    assert sol == Matrix([[w['tau0'] + 2*w['tau1'] + S.Half],
                         [-2*w['tau0'] - 3*w['tau1'] - Rational(1, 4)],
                         [w['tau0']], [w['tau1']]])
    assert params == Matrix([[w['tau0']], [w['tau1']]])
    # watch out for clashing symbols
    x0, x1, x2, _x0 = symbols('_tau0 _tau1 _tau2 tau1')
    M = Matrix([[0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, _x0]])
    A = M[:, :-1]
    b = M[:, -1:]
    sol, params = A.gauss_jordan_solve(b)
    assert params == Matrix(3, 1, [x0, x1, x2])
    assert sol == Matrix(5, 1, [x1, 0, x0, _x0, x2])

    # Rectangular, wide, reduced rank, no solution
    A = Matrix([[1, 2, 3, 4], [5, 6, 7, 8], [2, 4, 6, 8]])
    b = Matrix([1, 1, 1])
    raises(ValueError, lambda: A.gauss_jordan_solve(b))

    # Test for immutable matrix
    A = ImmutableMatrix([[1, 0], [0, 1]])
    B = ImmutableMatrix([1, 2])
    sol, params = A.gauss_jordan_solve(B)
    assert sol == ImmutableMatrix([1, 2])
    assert params == ImmutableMatrix(0, 1, [])
    assert sol.__class__ == ImmutableDenseMatrix
    assert params.__class__ == ImmutableDenseMatrix

def test_solve():
    A = Matrix([[1,2], [2,4]])
    b = Matrix([[3], [4]])
    raises(ValueError, lambda: A.solve(b)) #no solution
    b = Matrix([[ 4], [8]])
    raises(ValueError, lambda: A.solve(b)) #infinite solution

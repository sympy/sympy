from sympy import S
from sympy.core import symbols
from sympy.matrices import Matrix
from sympy.solvers.linear import solve_general_linear
from sympy.utilities.pytest import raises

t0, t1, t2 = symbols("t0 t1 t2")

def test_linear_solve_square():
    # Square, full rank, unique solution
    M = Matrix([[1,2,3], [4,5,6], [7,8,10]])
    v = Matrix([3,6,9])
    sol, params = solve_general_linear(M, v)
    assert sol == Matrix([[-1], [2], [0]])
    assert params == Matrix(0, 1, [])

    # Square, reduced rank, parametrised solution
    M = Matrix([[1,2,3], [4,5,6], [7,8,9]])
    v = Matrix([3,6,9])
    fs = (s for s in (t0,t1,t2))
    sol, params = solve_general_linear(M, v, fs)
    assert sol == Matrix([[t0 - 1], [-2*t0 + 2], [t0]])
    assert params == Matrix([[t0]])

    # Square, reduced rank, parametrised solution
    M = Matrix([[1,2,3], [2,4,6], [3,6,9]])
    v = Matrix([0,0,0])
    fs = (s for s in (t0,t1,t2))
    sol, params = solve_general_linear(M, v, fs)
    assert sol == Matrix([[-2*t0 - 3*t1], [t0], [t1]])
    assert params == Matrix([[t0], [t1]])

    # Square, reduced rank, parametrised solution
    M = Matrix([[0,0,0], [0,0,0], [0,0,0]])
    v = Matrix([0,0,0])
    fs = (s for s in (t0,t1,t2))
    sol, params = solve_general_linear(M, v, fs)
    assert sol == Matrix([[t0], [t1], [t2]])
    assert params == Matrix([[t0], [t1], [t2]])

    # Square, reduced rank, no solution
    M = Matrix([[1,2,3], [2,4,6], [3,6,9]])
    v = Matrix([0,0,1])
    raises(ValueError, lambda: solve_general_linear(M, v))


def test_linear_solve_rectangular_tall():
    # Rectangular, tall, full rank, unique solution
    M = Matrix([[1,5,3], [2,1,6], [1,7,9], [1,4,3]])
    v = Matrix([0,0,1,0])
    sol, params = solve_general_linear(M, v)
    assert sol == Matrix([[-S(1)/2], [0], [S(1)/6]])
    assert params == Matrix(0, 1, [])

    # Rectangular, tall, full rank, no solution
    M = Matrix([[1,5,3], [2,1,6], [1,7,9], [1,4,3]])
    v = Matrix([0,0,0,1])
    raises(ValueError, lambda: solve_general_linear(M, v))

    # Rectangular, tall, reduced rank, parametrised solution
    M = Matrix([[1,5,3], [2,10,6], [3,15,9], [1,4,3]])
    v = Matrix([0,0,0,1])
    fs = (s for s in (t0,t1,t2))
    sol, params = solve_general_linear(M, v, fs)
    assert sol == Matrix([[-3*t0 + 5], [-1], [t0]])
    assert params == Matrix([[t0]])

    # Rectangular, tall, reduced rank, no solution
    M = Matrix([[1,5,3], [2,10,6], [3,15,9], [1,4,3]])
    v = Matrix([0,0,1,1])
    raises(ValueError, lambda: solve_general_linear(M, v))


def test_linear_solve_rectangular_wide():
    # Rectangular, wide, full rank, parametrized solution
    M = Matrix([[1,2,3,4], [5,6,7,8], [9,10,1,12]])
    v = Matrix([1,1,1])
    fs = (s for s in (t0,t1,t2))
    sol, params = solve_general_linear(M, v, fs)
    assert sol == Matrix([[2*t0 - 1], [-3*t0 + 1], [0], [t0]])
    assert params == Matrix([[t0]])

    # Rectangular, wide, reduced rank, parametrized solution
    M = Matrix([[1,2,3,4], [5,6,7,8], [2,4,6,8]])
    v = Matrix([0,1,0])
    fs = (s for s in (t0,t1,t2))
    sol, params = solve_general_linear(M, v, fs)
    assert sol == Matrix([[t0 + 2*t1 + 1/S(2)], [-2*t0 - 3*t1 - 1/S(4)], [t0], [t1]])
    assert params == Matrix([[t0], [t1]])

    # Rectangular, wide, reduced rank, no solution
    M = Matrix([[1,2,3,4], [5,6,7,8], [2,4,6,8]])
    v = Matrix([1,1,1])
    raises(ValueError, lambda: solve_general_linear(M, v))

from sympy import (symbols, Symbol, diff, Function, Derivative, Matrix, Rational, S,
                   I, Eq, sqrt)
from sympy.core.containers import Tuple
from sympy.functions import exp, cos, sin, log
from sympy.matrices import dotprodsimp, NonSquareMatrixError
from sympy.solvers.ode import dsolve
from sympy.solvers.ode.ode import constant_renumber
from sympy.solvers.ode.subscheck import checksysodesol
from sympy.solvers.ode.systems import (neq_nth_linear_constant_coeff_match, linear_ode_to_matrix,
                                       ODEOrderError, ODENonlinearError, _simpsol, _solsimp,
                                       _is_commutative_anti_derivative, linodesolve,
                                       canonical_odes, dsolve_system, _component_division,
                                       _eqs2dict, _dict2graph)
from sympy.integrals.integrals import Integral
from sympy.testing.pytest import ON_TRAVIS, raises, slow, skip, XFAIL


C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 = symbols('C0:11')


def test_linear_ode_to_matrix():
    f, g, h = symbols("f, g, h", cls=Function)
    t = Symbol("t")
    funcs = [f(t), g(t), h(t)]
    f1 = f(t).diff(t)
    g1 = g(t).diff(t)
    h1 = h(t).diff(t)
    f2 = f(t).diff(t, 2)
    g2 = g(t).diff(t, 2)
    h2 = h(t).diff(t, 2)

    eqs_1 = [Eq(f1, g(t)), Eq(g1, f(t))]
    sol_1 = ([Matrix([[1, 0], [0, 1]]), Matrix([[ 0, 1], [1,  0]])], Matrix([[0],[0]]))
    assert linear_ode_to_matrix(eqs_1, funcs[:-1], t, 1) == sol_1

    eqs_2 = [Eq(f1, f(t) + 2*g(t)), Eq(g1, h(t)), Eq(h1, g(t) + h(t) + f(t))]
    sol_2 = ([Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), Matrix([[1, 2,  0], [ 0,  0, 1], [1, 1, 1]])],
             Matrix([[0], [0], [0]]))
    assert linear_ode_to_matrix(eqs_2, funcs, t, 1) == sol_2

    eqs_3 = [Eq(2*f1 + 3*h1, f(t) + g(t)), Eq(4*h1 + 5*g1, f(t) + h(t)), Eq(5*f1 + 4*g1, g(t) + h(t))]
    sol_3 = ([Matrix([[2, 0, 3], [0, 5, 4], [5, 4, 0]]), Matrix([[1, 1,  0], [1,  0, 1], [0, 1, 1]])],
             Matrix([[0], [0], [0]]))
    assert linear_ode_to_matrix(eqs_3, funcs, t, 1) == sol_3

    eqs_4 = [Eq(f2 + h(t), f1 + g(t)), Eq(2*h2 + g2 + g1 + g(t), 0), Eq(3*h1, 4)]
    sol_4 = ([Matrix([[1, 0, 0], [0, 1, 2], [0, 0, 0]]), Matrix([[1, 0, 0], [0, -1, 0], [0, 0, -3]]),
              Matrix([[0, 1, -1], [0,  -1, 0], [0, 0, 0]])], Matrix([[0], [0], [4]]))
    assert linear_ode_to_matrix(eqs_4, funcs, t, 2) == sol_4

    eqs_5 = [Eq(f2, g(t)), Eq(f1 + g1, f(t))]
    raises(ODEOrderError, lambda: linear_ode_to_matrix(eqs_5, funcs[:-1], t, 1))

    eqs_6 = [Eq(f1, f(t)**2), Eq(g1, f(t) + g(t))]
    raises(ODENonlinearError, lambda: linear_ode_to_matrix(eqs_6, funcs[:-1], t, 1))


def test_neq_nth_linear_constant_coeff_match():
    x, y, z, w = symbols('x, y, z, w', cls=Function)
    t = Symbol('t')
    x1 = diff(x(t), t)
    y1 = diff(y(t), t)
    z1 = diff(z(t), t)
    w1 = diff(w(t), t)
    x2 = diff(x(t), t, t)
    funcs = [x(t), y(t)]
    funcs_2 = funcs + [z(t), w(t)]

    eqs_1 = (5 * x1 + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * t * x(t) + 3 * y(t) + t))
    assert neq_nth_linear_constant_coeff_match(eqs_1, funcs, t) is None

    eqs_2 = (5 * (x1**2) + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * t * x(t) + 3 * y(t) + t))
    sol2 = {'is_implicit': True,
            'canon_eqs': [[Eq(Derivative(x(t), t), -sqrt(-12*x(t)/5 + 6*y(t)/5)),
            Eq(Derivative(y(t), t), 11*t*x(t)/2 - t/2 - 3*y(t)/2)],
            [Eq(Derivative(x(t), t), sqrt(-12*x(t)/5 + 6*y(t)/5)),
            Eq(Derivative(y(t), t), 11*t*x(t)/2 - t/2 - 3*y(t)/2)]]}
    assert neq_nth_linear_constant_coeff_match(eqs_2, funcs, t) == sol2

    eqs_2_1 = [Eq(Derivative(x(t), t), -sqrt(-12*x(t)/5 + 6*y(t)/5)),
            Eq(Derivative(y(t), t), 11*t*x(t)/2 - t/2 - 3*y(t)/2)]
    assert neq_nth_linear_constant_coeff_match(eqs_2_1, funcs, t) is None

    eqs_2_2 = [Eq(Derivative(x(t), t), sqrt(-12*x(t)/5 + 6*y(t)/5)),
            Eq(Derivative(y(t), t), 11*t*x(t)/2 - t/2 - 3*y(t)/2)]
    assert neq_nth_linear_constant_coeff_match(eqs_2_2, funcs, t) is None

    eqs_3 = (5 * x1 + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * x(t) + 3 * y(t)), (5 * w1 + z(t)), (z1 + w(t)))
    answer_3 = {'no_of_equation': 4,
     'eq': (12*x(t) - 6*y(t) + 5*Derivative(x(t), t),
      -11*x(t) + 3*y(t) + 2*Derivative(y(t), t),
      z(t) + 5*Derivative(w(t), t),
      w(t) + Derivative(z(t), t)),
     'func': [x(t), y(t), z(t), w(t)],
     'order': {x(t): 1, y(t): 1, z(t): 1, w(t): 1},
     'is_linear': True,
     'is_constant': True,
     'is_homogeneous': True,
     'func_coeff': -Matrix([
     [Rational(12, 5), Rational(-6, 5), 0, 0],
     [Rational(-11, 2),  Rational(3, 2),   0, 0],
     [0, 0, 0, 1],
     [0, 0, Rational(1, 5), 0]]),
     'type_of_equation': 'type1',
     'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eqs_3, funcs_2, t) == answer_3

    eqs_4 = (5 * x1 + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * x(t) + 3 * y(t)), (z1 - w(t)), (w1 - z(t)))
    answer_4 = {'no_of_equation': 4,
     'eq': (12 * x(t) - 6 * y(t) + 5 * Derivative(x(t), t),
            -11 * x(t) + 3 * y(t) + 2 * Derivative(y(t), t),
            -w(t) + Derivative(z(t), t),
            -z(t) + Derivative(w(t), t)),
     'func': [x(t), y(t), z(t), w(t)],
     'order': {x(t): 1, y(t): 1, z(t): 1, w(t): 1},
     'is_linear': True,
     'is_constant': True,
     'is_homogeneous': True,
     'func_coeff': -Matrix([
         [Rational(12, 5), Rational(-6, 5), 0, 0],
         [Rational(-11, 2), Rational(3, 2), 0, 0],
         [0, 0, 0, -1],
         [0, 0, -1, 0]]),
     'type_of_equation': 'type1',
     'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eqs_4, funcs_2, t) == answer_4

    eqs_5 = (5*x1 + 12*x(t) - 6*(y(t)) + x2, (2*y1 - 11*x(t) + 3*y(t)), (z1 - w(t)), (w1 - z(t)))
    assert neq_nth_linear_constant_coeff_match(eqs_5, funcs_2, t) is None

    eqs_6 = (Eq(x1, 3*y(t) - 11*z(t)), Eq(y1, 7*z(t) - 3*x(t)), Eq(z1, 11*x(t) - 7*y(t)))
    answer_6 = {'no_of_equation': 3, 'eq': (Eq(Derivative(x(t), t), 3*y(t) - 11*z(t)), Eq(Derivative(y(t), t), -3*x(t) + 7*z(t)),
            Eq(Derivative(z(t), t), 11*x(t) - 7*y(t))), 'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1},
            'is_linear': True, 'is_constant': True, 'is_homogeneous': True,
            'func_coeff': -Matrix([
                         [  0, -3, 11],
                         [  3,  0, -7],
                         [-11,  7,  0]]),
            'type_of_equation': 'type1', 'is_general': True}

    assert neq_nth_linear_constant_coeff_match(eqs_6, funcs_2[:-1], t) == answer_6

    eqs_7 = (Eq(x1, y(t)), Eq(y1, x(t)))
    answer_7 = {'no_of_equation': 2, 'eq': (Eq(Derivative(x(t), t), y(t)), Eq(Derivative(y(t), t), x(t))),
                'func': [x(t), y(t)], 'order': {x(t): 1, y(t): 1}, 'is_linear': True, 'is_constant': True,
                'is_homogeneous': True, 'func_coeff': -Matrix([
                                                        [ 0, -1],
                                                        [-1,  0]]),
                'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eqs_7, funcs, t) == answer_7

    eqs_8 = (Eq(x1, 21*x(t)), Eq(y1, 17*x(t) + 3*y(t)), Eq(z1, 5*x(t) + 7*y(t) + 9*z(t)))
    answer_8 = {'no_of_equation': 3, 'eq': (Eq(Derivative(x(t), t), 21*x(t)), Eq(Derivative(y(t), t), 17*x(t) + 3*y(t)),
            Eq(Derivative(z(t), t), 5*x(t) + 7*y(t) + 9*z(t))), 'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1},
            'is_linear': True, 'is_constant': True, 'is_homogeneous': True,
            'func_coeff': -Matrix([
                            [-21,  0,  0],
                            [-17, -3,  0],
                            [ -5, -7, -9]]),
            'type_of_equation': 'type1', 'is_general': True}

    assert neq_nth_linear_constant_coeff_match(eqs_8, funcs_2[:-1], t) == answer_8

    eqs_9 = (Eq(x1, 4*x(t) + 5*y(t) + 2*z(t)), Eq(y1, x(t) + 13*y(t) + 9*z(t)), Eq(z1, 32*x(t) + 41*y(t) + 11*z(t)))
    answer_9 = {'no_of_equation': 3, 'eq': (Eq(Derivative(x(t), t), 4*x(t) + 5*y(t) + 2*z(t)),
                Eq(Derivative(y(t), t), x(t) + 13*y(t) + 9*z(t)), Eq(Derivative(z(t), t), 32*x(t) + 41*y(t) + 11*z(t))),
                'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True,
                'is_constant': True, 'is_homogeneous': True,
                'func_coeff': -Matrix([
                            [ -4,  -5,  -2],
                            [ -1, -13,  -9],
                            [-32, -41, -11]]),
                'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eqs_9, funcs_2[:-1], t) == answer_9

    eqs_10 = (Eq(3*x1, 4*5*(y(t) - z(t))), Eq(4*y1, 3*5*(z(t) - x(t))), Eq(5*z1, 3*4*(x(t) - y(t))))
    answer_10 = {'no_of_equation': 3, 'eq': (Eq(3*Derivative(x(t), t), 20*y(t) - 20*z(t)),
                Eq(4*Derivative(y(t), t), -15*x(t) + 15*z(t)), Eq(5*Derivative(z(t), t), 12*x(t) - 12*y(t))),
                'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True,
                'is_constant': True, 'is_homogeneous': True,
                'func_coeff': -Matrix([
                                [  0, Rational(-20, 3),  Rational(20, 3)],
                                [Rational(15, 4),     0, Rational(-15, 4)],
                                [Rational(-12, 5), Rational(12, 5),  0]]),
                'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eqs_10, funcs_2[:-1], t) == answer_10

    eq11 = (Eq(x1, 3*y(t) - 11*z(t)), Eq(y1, 7*z(t) - 3*x(t)), Eq(z1, 11*x(t) - 7*y(t)))
    sol11 = {'no_of_equation': 3, 'eq': (Eq(Derivative(x(t), t), 3*y(t) - 11*z(t)), Eq(Derivative(y(t), t), -3*x(t) + 7*z(t)),
            Eq(Derivative(z(t), t), 11*x(t) - 7*y(t))), 'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1},
            'is_linear': True, 'is_constant': True, 'is_homogeneous': True, 'func_coeff': -Matrix([
            [  0, -3, 11], [  3,  0, -7], [-11,  7,  0]]), 'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eq11, funcs_2[:-1], t) == sol11

    eq12 = (Eq(Derivative(x(t), t), y(t)), Eq(Derivative(y(t), t), x(t)))
    sol12 = {'no_of_equation': 2, 'eq': (Eq(Derivative(x(t), t), y(t)), Eq(Derivative(y(t), t), x(t))),
             'func': [x(t), y(t)], 'order': {x(t): 1, y(t): 1}, 'is_linear': True, 'is_constant': True,
             'is_homogeneous': True, 'func_coeff': -Matrix([
            [0, -1],
            [-1, 0]]), 'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eq12, [x(t), y(t)], t) == sol12

    eq13 = (Eq(Derivative(x(t), t), 21*x(t)), Eq(Derivative(y(t), t), 17*x(t) + 3*y(t)),
            Eq(Derivative(z(t), t), 5*x(t) + 7*y(t) + 9*z(t)))
    sol13 = {'no_of_equation': 3, 'eq': (
    Eq(Derivative(x(t), t), 21 * x(t)), Eq(Derivative(y(t), t), 17 * x(t) + 3 * y(t)),
    Eq(Derivative(z(t), t), 5 * x(t) + 7 * y(t) + 9 * z(t))), 'func': [x(t), y(t), z(t)],
             'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True, 'is_constant': True, 'is_homogeneous': True,
             'func_coeff': -Matrix([
                 [-21, 0, 0],
                 [-17, -3, 0],
                 [-5, -7, -9]]), 'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eq13, [x(t), y(t), z(t)], t) == sol13

    eq14 = (
        Eq(Derivative(x(t), t), 4*x(t) + 5*y(t) + 2*z(t)), Eq(Derivative(y(t), t), x(t) + 13*y(t) + 9*z(t)),
        Eq(Derivative(z(t), t), 32*x(t) + 41*y(t) + 11*z(t)))
    sol14 = {'no_of_equation': 3, 'eq': (
    Eq(Derivative(x(t), t), 4 * x(t) + 5 * y(t) + 2 * z(t)), Eq(Derivative(y(t), t), x(t) + 13 * y(t) + 9 * z(t)),
    Eq(Derivative(z(t), t), 32 * x(t) + 41 * y(t) + 11 * z(t))), 'func': [x(t), y(t), z(t)],
             'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True, 'is_constant': True, 'is_homogeneous': True,
             'func_coeff': -Matrix([
                 [-4, -5, -2],
                 [-1, -13, -9],
                 [-32, -41, -11]]), 'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eq14, [x(t), y(t), z(t)], t) == sol14

    eq15 = (Eq(3*Derivative(x(t), t), 20*y(t) - 20*z(t)), Eq(4*Derivative(y(t), t), -15*x(t) + 15*z(t)),
            Eq(5*Derivative(z(t), t), 12*x(t) - 12*y(t)))
    sol15 = {'no_of_equation': 3, 'eq': (
    Eq(3 * Derivative(x(t), t), 20 * y(t) - 20 * z(t)), Eq(4 * Derivative(y(t), t), -15 * x(t) + 15 * z(t)),
    Eq(5 * Derivative(z(t), t), 12 * x(t) - 12 * y(t))), 'func': [x(t), y(t), z(t)],
             'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True, 'is_constant': True, 'is_homogeneous': True,
             'func_coeff': -Matrix([
                 [0, Rational(-20, 3), Rational(20, 3)],
                 [Rational(15, 4), 0, Rational(-15, 4)],
                 [Rational(-12, 5), Rational(12, 5), 0]]), 'type_of_equation': 'type1', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eq15, [x(t), y(t), z(t)], t) == sol15

    # Constant coefficient homogeneous ODEs
    eq1 = (Eq(diff(x(t), t), x(t) + y(t) + 9), Eq(diff(y(t), t), 2*x(t) + 5*y(t) + 23))
    sol1 = {'no_of_equation': 2, 'eq': (Eq(Derivative(x(t), t), x(t) + y(t) + 9),
        Eq(Derivative(y(t), t), 2*x(t) + 5*y(t) + 23)), 'func': [x(t), y(t)],
        'order': {x(t): 1, y(t): 1}, 'is_linear': True, 'is_constant': True, 'is_homogeneous': False, 'is_general': True,
        'func_coeff': -Matrix([[-1, -1], [-2, -5]]), 'rhs': Matrix([[ 9], [23]]), 'type_of_equation': 'type2'}
    assert neq_nth_linear_constant_coeff_match(eq1, funcs, t) == sol1

    # Non constant coefficient homogeneous ODEs
    eq1 = (Eq(diff(x(t), t), 5*t*x(t) + 2*y(t)), Eq(diff(y(t), t), 2*x(t) + 5*t*y(t)))
    sol1 = {'no_of_equation': 2, 'eq': (Eq(Derivative(x(t), t), 5*t*x(t) + 2*y(t)), Eq(Derivative(y(t), t), 5*t*y(t) + 2*x(t))),
            'func': [x(t), y(t)], 'order': {x(t): 1, y(t): 1}, 'is_linear': True, 'is_constant': False,
            'is_homogeneous': True, 'func_coeff': -Matrix([ [-5*t,   -2], [  -2, -5*t]]), 'commutative_antiderivative': Matrix([
            [5*t**2/2,      2*t], [     2*t, 5*t**2/2]]), 'type_of_equation': 'type3', 'is_general': True}
    assert neq_nth_linear_constant_coeff_match(eq1, funcs, t) == sol1

    # Non constant coefficient non-homogeneous ODEs
    eq1 = [Eq(x1, x(t) + t*y(t) + t), Eq(y1, t*x(t) + y(t))]
    sol1 = {'no_of_equation': 2, 'eq': [Eq(Derivative(x(t), t), t*y(t) + t + x(t)), Eq(Derivative(y(t), t),
            t*x(t) + y(t))], 'func': [x(t), y(t)], 'order': {x(t): 1, y(t): 1}, 'is_linear': True,
            'is_constant': False, 'is_homogeneous': False, 'is_general': True, 'func_coeff': -Matrix([ [-1, -t],
            [-t, -1]]), 'commutative_antiderivative': Matrix([ [     t, t**2/2], [t**2/2,      t]]), 'rhs':
            Matrix([ [t], [0]]), 'type_of_equation': 'type4'}
    assert neq_nth_linear_constant_coeff_match(eq1, funcs, t) == sol1

    eq2 = [Eq(x1, t*x(t) + t*y(t) + t), Eq(y1, t*x(t) + t*y(t) + cos(t))]
    sol2 = {'no_of_equation': 2, 'eq': [Eq(Derivative(x(t), t), t*x(t) + t*y(t) + t), Eq(Derivative(y(t), t),
            t*x(t) + t*y(t) + cos(t))], 'func': [x(t), y(t)], 'order': {x(t): 1, y(t): 1}, 'is_linear': True,
            'is_constant': False, 'is_homogeneous': False, 'is_general': True, 'func_coeff': -Matrix([ [-t, -t],
            [-t, -t]]), 'commutative_antiderivative': Matrix([ [t**2/2, t**2/2], [t**2/2, t**2/2]]), 'rhs':
            Matrix([ [     t], [cos(t)]]), 'type_of_equation': 'type4'}
    assert neq_nth_linear_constant_coeff_match(eq2, funcs, t) == sol2

    eq3 = [Eq(x1, t*(x(t) + y(t) + z(t) + 1)), Eq(y1, t*(x(t) + y(t) + z(t))), Eq(z1, t*(x(t) + y(t) + z(t)))]
    sol3 = {'no_of_equation': 3, 'eq': [Eq(Derivative(x(t), t), t*(x(t) + y(t) + z(t) + 1)),
            Eq(Derivative(y(t), t), t*(x(t) + y(t) + z(t))), Eq(Derivative(z(t), t), t*(x(t) + y(t) + z(t)))],
            'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True, 'is_constant':
            False, 'is_homogeneous': False, 'is_general': True, 'func_coeff': -Matrix([ [-t, -t, -t], [-t, -t,
            -t], [-t, -t, -t]]), 'commutative_antiderivative': Matrix([ [t**2/2, t**2/2, t**2/2], [t**2/2,
            t**2/2, t**2/2], [t**2/2, t**2/2, t**2/2]]), 'rhs': Matrix([ [t], [0], [0]]), 'type_of_equation':
            'type4'}
    assert neq_nth_linear_constant_coeff_match(eq3, funcs_2[:-1], t) == sol3

    eq4 = [Eq(x1, x(t) + y(t) + t*z(t) + 1), Eq(y1, x(t) + t*y(t) + z(t) + 10), Eq(z1, t*x(t) + y(t) + z(t) + t)]
    sol4 = {'no_of_equation': 3, 'eq': [Eq(Derivative(x(t), t), t*z(t) + x(t) + y(t) + 1), Eq(Derivative(y(t),
            t), t*y(t) + x(t) + z(t) + 10), Eq(Derivative(z(t), t), t*x(t) + t + y(t) + z(t))], 'func': [x(t),
            y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True, 'is_constant': False,
            'is_homogeneous': False, 'is_general': True, 'func_coeff': -Matrix([ [-1, -1, -t], [-1, -t, -1], [-t,
            -1, -1]]), 'commutative_antiderivative': Matrix([ [     t,      t, t**2/2], [     t, t**2/2,
            t], [t**2/2,      t,      t]]), 'rhs': Matrix([ [ 1], [10], [ t]]), 'type_of_equation': 'type4'}
    assert neq_nth_linear_constant_coeff_match(eq4, funcs_2[:-1], t) == sol4

    sum_terms = t*(x(t) + y(t) + z(t) + w(t))
    eq5 = [Eq(x1, sum_terms), Eq(y1, sum_terms), Eq(z1, sum_terms + 1), Eq(w1, sum_terms)]
    sol5 = {'no_of_equation': 4, 'eq': [Eq(Derivative(x(t), t), t*(w(t) + x(t) + y(t) + z(t))),
            Eq(Derivative(y(t), t), t*(w(t) + x(t) + y(t) + z(t))), Eq(Derivative(z(t), t), t*(w(t) + x(t) +
            y(t) + z(t)) + 1), Eq(Derivative(w(t), t), t*(w(t) + x(t) + y(t) + z(t)))], 'func': [x(t), y(t),
            z(t), w(t)], 'order': {x(t): 1, y(t): 1, z(t): 1, w(t): 1}, 'is_linear': True, 'is_constant': False,
            'is_homogeneous': False, 'is_general': True, 'func_coeff': -Matrix([ [-t, -t, -t, -t], [-t, -t, -t,
            -t], [-t, -t, -t, -t], [-t, -t, -t, -t]]), 'commutative_antiderivative': Matrix([ [t**2/2, t**2/2,
            t**2/2, t**2/2], [t**2/2, t**2/2, t**2/2, t**2/2], [t**2/2, t**2/2, t**2/2, t**2/2], [t**2/2,
            t**2/2, t**2/2, t**2/2]]), 'rhs': Matrix([ [0], [0], [1], [0]]), 'type_of_equation': 'type4'}
    assert neq_nth_linear_constant_coeff_match(eq5, funcs_2, t) == sol5

    # Multiple matchs

    f, g = symbols("f g", cls=Function)
    y = symbols("y")
    funcs = [f(t), g(t)]

    eq1 = [Eq(Derivative(f(t), t)**2 - 2*Derivative(f(t), t) + 1, 4),
            Eq(-y*f(t) + Derivative(g(t), t), 0)]
    sol1 = {'is_implicit': True,
             'canon_eqs': [[Eq(Derivative(f(t), t), -1), Eq(Derivative(g(t), t), y*f(t))],
              [Eq(Derivative(f(t), t), 3), Eq(Derivative(g(t), t), y*f(t))]]}
    assert neq_nth_linear_constant_coeff_match(eq1, funcs, t) == sol1

    raises(ValueError, lambda: neq_nth_linear_constant_coeff_match(eq1, funcs[:1], t))


def test_matrix_exp():
    from sympy.matrices.dense import Matrix, eye, zeros
    from sympy.solvers.ode.systems import matrix_exp
    t = Symbol('t')

    for n in range(1, 6+1):
        assert matrix_exp(zeros(n), t) == eye(n)

    for n in range(1, 6+1):
        A = eye(n)
        expAt = exp(t) * eye(n)
        assert matrix_exp(A, t) == expAt

    for n in range(1, 6+1):
        A = Matrix(n, n, lambda i,j: i+1 if i==j else 0)
        expAt = Matrix(n, n, lambda i,j: exp((i+1)*t) if i==j else 0)
        assert matrix_exp(A, t) == expAt

    A = Matrix([[0, 1], [-1, 0]])
    expAt = Matrix([[cos(t), sin(t)], [-sin(t), cos(t)]])
    assert matrix_exp(A, t) == expAt

    A = Matrix([[2, -5], [2, -4]])
    expAt = Matrix([
            [3*exp(-t)*sin(t) + exp(-t)*cos(t), -5*exp(-t)*sin(t)],
            [2*exp(-t)*sin(t), -3*exp(-t)*sin(t) + exp(-t)*cos(t)]
            ])
    assert matrix_exp(A, t) == expAt

    A = Matrix([[21, 17, 6], [-5, -1, -6], [4, 4, 16]])
    # TO update this.
    # expAt = Matrix([
    #     [(8*t*exp(12*t) + 5*exp(12*t) - 1)*exp(4*t)/4,
    #      (8*t*exp(12*t) + 5*exp(12*t) - 5)*exp(4*t)/4,
    #      (exp(12*t) - 1)*exp(4*t)/2],
    #     [(-8*t*exp(12*t) - exp(12*t) + 1)*exp(4*t)/4,
    #      (-8*t*exp(12*t) - exp(12*t) + 5)*exp(4*t)/4,
    #      (-exp(12*t) + 1)*exp(4*t)/2],
    #     [4*t*exp(16*t), 4*t*exp(16*t), exp(16*t)]])
    expAt = Matrix([
        [2*t*exp(16*t) + 5*exp(16*t)/4 - exp(4*t)/4, 2*t*exp(16*t) + 5*exp(16*t)/4 - 5*exp(4*t)/4,  exp(16*t)/2 - exp(4*t)/2],
        [ -2*t*exp(16*t) - exp(16*t)/4 + exp(4*t)/4,  -2*t*exp(16*t) - exp(16*t)/4 + 5*exp(4*t)/4, -exp(16*t)/2 + exp(4*t)/2],
        [                             4*t*exp(16*t),                                4*t*exp(16*t),                 exp(16*t)]
        ])
    assert matrix_exp(A, t) == expAt

    A = Matrix([[1, 1, 0, 0],
                [0, 1, 1, 0],
                [0, 0, 1, -S(1)/8],
                [0, 0, S(1)/2, S(1)/2]])
    expAt = Matrix([
        [exp(t), t*exp(t), 4*t*exp(3*t/4) + 8*t*exp(t) + 48*exp(3*t/4) - 48*exp(t),
                            -2*t*exp(3*t/4) - 2*t*exp(t) - 16*exp(3*t/4) + 16*exp(t)],
        [0, exp(t), -t*exp(3*t/4) - 8*exp(3*t/4) + 8*exp(t), t*exp(3*t/4)/2 + 2*exp(3*t/4) - 2*exp(t)],
        [0, 0, t*exp(3*t/4)/4 + exp(3*t/4), -t*exp(3*t/4)/8],
        [0, 0, t*exp(3*t/4)/2, -t*exp(3*t/4)/4 + exp(3*t/4)]
        ])
    assert matrix_exp(A, t) == expAt

    A = Matrix([
    [ 0, 1,  0, 0],
    [-1, 0,  0, 0],
    [ 0, 0,  0, 1],
    [ 0, 0, -1, 0]])

    expAt = Matrix([
    [ cos(t), sin(t),         0,        0],
    [-sin(t), cos(t),         0,        0],
    [      0,      0,    cos(t),   sin(t)],
    [      0,      0,   -sin(t),   cos(t)]])
    assert matrix_exp(A, t) == expAt

    A = Matrix([
    [ 0, 1,  1, 0],
    [-1, 0,  0, 1],
    [ 0, 0,  0, 1],
    [ 0, 0, -1, 0]])

    expAt = Matrix([
    [ cos(t), sin(t),  t*cos(t), t*sin(t)],
    [-sin(t), cos(t), -t*sin(t), t*cos(t)],
    [      0,      0,    cos(t),   sin(t)],
    [      0,      0,   -sin(t),   cos(t)]])
    assert matrix_exp(A, t) == expAt

    # This case is unacceptably slow right now but should be solvable...
    #a, b, c, d, e, f = symbols('a b c d e f')
    #A = Matrix([
    #[-a,  b,          c,  d],
    #[ a, -b,          e,  0],
    #[ 0,  0, -c - e - f,  0],
    #[ 0,  0,          f, -d]])

    A = Matrix([[0, I], [I, 0]])
    expAt = Matrix([
    [exp(I*t)/2 + exp(-I*t)/2, exp(I*t)/2 - exp(-I*t)/2],
    [exp(I*t)/2 - exp(-I*t)/2, exp(I*t)/2 + exp(-I*t)/2]])
    assert matrix_exp(A, t) == expAt

    # Testing Errors
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 7, 7]])
    M1 = Matrix([[t, 1], [1, 1]])

    raises(ValueError, lambda: matrix_exp(M[:, :2], t))
    raises(ValueError, lambda: matrix_exp(M[:2, :], t))
    raises(ValueError, lambda: matrix_exp(M1, t))
    raises(ValueError, lambda: matrix_exp(M1[:1, :1], t))


def test_canonical_odes():
    f, g, h = symbols('f g h', cls=Function)
    x = symbols('x')
    funcs = [f(x), g(x), h(x)]

    eqs1 = [Eq(f(x).diff(x, x), f(x) + 2*g(x)), Eq(g(x) + 1, g(x).diff(x) + f(x))]
    sol1 = [[Eq(Derivative(f(x), (x, 2)), f(x) + 2*g(x)), Eq(Derivative(g(x), x), -f(x) + g(x) + 1)]]
    assert canonical_odes(eqs1, funcs[:2], x) == sol1

    eqs2 = [Eq(f(x).diff(x), h(x).diff(x) + f(x)), Eq(g(x).diff(x)**2, f(x) + h(x)), Eq(h(x).diff(x), f(x))]
    sol2 = [[Eq(Derivative(f(x), x), 2*f(x)), Eq(Derivative(g(x), x), -sqrt(f(x) + h(x))), Eq(Derivative(h(x), x), f(x))],
            [Eq(Derivative(f(x), x), 2*f(x)), Eq(Derivative(g(x), x), sqrt(f(x) + h(x))), Eq(Derivative(h(x), x), f(x))]]
    assert canonical_odes(eqs2, funcs, x) == sol2


def test_sysode_linear_neq_order1_type1():

    f, g, x, y, h = symbols('f g x y h', cls=Function)
    a, b, c, t = symbols('a b c t')

    eq1 = [Eq(x(t).diff(t), x(t)), Eq(y(t).diff(t), y(t))]
    sol1 = [Eq(x(t), C1*exp(t)), Eq(y(t), C2*exp(t))]
    assert dsolve(eq1) == sol1
    assert checksysodesol(eq1, sol1) == (True, [0, 0])

    eq2 = [Eq(x(t).diff(t), 2*x(t)), Eq(y(t).diff(t), 3*y(t))]
    #sol2 = [Eq(x(t), C1*exp(2*t)), Eq(y(t), C2*exp(3*t))]
    sol2 = [Eq(x(t), C1*exp(2*t)), Eq(y(t), C2*exp(3*t))]
    assert dsolve(eq2) == sol2
    assert checksysodesol(eq2, sol2) == (True, [0, 0])

    eq3 = [Eq(x(t).diff(t), a*x(t)), Eq(y(t).diff(t), a*y(t))]
    sol3 = [Eq(x(t), C1*exp(a*t)), Eq(y(t), C2*exp(a*t))]
    assert dsolve(eq3) == sol3
    assert checksysodesol(eq3, sol3) == (True, [0, 0])

    # Regression test case for issue #15474
    # https://github.com/sympy/sympy/issues/15474
    eq4 = [Eq(x(t).diff(t), a*x(t)), Eq(y(t).diff(t), b*y(t))]
    sol4 = [Eq(x(t), C1*exp(a*t)), Eq(y(t), C2*exp(b*t))]
    assert dsolve(eq4) == sol4
    assert checksysodesol(eq4, sol4) == (True, [0, 0])

    eq5 = [Eq(x(t).diff(t), -y(t)), Eq(y(t).diff(t), x(t))]
    sol5 = [Eq(x(t), -C1*sin(t) - C2*cos(t)), Eq(y(t), C1*cos(t) - C2*sin(t))]
    assert dsolve(eq5) == sol5
    assert checksysodesol(eq5, sol5) == (True, [0, 0])

    eq6 = [Eq(x(t).diff(t), -2*y(t)), Eq(y(t).diff(t), 2*x(t))]
    sol6 = [Eq(x(t), -C1*sin(2*t) - C2*cos(2*t)), Eq(y(t), C1*cos(2*t) - C2*sin(2*t))]
    assert dsolve(eq6) == sol6
    assert checksysodesol(eq6, sol6) == (True, [0, 0])

    eq7 = [Eq(x(t).diff(t), I*y(t)), Eq(y(t).diff(t), I*x(t))]
    sol7 = [Eq(x(t), -C1*exp(-I*t) + C2*exp(I*t)), Eq(y(t), C1*exp(-I*t) + C2*exp(I*t))]
    assert dsolve(eq7) == sol7
    assert checksysodesol(eq7, sol7) == (True, [0, 0])

    eq8 = [Eq(x(t).diff(t), -a*y(t)), Eq(y(t).diff(t), a*x(t))]
    sol8 = [Eq(x(t), -I*C1*exp(-I*a*t) + I*C2*exp(I*a*t)), Eq(y(t), C1*exp(-I*a*t) + C2*exp(I*a*t))]
    assert dsolve(eq8) == sol8
    assert checksysodesol(eq8, sol8) == (True, [0, 0])

    eq9 = [Eq(x(t).diff(t), x(t) + y(t)), Eq(y(t).diff(t), x(t) - y(t))]
    sol9 = [Eq(x(t), C1*(1 - sqrt(2))*exp(-sqrt(2)*t) + C2*(1 + sqrt(2))*exp(sqrt(2)*t)),
           Eq(y(t), C1*exp(-sqrt(2)*t) + C2*exp(sqrt(2)*t))]
    assert dsolve(eq9) == sol9
    assert checksysodesol(eq9, sol9) == (True, [0, 0])

    eq10 = [Eq(x(t).diff(t), x(t) + y(t)), Eq(y(t).diff(t), x(t) + y(t))]
    sol10 = [Eq(x(t), -C1 + C2*exp(2*t)), Eq(y(t), C1 + C2*exp(2*t))]
    assert dsolve(eq10) == sol10
    assert checksysodesol(eq10, sol10) == (True, [0, 0])

    eq11 = [Eq(x(t).diff(t), 2*x(t) + y(t)), Eq(y(t).diff(t), -x(t) + 2*y(t))]
    sol11 = [Eq(x(t), (C1*sin(t) + C2*cos(t))*exp(2*t)),
            Eq(y(t), (C1*cos(t) - C2*sin(t))*exp(2*t))]
    assert dsolve(eq11) == sol11
    assert checksysodesol(eq11, sol11) == (True, [0, 0])

    eq12 = [Eq(x(t).diff(t), x(t) + 2*y(t)), Eq(y(t).diff(t), 2*x(t) + y(t))]
    sol12 = [Eq(x(t), -C1*exp(-t) + C2*exp(3*t)), Eq(y(t), C1*exp(-t) + C2*exp(3*t))]
    assert dsolve(eq12) == sol12
    assert checksysodesol(eq12, sol12) == (True, [0, 0])

    eq13 = [Eq(x(t).diff(t), 4*x(t) + y(t)), Eq(y(t).diff(t), -x(t) + 2*y(t))]
    sol13 = [Eq(x(t), (C1 + C2*t + C2)*exp(3*t)), Eq(y(t), (-C1 - C2*t)*exp(3*t))]
    assert dsolve(eq13) == sol13
    assert checksysodesol(eq13, sol13) == (True, [0, 0])

    eq14 = [Eq(x(t).diff(t), a*y(t)), Eq(y(t).diff(t), a*x(t))]
    sol14 = [Eq(x(t), -C1*exp(-a*t) + C2*exp(a*t)), Eq(y(t), C1*exp(-a*t) + C2*exp(a*t))]
    assert dsolve(eq14) == sol14
    assert checksysodesol(eq14, sol14) == (True, [0, 0])

    eq15 = [Eq(x(t).diff(t), a*y(t)), Eq(y(t).diff(t), b*x(t))]
    sol15 = [Eq(x(t), -C1*a*exp(-t*sqrt(a*b))/sqrt(a*b) + C2*a*exp(t*sqrt(a*b))/sqrt(a*b)),
            Eq(y(t), C1*exp(-t*sqrt(a*b)) + C2*exp(t*sqrt(a*b)))]
    assert dsolve(eq15) == sol15
    assert checksysodesol(eq15, sol15) == (True, [0, 0])

    eq16 = [Eq(x(t).diff(t), a*x(t) + b*y(t)), Eq(y(t).diff(t), c*x(t))]
    sol16 = [Eq(x(t), -2*C1*b*exp(t*(a/2 + sqrt(a**2 + 4*b*c)/2))/(a - sqrt(a**2 + 4*b*c)) - 2*C2*b*exp(t*(a/2 - sqrt(a**2 + 4*b*c)/2))/(a + sqrt(a**2 + 4*b*c))),
            Eq(y(t), C1*exp(t*(a/2 + sqrt(a**2 + 4*b*c)/2)) + C2*exp(t*(a/2 - sqrt(a**2 + 4*b*c)/2)))]
    assert dsolve(eq16) == sol16
    assert checksysodesol(eq16, sol16) == (True, [0, 0])

    # Regression test case for issue #18562
    # https://github.com/sympy/sympy/issues/18562
    eq17 = [Eq(x(t).diff(t), x(t) + a*y(t)), Eq(y(t).diff(t), x(t)*a - y(t))]
    sol17 = [Eq(x(t), C1*a*exp(t*sqrt(a**2 + 1))/(sqrt(a**2 + 1) - 1) - C2*a*exp(-t*sqrt(a**2 + 1))/(sqrt(a**2 + 1) + 1)),
            Eq(y(t), C1*exp(t*sqrt(a**2 + 1)) + C2*exp(-t*sqrt(a**2 + 1)))]
    assert dsolve(eq17) == sol17
    assert checksysodesol(eq17, sol17) == (True, [0, 0])

    eq18 = [Eq(x(t).diff(t), 0), Eq(y(t).diff(t), 0)]
    sol18 = [Eq(x(t), C1), Eq(y(t), C2)]
    assert dsolve(eq18) == sol18
    assert checksysodesol(eq18, sol18) == (True, [0, 0])

    eq19 = [Eq(x(t).diff(t), 2*x(t) - y(t)), Eq(y(t).diff(t), x(t))]
    sol19 = [Eq(x(t), (C1 + C2*t + C2)*exp(t)), Eq(y(t), (C1 + C2*t)*exp(t))]
    assert dsolve(eq19) == sol19
    assert checksysodesol(eq19, sol19) == (True, [0, 0])

    eq20 = [Eq(x(t).diff(t), x(t)), Eq(y(t).diff(t), x(t) + y(t))]
    sol20 = [Eq(x(t), C1*exp(t)), Eq(y(t), (C1*t + C2)*exp(t))]
    assert dsolve(eq20) == sol20
    assert checksysodesol(eq20, sol20) == (True, [0, 0])

    eq21 = [Eq(Derivative(x(t), t), 3*x(t)), Eq(Derivative(y(t), t), x(t) + y(t))]
    sol21 = [Eq(x(t), 2*C1*exp(3*t)), Eq(y(t), C1*exp(3*t) + C2*exp(t))]
    assert dsolve(eq21) == sol21
    assert checksysodesol(eq21, sol21) == (True, [0, 0])

    eq22 = [Eq(Derivative(x(t), t), 3*x(t)), Eq(Derivative(y(t), t), y(t))]
    sol22 = [Eq(x(t), C1*exp(3*t)), Eq(y(t), C2*exp(t))]
    assert dsolve(eq22) == sol22
    assert checksysodesol(eq22, sol22) == (True, [0, 0])

    Z0 = Function('Z0')
    Z1 = Function('Z1')
    Z2 = Function('Z2')
    Z3 = Function('Z3')

    k01, k10, k20, k21, k23, k30 = symbols('k01 k10 k20 k21 k23 k30')

    eq1 = (Eq(Derivative(Z0(t), t), -k01*Z0(t) + k10*Z1(t) + k20*Z2(t) + k30*Z3(t)),
           Eq(Derivative(Z1(t), t), k01*Z0(t) - k10*Z1(t) + k21*Z2(t)),
           Eq(Derivative(Z2(t), t), -(k20 + k21 + k23)*Z2(t)), Eq(Derivative(Z3(t), t), k23*Z2(t) - k30*Z3(t)))
    sol1 = [Eq(Z0(t), C1*k10/k01 + C2*(-k10 + k30)*exp(-k30*t)/(k01 + k10 - k30) + C3*(k10*k20 + k10*k21 - k10*k30 -
                k20**2 - k20*k21 - k20*k23 + k20*k30 + k23*k30)*exp(t*(-k20 - k21 - k23))/(k23*(k01 + k10 - k20 - k21 -
                k23)) - C4*exp(t*(-k01 - k10))),
            Eq(Z1(t), C1 - C2*k01*exp(-k30*t)/(k01 + k10 - k30) + C3*(k01*k20 + k01*k21 - k01*k30 - k20*k21 - k21**2 -
                k21*k23 + k21*k30)*exp(t*(-k20 - k21 - k23))/(k23*(k01 + k10 - k20 - k21 - k23)) + C4*exp(t*(-k01 - k10))),
            Eq(Z2(t), C3*(-k20 - k21 - k23 + k30)*exp(t*(-k20 - k21 - k23))/k23),
            Eq(Z3(t), C2*exp(-k30*t) + C3*exp(t*(-k20 - k21 - k23)))]

    assert dsolve(eq1, simplify=False) == sol1
    assert checksysodesol(eq1, sol1) == (True, [0, 0, 0, 0])

    x, y, z, u, v, w = symbols('x y z u v w', cls=Function)
    k2, k3 = symbols('k2 k3')
    a_b, a_c = symbols('a_b a_c', real=True)
    eq2 = (
        Eq(Derivative(z(t), t), k2*y(t)),
        Eq(Derivative(x(t), t), k3*y(t)),
        Eq(Derivative(y(t), t), (-k2 - k3)*y(t))
    )
    sol2 = [Eq(z(t), C1 - C2*k2*exp(t*(-k2 - k3))/(k2 + k3)),
            Eq(x(t), -C2*k3*exp(t*(-k2 - k3))/(k2 + k3) + C3),
            Eq(y(t), C2*exp(t*(-k2 - k3)))]
    assert dsolve(eq2) == sol2
    assert checksysodesol(eq2, sol2) == (True, [0, 0, 0])

    eq3 = [4*u(t) - v(t) - 2*w(t) + Derivative(u(t), t),
           2*u(t) + v(t) - 2*w(t) + Derivative(v(t), t),
           5*u(t) + v(t) - 3*w(t) + Derivative(w(t), t)]
    sol3 = [Eq(u(t), C1*cos(sqrt(3)*t)/2 - C2*sin(sqrt(3)*t)/2 + C3*exp(-2*t) + sqrt(3)*(C1*sin(sqrt(3)*t) + C2*cos(sqrt(3)*t))/6),
            Eq(v(t), C1*cos(sqrt(3)*t)/2 - C2*sin(sqrt(3)*t)/2 + sqrt(3)*(C1*sin(sqrt(3)*t) + C2*cos(sqrt(3)*t))/6),
            Eq(w(t), C1*cos(sqrt(3)*t) - C2*sin(sqrt(3)*t) + C3*exp(-2*t))]
    assert dsolve(eq3) == sol3
    assert checksysodesol(eq3, sol3) == (True, [0, 0, 0])

    tw = Rational(2, 9)
    eq4 = [Eq(x(t).diff(t), 2*x(t) + y(t) - tw*4*z(t) - tw*w(t)),
           Eq(y(t).diff(t), 2*y(t) + 8*tw*z(t) + 2*tw*w(t)),
           Eq(z(t).diff(t), Rational(37, 9)*z(t) - tw*w(t)),
           Eq(w(t).diff(t), 22*tw*w(t) - 2*tw*z(t))]

    sol4 = [Eq(x(t), (C1 + C2*t)*exp(2*t)),
            Eq(y(t), C2*exp(2*t) + 2*C3*exp(4*t)),
            Eq(z(t), 2*C3*exp(4*t) - C4*exp(5*t)/4),
            Eq(w(t), C3*exp(4*t) + C4*exp(5*t))]

    assert dsolve(eq4) == sol4
    assert checksysodesol(eq4, sol4) == (True, [0, 0, 0, 0])

    # Regression test case for issue #15574
    # https://github.com/sympy/sympy/issues/15574
    eq5 = [Eq(x(t).diff(t), x(t)), Eq(y(t).diff(t), y(t)), Eq(z(t).diff(t), z(t)), Eq(w(t).diff(t), w(t))]
    sol5 = [Eq(x(t), C1*exp(t)), Eq(y(t), C2*exp(t)), Eq(z(t), C3*exp(t)), Eq(w(t), C4*exp(t))]
    assert dsolve(eq5) == sol5
    assert checksysodesol(eq5, sol5) == (True, [0, 0, 0, 0])

    eq6 = [Eq(Derivative(x(t), t), x(t) + y(t)),
           Eq(Derivative(y(t), t), y(t) + z(t)),
           Eq(Derivative(z(t), t), -w(t)/8 + z(t)),
           Eq(Derivative(w(t), t), w(t)/2 + z(t)/2)]
    sol6 = [Eq(x(t), (C1 + C2*t)*exp(t) + (4*C3 + 4*C4*t + 48*C4)*exp(3*t/4)),
            Eq(y(t), C2*exp(t) + (-C3 - C4*t - 8*C4)*exp(3*t/4)),
            Eq(z(t), (C3/4 + C4*t/4 + C4)*exp(3*t/4)),
            Eq(w(t), (C3/2 + C4*t/2)*exp(3*t/4))]
    assert dsolve(eq6) == sol6
    assert checksysodesol(eq6, sol6) == (True, [0, 0, 0, 0])

    # Regression test case for issue #15574
    # https://github.com/sympy/sympy/issues/15574
    eq7 = [Eq(Derivative(x(t), t), x(t)), Eq(Derivative(y(t), t), y(t)), Eq(Derivative(z(t), t), z(t)),
           Eq(Derivative(w(t), t), w(t)), Eq(Derivative(u(t), t), u(t))]
    sol7 = [Eq(x(t), C1*exp(t)), Eq(y(t), C2*exp(t)), Eq(z(t), C3*exp(t)), Eq(w(t), C4*exp(t)),
            Eq(u(t), C5*exp(t))]
    assert dsolve(eq7) == sol7
    assert checksysodesol(eq7, sol7) == (True, [0, 0, 0, 0, 0])

    eq8 = [Eq(Derivative(x(t), t), 2*x(t) + y(t)), Eq(Derivative(y(t), t), 2*y(t)),
           Eq(Derivative(z(t), t), 4*z(t)),
           Eq(Derivative(w(t), t), u(t) + 5*w(t)), Eq(Derivative(u(t), t), 5*u(t))]
    sol8 = [Eq(x(t), (C1 + C2*t)*exp(2*t)), Eq(y(t), C2*exp(2*t)), Eq(z(t), C3*exp(4*t)),
            Eq(w(t), (C4 + C5*t)*exp(5*t)),
            Eq(u(t), C5*exp(5*t))]
    assert dsolve(eq8) == sol8
    assert checksysodesol(eq8, sol8) == (True, [0, 0, 0, 0, 0])

    # Regression test case for issue #15574
    # https://github.com/sympy/sympy/issues/15574
    eq9 = [Eq(Derivative(x(t), t), x(t)), Eq(Derivative(y(t), t), y(t)), Eq(Derivative(z(t), t), z(t))]
    sol9 = [Eq(x(t), C1*exp(t)), Eq(y(t), C2*exp(t)), Eq(z(t), C3*exp(t))]
    assert dsolve(eq9) == sol9
    assert checksysodesol(eq9, sol9) == (True, [0, 0, 0])

    # Regression test case for issue #15407
    # https://github.com/sympy/sympy/issues/15407
    eq10 = [Eq(Derivative(x(t), t), (-a_b - a_c)*x(t)), Eq(Derivative(y(t), t), a_b*y(t)), Eq(Derivative(z(t), t),
                                                                                              a_c*x(t))]
    sol10 = [Eq(x(t), -C1*(a_b + a_c)*exp(t*(-a_b - a_c))/a_c),
             Eq(y(t), C2*exp(a_b*t)),
             Eq(z(t), C1*exp(t*(-a_b - a_c)) + C3)]
    assert dsolve(eq10) == sol10
    assert checksysodesol(eq10, sol10) == (True, [0, 0, 0])

    # Regression test case for issue #14312
    # https://github.com/sympy/sympy/issues/14312
    eq11 = (Eq(Derivative(x(t), t), k3*y(t)), Eq(Derivative(y(t), t), (-k2 - k3)*y(t)), Eq(Derivative(z(t), t), k2*y(t)))
    sol11 = [Eq(x(t), C1 + C2*k3*exp(t*(-k2 - k3))/k2),
             Eq(y(t), -C2*(k2 + k3)*exp(t*(-k2 - k3))/k2),
             Eq(z(t), C2*exp(t*(-k2 - k3)) + C3)]
    assert dsolve(eq11) == sol11
    assert checksysodesol(eq11, sol11) == (True, [0, 0, 0])

    # Regression test case for issue #14312
    # https://github.com/sympy/sympy/issues/14312
    eq12 = (Eq(Derivative(z(t), t), k2*y(t)), Eq(Derivative(x(t), t), k3*y(t)), Eq(Derivative(y(t), t), (-k2 - k3)*y(t)))
    sol12 = [Eq(z(t), C1 - C2*k2*exp(t*(-k2 - k3))/(k2 + k3)),
             Eq(x(t), -C2*k3*exp(t*(-k2 - k3))/(k2 + k3) + C3),
             Eq(y(t), C2*exp(t*(-k2 - k3)))]
    assert dsolve(eq12) == sol12
    assert checksysodesol(eq12, sol12) == (True, [0, 0, 0])

    # Regression test case for issue #15474
    # https://github.com/sympy/sympy/issues/15474
    eq13 = [Eq(Derivative(f(t), t), 2*f(t) + g(t)), Eq(Derivative(g(t), t), a*f(t))]
    sol13 = [Eq(f(t), C1*exp(t*(sqrt(a + 1) + 1))/(sqrt(a + 1) - 1) - C2*exp(t*(1 - sqrt(a + 1)))/(sqrt(a + 1) + 1)),
             Eq(g(t), C1*exp(t*(sqrt(a + 1) + 1)) + C2*exp(t*(1 - sqrt(a + 1))))]
    assert dsolve(eq13) == sol13
    assert checksysodesol(eq13, sol13) == (True, [0, 0])

    eq14 = [Eq(f(t).diff(t), 2*g(t) - 3*h(t)),
            Eq(g(t).diff(t), 4*h(t) - 2*f(t)),
            Eq(h(t).diff(t), 3*f(t) - 4*g(t))]
    sol14 = [Eq(f(t), 2*C1 - 8*C2*cos(sqrt(29)*t)/25 + 8*C3*sin(sqrt(29)*t)/25 - 3*sqrt(29)*(C2*sin(sqrt(29)*t)
                + C3*cos(sqrt(29)*t))/25),
             Eq(g(t), 3*C1/2 - 6*C2*cos(sqrt(29)*t)/25 + 6*C3*sin(sqrt(29)*t)/25
                + 4*sqrt(29)*(C2*sin(sqrt(29)*t) + C3*cos(sqrt(29)*t))/25),
             Eq(h(t), C1 + C2*cos(sqrt(29)*t)
                - C3*sin(sqrt(29)*t))]

    assert dsolve(eq14) == sol14
    assert checksysodesol(eq14, sol14) == (True, [0, 0, 0])

    eq15 = [Eq(2*f(t).diff(t), 3*4*(g(t) - h(t))),
            Eq(3*g(t).diff(t), 2*4*(h(t) - f(t))),
            Eq(4*h(t).diff(t), 2*3*(f(t) - g(t)))]
    sol15 = [Eq(f(t), C1 - 16*C2*cos(sqrt(29)*t)/13 + 16*C3*sin(sqrt(29)*t)/13 - 6*sqrt(29)*(
                C2*sin(sqrt(29)*t) + C3*cos(sqrt(29)*t))/13),
            Eq(g(t), C1 - 16*C2*cos(sqrt(29)*t)/13 + 16*C3*sin(sqrt(29)*t)/13 + 8*sqrt(29)*(C2*sin(sqrt(29)*t) +
                C3*cos(sqrt(29)*t))/39),
            Eq(h(t), C1 + C2*cos(sqrt(29)*t) - C3*sin(sqrt(29)*t))]

    assert dsolve(eq15) == sol15
    assert checksysodesol(eq15, sol15) == (True, [0, 0, 0])

    eq16 = (Eq(diff(x(t), t), 21*x(t)), Eq(diff(y(t), t), 17*x(t) + 3*y(t)),
            Eq(diff(z(t), t), 5*x(t) + 7*y(t) + 9*z(t)))
    sol16 = [Eq(x(t), 216*C1*exp(21*t)/209),
             Eq(y(t), 204*C1*exp(21*t)/209 - 6*C2*exp(3*t)/7),
             Eq(z(t), C1*exp(21*t) + C2*exp(3*t) + C3*exp(9*t))]
    assert dsolve(eq16) == sol16
    assert checksysodesol(eq16, sol16) == (True, [0, 0, 0])

    eq17 = (Eq(diff(x(t), t), 3*y(t) - 11*z(t)), Eq(diff(y(t), t), 7*z(t) - 3*x(t)), Eq(diff(z(t), t), 11*x(t) - 7*y(t)))
    sol17 = [Eq(x(t), 7*C1/3 - 21*C2*cos(sqrt(179)*t)/170 + 21*C3*sin(sqrt(179)*t)/170 - 11*sqrt(179)*(
                C2*sin(sqrt(179)*t) + C3*cos(sqrt(179)*t))/170),
             Eq(y(t), 11*C1/3 - 33*C2*cos(sqrt(179)*t)/170 + 33*C3*sin(sqrt(179)*t)/170 + 7*sqrt(179)*(
                 C2*sin(sqrt(179)*t) + C3*cos(sqrt(179)*t))/170),
             Eq(z(t), C1 + C2*cos(sqrt(179)*t) - C3*sin(sqrt(179)*t))]

    assert dsolve(eq17) == sol17
    assert checksysodesol(eq17, sol17) == (True, [0, 0, 0])

    eq18 = (Eq(3*diff(x(t), t), 4*5*(y(t) - z(t))), Eq(4*diff(y(t), t), 3*5*(z(t) - x(t))),
            Eq(5*diff(z(t), t), 3*4*(x(t) - y(t))))
    sol18 = [Eq(x(t), C1 - C2*cos(5*sqrt(2)*t) + C3*sin(5*sqrt(2)*t) - 4*sqrt(2)*(C2*sin(5*sqrt(2)*t) + C3*cos(5*sqrt(2)*t))/3),
             Eq(y(t), C1 - C2*cos(5*sqrt(2)*t) + C3*sin(5*sqrt(2)*t) + 3*sqrt(2)*(C2*sin(5*sqrt(2)*t) + C3*cos(5*sqrt(2)*t))/4),
             Eq(z(t), C1 + C2*cos(5*sqrt(2)*t) - C3*sin(5*sqrt(2)*t))]

    assert dsolve(eq18) == sol18
    assert checksysodesol(eq18, sol18) == (True, [0, 0, 0])

    eq19 = (Eq(diff(x(t), t), 4*x(t) - z(t)), Eq(diff(y(t), t), 2*x(t) + 2*y(t) - z(t)), Eq(diff(z(t), t), 3*x(t) + y(t)))
    sol19 = [Eq(x(t), (C1 + C2*t**2/2 + 2*C2*t + C2 + C3*t + 2*C3)*exp(2*t)),
             Eq(y(t), (C1 + C2*t**2/2 + 2*C2*t + C3*t + 2*C3)*exp(2*t)),
             Eq(z(t), (2*C1 + C2*t**2 + 3*C2*t + 2*C3*t + 3*C3)*exp(2*t))]
    assert dsolve(eq19) == sol19
    assert checksysodesol(eq19, sol19) == (True, [0, 0, 0])

    eq20 = (Eq(diff(x(t), t), 4*x(t) - y(t) - 2*z(t)), Eq(diff(y(t), t), 2*x(t) + y(t) - 2*z(t)),
            Eq(diff(z(t), t), 5*x(t) - 3*z(t)))
    sol20 = [Eq(x(t), C1*exp(2*t) - 3*C2*sin(t)/5 - C2*cos(t)/5 - C3*sin(t)/5 + 3*C3*cos(t)/5),
             Eq(y(t), -3*C2*sin(t)/5 - C2*cos(t)/5 - C3*sin(t)/5 + 3*C3*cos(t)/5),
             Eq(z(t), C1*exp(2*t) - C2*sin(t) + C3*cos(t))]
    assert dsolve(eq20) == sol20
    assert checksysodesol(eq20, sol20) == (True, [0, 0, 0])

    eq21 = (Eq(diff(x(t), t), 9*y(t)), Eq(diff(y(t), t), 12*x(t)))
    sol21 = [Eq(x(t), -sqrt(3)*C1*exp(-6*sqrt(3)*t)/2 + sqrt(3)*C2*exp(6*sqrt(3)*t)/2),
             Eq(y(t), C1*exp(-6*sqrt(3)*t) + C2*exp(6*sqrt(3)*t))]

    assert dsolve(eq21) == sol21
    assert checksysodesol(eq21, sol21) == (True, [0, 0])

    eq22 = (Eq(Derivative(x(t), t), 2*x(t) + 4*y(t)), Eq(Derivative(y(t), t), 12*x(t) + 41*y(t)))
    sol22 = [Eq(x(t), C1*(-Rational(13, 8) + sqrt(1713)/24)*exp(t*(sqrt(1713)/2 + Rational(43, 2))) + C2*(
                - sqrt(1713)/24 - Rational(13, 8))*exp(t*(Rational(43, 2) - sqrt(1713)/2))),
             Eq(y(t), C1*exp(t*(sqrt(1713)/2 + Rational(43, 2))) + C2*exp(t*(Rational(43, 2) - sqrt(1713)/2)))]
    assert dsolve(eq22) == sol22
    assert checksysodesol(eq22, sol22) == (True, [0, 0])

    eq23 = (Eq(Derivative(x(t), t), x(t) + y(t)), Eq(Derivative(y(t), t), -2*x(t) + 2*y(t)))
    sol23 = [
        Eq(x(t), (C1*cos(sqrt(7)*t/2)/4 - C2*sin(sqrt(7)*t/2)/4 + sqrt(7)*(C1*sin(sqrt(7)*t/2) + C2*cos(
            sqrt(7)*t/2))/4)*exp(3*t/2)),
        Eq(y(t), (C1*cos(sqrt(7)*t/2) - C2*sin(sqrt(7)*t/2))*exp(3*t/2))
    ]
    assert dsolve(eq23) == sol23
    assert checksysodesol(eq23, sol23) == (True, [0, 0])

    # Regression test case for issue #15474
    # https://github.com/sympy/sympy/issues/15474
    a = Symbol("a", real=True)
    eq24 = [x(t).diff(t) - a*y(t), y(t).diff(t) + a*x(t)]
    sol24 = [Eq(x(t), C1*sin(a*t) + C2*cos(a*t)), Eq(y(t), C1*cos(a*t) - C2*sin(a*t))]
    assert dsolve(eq24) == sol24
    assert checksysodesol(eq24, sol24) == (True, [0, 0])

    # Regression test case for issue #19150
    # https://github.com/sympy/sympy/issues/19150
    eq25 = [Eq(Derivative(f(t), t), 0),
            Eq(Derivative(g(t), t), 1/(c*b)*(-2*g(t) + x(t) + f(t))),
            Eq(Derivative(x(t), t), 1/(c*b)*(-2*x(t) + g(t) + y(t))),
            Eq(Derivative(y(t), t), 1/(c*b)*(-2*y(t) + x(t) + h(t))),
            Eq(Derivative(h(t), t), 0)]

    sol25 = [Eq(f(t), -3*C1 + 4*C2),
             Eq(g(t), -2*C1 + 3*C2 - C3*exp(-2*t/(b*c)) + C4*exp(t*(-2 - sqrt(2))/(b*c)) +
                C5*exp(t*(-2 + sqrt(2))/(b*c))),
             Eq(x(t), -C1 + 2*C2 - sqrt(2)*C4*exp(t*(-2 - sqrt(2))/(b*c)) + sqrt(2)*C5*exp(t*(-2 + sqrt(2))/(b*c))),
             Eq(y(t), C2 + C3*exp(-2*t/(b*c)) + C4*exp(t*(-2 - sqrt(2))/(b*c)) + C5*exp(
                 t*(-2 + sqrt(2))/(b*c))),
             Eq(h(t), C1)]

    assert dsolve(eq25) == sol25
    assert checksysodesol(eq25, sol25) == (True, [0, 0, 0, 0, 0])

    eq26 = [Eq(Derivative(f(t), t), 2*f(t)), Eq(Derivative(g(t), t), 3*f(t) + 7*g(t))]
    sol26 = [Eq(f(t), -5*C1*exp(2*t)/3), Eq(g(t), C1*exp(2*t) + C2*exp(7*t))]
    assert dsolve(eq26) == sol26
    assert checksysodesol(eq26, sol26) == (True, [0, 0])

    eq27 = [Eq(Derivative(f(t), t), -9*I*f(t) - 4*g(t)), Eq(Derivative(g(t), t), -4*I*g(t))]
    sol27 = [Eq(f(t), 4*I*C1*exp(-4*I*t)/5 + C2*exp(-9*I*t)), Eq(g(t), C1*exp(-4*I*t))]
    assert dsolve(eq27) == sol27
    assert checksysodesol(eq27, sol27) == (True, [0, 0])

    eq28 = [Eq(Derivative(f(t), t), -9*I*f(t)), Eq(Derivative(g(t), t), -4*I*g(t))]
    sol28 = [Eq(f(t), C1*exp(-9*I*t)), Eq(g(t), C2*exp(-4*I*t))]
    assert dsolve(eq28) == sol28
    assert checksysodesol(eq28, sol28) == (True, [0, 0])

    eq29 = [Eq(Derivative(f(t), t), 0), Eq(Derivative(g(t), t), 0)]
    sol29 = [Eq(f(t), C1), Eq(g(t), C2)]
    assert dsolve(eq29) == sol29
    assert checksysodesol(eq29, sol29) == (True, [0, 0])

    eq30 = [Eq(Derivative(f(t), t), f(t)), Eq(Derivative(g(t), t), 0)]
    sol30 = [Eq(f(t), C1*exp(t)), Eq(g(t), C2)]
    assert dsolve(eq30) == sol30
    assert checksysodesol(eq30, sol30) == (True, [0, 0])

    eq31 = [Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), 0)]
    sol31 = [Eq(f(t), C1 + C2*t), Eq(g(t), C2)]
    assert dsolve(eq31) == sol31
    assert checksysodesol(eq31, sol31) == (True, [0, 0])

    eq32 = [Eq(Derivative(f(t), t), 0), Eq(Derivative(g(t), t), f(t))]
    sol32 = [Eq(f(t), C1), Eq(g(t), C1*t + C2)]
    assert dsolve(eq32) == sol32
    assert checksysodesol(eq32, sol32) == (True, [0, 0])

    eq33 = [Eq(Derivative(f(t), t), 0), Eq(Derivative(g(t), t), g(t))]
    sol33 = [Eq(f(t), C1), Eq(g(t), C2*exp(t))]
    assert dsolve(eq33) == sol33
    assert checksysodesol(eq33, sol33) == (True, [0, 0])

    eq34 = [Eq(Derivative(f(t), t), f(t)), Eq(Derivative(g(t), t), I*g(t))]
    sol34 = [Eq(f(t), C1*exp(t)), Eq(g(t), C2*exp(I*t))]
    assert dsolve(eq34) == sol34
    assert checksysodesol(eq34, sol34) == (True, [0, 0])

    eq35 = [Eq(Derivative(f(t), t), I*f(t)), Eq(Derivative(g(t), t), -I*g(t))]
    sol35 = [Eq(f(t), C1*exp(I*t)), Eq(g(t), C2*exp(-I*t))]
    assert dsolve(eq35) == sol35
    assert checksysodesol(eq35, sol35) == (True, [0, 0])

    eq36 = [Eq(Derivative(f(t), t), I*g(t)), Eq(Derivative(g(t), t), 0)]
    sol36 = [Eq(f(t), I*(C1 + C2*t)), Eq(g(t), C2)]
    assert dsolve(eq36) == sol36
    assert checksysodesol(eq36, sol36) == (True, [0, 0])

    eq37 = [Eq(Derivative(f(t), t), I*g(t)), Eq(Derivative(g(t), t), I*f(t))]
    sol37 = [Eq(f(t), -C1*exp(-I*t) + C2*exp(I*t)), Eq(g(t), C1*exp(-I*t) + C2*exp(I*t))]
    assert dsolve(eq37) == sol37
    assert checksysodesol(eq37, sol37) == (True, [0, 0])

    # Multiple systems
    eq1 = [Eq(Derivative(f(t), t)**2, g(t)**2), Eq(-f(t) + Derivative(g(t), t), 0)]
    sol1 = [[Eq(f(t), -C1*sin(t) - C2*cos(t)),
             Eq(g(t), C1*cos(t) - C2*sin(t))],
            [Eq(f(t), -C1*exp(-t) + C2*exp(t)),
             Eq(g(t), C1*exp(-t) + C2*exp(t))]]
    assert dsolve(eq1) == sol1
    for sol in sol1:
        assert checksysodesol(eq1, sol) == (True, [0, 0])


def test_sysode_linear_neq_order1_type2():

    f, g, h, k = symbols('f g h k', cls=Function)
    x, t, a, b, c, d, y = symbols('x t a b c d y')

    eq1 = [Eq(diff(f(x), x),  f(x) + g(x) + 5),
           Eq(diff(g(x), x), -f(x) - g(x) + 7)]
    sol1 = [Eq(f(x), C1 + C2*x + C2 + x*Integral(12, x) + Integral(12, x) + Integral(-12*x - 7, x)),
            Eq(g(x), -C1 - C2*x - x*Integral(12, x) - Integral(-12*x - 7, x))]
    assert dsolve(eq1) == sol1
    assert checksysodesol(eq1, sol1) == (True, [0, 0])

    eq2 = [Eq(diff(f(x), x), f(x) + g(x) + 5),
           Eq(diff(g(x), x), f(x) + g(x) + 7)]
    sol2 = [Eq(f(x), -C1 + C2*exp(2*x) + exp(2*x)*Integral(6*exp(-2*x), x) - Integral(1, x)),
            Eq(g(x), C1 + C2*exp(2*x) + exp(2*x)*Integral(6*exp(-2*x), x) + Integral(1, x))]
    assert dsolve(eq2) == sol2
    assert checksysodesol(eq2, sol2) == (True, [0, 0])

    eq3 = [Eq(diff(f(x), x), f(x) + 5), Eq(diff(g(x), x), f(x) + 7)]
    sol3 = [Eq(f(x), C1*exp(x) + exp(x)*Integral(5*exp(-x), x)),
            Eq(g(x), C1*exp(x) + C2 + exp(x)*Integral(5*exp(-x), x) + Integral(2, x))]
    assert dsolve(eq3) == sol3
    assert checksysodesol(eq3, sol3) == (True, [0, 0])

    eq4 = [Eq(diff(f(x), x), f(x) + exp(x)), Eq(diff(g(x), x), f(x) + g(x) + x*exp(x))]
    sol4 = [Eq(f(x), (C1 + Integral(1, x))*exp(x)), Eq(g(x), (C1*x + C2 + x*Integral(1, x) + Integral(0, x))*exp(x))]
    assert dsolve(eq4) == sol4
    assert checksysodesol(eq4, sol4) == (True, [0, 0])

    eq5 = [Eq(diff(f(x), x), f(x) + g(x) + 5*x), Eq(diff(g(x), x), f(x) - g(x))]
    sol5 = [Eq(f(x), C1*exp(sqrt(2)*x) + sqrt(2)*C1*exp(sqrt(2)*x) + (-sqrt(2)*C2 + C2 - sqrt(2)*Integral(
                -5*sqrt(2)*x*exp(sqrt(2)*x)/4, x) + Integral(-5*sqrt(2)*x*exp(sqrt(2)*x)/4, x))*exp(-sqrt(2)*x) +
                exp(sqrt(2)*x)*Integral(5*sqrt(2)*x*exp(-sqrt(2)*x)/4, x) + sqrt(2)*exp(sqrt(2)*x)*Integral(
                5*sqrt(2)*x*exp(-sqrt(2)*x)/4, x)),
            Eq(g(x), C1*exp(sqrt(2)*x) + (C2 + Integral(-5*sqrt(2)*x*exp(sqrt(2)*x)/4, x))*exp(-sqrt(2)*x) + exp(
                sqrt(2)*x)*Integral(5*sqrt(2)*x*exp(-sqrt(2)*x)/4, x))]
    assert dsolve(eq5) == sol5
    assert checksysodesol(eq5, sol5) == (True, [0, 0])

    eq6 = [Eq(diff(f(x), x), -9*f(x) - 4*g(x)),
           Eq(diff(g(x), x), -4*g(x)),
           Eq(diff(h(x), x), h(x) + exp(x))]
    sol6 = [Eq(f(x), (-4*C1/5 - 4*Integral(0, x)/5)*exp(-4*x) + (C2 + Integral(0, x))*exp(-9*x)),
            Eq(g(x), (C1 + Integral(0, x))*exp(-4*x)), Eq(h(x), (C3 + Integral(1, x))*exp(x))]
    assert dsolve(eq6) == sol6
    assert checksysodesol(eq6, sol6) == (True, [0, 0, 0])

    # Regression test case for issue #8859
    # https://github.com/sympy/sympy/issues/8859
    eq7 = [Eq(diff(f(t), t), f(t) + 3*t), Eq(diff(g(t), t), g(t))]
    sol7 = [Eq(f(t), C1*exp(t) + exp(t)*Integral(3*t*exp(-t), t)), Eq(g(t), C2*exp(t) + exp(t)*Integral(0, t))]
    assert dsolve(eq7) == sol7
    assert checksysodesol(eq7, sol7) == (True, [0, 0])

    # Regression test case for issue #8567
    # https://github.com/sympy/sympy/issues/8567
    eq8 = [Eq(f(t).diff(t), f(t) + 2*g(t)), Eq(g(t).diff(t), -2*f(t) + g(t) + 2*exp(t))]
    sol8 = [Eq(f(t), (C1*sin(2*t) + C2*cos(2*t) + sin(2*t)*Integral(-2*sin(2*t)**2/cos(2*t) + 2/cos(2*t), t) +
                cos(2*t)*Integral(-2*sin(2*t), t))*exp(t)),
            Eq(g(t), (C1*cos(2*t) - C2*sin(2*t) - sin(2*t)*Integral(-2*sin(2*t), t) +
                cos(2*t)*Integral(-2*sin(2*t)**2/cos(2*t) + 2/cos(2*t), t))*exp(t))]
    assert dsolve(eq8) == sol8
    assert checksysodesol(eq8, sol8) == (True, [0, 0])

    # Regression test case for issue #19150
    # https://github.com/sympy/sympy/issues/19150
    eq9 = [Eq(Derivative(f(t), t), 1/(a*b)*(-2*f(t) + g(t) + c)),
           Eq(Derivative(g(t), t), 1/(a*b)*(-2*g(t) + f(t) + h(t))),
           Eq(Derivative(h(t), t), 1/(a*b)*(-2*h(t) + g(t) + d))]
    sol9 = [Eq(f(t), (-C1 + C2*exp(-sqrt(2)*t/(a*b)) + C3*exp(sqrt(2)*t/(a*b)) + exp(sqrt(2)*t/(a*b))*Integral(
                c*exp(-sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b) + d*exp(-sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b), t) - Integral(
                -c*exp(2*t/(a*b))/(2*a*b) + d*exp(2*t/(a*b))/(2*a*b), t) + exp(-sqrt(2)*t/(a*b))*Integral(
                c*exp(sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b) + d*exp(sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b), t))*exp(
                -2*t/(a*b))),
            Eq(g(t), (-sqrt(2)*C2*exp(-sqrt(2)*t/(a*b)) + sqrt(2)*C3*exp(sqrt(2)*t/(a*b)) +
                sqrt(2)*exp(sqrt(2)*t/(a*b))*Integral(c*exp(-sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b) + d*exp(-sqrt(2)*t/(a*b) +
                2*t/(a*b))/(4*a*b), t) - sqrt(2)*exp(-sqrt(2)*t/(a*b))*Integral(c*exp(sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b) +
                d*exp(sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b), t))*exp(-2*t/(a*b))),
            Eq(h(t), (C1 + C2*exp(-sqrt(2)*t/(a*b)) +
                C3*exp(sqrt(2)*t/(a*b)) + exp(sqrt(2)*t/(a*b))*Integral(c*exp(-sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b) +
                d*exp(-sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b), t) + Integral(-c*exp(2*t/(a*b))/(2*a*b) + d*exp(2*t/(a*b))/(
                2*a*b), t) + exp(-sqrt(2)*t/(a*b))*Integral(c*exp(sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b) + d*exp(
                sqrt(2)*t/(a*b) + 2*t/(a*b))/(4*a*b), t))*exp(-2*t/(a*b)))]
    assert dsolve(eq9) == sol9
    assert checksysodesol(eq9, sol9) == (True, [0, 0, 0])

    # Simpsol and Solsimp testing
    _x1 = sqrt(2)
    _x2 = exp(2*_x1*t/(a*b))
    _x3 = 4*C3*_x2*a*b
    _x4 = exp(2*t/(a*b))
    _x5 = exp(_x1*t/(a*b))
    _x6 = Integral(_x4*_x5, t)
    _x7 = exp(-2*t/(a*b))
    _x8 = exp(-_x1*t/(a*b))
    _x9 = Integral(_x4*_x8, t)
    _x10 = _x2*_x9*d
    _x11 = 4*C2*a*b
    sol9_simpsol = [
        Eq(f(t), _x7*_x8*(-4*C1*_x5*a*b + _x10 + _x11 + _x2*_x9*c + _x3 + 2*_x5*c*Integral(_x4,
                                                                                           t) - 2*_x5*d*Integral(
            _x4, t) + _x6*c + _x6*d)/(4*a*b)),
        Eq(g(t), _x7*_x8*(
            -4*C2*_x1*a*b + 4*C3*_x1*_x2*a*b + _x1*_x2*_x9*c + _x1*_x2*_x9*d - _x1*_x6*c - _x1*_x6*d)/(
               4*a*b)),
        Eq(h(t), _x7*_x8*(4*C1*_x5*a*b + _x10 + _x11 + _x2*_x9*c + _x3 - 2*_x5*c*Integral(_x4,
                                                                                          t) + 2*_x5*d*Integral(
            _x4, t) + _x6*c + _x6*d)/(4*a*b)),
    ]
    assert [_simpsol(s) for s in sol9] == sol9_simpsol

    _x1 = sqrt(2)
    _x2 = 1/b
    _x3 = exp(_x2*t*(2 - _x1)/a)
    _x4 = exp(_x2*t*(_x1 + 2)/a)
    _x5 = exp(-2*_x2*t/a)
    _x6 = exp(2*_x2*t/a)
    _x7 = Integral(-_x2*_x6*c/(2*a) + _x2*_x6*d/(2*a), t)
    _x8 = exp(_x2*t*(_x1 - 2)/a)
    _x9 = exp(-_x2*t*(_x1 + 2)/a)
    _x10 = Integral(_x2*_x3*c/(4*a) + _x2*_x3*d/(4*a), t)
    _x11 = Integral(_x2*_x4*c/(4*a) + _x2*_x4*d/(4*a), t)
    sol9_solsimp = [
        Eq(f(t), -C1*_x5 + C2*_x9 + C3*_x8 + _x10*_x8 + _x11*_x9 - _x5*_x7),
        Eq(g(t), -C2*_x1*_x9 + C3*_x1*_x8 + _x1*_x10*_x8 - _x1*_x11*_x9),
        Eq(h(t), C1*_x5 + C2*_x9 + C3*_x8 + _x10*_x8 + _x11*_x9 + _x5*_x7),
    ]
    assert [Eq(s.lhs, _solsimp(s.rhs, t)) for s in sol9] == sol9_solsimp

    # Regression test case for issue #16635
    # https://github.com/sympy/sympy/issues/16635
    eq10 = [Eq(f(t).diff(t), f(t) - g(t) + 15*t - 10), Eq(g(t).diff(t), f(t) - g(t) - 15*t - 5)]
    sol10 = [Eq(f(t), C1 + C2*t + C2 + t*Integral(30*t - 5, t) + Integral(30*t - 5, t) + Integral(-30*t**2 - 10*t - 5, t)),
             Eq(g(t), C1 + C2*t + t*Integral(30*t - 5, t) + Integral(-30*t**2 - 10*t - 5, t))]
    assert dsolve(eq10) == sol10
    assert checksysodesol(eq10, sol10) == (True, [0, 0])

    # Multiple equations
    eq1 = [Eq(Derivative(f(t), t)**2 - 2*Derivative(f(t), t) + 1, 4),
            Eq(-y*f(t) + Derivative(g(t), t), 0)]
    sol1 = [[Eq(f(t), C1 + Integral(-1, t)),
             Eq(g(t), C1*t*y + C2*y + t*y*Integral(-1, t) + y*Integral(t, t))],
            [Eq(f(t), C1 + Integral(3, t)),
             Eq(g(t), C1*t*y + C2*y + t*y*Integral(3, t) + y*Integral(-3*t, t))]]
    assert dsolve(eq1) == sol1
    for sol in sol1:
        assert checksysodesol(eq1, sol) == (True, [0, 0])


def test_sysode_linear_neq_order1_type3():

    f, g, h, k = symbols('f g h k', cls=Function)
    x, t, a = symbols('x t a')
    r = symbols('r', real=True)
    eqs1 = [Eq(diff(f(r), r),  f(r) + r*g(r)),
           Eq(diff(g(r), r),-r*f(r) + g(r))]
    sol1 = [Eq(f(r), (C1*sin(r**2/2) + C2*cos(r**2/2))*exp(r)), Eq(g(r), (C1*cos(r**2/2) - C2*sin(r**2/2))*exp(r))]
    assert dsolve(eqs1) == sol1
    assert checksysodesol(eqs1, sol1) == (True, [0, 0])

    eqs2 = [Eq(diff(f(x), x),  x*f(x) + x**2*g(x)),
           Eq(diff(g(x), x),  2*x**2*f(x) + (x + 3*x**2)*g(x))]
    sol2 = [Eq(f(x), (3*C1/17 + 4*C2/17 + (8*C1 - 12*C2)/(51 + 17*sqrt(17)))*exp(x**2*(3*x + sqrt(17)*x + 3)/6) + (13*C1/51 +
                2*C2/17 + (16*C1 - 24*C2)/(-663 + 153*sqrt(17)))*exp(x**2*(-sqrt(17)*x + 3*x + 3)/6)),
            Eq(g(x), (4*C1/17 - 6*C2/17 + (12*C1 + 16*C2)/(-51 + 17*sqrt(17)))*exp(x**2*(3*x + sqrt(17)*x + 3)/6) + (
                13*C1/17 + 6*C2/17 + (-12*C1 - 16*C2)/(-51 + 17*sqrt(17)))*exp(x**2*(-sqrt(17)*x + 3*x + 3)/6))]

    assert [_simpsol(s) for s in dsolve(eqs2)] == sol2
    assert checksysodesol(eqs2, sol2) == (True, [0, 0])

    eqs3 = [Eq(f(x).diff(x), x*f(x) + g(x)), Eq(g(x).diff(x), -f(x) + x*g(x))]
    sol3 = [Eq(f(x), (C1/2 - I*C2/2)*exp(x**2/2 + I*x) + (C1/2 + I*C2/2)*exp(x**2/2 - I*x)),
            Eq(g(x), (-I*C1/2 + C2/2)*exp(x**2/2 - I*x) + (I*C1/2 + C2/2)*exp(x**2/2 + I*x))]
    assert dsolve(eqs3) == sol3
    assert checksysodesol(eqs3, sol3) == (True, [0, 0])

    eqs4 = [Eq(f(x).diff(x), x*(f(x) + g(x) + h(x))), Eq(g(x).diff(x), x*(f(x) + g(x) + h(x))),
            Eq(h(x).diff(x), x*(f(x) + g(x) + h(x)))]
    sol4 = [Eq(f(x), -C1/3 - C2/3 + 2*C3/3 + (C1/3 + C2/3 + C3/3)*exp(3*x**2/2)),
            Eq(g(x), 2*C1/3 - C2/3 - C3/3 + (C1/3 + C2/3 + C3/3)*exp(3*x**2/2)),
            Eq(h(x), -C1/3 + 2*C2/3 - C3/3 + (C1/3 + C2/3 + C3/3)*exp(3*x**2/2))]
    assert dsolve(eqs4) == sol4
    assert checksysodesol(eqs4, sol4) == (True, [0, 0, 0])

    eqs5 = [Eq(f(x).diff(x), x**2*(f(x) + g(x) + h(x))), Eq(g(x).diff(x), x**2*(f(x) + g(x) + h(x))),
            Eq(h(x).diff(x), x**2*(f(x) + g(x) + h(x)))]
    sol5 = [Eq(f(x), -C1/3 - C2/3 + 2*C3/3 + (C1/3 + C2/3 + C3/3)*exp(x**3)),
            Eq(g(x), 2*C1/3 - C2/3 - C3/3 + (C1/3 + C2/3 + C3/3)*exp(x**3)),
            Eq(h(x), -C1/3 + 2*C2/3 - C3/3 + (C1/3 + C2/3 + C3/3)*exp(x**3))]
    assert dsolve(eqs5) == sol5
    assert checksysodesol(eqs5, sol5) == (True, [0, 0, 0])

    eqs6 = [Eq(Derivative(f(x), x), x*(f(x) + g(x) + h(x) + k(x))),
            Eq(Derivative(g(x), x), x*(f(x) + g(x) + h(x) + k(x))),
            Eq(Derivative(h(x), x), x*(f(x) + g(x) + h(x) + k(x))),
            Eq(Derivative(k(x), x), x*(f(x) + g(x) + h(x) + k(x)))]
    sol6 = [Eq(f(x), -C1/4 - C2/4 - C3/4 + 3*C4/4 + (C1/4 + C2/4 + C3/4 + C4/4)*exp(2*x**2)),
            Eq(g(x), 3*C1/4 - C2/4 - C3/4 - C4/4 + (C1/4 + C2/4 + C3/4 + C4/4)*exp(2*x**2)),
            Eq(h(x), -C1/4 + 3*C2/4 - C3/4 - C4/4 + (C1/4 + C2/4 + C3/4 + C4/4)*exp(2*x**2)),
            Eq(k(x), -C1/4 - C2/4 + 3*C3/4 - C4/4 + (C1/4 + C2/4 + C3/4 + C4/4)*exp(2*x**2))]
    assert dsolve(eqs6) == sol6
    assert checksysodesol(eqs6, sol6) == (True, [0, 0, 0, 0])

    y = symbols("y", real=True)

    eqs7 = [Eq(Derivative(f(y), y), y*f(y) + g(y)), Eq(Derivative(g(y), y), y*g(y) - f(y))]
    sol7 = [Eq(f(y), (C1*sin(y) + C2*cos(y))*exp(y**2/2)), Eq(g(y), (C1*cos(y) - C2*sin(y))*exp(y**2/2))]
    assert dsolve(eqs7) == sol7
    assert checksysodesol(eqs7, sol7) == (True, [0, 0])


def test_sysode_linear_neq_order1_type4():

    f, g, h, k = symbols('f g h k', cls=Function)
    x, t, a = symbols('x t a')
    r = symbols('r', real=True)

    eqs1 = [Eq(diff(f(r), r), f(r) + r*g(r) + r**2), Eq(diff(g(r), r), -r*f(r) + g(r) + r)]
    sol1 = [Eq(f(r), C1*exp(r)*sin(r**2/2) + C2*exp(r)*cos(r**2/2) + exp(r)*sin(r**2/2)*Integral(r**2*exp(-r)*sin(r**2/2) +
                r*exp(-r)*cos(r**2/2), r) + exp(r)*cos(r**2/2)*Integral(r**2*exp(-r)*cos(r**2/2) - r*exp(-r)*sin(r**2/2), r)),
            Eq(g(r), C1*exp(r)*cos(r**2/2) - C2*exp(r)*sin(r**2/2) - exp(r)*sin(r**2/2)*Integral(r**2*exp(-r)*cos(r**2/2) -
                r*exp(-r)*sin(r**2/2), r) + exp(r)*cos(r**2/2)*Integral(r**2*exp(-r)*sin(r**2/2) + r*exp(-r)*cos(r**2/2), r))]
    assert dsolve(eqs1) == sol1
    assert checksysodesol(eqs1, sol1) == (True, [0, 0])

    eqs2 = [Eq(diff(f(r), r), f(r) + r*g(r) + r), Eq(diff(g(r), r), -r*f(r) + g(r) + log(r))]
    sol2 = [Eq(f(r), C1*exp(r)*sin(r**2/2) + C2*exp(r)*cos(r**2/2) + exp(r)*sin(r**2/2)*Integral(r*exp(-r)*sin(r**2/2) +
                exp(-r)*log(r)*cos(r**2/2), r) + exp(r)*cos(r**2/2)*Integral(r*exp(-r)*cos(r**2/2) - exp(-r)*log(r)*sin(
                r**2/2), r)),
            Eq(g(r), C1*exp(r)*cos(r**2/2) - C2*exp(r)*sin(r**2/2) - exp(r)*sin(r**2/2)*Integral(r*exp(-r)*cos(r**2/2) -
                exp(-r)*log(r)*sin(r**2/2), r) + exp(r)*cos(r**2/2)*Integral(r*exp(-r)*sin(r**2/2) + exp(-r)*log(r)*cos(
                r**2/2), r))]
    assert dsolve(eqs2) == sol2
    assert checksysodesol(eqs2, sol2) == (True, [0, 0])

    eqs3 = [Eq(f(x).diff(x), x*(f(x) + g(x) + h(x)) + x), Eq(g(x).diff(x), x*(f(x) + g(x) + h(x)) + x),
            Eq(h(x).diff(x), x*(f(x) + g(x) + h(x)) + 1)]
    _x1 = exp(-3*x**2/2)
    _x2 = exp(3*x**2/2)
    _x3 = Integral(2*_x1*x/3 + _x1/3 + x/3 - Rational(1, 3), x)
    _x4 = Integral(2*_x1*x/3 + _x1/3 - 2*x/3 + Rational(2, 3), x)
    sol3 = [
        Eq(f(x),
           C1*_x2/3 - C1/3 + C2*_x2/3 - C2/3 + C3*_x2/3 + 2*C3/3 + 2*_x2*_x3/3 + _x2*_x4/3 + _x3/3 - _x4/3),
        Eq(g(x),
           C1*_x2/3 + 2*C1/3 + C2*_x2/3 - C2/3 + C3*_x2/3 - C3/3 + 2*_x2*_x3/3 + _x2*_x4/3 + _x3/3 - _x4/3),
        Eq(h(x),
           C1*_x2/3 - C1/3 + C2*_x2/3 + 2*C2/3 + C3*_x2/3 - C3/3 + 2*_x2*_x3/3 + _x2*_x4/3 - 2*_x3/3 + 2*_x4/3),
    ]
    assert dsolve(eqs3) == sol3
    assert checksysodesol(eqs3, sol3) == (True, [0, 0, 0])

    eqs4 = [Eq(f(x).diff(x), x*(f(x) + g(x) + h(x)) + sin(x)), Eq(g(x).diff(x), x*(f(x) + g(x) + h(x)) + sin(x)),
            Eq(h(x).diff(x), x*(f(x) + g(x) + h(x)) + sin(x))]
    sol4 = [Eq(f(x), C1*exp(3*x**2/2)/3 - C1/3 + C2*exp(3*x**2/2)/3 - C2/3 + C3*exp(3*x**2/2)/3 + 2*C3/3 +
                exp(3*x**2/2)*Integral(exp(-3*x**2/2)*sin(x), x)),
            Eq(g(x), C1*exp(3*x**2/2)/3 + 2*C1/3 + C2*exp(3*x**2/2)/3 - C2/3 +
                C3*exp(3*x**2/2)/3 - C3/3 + exp(3*x**2/2)*Integral(exp(-3*x**2/2)*sin(x), x)),
            Eq(h(x), C1*exp(3*x**2/2)/3 - C1/3 + C2*exp(3*x**2/2)/3 + 2*C2/3 + C3*exp(3*x**2/2)/3 - C3/3 + exp(
                3*x**2/2)*Integral(exp(-3*x**2/2)*sin(x), x))]
    assert dsolve(eqs4) == sol4
    assert checksysodesol(eqs4, sol4) == (True, [0, 0, 0])

    eqs5 = [Eq(Derivative(f(x), x), x*(f(x) + g(x) + h(x) + k(x) + 1)), Eq(Derivative(g(x), x), x*(f(x) + g(x)
            + h(x) + k(x) + 1)), Eq(Derivative(h(x), x), x*(f(x) + g(x) + h(x) + k(x) + 1)), Eq(Derivative(k(x),
            x), x*(f(x) + g(x) + h(x) + k(x) + 1))]
    sol5 = [Eq(f(x), C1*exp(2*x**2)/4 - C1/4 + C2*exp(2*x**2)/4 - C2/4 + C3*exp(2*x**2)/4 - C3/4 + C4*exp(2*x**2)/4 + 3*C4/4 +
                exp(2*x**2)*Integral(x*exp(-2*x**2), x)),
            Eq(g(x), C1*exp(2*x**2)/4 + 3*C1/4 + C2*exp(2*x**2)/4 - C2/4 +
                C3*exp(2*x**2)/4 - C3/4 + C4*exp(2*x**2)/4 - C4/4 + exp(2*x**2)*Integral(x*exp(-2*x**2), x)),
            Eq(h(x), C1*exp(2*x**2)/4 - C1/4 + C2*exp(2*x**2)/4 + 3*C2/4 + C3*exp(2*x**2)/4 - C3/4 + C4*exp(2*x**2)/4 - C4/4 +
                exp(2*x**2)*Integral(x*exp(-2*x**2), x)),
            Eq(k(x), C1*exp(2*x**2)/4 - C1/4 + C2*exp(2*x**2)/4 - C2/4 + C3*exp(2*x**2)/4
                + 3*C3/4 + C4*exp(2*x**2)/4 - C4/4 + exp(2*x**2)*Integral(x*exp(-2*x**2), x))]
    assert dsolve(eqs5) == sol5
    assert checksysodesol(eqs5, sol5) == (True, [0, 0, 0, 0])

    eqs6 = [Eq(Derivative(f(x), x), x**2*(f(x) + g(x) + h(x) + k(x) + 1)), Eq(Derivative(g(x), x), x**2*(f(x) +
            g(x) + h(x) + k(x) + 1)), Eq(Derivative(h(x), x), x**2*(f(x) + g(x) + h(x) + k(x) + 1)),
            Eq(Derivative(k(x), x), x**2*(f(x) + g(x) + h(x) + k(x) + 1))]
    sol6 = [Eq(f(x), C1*exp(4*x**3/3)/4 - C1/4 + C2*exp(4*x**3/3)/4 - C2/4 + C3*exp(4*x**3/3)/4 - C3/4 + C4*exp(4*x**3/3)/4 +
                3*C4/4 + exp(4*x**3/3)*Integral(x**2*exp(-4*x**3/3), x)),
            Eq(g(x), C1*exp(4*x**3/3)/4 + 3*C1/4 + C2*exp(4*x**3/3)/4 - C2/4 + C3*exp(4*x**3/3)/4 - C3/4 + C4*exp(
                4*x**3/3)/4 - C4/4 + exp(4*x**3/3)*Integral(x**2*exp(-4*x**3/3), x)),
            Eq(h(x), C1*exp(4*x**3/3)/4 - C1/4 + C2*exp(4*x**3/3)/4 + 3*C2/4 + C3*exp(4*x**3/3)/4 - C3/4 + C4*exp(
                4*x**3/3)/4 - C4/4 + exp(4*x**3/3)*Integral(x**2*exp(-4*x**3/3), x)),
            Eq(k(x), C1*exp(4*x**3/3)/4 - C1/4 + C2*exp(4*x**3/3)/4 - C2/4 +
                C3*exp(4*x**3/3)/4 + 3*C3/4 + C4*exp(4*x**3/3)/4 - C4/4 + exp(4*x**3/3)*Integral(x**2*exp(-4*x**3/3), x))]
    assert dsolve(eqs6) == sol6
    assert checksysodesol(eqs6, sol6) == (True, [0, 0, 0, 0])

    eqs7 = [Eq(Derivative(f(x), x), (f(x) + g(x) + h(x))*log(x) + sin(x)), Eq(Derivative(g(x), x), (f(x) + g(x)
            + h(x))*log(x) + sin(x)), Eq(Derivative(h(x), x), (f(x) + g(x) + h(x))*log(x) + sin(x))]
    sol7 = [Eq(f(x), -C1/3 - C2/3 + 2*C3/3 + (C1/3 + C2/3 + C3/3 + Integral(exp(-3*x*log(x) + 3*x)*sin(x), x))*exp(3*x*log(x) -
                3*x)),
            Eq(g(x), 2*C1/3 - C2/3 - C3/3 + (C1/3 + C2/3 + C3/3 + Integral(exp(-3*x*log(x) + 3*x)*sin(x), x))*exp(3*x*log(x)
                - 3*x)),
            Eq(h(x), -C1/3 + 2*C2/3 - C3/3 + (C1/3 + C2/3 + C3/3 + Integral(exp(-3*x*log(x) + 3*x)*sin(x),
                x))*exp(3*x*log(x) - 3*x))]
    with dotprodsimp(True):
        assert dsolve(eqs7) == sol7
    assert checksysodesol(eqs7, sol7) == (True, [0, 0, 0])

    eqs8 = [Eq(Derivative(f(x), x), (f(x) + g(x) + h(x) + k(x))*log(x) + sin(x)), Eq(Derivative(g(x), x), (f(x)
            + g(x) + h(x) + k(x))*log(x) + sin(x)), Eq(Derivative(h(x), x), (f(x) + g(x) + h(x) + k(x))*log(x) +
            sin(x)), Eq(Derivative(k(x), x), (f(x) + g(x) + h(x) + k(x))*log(x) + sin(x))]
    sol8 = [Eq(f(x), -C1/4 - C2/4 - C3/4 + 3*C4/4 + (C1/4 + C2/4 + C3/4 + C4/4 + Integral(exp(-4*x*log(x) + 4*x)*sin(x),
                x))*exp(4*x*log(x) - 4*x)),
            Eq(g(x), 3*C1/4 - C2/4 - C3/4 - C4/4 + (C1/4 + C2/4 + C3/4 + C4/4 + Integral(exp(-4*x*log(x)
            + 4*x)*sin(x), x))*exp(4*x*log(x) - 4*x)),
            Eq(h(x), -C1/4 + 3*C2/4 - C3/4 - C4/4 + (C1/4 + C2/4 + C3/4 + C4/4 +
                Integral(exp(-4*x*log(x) + 4*x)*sin(x), x))*exp(4*x*log(x) - 4*x)),
            Eq(k(x), -C1/4 - C2/4 + 3*C3/4 - C4/4 + (C1/4 + C2/4
                + C3/4 + C4/4 + Integral(exp(-4*x*log(x) + 4*x)*sin(x), x))*exp(4*x*log(x) - 4*x))]
    with dotprodsimp(True):
        assert dsolve(eqs8) == sol8
    assert checksysodesol(eqs8, sol8) == (True, [0, 0, 0, 0])


def test_higher_order_to_first_order():
    f, g = symbols('f g', cls=Function)
    x = symbols('x')

    eqs1 = [Eq(f(x).diff(x, 2), 2*f(x) + g(x)), Eq(g(x).diff(x, 2), -f(x))]
    sol1 = [Eq(f(x), C3*x*exp(x) + C3*exp(x) + C4*exp(x) + (-C1 - C2*x + C2)*exp(-x)),
            Eq(g(x), -C3*x*exp(x) + C3*exp(x) - C4*exp(x) + (C1 + C2*x + C2)*exp(-x))]
    assert dsolve(eqs1) == sol1
    assert checksysodesol(eqs1, sol1) == (True, [0, 0])

    eqs2 = [Eq(f(x).diff(x, 2), 0), Eq(g(x).diff(x, 2), f(x))]
    sol2 = [Eq(f(x), C1 + C2*x), Eq(g(x), C1*x**2/2 + C2*x**3/6 + C3 + C4*x)]
    assert dsolve(eqs2) == sol2
    assert checksysodesol(eqs2, sol2) == (True, [0, 0])

    eqs3 = [Eq(f(x).diff(x, 2), 2*f(x)), Eq(g(x).diff(x, 2), -f(x) + 2*g(x))]
    sol3 = [Eq(f(x), 4*C1*exp(-sqrt(2)*x) + 4*C2*exp(sqrt(2)*x)),
            Eq(g(x), C2*exp(sqrt(2)*x) + (C1 + sqrt(2)*(C1*x + C4))*exp(-sqrt(2)*x) - sqrt(2)*(C2*x*exp(sqrt(2)*x) +
                C3*exp(sqrt(2)*x)))]
    assert dsolve(eqs3) == sol3
    assert checksysodesol(eqs3, sol3) == (True, [0, 0])

    eqs4 = [Eq(f(x).diff(x, 2), 2*f(x) + g(x)), Eq(g(x).diff(x, 2), + 2*g(x))]
    sol4 = [Eq(f(x), C1*x*exp(sqrt(2)*x)/4 - sqrt(2)*C1*exp(sqrt(2)*x)/8 + C4*exp(sqrt(2)*x)/4 + (C2/4 + C3*x/4 +
                sqrt(2)*C3/8)*exp(-sqrt(2)*x)),
            Eq(g(x), sqrt(2)*C1*exp(sqrt(2)*x)/2 - sqrt(2)*C3*exp(-sqrt(2)*x)/2)]
    assert dsolve(eqs4) == sol4
    assert checksysodesol(eqs4, sol4) == (True, [0, 0])

    eqs5 = [Eq(f(x).diff(x, 2), f(x)), Eq(g(x).diff(x, 2), f(x))]
    sol5 = [Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), -C1*exp(-x) + C2*exp(x) + C3 + C4*x)]
    assert dsolve(eqs5) == sol5
    assert checksysodesol(eqs5, sol5) == (True, [0, 0])

    eqs6 = [Eq(f(x).diff(x, 2), f(x) + g(x)), Eq(g(x).diff(x, 2), -f(x) - g(x))]
    sol6 = [Eq(f(x), C1 + C2*x**2/2 + C2 + C3*x + C4*x**3/6 + C4*x), Eq(g(x), -C1 - C2*x**2/2 - C3*x - C4*x**3/6)]
    assert dsolve(eqs6) == sol6
    assert checksysodesol(eqs6, sol6) == (True, [0, 0])

    eqs7 = [Eq(f(x).diff(x, 2), f(x) + g(x) + 1), Eq(g(x).diff(x, 2), f(x) + g(x) + 1)]
    sol7 = [Eq(f(x), -C1 - C2*x + sqrt(2)*C3*exp(sqrt(2)*x)/2 - x*Integral(0, x) + (-sqrt(2)*C4/2 -
                sqrt(2)*Integral(exp(sqrt(2)*x)/2, x)/2)*exp(-sqrt(2)*x) + sqrt(2)*exp(sqrt(2)*x)*
                Integral(exp(-sqrt(2)*x)/2, x)/2 - Integral(0, x)),
            Eq(g(x), C1 + C2*x + sqrt(2)*C3*exp(sqrt(2)*x)/2 + x*Integral(0, x) + (-sqrt(2)*C4/2 - sqrt(2)*
                Integral(exp(sqrt(2)*x)/2, x)/2)*exp(-sqrt(2)*x) + sqrt(2)*exp(sqrt(2)*x)*Integral(
                exp(-sqrt(2)*x)/2, x)/2 + Integral(0, x))]
    assert dsolve(eqs7) == sol7
    assert checksysodesol(eqs7, sol7) == (True, [0, 0])

    eqs8 = [Eq(f(x).diff(x, 2), f(x) + g(x) + 1), Eq(g(x).diff(x, 2), -f(x) - g(x) + 1)]
    sol8 = [Eq(f(x), C1 + C2*x**2/2 + C2 + C3*x + C4*x**3/6 + C4*x + x**3*Integral(2, x)/6 +
                x**2*Integral(-2*x, x)/2 + x*Integral(2, x) + x*Integral(x**2 - 1, x) + Integral(-2*x, x) +
                Integral(-x**3/3 + x, x)),
            Eq(g(x), -C1 - C2*x**2/2 - C3*x - C4*x**3/6 - x**3*Integral(2, x)/6 -
                x**2*Integral(-2*x, x)/2 - x*Integral(x**2 - 1, x) - Integral(-x**3/3 + x, x))]
    assert dsolve(eqs8) == sol8
    assert checksysodesol(eqs8, sol8) == (True, [0, 0])

    eqs9 = [f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) + g(x) - 2*exp(I*x),
            g(x).diff(x, 2) + 2*g(x).diff(x) + f(x) + g(x) - 2*exp(I*x)]
    sol9 = [Eq(f(x), -C1 + (C2/2 + Integral(0, x)/2)*exp(-2*x) + (C3*sin(x)/2 - C3*cos(x)/2 + C4*sin(x)/2 + C4*cos(x)/2 +
                sin(x)*Integral(-2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)*cos(x), x)/2 + sin(x)*Integral(-2*exp(x)*
                exp(I*x)*sin(x)**2/cos(x) + 2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)/cos(x), x)/2 + cos(x)*
                Integral(-2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)*cos(x), x)/2 - cos(x)*Integral(-2*exp(x)*exp(I*x)*
                sin(x)**2/cos(x) + 2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)/cos(x), x)/2)*exp(-x) - Integral(0, x)),
            Eq(g(x), C1 + (-C2/2 - Integral(0, x)/2)*exp(-2*x) + (C3*sin(x)/2 - C3*cos(x)/2 + C4*sin(x)/2 + C4*cos(x)/2 +
                sin(x)*Integral(-2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)*cos(x), x)/2 + sin(x)*Integral(-2*exp(x)*
                exp(I*x)*sin(x)**2/cos(x) + 2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)/cos(x), x)/2 + cos(x)*
                Integral(-2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)*cos(x), x)/2 - cos(x)*Integral(-2*exp(x)*exp(I*x)*
                sin(x)**2/cos(x) + 2*exp(x)*exp(I*x)*sin(x) + 2*exp(x)*exp(I*x)/cos(x), x)/2)*exp(-x) + Integral(0, x))]
    assert dsolve(eqs9) == sol9
    assert checksysodesol(eqs9, sol9) == (True, [0, 0])

    x, y = symbols('x, y', cls=Function)
    t, l = symbols('t, l')

    eq1 = (Eq(diff(x(t),t,t), 5*x(t) + 43*y(t)), Eq(diff(y(t),t,t), x(t) + 9*y(t)))
    sol1 = [Eq(x(t), 43*C1*exp(-t*sqrt(7 - sqrt(47)))/((-2 + sqrt(47))*sqrt(7 - sqrt(47))) + 43*C2*exp(t*sqrt(7
                - sqrt(47)))/((2 - sqrt(47))*sqrt(7 - sqrt(47))) - 43*C3*exp(-t*sqrt(sqrt(47) + 7))/((2 +
                sqrt(47))*sqrt(sqrt(47) + 7)) + 43*C4*exp(t*sqrt(sqrt(47) + 7))/((2 + sqrt(47))*sqrt(sqrt(47) + 7))),
            Eq(y(t), -C1*exp(-t*sqrt(7 - sqrt(47)))/sqrt(7 - sqrt(47)) + C2*exp(t*sqrt(7 - sqrt(47)))/sqrt(7 -
                sqrt(47)) - C3*exp(-t*sqrt(sqrt(47) + 7))/sqrt(sqrt(47) + 7) + C4*exp(t*sqrt(sqrt(47) +
                7))/sqrt(sqrt(47) + 7))]
    assert dsolve(eq1) == sol1
    assert checksysodesol(eq1, sol1) == (True, [0, 0])

    # Note: Solution to be updated after _extract_funcs is fixed
    eq3 = (Eq(diff(x(t),t,t), t*(4*diff(x(t),t) + 9*diff(y(t),t))), Eq(diff(y(t),t,t), t*(12*diff(x(t),t) - 6*diff(y(t),t))))
    # sol3 = [Eq(x(t), C1 + Integral(9*sqrt(133)*C2*exp(-t**2/2 + sqrt(133)*t**2/2)/266 - 9*sqrt(133)*C2*exp(-sqrt(133)*t**2/2
    #             - t**2/2)/266 + 27*sqrt(133)*C3*t**2*exp(-t**2/2 + sqrt(133)*t**2/2)/(133*(-sqrt(133)*t**2/2 - 5*t**2/2))
    #             - 27*sqrt(133)*C3*t**2*exp(-sqrt(133)*t**2/2 - t**2/2)/(133*(-sqrt(133)*t**2/2 - 5*t**2/2)) +
    #             C3*exp(-t**2/2 + sqrt(133)*t**2/2), t)),
    #         Eq(y(t), C4 + Integral(-27*sqrt(133)*C2*t**2*exp(-t**2/2 + sqrt(133)*t**2/2)/(133*(-sqrt(133)*t**2/2 -
    #             5*t**2/2)) + 27*sqrt(133)*C2*t**2*exp(-sqrt(133)*t**2/2 - t**2/2)/(133*(-5*t**2/2 + sqrt(133)*t**2/2)) -
    #             162*sqrt(133)*C3*t**4*exp(-t**2/2 + sqrt(133)*t**2/2)/(133*(-sqrt(133)*t**2/2 - 5*t**2/2)**2) -
    #             6*C3*t**2*exp(-t**2/2 + sqrt(133)*t**2/2)/(-sqrt(133)*t**2/2 - 5*t**2/2) -
    #             6*sqrt(133)*C3*exp(-sqrt(133)*t**2/2 - t**2/2)/133, t))]
    # assert dsolve(eq3) == sol3
    sol3 = dsolve(eq3)
    assert checksysodesol(eq3, sol3) == (True, [0, 0])


def test_component_division():
    f, g, h, k = symbols('f g h k', cls=Function)
    x = symbols("x")
    funcs = [f(x), g(x), h(x), k(x)]

    eqs1 = [Eq(Derivative(f(x), x), 2*f(x)),
            Eq(Derivative(g(x), x), f(x)),
            Eq(Derivative(h(x), x), h(x)),
            Eq(Derivative(k(x), x), h(x)**4 + k(x))]
    sol1 = [Eq(f(x), 2*C1*exp(2*x)),
            Eq(g(x), C1*exp(2*x) + C2),
            Eq(h(x), C3*exp(x)),
            Eq(k(x), (C4 + Integral(C3**4*exp(3*x), x))*exp(x))]
    components1 = {((Eq(Derivative(f(x), x), 2*f(x)),), (Eq(Derivative(g(x), x), f(x)),)),
                   ((Eq(Derivative(h(x), x), h(x)),), (Eq(Derivative(k(x), x), h(x)**4 + k(x)),))}
    eqsdict1 = ({f(x): set(), g(x): {f(x)}, h(x): set(), k(x): {h(x)}},
                {f(x): Eq(Derivative(f(x), x), 2*f(x)),
                g(x): Eq(Derivative(g(x), x), f(x)),
                h(x): Eq(Derivative(h(x), x), h(x)),
                k(x): Eq(Derivative(k(x), x), h(x)**4 + k(x))})
    graph1 = [{f(x), g(x), h(x), k(x)}, {(g(x), f(x)), (k(x), h(x))}]
    assert {tuple(tuple(scc) for scc in wcc) for wcc in _component_division(eqs1, funcs, x)} == components1
    assert _eqs2dict(eqs1, funcs) == eqsdict1
    assert [set(element) for element in _dict2graph(eqsdict1[0])] == graph1
    assert dsolve(eqs1) == sol1
    assert checksysodesol(eqs1, sol1) == (True, [0, 0, 0, 0])

    eqs2 = [Eq(Derivative(f(x), x), 2*f(x)),
            Eq(Derivative(g(x), x), f(x)),
            Eq(Derivative(h(x), x), h(x)),
            Eq(Derivative(k(x), x), f(x)**4 + k(x))]
    sol2 = [Eq(f(x), C1*exp(2*x)),
            Eq(g(x), C2 + Integral(C1*exp(2*x), x)),
            Eq(h(x), C3*exp(x)),
            Eq(k(x), (C4 + Integral(C1**4*exp(7*x), x))*exp(x))]
    components2 = {frozenset([(Eq(Derivative(f(x), x), 2*f(x)),),
                    (Eq(Derivative(g(x), x), f(x)),),
                    (Eq(Derivative(k(x), x), f(x)**4 + k(x)),)]),
                   frozenset([(Eq(Derivative(h(x), x), h(x)),)])}
    eqsdict2 = ({f(x): set(), g(x): {f(x)}, h(x): set(), k(x): {f(x)}},
                 {f(x): Eq(Derivative(f(x), x), 2*f(x)),
                  g(x): Eq(Derivative(g(x), x), f(x)),
                  h(x): Eq(Derivative(h(x), x), h(x)),
                  k(x): Eq(Derivative(k(x), x), f(x)**4 + k(x))})
    graph2 = [{f(x), g(x), h(x), k(x)}, {(g(x), f(x)), (k(x), f(x))}]
    assert {frozenset(tuple(scc) for scc in wcc) for wcc in _component_division(eqs2, funcs, x)} == components2
    assert _eqs2dict(eqs2, funcs) == eqsdict2
    assert [set(element) for element in _dict2graph(eqsdict2[0])] == graph2
    assert dsolve(eqs2) == sol2
    assert checksysodesol(eqs2, sol2) == (True, [0, 0, 0, 0])

    eqs3 = [Eq(Derivative(f(x), x), 2*f(x)),
            Eq(Derivative(g(x), x), f(x) + x),
            Eq(Derivative(h(x), x), h(x)),
            Eq(Derivative(k(x), x), f(x)**4 + k(x))]
    sol3 = [Eq(f(x), C1*exp(2*x)),
            Eq(g(x), C2 + Integral(C1*exp(2*x) + x, x)),
            Eq(h(x), C3*exp(x)),
            Eq(k(x), (C4 + Integral(C1**4*exp(7*x), x))*exp(x))]
    components3 = {frozenset([(Eq(Derivative(f(x), x), 2*f(x)),),
                    (Eq(Derivative(g(x), x), x + f(x)),),
                    (Eq(Derivative(k(x), x), f(x)**4 + k(x)),)]),
                    frozenset([(Eq(Derivative(h(x), x), h(x)),),])}
    eqsdict3 = ({f(x): set(), g(x): {f(x)}, h(x): set(), k(x): {f(x)}},
                {f(x): Eq(Derivative(f(x), x), 2*f(x)),
                g(x): Eq(Derivative(g(x), x), x + f(x)),
                h(x): Eq(Derivative(h(x), x), h(x)),
                k(x): Eq(Derivative(k(x), x), f(x)**4 + k(x))})
    graph3 = [{f(x), g(x), h(x), k(x)}, {(g(x), f(x)), (k(x), f(x))}]
    assert {frozenset(tuple(scc) for scc in wcc) for wcc in _component_division(eqs3, funcs, x)} == components3
    assert _eqs2dict(eqs3, funcs) == eqsdict3
    assert [set(l) for l in _dict2graph(eqsdict3[0])] == graph3
    assert dsolve(eqs3) == sol3
    assert checksysodesol(eqs3, sol3) == (True, [0, 0, 0, 0])

    eqs4 = [Eq(Derivative(f(x), x), x*f(x) + 2*g(x)),
            Eq(Derivative(g(x), x), f(x) + x*g(x) + x),
            Eq(Derivative(h(x), x), h(x)),
            Eq(Derivative(k(x), x), f(x)**4 + k(x))]
    sol4 = [Eq(f(x), (C1/2 - sqrt(2)*C2/2 - sqrt(2)*Integral(x*exp(-x**2/2 - sqrt(2)*x)/2 + x*exp(-x**2/2 +\
                sqrt(2)*x)/2, x)/2 + Integral(sqrt(2)*x*exp(-x**2/2 - sqrt(2)*x)/2 - sqrt(2)*x*exp(-x**2/2 +\
                sqrt(2)*x)/2, x)/2)*exp(x**2/2 - sqrt(2)*x) + (C1/2 + sqrt(2)*C2/2 + sqrt(2)*Integral(x*exp(-x**2/2
                - sqrt(2)*x)/2 + x*exp(-x**2/2 + sqrt(2)*x)/2, x)/2 + Integral(sqrt(2)*x*exp(-x**2/2 - sqrt(2)*x)/2
                - sqrt(2)*x*exp(-x**2/2 + sqrt(2)*x)/2, x)/2)*exp(x**2/2 + sqrt(2)*x)),
            Eq(g(x), (-sqrt(2)*C1/4 + C2/2 + Integral(x*exp(-x**2/2 - sqrt(2)*x)/2 + x*exp(-x**2/2 + sqrt(2)*x)/2, x)/2 -\
                sqrt(2)*Integral(sqrt(2)*x*exp(-x**2/2 - sqrt(2)*x)/2 - sqrt(2)*x*exp(-x**2/2 + sqrt(2)*x)/2,
                x)/4)*exp(x**2/2 - sqrt(2)*x) + (sqrt(2)*C1/4 + C2/2 + Integral(x*exp(-x**2/2 - sqrt(2)*x)/2 +
                x*exp(-x**2/2 + sqrt(2)*x)/2, x)/2 + sqrt(2)*Integral(sqrt(2)*x*exp(-x**2/2 - sqrt(2)*x)/2 -
                sqrt(2)*x*exp(-x**2/2 + sqrt(2)*x)/2, x)/4)*exp(x**2/2 + sqrt(2)*x)),
            Eq(h(x), C3*exp(x)),
            Eq(k(x), C4*exp(x) + exp(x)*Integral((C1*exp(x**2/2 - sqrt(2)*x)/2 + C1*exp(x**2/2 + sqrt(2)*x)/2 -
                sqrt(2)*C2*exp(x**2/2 - sqrt(2)*x)/2 + sqrt(2)*C2*exp(x**2/2 + sqrt(2)*x)/2 - sqrt(2)*exp(x**2/2 -
                sqrt(2)*x)*Integral(x*exp(-x**2/2 - sqrt(2)*x)/2 + x*exp(-x**2/2 + sqrt(2)*x)/2, x)/2 + exp(x**2/2 -
                sqrt(2)*x)*Integral(sqrt(2)*x*exp(-x**2/2 - sqrt(2)*x)/2 - sqrt(2)*x*exp(-x**2/2 + sqrt(2)*x)/2,
                x)/2 + sqrt(2)*exp(x**2/2 + sqrt(2)*x)*Integral(x*exp(-x**2/2 - sqrt(2)*x)/2 + x*exp(-x**2/2 +
                sqrt(2)*x)/2, x)/2 + exp(x**2/2 + sqrt(2)*x)*Integral(sqrt(2)*x*exp(-x**2/2 - sqrt(2)*x)/2 -
                sqrt(2)*x*exp(-x**2/2 + sqrt(2)*x)/2, x)/2)**4*exp(-x), x))]
    components4 = {(frozenset([Eq(Derivative(f(x), x), x*f(x) + 2*g(x)),
                    Eq(Derivative(g(x), x), x*g(x) + x + f(x))]),
                    frozenset([Eq(Derivative(k(x), x), f(x)**4 + k(x)),])),
                    (frozenset([Eq(Derivative(h(x), x), h(x)),]),)}
    eqsdict4 = ({f(x): {g(x)}, g(x): {f(x)}, h(x): set(), k(x): {f(x)}},
                {f(x): Eq(Derivative(f(x), x), x*f(x) + 2*g(x)),
                g(x): Eq(Derivative(g(x), x), x*g(x) + x + f(x)),
                h(x): Eq(Derivative(h(x), x), h(x)),
                k(x): Eq(Derivative(k(x), x), f(x)**4 + k(x))})
    graph4 = [{f(x), g(x), h(x), k(x)}, {(f(x), g(x)), (g(x), f(x)), (k(x), f(x))}]
    assert {tuple(frozenset(scc) for scc in wcc) for wcc in _component_division(eqs4, funcs, x)} == components4
    assert _eqs2dict(eqs4, funcs) == eqsdict4
    assert [set(element) for element in _dict2graph(eqsdict4[0])] == graph4
    assert dsolve(eqs4) == sol4
    assert checksysodesol(eqs4, sol4) == (True, [0, 0, 0, 0])

    eqs5 = [Eq(Derivative(f(x), x), x*f(x) + 2*g(x)),
            Eq(Derivative(g(x), x), x*g(x) + f(x)),
            Eq(Derivative(h(x), x), h(x)),
            Eq(Derivative(k(x), x), f(x)**4 + k(x))]
    sol5 = [Eq(f(x), (C1/2 - sqrt(2)*C2/2)*exp(x**2/2 - sqrt(2)*x) + (C1/2 + sqrt(2)*C2/2)*exp(x**2/2 + sqrt(2)*x)),
            Eq(g(x), (-sqrt(2)*C1/4 + C2/2)*exp(x**2/2 - sqrt(2)*x) + (sqrt(2)*C1/4 + C2/2)*exp(x**2/2 + sqrt(2)*x)),
            Eq(h(x), C3*exp(x)),
            Eq(k(x), C4*exp(x) + exp(x)*Integral((C1*exp(x**2/2 - sqrt(2)*x)/2 + C1*exp(x**2/2 + sqrt(2)*x)/2 -
                sqrt(2)*C2*exp(x**2/2 - sqrt(2)*x)/2 + sqrt(2)*C2*exp(x**2/2 + sqrt(2)*x)/2)**4*exp(-x), x))]
    components5 = {(frozenset([Eq(Derivative(f(x), x), x*f(x) + 2*g(x)),
                    Eq(Derivative(g(x), x), x*g(x) + f(x))]),
                    frozenset([Eq(Derivative(k(x), x), f(x)**4 + k(x)),])),
                    (frozenset([Eq(Derivative(h(x), x), h(x)),]),)}
    eqsdict5 = ({f(x): {g(x)}, g(x): {f(x)}, h(x): set(), k(x): {f(x)}},
                {f(x): Eq(Derivative(f(x), x), x*f(x) + 2*g(x)),
                g(x): Eq(Derivative(g(x), x), x*g(x) + f(x)),
                h(x): Eq(Derivative(h(x), x), h(x)),
                k(x): Eq(Derivative(k(x), x), f(x)**4 + k(x))})
    graph5 = [{f(x), g(x), h(x), k(x)}, {(f(x), g(x)), (g(x), f(x)), (k(x), f(x))}]
    assert {tuple(frozenset(scc) for scc in wcc) for wcc in _component_division(eqs5, funcs, x)} == components5
    assert _eqs2dict(eqs5, funcs) == eqsdict5
    assert [set(element) for element in _dict2graph(eqsdict5[0])] == graph5
    assert dsolve(eqs5) == sol5
    assert checksysodesol(eqs5, sol5) == (True, [0, 0, 0, 0])


def test_linodesolve():
    t, x, a = symbols("t x a")
    f, g, h = symbols("f g h", cls=Function)

    # Testing the Errors
    raises(ValueError, lambda: linodesolve(1, t))
    raises(ValueError, lambda: linodesolve(a, t))

    A1 = Matrix([[1, 2], [2, 4], [4, 6]])
    raises(NonSquareMatrixError, lambda: linodesolve(A1, t))

    A2 = Matrix([[1, 2, 1], [3, 1, 2]])
    raises(NonSquareMatrixError, lambda: linodesolve(A2, t))

    # Testing auto functionality
    func = [f(t), g(t)]
    eq = [Eq(f(t).diff(t) + g(t).diff(t), g(t)), Eq(g(t).diff(t), f(t))]
    ceq = canonical_odes(eq, func, t)
    (A1, A0), b = linear_ode_to_matrix(ceq[0], func, t, 1)
    A = A0
    sol = [C1*(-Rational(1, 2) + sqrt(5)/2)*exp(t*(-Rational(1, 2) + sqrt(5)/2)) + C2*(-sqrt(5)/2 - Rational(1, 2))*
           exp(t*(-sqrt(5)/2 - Rational(1, 2))),
           C1*exp(t*(-Rational(1, 2) + sqrt(5)/2)) + C2*exp(t*(-sqrt(5)/2 - Rational(1, 2)))]
    assert constant_renumber(linodesolve(A, t), variables=Tuple(*eq).free_symbols) == sol

    # Testing the Errors
    raises(ValueError, lambda: linodesolve(1, t, b=Matrix([t+1])))
    raises(ValueError, lambda: linodesolve(a, t, b=Matrix([log(t) + sin(t)])))

    raises(ValueError, lambda: linodesolve(Matrix([7]), t, b=t**2))
    raises(ValueError, lambda: linodesolve(Matrix([a+10]), t, b=log(t)*cos(t)))

    raises(ValueError, lambda: linodesolve(7, t, b=t**2))
    raises(ValueError, lambda: linodesolve(a, t, b=log(t) + sin(t)))

    A1 = Matrix([[1, 2], [2, 4], [4, 6]])
    b1 = Matrix([t, 1, t**2])
    raises(NonSquareMatrixError, lambda: linodesolve(A1, t, b=b1))

    A2 = Matrix([[1, 2, 1], [3, 1, 2]])
    b2 = Matrix([t, t**2])
    raises(NonSquareMatrixError, lambda: linodesolve(A2, t, b=b2))

    raises(ValueError, lambda: linodesolve(A1[:2, :], t, b=b1))
    raises(ValueError, lambda: linodesolve(A1[:2, :], t, b=b1[:1]))

    # DOIT check
    A1 = Matrix([[1, -1], [1, -1]])
    b1 = Matrix([15*t - 10, -15*t - 5])
    sol1 = [C1 + C2*t + C2 - 10*t**3 + 10*t**2 + t*(15*t**2 - 5*t) - 10*t,
            C1 + C2*t - 10*t**3 - 5*t**2 + t*(15*t**2 - 5*t) - 5*t]
    assert constant_renumber(linodesolve(A1, t, b=b1, type="type2", doit=True),
                             variables=[t]) == sol1

    # Testing auto functionality
    func = [f(t), g(t)]
    eq = [Eq(f(t).diff(t) + g(t).diff(t), g(t) + t), Eq(g(t).diff(t), f(t))]
    ceq = canonical_odes(eq, func, t)
    (A1, A0), b = linear_ode_to_matrix(ceq[0], func, t, 1)
    A = A0
    sol = [-C1*exp(-t/2 + sqrt(5)*t/2)/2 + sqrt(5)*C1*exp(-t/2 + sqrt(5)*t/2)/2 - sqrt(5)*C2*exp(-sqrt(5)*t/2
                - t/2)/2 - C2*exp(-sqrt(5)*t/2 - t/2)/2 - exp(-t/2 +
                sqrt(5)*t/2)*Integral(sqrt(5)*t*exp(-sqrt(5)*t/2 + t/2)/5, t)/2 + sqrt(5)*exp(-t/2 +
                sqrt(5)*t/2)*Integral(sqrt(5)*t*exp(-sqrt(5)*t/2 + t/2)/5, t)/2 - sqrt(5)*exp(-sqrt(5)*t/2 -
                t/2)*Integral(-sqrt(5)*t*exp(t/2 + sqrt(5)*t/2)/5, t)/2 - exp(-sqrt(5)*t/2 -
                t/2)*Integral(-sqrt(5)*t*exp(t/2 + sqrt(5)*t/2)/5, t)/2,
            C1*exp(-t/2 + sqrt(5)*t/2) +
                C2*exp(-sqrt(5)*t/2 - t/2) + exp(-t/2 + sqrt(5)*t/2)*Integral(sqrt(5)*t*exp(-sqrt(5)*t/2 + t/2)/5,
                t) + exp(-sqrt(5)*t/2 - t/2)*Integral(-sqrt(5)*t*exp(t/2 + sqrt(5)*t/2)/5, t)]
    assert constant_renumber(linodesolve(A, t, b=b), variables=[t]) == sol

    # non-homogeneous term assumed to be 0
    sol1 = [-C1*exp(-t/2 + sqrt(5)*t/2)/2 + sqrt(5)*C1*exp(-t/2 + sqrt(5)*t/2)/2 - sqrt(5)*C2*exp(-sqrt(5)*t/2
                - t/2)/2 - C2*exp(-sqrt(5)*t/2 - t/2)/2 - exp(-t/2 + sqrt(5)*t/2)*Integral(0, t)/2 +
                sqrt(5)*exp(-t/2 + sqrt(5)*t/2)*Integral(0, t)/2 - sqrt(5)*exp(-sqrt(5)*t/2 - t/2)*Integral(0, t)/2
                - exp(-sqrt(5)*t/2 - t/2)*Integral(0, t)/2,
            C1*exp(-t/2 + sqrt(5)*t/2) + C2*exp(-sqrt(5)*t/2 - t/2)
                + exp(-t/2 + sqrt(5)*t/2)*Integral(0, t) + exp(-sqrt(5)*t/2 - t/2)*Integral(0, t)]
    assert constant_renumber(linodesolve(A, t, type="type2"), variables=[t]) == sol1

    # Testing the Errors
    raises(ValueError, lambda: linodesolve(t+10, t))
    raises(ValueError, lambda: linodesolve(a*t, t))

    A1 = Matrix([[1, t], [-t, 1]])
    B1, _ = _is_commutative_anti_derivative(A1, t)
    raises(NonSquareMatrixError, lambda: linodesolve(A1[:, :1], t, B=B1))
    raises(ValueError, lambda: linodesolve(A1, t, B=1))

    A2 = Matrix([[t, t, t], [t, t, t], [t, t, t]])
    B2, _ = _is_commutative_anti_derivative(A2, t)
    raises(NonSquareMatrixError, lambda: linodesolve(A2, t, B=B2[:2, :]))
    raises(ValueError, lambda: linodesolve(A2, t, B=2))
    raises(ValueError, lambda: linodesolve(A2, t, B=B2, type="type31"))

    raises(ValueError, lambda: linodesolve(A1, t, B=B2))
    raises(ValueError, lambda: linodesolve(A2, t, B=B1))

    # Testing auto functionality
    func = [f(t), g(t)]
    eq = [Eq(f(t).diff(t), f(t) + t*g(t)), Eq(g(t).diff(t), -t*f(t) + g(t))]
    ceq = canonical_odes(eq, func, t)
    (A1, A0), b = linear_ode_to_matrix(ceq[0], func, t, 1)
    A = A0
    sol = [(C1/2 - I*C2/2)*exp(I*t**2/2 + t) + (C1/2 + I*C2/2)*exp(-I*t**2/2 + t),
           (-I*C1/2 + C2/2)*exp(-I*t**2/2 + t) + (I*C1/2 + C2/2)*exp(I*t**2/2 + t)]
    assert constant_renumber(linodesolve(A, t), variables=Tuple(*eq).free_symbols) == sol
    assert constant_renumber(linodesolve(A, t, type="type3"), variables=Tuple(*eq).free_symbols) == sol

    A1 = Matrix([[t, 1], [t, -1]])
    raises(NotImplementedError, lambda: linodesolve(A1, t))

    # Testing the Errors
    raises(ValueError, lambda: linodesolve(t+10, t, b=Matrix([t+1])))
    raises(ValueError, lambda: linodesolve(a*t, t, b=Matrix([log(t) + sin(t)])))

    raises(ValueError, lambda: linodesolve(Matrix([7*t]), t, b=t**2))
    raises(ValueError, lambda: linodesolve(Matrix([a + 10*log(t)]), t, b=log(t)*cos(t)))

    raises(ValueError, lambda: linodesolve(7*t, t, b=t**2))
    raises(ValueError, lambda: linodesolve(a*t**2, t, b=log(t) + sin(t)))

    A1 = Matrix([[1, t], [-t, 1]])
    b1 = Matrix([t, t ** 2])
    B1, _ = _is_commutative_anti_derivative(A1, t)
    raises(NonSquareMatrixError, lambda: linodesolve(A1[:, :1], t, b=b1))

    A2 = Matrix([[t, t, t], [t, t, t], [t, t, t]])
    b2 = Matrix([t, 1, t**2])
    B2, _ = _is_commutative_anti_derivative(A2, t)
    raises(NonSquareMatrixError, lambda: linodesolve(A2[:2, :], t, b=b2))

    raises(ValueError, lambda: linodesolve(A1, t, b=b2))
    raises(ValueError, lambda: linodesolve(A2, t, b=b1))

    raises(ValueError, lambda: linodesolve(A1, t, b=b1, B=B2))
    raises(ValueError, lambda: linodesolve(A2, t, b=b2, B=B1))

    # Testing auto functionality
    func = [f(x), g(x), h(x)]
    eq = [Eq(f(x).diff(x), x*(f(x) + g(x) + h(x)) + x),
          Eq(g(x).diff(x), x*(f(x) + g(x) + h(x)) + x),
          Eq(h(x).diff(x), x*(f(x) + g(x) + h(x)) + 1)]
    ceq = canonical_odes(eq, func, x)
    (A1, A0), b = linear_ode_to_matrix(ceq[0], func, x, 1)
    A = A0
    _x1 = exp(-3*x**2/2)
    _x2 = exp(3*x**2/2)
    _x3 = Integral(2*_x1*x/3 + _x1/3 + x/3 - Rational(1, 3), x)
    _x4 = 2*_x2*_x3/3
    _x5 = Integral(2*_x1*x/3 + _x1/3 - 2*x/3 + Rational(2, 3), x)
    sol = [
        C1*_x2/3 - C1/3 + C2*_x2/3 - C2/3 + C3*_x2/3 + 2*C3/3 + _x2*_x5/3 + _x3/3 + _x4 - _x5/3,
        C1*_x2/3 + 2*C1/3 + C2*_x2/3 - C2/3 + C3*_x2/3 - C3/3 + _x2*_x5/3 + _x3/3 + _x4 - _x5/3,
        C1*_x2/3 - C1/3 + C2*_x2/3 + 2*C2/3 + C3*_x2/3 - C3/3 + _x2*_x5/3 - 2*_x3/3 + _x4 + 2*_x5/3,
    ]
    assert constant_renumber(linodesolve(A, x, b=b), variables=Tuple(*eq).free_symbols) == sol
    assert constant_renumber(linodesolve(A, x, b=b, type="type4"),
                             variables=Tuple(*eq).free_symbols) == sol

    A1 = Matrix([[t, 1], [t, -1]])
    raises(NotImplementedError, lambda: linodesolve(A1, t, b=b1))

    # non-homogeneous term not passed
    sol1 = [-C1/3 - C2/3 + 2*C3/3 + (C1/3 + C2/3 + C3/3)*exp(3*x**2/2), 2*C1/3 - C2/3 - C3/3 + (C1/3 + C2/3 + C3/3)*exp(3*x**2/2),
            -C1/3 + 2*C2/3 - C3/3 + (C1/3 + C2/3 + C3/3)*exp(3*x**2/2)]
    assert constant_renumber(linodesolve(A, x, type="type4", doit=True), variables=Tuple(*eq).free_symbols) == sol1


@slow
def test_linear_3eq_order1_type4_slow():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')

    f = t ** 3 + log(t)
    g = t ** 2 + sin(t)
    eq1 = (Eq(diff(x(t), t), (4 * f + g) * x(t) - f * y(t) - 2 * f * z(t)),
                Eq(diff(y(t), t), 2 * f * x(t) + (f + g) * y(t) - 2 * f * z(t)), Eq(diff(z(t), t), 5 * f * x(t) + f * y(
        t) + (-3 * f + g) * z(t)))
    with dotprodsimp(True):
        dsolve(eq1)


@slow
def test_linear_neq_order1_type2_slow1():
    i, r1, c1, r2, c2, t = symbols('i, r1, c1, r2, c2, t')
    x1 = Function('x1')
    x2 = Function('x2')

    eq1 = r1*c1*Derivative(x1(t), t) + x1(t) - x2(t) - r1*i
    eq2 = r2*c1*Derivative(x1(t), t) + r2*c2*Derivative(x2(t), t) + x2(t) - r2*i
    eq = [eq1, eq2]

    _x1 = sqrt(
        c1**2*r1**2 + 2*c1**2*r1*r2 + c1**2*r2**2 - 2*c1*c2*r1*r2 + 2*c1*c2*r2**2 + c2**2*r2**2)
    _x2 = -_x1*t/(2*c1*c2*r1*r2)
    _x3 = 1/(_x1 - c1*r1 - c1*r2 + c2*r2)
    _x4 = 1/(_x1 + c1*r1 + c1*r2 - c2*r2)
    _x5 = 1/(2*_x3*c1*c2*r2 + 2*_x4*c1*c2*r2)
    _x6 = Integral(_x5*i*exp(_x2 + t/(2*c2*r2) + t/(2*c2*r1) + t/(2*c1*r1)), t)
    _x7 = Integral(
        -_x5*i*exp(_x1*t/(2*c1*c2*r1*r2) + t/(2*c2*r2) + t/(2*c2*r1) + t/(2*c1*r1)),
        t)
    _x8 = exp(_x1*t/(2*c1*c2*r1*r2) - t/(2*c2*r2) - t/(2*c2*r1) - t/(2*c1*r1))
    _x9 = exp(_x2 - t/(2*c2*r2) - t/(2*c2*r1) - t/(2*c1*r1))
    sol = [
        Eq(x1(t),
           2*C1*_x3*_x8*c2*r2 - 2*C2*_x4*_x9*c2*r2 + 2*_x3*_x6*_x8*c2*r2 - 2*_x4*_x7*_x9*c2*r2),
        Eq(x2(t), C1*_x8 + C2*_x9 + _x6*_x8 + _x7*_x9),
    ]

    assert dsolve(eq) == sol
    assert checksysodesol(eq, sol) == (True, [0, 0])


def _de_lorentz_solution():
    m = Symbol("m", real=True)
    q = Symbol("q", real=True)
    t = Symbol("t", real=True)

    e1, e2, e3 = symbols("e1:4", real=True)
    b1, b2, b3 = symbols("b1:4", real=True)
    v1, v2, v3 = symbols("v1:4", cls=Function, real=True)

    eqs = [
        -e1*q + m*Derivative(v1(t), t) - q*(-b2*v3(t) + b3*v2(t)),
        -e2*q + m*Derivative(v2(t), t) - q*(b1*v3(t) - b3*v1(t)),
        -e3*q + m*Derivative(v3(t), t) - q*(-b1*v2(t) + b2*v1(t))
    ]

    # The code for the solution here is made using
    # printsol from https://github.com/sympy/sympy/issues/19574
    _x1 = 1/m
    _x2 = b1**2
    _x3 = b2**2
    _x4 = b3**2
    _x5 = sqrt(-_x2 - _x3 - _x4)
    _x6 = exp(2*_x1*_x5*q*t)
    _x7 = exp(_x1*_x5*q*t)
    _x8 = 1/(-_x2*_x7*b3 - _x3*_x7*b3)
    _x9 = sqrt(_x2 + _x3 + _x4)
    _x10 = 1/(_x2*m + _x3*m + _x4*m)
    _x11 = Integral(_x10*_x4*e3*q + _x10*b1*b3*e1*q + _x10*b2*b3*e2*q, t)
    _x12 = b1**3
    _x13 = b2**3
    _x14 = 1/(-2*I*_x2*_x7*_x9*m - 2*I*_x3*_x7*_x9*m - 2*I*_x4*_x7*_x9*m)
    _x15 = Integral(
        _x12*_x14*e2*q - _x13*_x14*e1*q - I*_x14*_x2*_x9*e3*q - _x14*_x2*b2*e1*q - I*_x14*_x3*_x9*e3*q + _x14*_x3*b1*e2*q + _x14*_x4*b1*e2*q - _x14*_x4*b2*e1*q + I*_x14*_x9*b1*b3*e1*q + I*_x14*_x9*b2*b3*e2*q,
        t)
    _x16 = 1/(
        -2*_x12*b3*m - 2*I*_x13*_x9*m - 2*I*_x2*_x9*b2*m - 2*_x3*b1*b3*m - 2*I*_x4*_x9*b2*m - 2*b1*b3**3*m)
    _x17 = Integral(
        -_x12*_x16*_x7*b2*e2*q - _x12*_x16*_x7*b3*e3*q - I*_x13*_x16*_x7*_x9*e3*q - _x13*_x16*_x7*b1*e2*q + _x16*_x2*_x3*_x7*e1*q + _x16*_x2*_x4*_x7*e1*q - I*_x16*_x2*_x7*_x9*b2*e3*q + I*_x16*_x2*_x7*_x9*b3*e2*q + _x16*_x3*_x4*_x7*e1*q + I*_x16*_x3*_x7*_x9*b3*e2*q - _x16*_x3*_x7*b1*b3*e3*q + _x16*_x7*b2**4*e1*q,
        t)
    _x18 = 1/(_x2*_x7*b3 + _x3*_x7*b3)
    sol = [
        Eq(v1(t),
           -C1*_x18*_x4*b1 - I*C1*_x18*_x9*b2*b3 + C2*_x12*_x18*_x7 + C2*_x18*_x3*_x7*b1 - C3*_x18*_x4*_x6*b1 + I*C3*_x18*_x6*_x9*b2*b3 + _x11*_x12*_x18*_x7 + _x11*_x18*_x3*_x7*b1 - _x15*_x18*_x4*_x6*b1 + I*_x15*_x18*_x6*_x9*b2*b3 - _x17*_x18*_x4*b1 - I*_x17*_x18*_x9*b2*b3),
        Eq(v2(t),
           C1*_x4*_x8*b2 - I*C1*_x8*_x9*b1*b3 - C2*_x13*_x7*_x8 - C2*_x2*_x7*_x8*b2 + C3*_x4*_x6*_x8*b2 + I*C3*_x6*_x8*_x9*b1*b3 - _x11*_x13*_x7*_x8 - _x11*_x2*_x7*_x8*b2 + _x15*_x4*_x6*_x8*b2 + I*_x15*_x6*_x8*_x9*b1*b3 + _x17*_x4*_x8*b2 - I*_x17*_x8*_x9*b1*b3),
        Eq(v3(t), C2 + C3*_x7 + _x11 + _x15*_x7 + (C1 + _x17)*exp(-_x1*_x5*q*t)),
    ]

    return eqs, sol

# Regression test case for issue #9204
# https://github.com/sympy/sympy/issues/9204

# A very big solution is obtained for this
# test case. To be simplified in future.
@slow
def test_linear_new_order1_type2_de_lorentz():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eqs, sol = _de_lorentz_solution()

    with dotprodsimp(True):
        assert dsolve(eqs) == sol


@slow
def test_linear_new_order1_type2_de_lorentz_slow_check():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eqs, sol = _de_lorentz_solution()
    assert checksysodesol(eqs, sol) == (True, [0, 0, 0])


def _neq_order1_type2_slow():
    RC, t, C, Vs, L, R1, V0, I0 = symbols("RC t C Vs L R1 V0 I0")
    V = Function("V")
    I = Function("I")
    system = [Eq(V(t).diff(t), -1/RC*V(t) + I(t)/C), Eq(I(t).diff(t), -R1/L*I(t) - 1/L*V(t) + Vs/L)]

    z1 = sqrt(C**2*L**2 - 2*C**2*L*R1*RC + C**2*R1**2*RC**2 - 4*C*L*RC**2)
    z2 = 1/(C*L - C*R1*RC - z1)
    z3 = 1/(C*L - C*R1*RC + z1)
    z4 = exp(-t/(2*RC) - R1*t/(2*L) + t*z1/(2*C*L*RC))
    z5 = exp(-t/(2*RC) - R1*t/(2*L) - t*z1/(2*C*L*RC))
    z6 = Integral(2*RC*Vs*exp(t/(2*RC) + R1*t/(2*L) +
        t*z1/(2*C*L*RC))/(-2*C*L**2*RC*z2 + 2*C*L**2*RC*z3 + 2*C*L*R1*RC**2*z2
            - 2*C*L*R1*RC**2*z3 - 2*L*RC*z1*z2 + 2*L*RC*z1*z3), t)
    z7 = Integral(-2*RC*Vs*exp(t/(2*RC) + R1*t/(2*L) -
        t*z1/(2*C*L*RC))/(-2*C*L**2*RC*z2 + 2*C*L**2*RC*z3 + 2*C*L*R1*RC**2*z2
            - 2*C*L*R1*RC**2*z3 + 2*L*RC*z1*z2 - 2*L*RC*z1*z3), t)
    sol = [
        Eq(V(t), 2*C1*L*RC*z2*z5 + 2*C2*L*RC*z3*z4 + 2*L*RC*z2*z5*z6 + 2*L*RC*z3*z4*z7),
        Eq(I(t), C1*z5 + C2*z4 + z4*z7 + z5*z6),
    ]

    return system, sol


# A very big solution is obtained for this
# test case. To be simplified in future.
def test_linear_neq_order1_type2_slow():
    system, sol = _neq_order1_type2_slow()

    assert dsolve(system) == sol

# Regression test case for issue #14001
# https://github.com/sympy/sympy/issues/14001

@slow
def test_linear_neq_order1_type2_slow_check():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    system, sol = _neq_order1_type2_slow()

    assert checksysodesol(system, sol) == (True, [0, 0])


def _linear_3eq_order1_type4_long():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')

    f = t ** 3 + log(t)
    g = t ** 2 + sin(t)

    eq1 = (Eq(diff(x(t), t), (4*f + g)*x(t) - f*y(t) - 2*f*z(t)),
           Eq(diff(y(t), t), 2*f*x(t) + (f + g)*y(t) - 2*f*z(t)), Eq(diff(z(t), t), 5*f*x(t) + f*y(
        t) + (-3*f + g)*z(t)))

    dsolve_sol = dsolve(eq1)
    dsolve_sol1 = [_simpsol(sol) for sol in dsolve_sol]

    x_1 = sqrt(-t**6 - 8*t**3*log(t) + 8*t**3 - 16*log(t)**2 + 32*log(t) - 16)
    x_2 = sqrt(3)
    x_3 = 8324372644*C1*x_1*x_2 + 4162186322*C2*x_1*x_2 - 8324372644*C3*x_1*x_2
    x_4 = 1 / (1903457163*t**3 + 3825881643*x_1*x_2 + 7613828652*log(t) - 7613828652)
    x_5 = exp(t**3/3 + t*x_1*x_2/4 - cos(t))
    x_6 = exp(t**3/3 - t*x_1*x_2/4 - cos(t))
    x_7 = exp(t**4/2 + t**3/3 + 2*t*log(t) - 2*t - cos(t))
    x_8 = 91238*C1*x_1*x_2 + 91238*C2*x_1*x_2 - 91238*C3*x_1*x_2
    x_9 = 1 / (66049*t**3 - 50629*x_1*x_2 + 264196*log(t) - 264196)
    x_10 = 50629 * C1 / 25189 + 37909*C2/25189 - 50629*C3/25189 - x_3*x_4
    x_11 = -50629*C1/25189 - 12720*C2/25189 + 50629*C3/25189 + x_3*x_4
    sol = [Eq(x(t), x_10*x_5 + x_11*x_6 + x_7*(C1 - C2)), Eq(y(t), x_10*x_5 + x_11*x_6), Eq(z(t), x_5*(
            -424*C1/257 - 167*C2/257 + 424*C3/257 - x_8*x_9) + x_6*(167*C1/257 + 424*C2/257 -
            167*C3/257 + x_8*x_9) + x_7*(C1 - C2))]

    assert dsolve_sol1 == sol
    assert checksysodesol(eq1, dsolve_sol1) == (True, [0, 0, 0])


def _neq_order1_type4_slow1():
    f, g = symbols("f, g", cls=Function)
    x = Symbol("x")

    eqs = [Eq(diff(f(x), x), x*f(x) + x**2*g(x) + x),
           Eq(diff(g(x), x), 2*x**2*f(x) + (x + 3*x**2)*g(x) + 1)]
    _x1 = sqrt(17)
    _x2 = 1/(45*_x1 + 187)
    _x3 = 1/(161*_x1 + 663)
    _x4 = 1/(51*_x1 - 221)
    _x5 = 1/(161*_x1 - 663)
    _x6 = -x**3/2
    _x7 = _x1*x**3/6
    _x8 = exp(_x6 + _x7 - x**2/2)
    _x9 = -_x1*x**3/6
    _x10 = exp(_x6 + _x9 - x**2/2)
    _x11 = Integral(
        22*_x1*_x10*_x3*x + 39*_x1*_x10*_x3 + 139*_x1*_x3*_x8*x - 39*_x1*_x3*_x8 + 90*_x10*_x3*x + 161*_x10*_x3 + 573*_x3*_x8*x - 161*_x3*_x8,
        x)
    _x12 = exp(_x9 + x**3/2 + x**2/2)
    _x13 = exp(_x7 + x**3/2 + x**2/2)
    _x14 = Integral(
        22*_x1*_x10*_x2*x + 39*_x1*_x10*_x2 - 22*_x1*_x2*_x8*x + 6*_x1*_x2*_x8 + 90*_x10*_x2*x + 161*_x10*_x2 - 90*_x2*_x8*x + 26*_x2*_x8,
        x)
    sol = [
        Eq(f(x), _x12*(
            39*C1*_x1*_x5 - 161*C1*_x5 + 22*C2*_x1*_x5 - 90*C2*_x5 + 22*_x1*_x11*_x5 + 39*_x1*_x14*_x5 - 90*_x11*_x5 - 161*_x14*_x5) + _x13*(
               -39*C1*_x1*_x5 + 161*C1*_x5 + 139*C2*_x1*_x5 - 573*C2*_x5 + 139*_x1*_x11*_x5 - 39*_x1*_x14*_x5 - 573*_x11*_x5 + 161*_x14*_x5)),
        Eq(g(x), _x12*(
            45*C1*_x1*_x4 - 187*C1*_x4 + 26*C2*_x1*_x4 - 102*C2*_x4 + 26*_x1*_x11*_x4 + 45*_x1*_x14*_x4 - 102*_x11*_x4 - 187*_x14*_x4) + _x13*(
               6*C1*_x1*_x4 - 34*C1*_x4 - 26*C2*_x1*_x4 + 102*C2*_x4 - 26*_x1*_x11*_x4 + 6*_x1*_x14*_x4 + 102*_x11*_x4 - 34*_x14*_x4)),
    ]

    return eqs, sol


def test_neq_order1_type4_slow1():
    eqs, sol = _neq_order1_type4_slow1()
    with dotprodsimp(True):
        assert dsolve(eqs) == sol


@slow
def test_neq_order1_type4_slow_check1():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eqs, sol = _neq_order1_type4_slow1()
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def _neq_order1_type4_slow2():
    f, g, h = symbols("f, g, h", cls=Function)
    x = Symbol("x")

    eqs = [Eq(Derivative(f(x), x), x*h(x) + f(x) + g(x) + 1), Eq(Derivative(g(x), x), x*g(x) + f(x) + h(x) +
            10), Eq(Derivative(h(x), x), x*f(x) + x + g(x) + h(x))]
    _x1 = -x**2/2
    _x2 = Integral(x*exp(_x1 - 2*x)/3 - x*exp(_x1 + x)/3 + 11*exp(_x1 - 2*x)/3 + 19*exp(_x1 + x)/3,
                   x)
    _x3 = Integral(x*exp(_x1 - 2*x)/3 + x*exp(_x1 + x)/6 - x*exp(x**2/2 - x)/2 + 11*exp(
        _x1 - 2*x)/3 - 19*exp(_x1 + x)/6 + exp(x**2/2 - x)/2, x)
    _x4 = Integral(x*exp(_x1 - 2*x)/3 + x*exp(_x1 + x)/6 + x*exp(x**2/2 - x)/2 + 11*exp(
        _x1 - 2*x)/3 - 19*exp(_x1 + x)/6 - exp(x**2/2 - x)/2, x)
    _x5 = (C1/3 + C2/3 + C3/3 + _x2/3 + _x3/3 + _x4/3)*exp(x**2/2 + 2*x)
    sol = [
        Eq(f(x), _x5 + (-C1/2 + C2/2 + _x3/2 - _x4/2)*exp(_x1 + x) + (
            C1/6 + C2/6 - C3/3 - _x2/3 + _x3/6 + _x4/6)*exp(x**2/2 - x)),
        Eq(g(x), _x5 + (-C1/3 - C2/3 + 2*C3/3 + 2*_x2/3 - _x3/3 - _x4/3)*exp(x**2/2 - x)),
        Eq(h(x), _x5 + (C1/2 - C2/2 - _x3/2 + _x4/2)*exp(_x1 + x) + (
            C1/6 + C2/6 - C3/3 - _x2/3 + _x3/6 + _x4/6)*exp(x**2/2 - x)),
    ]

    return eqs, sol


def test_neq_order1_type4_slow2():
    eqs, sol = _neq_order1_type4_slow2()
    with dotprodsimp(True):
        assert dsolve(eqs) == sol


@slow
def test_neq_order1_type4_slow_check2():
    eqs, sol = _neq_order1_type4_slow2()
    assert checksysodesol(eqs, sol) == (True, [0, 0, 0])


def _neq_order1_type4_slow3():
    f, g = symbols("f g", cls=Function)
    x = symbols("x")

    eqs = [Eq(Derivative(f(x), x), x*f(x) + g(x) + sin(x)), Eq(Derivative(g(x), x), x**2 + x*g(x) - f(x))]
    _x1 = exp(x**2/2 + I*x)
    _x2 = exp(-x**2/2 + I*x)
    _x3 = exp(-x**2/2 - I*x)
    _x4 = Integral(_x2*x**2/2 - I*_x2*sin(x)/2 + _x3*x**2/2 + I*_x3*sin(x)/2, x)
    _x5 = Integral(I*_x2*x**2/2 + _x2*sin(x)/2 - I*_x3*x**2/2 + _x3*sin(x)/2, x)
    sol = [
        Eq(f(x),
           _x1*(C1/2 - I*C2/2 - I*_x4/2 + _x5/2) + (C1/2 + I*C2/2 + I*_x4/2 + _x5/2)*exp(
               x**2/2 - I*x)),
        Eq(g(x),
           _x1*(I*C1/2 + C2/2 + _x4/2 + I*_x5/2) + (-I*C1/2 + C2/2 + _x4/2 - I*_x5/2)*exp(
               x**2/2 - I*x)),
    ]

    return eqs, sol


def test_neq_order1_type4_slow3():
    eqs, sol = _neq_order1_type4_slow3()
    assert dsolve(eqs) == sol


@slow
def test_neq_order1_type4_slow_check3():
    eqs, sol = _neq_order1_type4_slow3()
    assert checksysodesol(eqs, sol) == (True, [0, 0])


@XFAIL
@slow
def test_linear_3eq_order1_type4_long_dsolve_slow_xfail():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eq, sol = _linear_3eq_order1_type4_long()

    dsolve_sol = dsolve(eq)
    dsolve_sol1 = [_simpsol(sol) for sol in dsolve_sol]

    assert dsolve_sol1 == sol


@slow
def test_linear_3eq_order1_type4_long_dsolve_dotprodsimp():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eq, sol = _linear_3eq_order1_type4_long()

    # XXX: Only works with dotprodsimp see
    # test_linear_3eq_order1_type4_long_dsolve_slow_xfail which is too slow
    with dotprodsimp(True):
        dsolve_sol = dsolve(eq)

    dsolve_sol1 = [_simpsol(sol) for sol in dsolve_sol]
    assert dsolve_sol1 == sol


@slow
def test_linear_3eq_order1_type4_long_check():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eq, sol = _linear_3eq_order1_type4_long()
    assert checksysodesol(eq, sol) == (True, [0, 0, 0])


def test_dsolve_system():
    f, g = symbols("f g", cls=Function)
    x = symbols("x")
    eqs = [Eq(f(x).diff(x), f(x) + g(x)), Eq(g(x).diff(x), f(x) + g(x))]
    funcs = [f(x), g(x)]

    sol = [[Eq(f(x), -C1 + C2*exp(2*x)), Eq(g(x), C1 + C2*exp(2*x))]]
    assert dsolve_system(eqs, funcs=funcs, t=x, doit=True) == sol

    raises(ValueError, lambda: dsolve_system(1))
    raises(ValueError, lambda: dsolve_system(eqs, 1))
    raises(ValueError, lambda: dsolve_system(eqs, funcs, 1))
    raises(ValueError, lambda: dsolve_system(eqs, funcs[:1], x))

    eq = (Eq(f(x).diff(x), 12 * f(x) - 6 * g(x)), Eq(g(x).diff(x) ** 2, 11 * f(x) + 3 * g(x)))
    raises(NotImplementedError, lambda: dsolve_system(eq) == ([], []))

    raises(NotImplementedError, lambda: dsolve_system(eq, funcs=[f(x), g(x)]) == ([], []))
    raises(NotImplementedError, lambda: dsolve_system(eq, funcs=[f(x), g(x)], t=x) == ([], []))
    raises(NotImplementedError, lambda: dsolve_system(eq, funcs=[f(x), g(x)], t=x, ics={f(0): 1, g(0): 1}) == ([], []))
    raises(NotImplementedError, lambda: dsolve_system(eq, t=x, ics={f(0): 1, g(0): 1}) == ([], []))
    raises(NotImplementedError, lambda: dsolve_system(eq, ics={f(0): 1, g(0): 1}) == ([], []))
    raises(NotImplementedError, lambda: dsolve_system(eq, funcs=[f(x), g(x)], ics={f(0): 1, g(0): 1}) == ([], []))

def test_dsolve():

    f, g = symbols('f g', cls=Function)
    x, y = symbols('x y')

    eqs = [f(x).diff(x) - x, f(x).diff(x) + x]
    with raises(ValueError):
        dsolve(eqs)

    eqs = [f(x, y).diff(x)]
    with raises(ValueError):
        dsolve(eqs)

    eqs = [f(x, y).diff(x)+g(x).diff(x), g(x).diff(x)]
    with raises(ValueError):
        dsolve(eqs)


def _higher_order_slow1():
    x, y = symbols("x y", cls=Function)
    t = symbols("t")

    eq = (Eq(diff(x(t),t,t), (log(t)+t**2)*diff(x(t),t)+(log(t)+t**2)*3*diff(y(t),t)), Eq(diff(y(t),t,t), \
    (log(t)+t**2)*2*diff(x(t),t)+(log(t)+t**2)*9*diff(y(t),t)))
    _x1 = log(t)**2
    _x2 = sqrt(22)
    _x3 = sqrt(9*_x1 + t**4 + 6*t**2*log(t) - 6*t**2 - 18*log(t) + 9)
    _x4 = 1/(-_x2*_x3*t/3 - 4*t**3/3 - 4*t*log(t) + 4*t)
    _x5 = 1/(_x2*_x3*t/3 - 4*t**3/3 - 4*t*log(t) + 4*t)
    _x6 = 1/(-_x4*t**3 - 3*_x4*t*log(t) + 3*_x4*t + _x5*t**3 + 3*_x5*t*log(t) - 3*_x5*t)
    _x7 = exp(-_x2*_x3*t/3 + 5*t**3/3 + 5*t*log(t) - 5*t)
    _x8 = exp(_x2*_x3*t/3 + 5*t**3/3 + 5*t*log(t) - 5*t)
    _x9 = 1/(18*_x1*_x4*t**5 - 54*_x1*_x4*t**3 - 18*_x1*_x5*t**5 + 54*_x1*_x5*t**3 + 2*_x4*t**9/3 + 6*_x4*t**7*log(t) -
             6*_x4*t**7 - 36*_x4*t**5*log(t) + 18*_x4*t**5 + 18*_x4*t**3*log(t)**3 + 54*_x4*t**3*log(t) - 18*_x4*t**3 -
             2*_x5*t**9/3 - 6*_x5*t**7*log(t) + 6*_x5*t**7 + 36*_x5*t**5*log(t) - 18*_x5*t**5 - 18*_x5*t**3*log(t)**3 -
             54*_x5*t**3*log(t) + 18*_x5*t**3)
    _x10 = -12*_x1*_x5*t**2
    _x11 = 24*_x5*t**2*log(t)
    _x12 = 8*_x4*t**4*log(t)
    _x13 = -24*_x4*t**2*log(t)
    _x14 = -8*_x5*t**4*log(t)
    _x15 = 1/(12*_x1*_x4*t**2 + _x10 + _x11 + _x12 + _x13 + _x14 + _x2*_x3*_x4*t**4/3 + _x2*_x3*_x4*t**2*log(t) -
              _x2*_x3*_x4*t**2 - _x2*_x3*_x5*t**4/3 - _x2*_x3*_x5*t**2*log(t) + _x2*_x3*_x5*t**2 + 4*_x4*t**6/3 -
              8*_x4*t**4 + 12*_x4*t**2 - 4*_x5*t**6/3 + 8*_x5*t**4 - 12*_x5*t**2)
    _x16 = 1/(12*_x1*_x4*t**2 + _x10 + _x11 + _x12 + _x13 + _x14 - _x2*_x3*_x4*t**4/3 - _x2*_x3*_x4*t**2*log(t) +
              _x2*_x3*_x4*t**2 + _x2*_x3*_x5*t**4/3 + _x2*_x3*_x5*t**2*log(t) - _x2*_x3*_x5*t**2 + 4*_x4*t**6/3 -
              8*_x4*t**4 + 12*_x4*t**2 - 4*_x5*t**6/3 + 8*_x5*t**4 - 12*_x5*t**2)
    _x17 = (t**3 + 3*t*log(t) - 3*t)**2
    sol = [
        Eq(x(t), C1 + Integral(-C2*_x15*_x8*t**3 - 3*C2*_x15*_x8*t*log(t) + 3*C2*_x15*_x8*t + C2*_x16*_x7*t**3 +
            3*C2*_x16*_x7*t*log(t) - 3*C2*_x16*_x7*t + C3*_x17*_x7*_x9 - C3*_x17*_x8*_x9, t)),
        Eq(y(t), C4 + Integral(-C2*_x6*_x7 + C2*_x6*_x8 - C3*_x15*_x7*t**3 - 3*C3*_x15*_x7*t*log(t) + 3*C3*_x15*_x7*t +
            C3*_x16*_x8*t**3 + 3*C3*_x16*_x8*t*log(t) - 3*C3*_x16*_x8*t, t)),
    ]

    return eq, sol


@slow
def test_higher_order_slow1():
    eq, sol = _higher_order_slow1()

    assert dsolve(eq) == sol


@slow
def test_higher_order1_slow1_check():
    if ON_TRAVIS:
        skip("Too slow for travis.")

    eq, sol = _higher_order_slow1()
    assert checksysodesol(eq, sol) == (True, [0, 0])


# Note: Test case fails not because of higher order reduction
# but due to type 1 solver's incorrect answer.
# Given a different order, the solution works.
# A specific test case will be added soon if it isn't
# fixed in the current PR.
@XFAIL
def test_higher_order_type1_fail():
    x, y = symbols("x y", cls=Function)
    t = symbols('t')

    eq2 = (Eq(diff(x(t),t,t) - 9*diff(y(t),t) + 7*x(t),0), Eq(diff(y(t),t,t) + 9*diff(x(t),t) + 7*y(t),0))
    sol2 = [Eq(x(t), (Rational(9, 14) - sqrt(109)/14)*(-C1*sin(t*sqrt(9*sqrt(109)/2 + Rational(95, 2))) + C2*
               cos(t*sqrt(9*sqrt(109)/2 + Rational(95, 2)))) + (Rational(9, 14) + sqrt(109)/14)*(-C3*sin(t*sqrt(Rational(95, 2)
               - 9*sqrt(109)/2)) + C4*cos(t*sqrt(Rational(95, 2) - 9*sqrt(109)/2)))),
            Eq(y(t), sqrt(2)*(C1*cos(t*sqrt(9*sqrt(109)/2 + Rational(95, 2))) + C2*sin(t*sqrt(9*sqrt(109)/2 + Rational(95, 2))))/sqrt(9*
                sqrt(109) + 95) + sqrt(2)*(C3*cos(t*sqrt(Rational(95, 2) - 9*sqrt(109)/2)) + C4*sin(t*sqrt(Rational(95, 2)
                - 9*sqrt(109)/2)))/sqrt(95 - 9*sqrt(109)))]
    assert dsolve(eq2) == sol2
    assert checksysodesol(eq2, sol2) == (True, [0, 0])

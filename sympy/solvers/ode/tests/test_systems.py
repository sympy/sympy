from sympy import (symbols, Symbol, diff, Function, Derivative, Matrix, Rational)
from sympy.solvers.ode.systems import neq_nth_linear_constant_coeff_match

def test_neq_nth_linear_constant_coeff_match():
    x, y = symbols('x, y', cls=Function)
    t = Symbol('t')
    x1 = diff(x(t), t)
    y1 = diff(y(t), t)
    x2 = diff(x(t), t, t)
    funcs = [x(t), y(t)]

    eqs_1 = (5 * x1 + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * t * x(t) + 3 * y(t) + t))
    answer_1 = {'no_of_equation': 2,
     'eq': (12*x(t) - 6*y(t) + 5*Derivative(x(t), t),
      -11*t*x(t) + t + 3*y(t) + 2*Derivative(y(t), t)),
     'func': [x(t), y(t)],
     'order': {x(t): 1, y(t): 1},
     'is_linear': True,
     'is_constant': False,
     'is_homogeneous': False}
    assert neq_nth_linear_constant_coeff_match(eqs_1, funcs, t) == answer_1

    eqs_2 = (5 * (x1**2) + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * t * x(t) + 3 * y(t) + t))
    answer_2 = None
    assert neq_nth_linear_constant_coeff_match(eqs_2, funcs, t) is answer_2

    eqs_3 = (5 * x1 + 12 * x(t) - 6 * (y(t)), (2 * y1 - 11 * x(t) + 3 * y(t)))
    answer_3 = {'no_of_equation': 2,
     'eq': (12 * x(t) - 6 * y(t) + 5 * Derivative(x(t), t),
            -11 * x(t) + 3 * y(t) + 2 * Derivative(y(t), t)),
     'func': [x(t), y(t)],
     'order': {x(t): 1, y(t): 1},
     'is_linear': True,
     'is_constant': True,
     'is_homogeneous': True,
     'func_coeff': Matrix([
         [Rational(12, 5), Rational(-6, 5)],
         [Rational(-11, 2), Rational(3, 2)]])}
    assert neq_nth_linear_constant_coeff_match(eqs_3, funcs, t) == answer_3

    eqs_4 = (5 * x1 + 12 * x(t) - 6 * (y(t)) + x2, (2 * y1 - 11 * x(t) + 3 * y(t)))
    answer_4 = {'no_of_equation': 2,
     'eq': (12*x(t) - 6*y(t) + 5*Derivative(x(t), t) + Derivative(x(t), (t, 2)),
      -11*x(t) + 3*y(t) + 2*Derivative(y(t), t)),
     'func': [x(t), y(t)],
     'order': {x(t): 2, y(t): 1},
     'is_linear': True,
     'is_constant': True,
     'is_homogeneous': True}
    assert neq_nth_linear_constant_coeff_match(eqs_4, funcs, t) == answer_4

from sympy import symbols, exp, sin, cos, pi, simplify
from sympy.simplify.precision_optimize import algebraic_rewrite

a, b, x = symbols('a b x')


def test_difference_of_squares_symbolic():
    expr = a**2 - b**2
    rewritten = algebraic_rewrite(expr)
    # Should simplify to (a - b)*(a + b)
    assert simplify(rewritten - (a - b)*(a + b)) == 0


def test_difference_of_squares_numeric():
    expr = a**2 - b**2
    val_plain = expr.subs({a: 1.0001, b: 1.0}).evalf()
    val_opt = expr.evalf(optimize_for_precision=True, subs={a: 1.0001, b: 1.0})
    # Both should be close and not throw
    assert abs(val_plain - val_opt) < 1e-12


def test_exponential_cancellation():
    expr = exp(x) - 1
    normal = expr.subs({x: 1e-10}).evalf()
    optimized = expr.evalf(optimize_for_precision=True, subs={x: 1e-10})
    # Optimized result should be closer to true value
    true_value = 1e-10 * exp(0)  # ≈ 1e-10
    assert abs(optimized - true_value) < abs(normal - true_value)


def test_trig_identity():
    expr = sin(pi/4)**2 + cos(pi/4)**2
    val_plain = expr.evalf()
    val_opt = expr.evalf(optimize_for_precision=True)
    assert abs(val_opt - 1) < 1e-12
    assert abs(val_opt - val_plain) < 1e-12


def test_nested_expression():
    expr = (a**2 - b**2) + (exp(x) - 1)
    val_opt = expr.evalf(optimize_for_precision=True, subs={a: 1.0001, b: 1.0, x: 1e-10})
    val_plain = expr.subs({a: 1.0001, b: 1.0, x: 1e-10}).evalf()
    assert abs(val_opt - val_plain) < 1e-8  # should be numerically stable


def test_no_error_on_non_numeric():
    expr = a + b
    # Should run without raising an exception
    expr.evalf(optimize_for_precision=True)


def test_invalid_kwarg_raises():
    # Passing an unexpected keyword should raise TypeError
    import pytest
    from sympy import pi
    with pytest.raises(TypeError):
        pi.evalf(random_kwarg=True)

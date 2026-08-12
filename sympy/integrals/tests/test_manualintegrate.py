from sympy import Symbol, symbols, sin, cos, tan, cot, sec, csc, Abs, sign, log, sqrt, I
from sympy.integrals.manualintegrate import manualintegrate

def test_manualintegrate_abs_trig_direct():
    x = Symbol('x')

    # Direct Abs(...) tests
    assert manualintegrate(Abs(sin(x)), x) == -sign(sin(x)) * cos(x)
    assert manualintegrate(Abs(cos(x)), x) == sign(cos(x)) * sin(x)
    assert manualintegrate(Abs(tan(x)), x) == -sign(tan(x)) * log(Abs(cos(x)))
    assert manualintegrate(Abs(cot(x)), x) == sign(cot(x)) * log(Abs(sin(x)))
    assert manualintegrate(Abs(sec(x)), x) == sign(cos(x)) * log(Abs(sec(x) + tan(x)))
    assert manualintegrate(Abs(csc(x)), x) == sign(sin(x)) * log(Abs(tan(x / 2)))

def test_manualintegrate_abs_trig_reciprocal():
    x = Symbol('x')

    # Reciprocal 1/Abs(...) tests
    assert manualintegrate(1 / Abs(sin(x)), x) == sign(sin(x)) * log(Abs(tan(x / 2)))
    assert manualintegrate(1 / Abs(cos(x)), x) == sign(cos(x)) * log(Abs(sec(x) + tan(x)))
    assert manualintegrate(1 / Abs(tan(x)), x) == sign(cot(x)) * log(Abs(sin(x)))
    assert manualintegrate(1 / Abs(cot(x)), x) == -sign(tan(x)) * log(Abs(cos(x)))

def test_manualintegrate_abs_trig_linear_args():
    x = Symbol('x')

    # Linear argument tests (a*x + b) across different functions
    assert manualintegrate(Abs(sin(2*x + 3)), x) == -sign(sin(2*x + 3)) * cos(2*x + 3) / 2
    assert manualintegrate(Abs(cos(3*x - 1)), x) == sign(cos(3*x - 1)) * sin(3*x - 1) / 3
    assert manualintegrate(1 / Abs(sin(2*x)), x) == sign(sin(2*x)) * log(Abs(tan(x))) / 2
    assert manualintegrate(1 / Abs(sin(x + 1)), x) == sign(sin(x + 1)) * log(Abs(tan((x + 1) / 2)))

def test_issue_30163_abs_trig_sqrt():
    # Test the radical normalization with positive assumptions (Original bug)
    x_pos, z_pos = symbols('x z', positive=True)

    result = manualintegrate(1 / sqrt(-z_pos * sin(x_pos)**2), x_pos)
    expected = -I * sign(sin(x_pos)) * log(Abs(tan(x_pos / 2))) / sqrt(z_pos)

    assert result == expected

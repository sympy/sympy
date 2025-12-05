from sympy import Symbol, simplify, Abs

def test_issue_20108():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)

    # 1. The main fix: x^2 / |x| should be |x|
    assert simplify(x**2 / Abs(x)) == Abs(x)

    # 2. Check it works with other variable names
    assert simplify(y**2 / Abs(y)) == Abs(y)

    # 3. Check it does NOT break complex numbers (Safety check)
    # (Complex numbers cannot simplify this way)
    z = Symbol('z') # z is complex by default
    assert simplify(z**2 / Abs(z)) != Abs(z)
    
from sympy import symbols, rf, Symbol, factorial, Factorial, ff, nan, oo

def test_rf_eval_apply():

    x, y = symbols('xy')

    assert rf(nan, y) == nan

    assert rf(x, y) == rf(x, y)

    assert rf(oo, 0) == 1
    assert rf(-oo, 0) == 1

    assert rf(oo, 6) == oo
    assert rf(-oo, 7) == -oo

    assert rf(oo, -6) == oo
    assert rf(-oo, -7) == oo

    assert rf(x, 0) == 1
    assert rf(x, 1) == x
    assert rf(x, 2) == x*(x+1)
    assert rf(x, 3) == x*(x+1)*(x+2)
    assert rf(x, 5) == x*(x+1)*(x+2)*(x+3)*(x+4)

    assert rf(x, -1) == 1/(x-1)
    assert rf(x, -2) == 1/((x-1)*(x-2))
    assert rf(x, -3) == 1/((x-1)*(x-2)*(x-3))

    assert rf(1, 100) == factorial(100)

def test_ff_eval_apply():

    x, y = symbols('xy')

    assert ff(nan, y) == nan

    assert ff(x, y) == ff(x, y)

    assert ff(oo, 0) == 1
    assert ff(-oo, 0) == 1

    assert ff(oo, 6) == oo
    assert ff(-oo, 7) == -oo

    assert ff(oo, -6) == oo
    assert ff(-oo, -7) == oo

    assert ff(x, 0) == 1
    assert ff(x, 1) == x
    assert ff(x, 2) == x*(x-1)
    assert ff(x, 3) == x*(x-1)*(x-2)
    assert ff(x, 5) == x*(x-1)*(x-2)*(x-3)*(x-4)

    assert ff(x, -1) == 1/(x+1)
    assert ff(x, -2) == 1/((x+1)*(x+2))
    assert ff(x, -3) == 1/((x+1)*(x+2)*(x+3))

    assert ff(100, 100) == factorial(100)

def test_factorials():
    n = Symbol('n', integer=True)
    assert factorial(-2) == 0
    assert factorial(0) == 1
    assert factorial(7) == 5040
    assert factorial(n).func == Factorial
    assert factorial(2*n).func == Factorial

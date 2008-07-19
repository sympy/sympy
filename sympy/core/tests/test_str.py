from sympy import Symbol, Wild, Rational, Derivative, I, log, sqrt, exp, sin

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
w = Symbol('w')


def test_pow():
    assert str(1/x) == "1/x"

def test_poly_str():
    #if any of these tests fails, it can still be correct, just the terms can
    #be in a different order. That happens for example when we change the
    #hash algorithm. If it is correct, just add another item in the list [] of
    #correct results.
    assert str((2*x-(7*x**2 - 2) + 3*y)) in ["2 - 7*x**2 + 2*x + 3*y",
            "2 + 3*y + 2*x - 7*x**2", "2 + 3*y - 7*x**2 + 2*x",
            "3*y + 2*x + 2 - 7*x**2", "2 + 2*x + 3*y - 7*x**2"]
    assert str(x-y) in ["x - y", "-y + x"]
    assert str(2+-x) in ["2 - x", "-x + 2"]
    assert str(x-2) in ["x - 2", "(-2) + x", "-2 + x"]
    assert str(x-y-z-w) in ["x - y - z - w","-w - y - z + x","x - w - y - z",
                            "-w + x - y - z","-z - w - y + x","-y + x - w - z",
                            "-y - z - w + x"]
    assert str(x-y-z-w) in [
            "-w - y - z + x","x - w - y - z","-w + x - z - y","-y - w - z + x",
            "-y + x - z - w","-y + x - w - z","-w + x - y - z","-z - w - y +x",
            "-y - z - w + x"]
    assert str(x-z*y**2*z*w) in ["-z**2*y**2*w + x", "x - w*y**2*z**2",
            "-y**2*z**2*w + x","x - w*z**2*y**2","x - y**2*z**2*w",
            "x - y**2*w*z**2","x - z**2*y**2*w","-w*z**2*y**2 + x",
            "-w*y**2*z**2 + x","x - z**2*w*y**2"]

def test_bug1():
    assert str(x-1*y*x*y) in ["x - x*y**2", "-x*y**2 + x"]

def test_bug2():
    e = x-y
    a = str(e)
    b = str(e)
    assert a == b

def test_bug3():
    x = Symbol("x")
    e = sqrt(x)
    assert str(e) == "x**(1/2)"

def test_bug4():
    x = Symbol("x")
    w = Symbol("w")
    e = -2*sqrt(x)-w/sqrt(x)/2
    assert str(e) not in ["(-2)*x**1/2(-1/2)*x**(-1/2)*w",
            "-2*x**1/2(-1/2)*x**(-1/2)*w","-2*x**1/2-1/2*x**-1/2*w"]
    assert str(e) in ["-2*x**(1/2) - 1/2*x**(-1/2)*w", "-2*x**(1/2) - 1/2*w*x**(-1/2)",
                      "-1/2*x**(-1/2)*w - 2*x**(1/2)", "-1/2*w*x**(-1/2) - 2*x**(1/2)"]

def test_Derivative():
    x = Symbol("x")
    y = Symbol("y")

    e = Derivative(x**2, x, evaluate=False)
    assert str(e) == "D(x**2, x)"

    e = Derivative(x**2/y, x, y, evaluate=False)
    assert str(e) == "D(x**2/y, x, y)"

def test_x_div_y():
    x = Symbol("x")
    y = Symbol("y")
    assert str(x/y) == "x/y"
    assert str(y/x) == "y/x"

def test_ordering():
    x = Symbol("x")
    assert str(sin(x).series(x, 0, 15)) == "x - 1/6*x**3 + (1/120)*x**5 - 1/5040*x**7 + (1/362880)*x**9 - 1/39916800*x**11 + (1/6227020800)*x**13 + O(x**15)"

def test_wild_str():
    # Check expressions containing Wild not causing infinite recursion
    a1 = Wild('a')
    a2 = Symbol('a')
    assert str(a1 + 1) == str(a2 + 1)
    assert str(exp(2**a1) + 5) == str(exp(2**a2) + 5)
    assert str(3*a1 + 1) == str(3*a2 + 1)
    assert str(1/a1 + 1) == str(1/a2 + 1)
    assert str(a1**2 + 1) == str(a2**2 + 1)
    assert str(1/(1-a1)) == str(1/(1-a2))

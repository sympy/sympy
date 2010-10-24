from sympy import residue, Symbol, Function, sin, S, I, pi, exp

def test_basic1():
    x = Symbol("x")
    assert residue(1/x, x, 0) == 1
    assert residue(-2/x, x, 0) == -2
    assert residue(81/x, x, 0) == 81
    assert residue(1/x**2, x, 0) == 0
    assert residue(0, x, 0) == 0
    assert residue(5, x, 0) == 0
    assert residue(x, x, 0) == 0
    assert residue(x**2, x, 0) == 0

def test_basic2():
    x = Symbol("x")
    assert residue(1/x, x, 1) == 0
    assert residue(-2/x, x, 1) == 0
    assert residue(81/x, x, -1) == 0
    assert residue(1/x**2, x, 1) == 0
    assert residue(0, x, 1) == 0
    assert residue(5, x, 1) == 0
    assert residue(x, x, 1) == 0
    assert residue(x**2, x, 5) == 0

def _test_f():
    # FIXME: we get infinite recursion here:
    x = Symbol("x")
    f = Function("f")
    assert residue(f(x)/x**5, x, 0) == f.diff(x, 4)/24

def test_functions():
    x = Symbol("x")
    assert residue(1/sin(x), x, 0) == 1
    assert residue(2/sin(x), x, 0) == 2
    assert residue(1/sin(x)**2, x, 0) == 0
    # FIXME: the series expansion fails to return the right answer:
    #assert residue(1/sin(x)**5, x, 0) == S(3)/8

def test_expressions():
    x = Symbol("x")
    assert residue(1/(x+1), x, 0) == 0
    assert residue(1/(x+1), x, -1) == 1
    assert residue(1/(x**2+1), x, -1) == 0
    assert residue(1/(x**2+1), x, I) == -I/2
    assert residue(1/(x**2+1), x, -I) == I/2
    assert residue(1/(x**4+1), x, 0) == 0
    # FIXME: this fails:
    #assert residue(1/(x**4+1), x, exp(I*pi/4)) == -(S(1)/4+I/4)/sqrt(2)

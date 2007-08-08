from sympy import symbols, lambdify, sqrt, sin, cos, Pi

x,y,z = symbols('xyz')

def test_no_vargs():
    f = lambdify(1, [])
    try:
        f(-1)
        raise Exception()
    except TypeError:
        pass
    assert f() == 1

def test_exponentiation():
    f = lambdify(x**2, [x])
    assert f(-1) == 1
    assert f(0) == 0
    assert f(1) == 1
    assert f(-2) == 4
    assert f(2) == 4
    assert f(2.5) == 6.25

def test_sqrt():
    f = lambdify(sqrt(x), [x])
    assert f(0) == 0.0
    assert f(1) == 1.0
    assert f(4) == 2.0
    assert abs(f(2) - 1.414) < 0.001
    assert f(6.25) == 2.5
    try:
        f(-1)
        raise Exception()
    except ValueError: pass

def test_vector_simple():
    f = lambdify([z,y,x], [x,y,z])
    assert f(3,2,1) == (1,2,3)
    assert f(1.0,2.0,3.0) == (3.0,2.0,1.0)
    # make sure correct number of args required
    try:
        f(0)
        raise Exception()
    except TypeError: pass

def test_vector_discontinuous():
    f = lambdify([-1/x, 1/x], [x])
    try:
        f(0)
        raise Exception()
    except ZeroDivisionError: pass
    assert f(1) == (-1.0, 1.0)
    assert f(2) == (-0.5, 0.5)
    assert f(-2) == (0.5, -0.5)

def test_trig_symbolic():
    f = lambdify([cos(x),sin(x)], [x])
    d = f(Pi)
    assert abs(d[0]+1) < 0.0001
    assert abs(d[1]-0) < 0.0001

def test_trig_float():
    f = lambdify([cos(x),sin(x)], [x])
    d = f(3.14159)
    print d
    assert abs(d[0]+1) < 0.0001
    assert abs(d[1]-0) < 0.0001

def test_bad_args():
    try:
        # no vargs given
        f = lambdify(1)
        raise Exception()
    except AssertionError: pass
    try:
        # same with vector exprs
        f = lambdify([1,2])
        raise Exception()
    except AssertionError: pass
    try:
        # vargs must be a list
        f = lambdify(1, x)
        raise Exception()
    except ValueError: pass
    try:
        # y is not in vargs
        f = lambdify(y, [x])
        raise Exception()
    except AssertionError: pass

def test_str_args():
    f = lambdify('z,y,x', 'x,y,z')
    assert f(3,2,1) == (1,2,3)
    assert f(1.0,2.0,3.0) == (3.0,2.0,1.0)
    # make sure correct number of args required
    try:
        f(0)
        raise Exception()
    except TypeError: pass

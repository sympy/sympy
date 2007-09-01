
from sympy import *

def test_re():

    x, y = symbols('xy')

    r = Symbol('r', real=True)

    assert re(nan) == nan

    assert re(oo) == oo
    assert re(-oo) == -oo

    assert re(0) == 0

    assert re(1) == 1
    assert re(-1) == -1

    assert re(E) == E
    assert re(-E) == -E

    assert re(x) == re(x)
    assert re(x*I) == -im(x)
    assert re(r*I) == 0
    assert re(r) == r

    assert re(x + y) == re(x + y)
    assert re(x + r) == re(x) + r

    assert re(re(x)) == re(x)

    assert re(2 + I) == 2
    assert re(x + I) == re(x)

    assert re(x + y*I) == re(x) - im(y)
    assert re(x + r*I) == re(x)

    assert re(log(2*I)) == log(2)

def test_im():

    x, y = symbols('xy')

    r = Symbol('r', real=True)

    assert im(nan) == nan

    assert im(oo*I) == oo
    assert im(-oo*I) == -oo

    assert im(0) == 0

    assert im(1) == 0
    assert im(-1) == 0

    assert im(E*I) == E
    assert im(-E*I) == -E

    assert im(x) == im(x)
    assert im(x*I) == re(x)
    assert im(r*I) == r
    assert im(r) == 0

    assert im(x + y) == im(x + y)
    assert im(x + r) == im(x)
    assert im(x + r*I) == im(x) + r

    assert im(im(x)*I) == im(x)

    assert im(2 + I) == 1
    assert im(x + I) == im(x) + 1

    assert im(x + y*I) == im(x) + re(y)
    assert im(x + r*I) == im(x) + r

    assert im(log(2*I)) == pi/2

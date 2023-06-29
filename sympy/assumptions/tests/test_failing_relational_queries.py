from sympy.testing.pytest import XFAIL
from sympy.assumptions.ask import ask, Q
from sympy.assumptions.assume import assuming
from sympy.assumptions.refine import refine
from sympy.logic.boolalg import Equivalent, Implies, And, Or
from sympy.core.relational import Gt
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.special.delta_functions import Heaviside
from sympy.abc import a, b, c, w, x, y, z


@XFAIL
def test_inequalities():
    # test basic defintions
    assert ask(Equivalent(Q.gt(x, 0), Q.positive(x))) is True
    assert ask(Equivalent(Q.eq(x, 0), Q.zero(x))) is True
    assert ask(Equivalent(Q.lt(x, 0), Q.negative(x))) is True
    assert ask(Equivalent(Q.ge(x, 0), Or(Q.positive(x), Q.zero(x)))) is True
    assert ask(Equivalent(Q.le(x, 0), Or(Q.negative(x), Q.zero(x)))) is True
    assert ask(Implies(And(Q.gt(x, y), Q.positive(y)), Q.positive(x))) is True
    assert ask(Implies(And(Q.lt(x, y), Q.negative(y)), Q.negative(x))) is True

    # test more complex problems
    assert ask(x > z, (x > y) & (y > z))


@XFAIL
def test_equality():
    # test substitution property of equality
    assert ask(Q.prime(x), Q.eq(x, y) & Q.prime(y)) is True
    assert ask(Q.real(x), Q.eq(x, y) & Q.real(y)) is True
    assert ask(Q.imaginary(x), Q.eq(x, y) & Q.imaginary(y)) is True


@XFAIL
def test_mixed():
    # test assumptions that mix inequalities and assumptions that the sat solver handles
    assert ask(x>0, Q.positive(x) & Q.prime(y))


@XFAIL
def test_relational():
    with assuming(x > 1):
        assert x > 1 == True
        assert x > 0 == True
        assert x < 0 == False


    #https://stackoverflow.com/questions/21958031/simplify-a-simple-inequity-with-sympy/21978199#21978199
    with assuming((x > y) & (x > 0) & (y > 0)):
        assert x+y < 2*x == True
        assert type(x > 2*y) == Gt


@XFAIL
def test_refine():
    # https://groups.google.com/g/sympy/c/tVo7iZx1ts0/m/qxRqBX0GAwAJ
    assert refine(z**2 + w**2 > 0, Q.positive(z) & Q.positive(w)) is True
    # inspired from https://stackoverflow.com/questions/19553652/sympy-limit-symbol-variable-to-interval/19579453#19579453
    assert refine(sqrt((x - 1) ** 2), x > 1) == x-1
    # https://stackoverflow.com/questions/67217022/simplify-expression-with-assumptions-involving-relations-between-variables-in-sy?rq=3
    assert refine(a+b-c, Q.eq(a+c, c)) == 0


@XFAIL
def test_piecewise():
    with assuming(x > 2):
        assert Heaviside(x) == 1
        assert Heaviside(x-1) == 1

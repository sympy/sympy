from sympy.testing.pytest import XFAIL
from sympy.assumptions.ask import ask, Q
from sympy.assumptions.assume import assuming
from sympy.assumptions.refine import refine
from sympy.logic.boolalg import Equivalent
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.special.delta_functions import Heaviside
from sympy.abc import a, b, c, w, x, y, z


@XFAIL
def test_inequalities():

    # inequalities only make sense over real numbers
    assert ask(Q.extended_real(x) & Q.extended_real(y), Q.gt(x, y)) is True

    # Test finite predicate 
    assert ask(Q.finite(x), Q.gt(x,-10) & Q.lt(x, 10)) is True
    assert ask(Q.finite(x), Q.lt(x 10)) is None
    
    # Test equivalent ways of writing Q.extended_positive, Q.zero, Q.extended_negative
    assert ask(Equivalent(Q.gt(x, 0), Q.extended_positive(x))) is True
    assert ask(Equivalent(Q.lt(x, 0), Q.extended_negative(x))) is True
    assert ask(Equivalent(Q.eq(x, 0), Q.zero(x))) is True
    assert ask(Equivalent(Q.ge(x, 0), Q.extended_positive(x) | Q.zero(x))) is True
    assert ask(Equivalent(Q.le(x, 0), Q.extended_negative(x) | Q.zero(x))) is True

    # test more complex problems
    assert ask(x > z, (x > y) & (y > z)) is True
    assert ask(x > z,  (x > w) & (w > y) & (y > z)) is True
    assert ask(x > z, ((x > y) & (y > z)) | ((x > w) & (w > y) & (y > z))) is True

    # test assumptions that mix inequalities and non-inequality unary assumptions
    assert ask(x > 0, Q.extended_positive(x) & Q.prime(y))

    #https://stackoverflow.com/questions/21958031/simplify-a-simple-inequity-with-sympy/21978199#21978199
    with assuming((x > y) & (x > 0) & (y > 0)):
        assert ask(x+y < 2*x) is True
        assert ask(x > 2*y) is None

@XFAIL
def test_number_line_properties():
    # From:
    # https://en.wikipedia.org/wiki/Inequality_(mathematics)#Properties_on_the_number_line

    # Converse property is currently supported

    # Transitivity
    # If a <= b and b <= c, then a <= c.
    assert ask(a <= c, (a <= b) & (b <= c)) is True
    # If a <= b and b < c, then a < c.
    assert ask(a < c, (a <= b) & (b < c)) is True
    # If a < b and b <= c, then a < c.
    assert ask(a < c, (a < b) & (b <= c)) is True

    # Addition and subtraction
    # If a <= b, then a + c <= b + c and a - c <= b - c.
    assert ask(a + c <= b + c, a <= b) is True
    assert ask(a - c <= b - c, a <= b) is True

    # Multiplication and division
    # If a <= b and c > 0, then ac <= bc and a/c <= b/c. (True for non-zero c)
    assert ask(a*c <= b*c, (a <= b) & (c > 0) & ~ Q.zero(c)) is True
    assert ask(a/c <= b/c, (a <= b) & (c > 0) & ~ Q.zero(c)) is True
    # If a <= b and c < 0, then ac >= bc and a/c >= b/c. (True for non-zero c)
    assert ask(a*c >= b*c, (a <= b) & (c < 0) & ~ Q.zero(c)) is True
    assert ask(a/c >= b/c, (a <= b) & (c < 0) & ~ Q.zero(c)) is True

    # Additive inverse
    # If a <= b, then -a >= -b.
    assert ask(-a >= -b, a <= b) is True

    # Multiplicative inverse
    # For a, b that are both negative or both positive:
    # If a <= b, then 1/a >= 1/b .
    assert ask(1/a >= 1/b, (a <= b) & Q.positive(x) & Q.positive(b)) is True
    assert ask(1/a >= 1/b, (a <= b) & Q.negative(x) & Q.negative(b)) is True


@XFAIL
def test_equality():
    # Currently reflexivity and symmetry are supported, but transitivity and substitution are not 
    
    # test transitive proprety
    assert ask(Q.eq(x,z), Q.eq(x,y) & Q.eq(y,z)) is True
    
    # test substitution property
    assert ask(Q.prime(x), Q.eq(x, y) & Q.prime(y)) is True
    assert ask(Q.real(x), Q.eq(x, y) & Q.real(y)) is True
    assert ask(Q.imaginary(x), Q.eq(x, y) & Q.imaginary(y)) is True


@XFAIL
def test_refine():
    # https://groups.google.com/g/sympy/c/tVo7iZx1ts0/m/qxRqBX0GAwAJ
    assert refine(z**2 + w**2 > 0, Q.positive(z) & Q.positive(w)) is True
    # inspired from https://stackoverflow.com/questions/19553652/sympy-limit-symbol-variable-to-interval/19579453#19579453
    assert refine(sqrt((x - 1) ** 2), x > 1) == x-1
    # https://stackoverflow.com/questions/67217022/simplify-expression-with-assumptions-involving-relations-between-variables-in-sy?rq=3
    assert refine(a+b-c, Q.eq(a+c, c)) == 0

    # test piecewise
    assert refine(Heaviside(x), x > 2) == 1
    assert refine(Heaviside(x-1), x > 2) == 1

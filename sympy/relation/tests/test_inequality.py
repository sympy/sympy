from sympy import Q, symbols
from sympy.printing import sstr, pretty, latex
from sympy.testing.pytest import raises

a,b,c,d,e = symbols('a b c d e')

eq = Q.eq(a,1)
gt = Q.gt(b,1)
ge = Q.ge(c,1)
lt = Q.lt(d,1)
le = Q.le(e,1)

def test_printing():
    assert sstr(gt) == "b > 1"
    assert pretty(gt) == "b > 1"
    assert latex(gt) == "b > 1"

    assert sstr(ge) == "c >= 1"
    assert pretty(ge) == "c >= 1"
    assert latex(ge) == r"c \geq 1"

    assert sstr(lt) == "d < 1"
    assert pretty(lt) == "d < 1"
    assert latex(lt) == "d < 1"

    assert sstr(le) == "e <= 1"
    assert pretty(le) == "e <= 1"
    assert latex(le) == r"e \leq 1"


def test_add():
    assert gt + 2 == Q.gt(b+2, 3)
    assert 3 + gt == Q.gt(3+b, 4)
    assert eq + gt == Q.gt(a+b, 2)
    assert gt + eq == Q.gt(b+a, 2)
    assert gt + gt == Q.gt(b+b, 2)

    assert ge + 2 == Q.ge(c+2, 3)
    assert 3 + ge == Q.ge(3+c, 4)
    assert eq + ge == Q.ge(a+c, 2)
    assert ge + eq == Q.ge(c+a, 2)
    assert ge  + ge == Q.ge(c+c, 2)

    assert gt + ge == Q.gt(b+c, 2)
    assert ge + gt == Q.gt(c+b, 2)

    assert lt + 2 == Q.lt(d+2, 3)
    assert 3 + lt == Q.lt(3+d, 4)
    assert eq + lt == Q.lt(a+d, 2)
    assert lt + eq == Q.lt(d+a, 2)
    assert lt + lt == Q.lt(d+d, 2)

    assert le + 2 == Q.le(e+2, 3)
    assert 3 + le == Q.le(3+e, 4)
    assert eq + le == Q.le(a+e, 2)
    assert le + eq == Q.le(e+a, 2)
    assert le  + le == Q.le(e+e, 2)

    assert lt + le == Q.lt(d+e, 2)
    assert le + lt == Q.lt(e+d, 2)

    with raises(TypeError):
        le + ge
        lt + ge
        le + gt
        lt + gt
        ge + le
        gt + le
        ge + lt
        gt + lt

def test_mul():
    assert -gt == Q.lt(-b, -1)
    assert gt*-2 == Q.lt(b*-2, -2)
    assert gt*3 == Q.gt(b*3, 3)
    assert -ge == Q.le(-c, -1)
    assert ge*-2 == Q.le(c*-2, -2)
    assert ge*3 == Q.ge(c*3, 3)
    assert -lt == Q.gt(-d, -1)
    assert lt*-2 == Q.gt(d*-2, -2)
    assert lt*3 == Q.lt(d*3, 3)
    assert -le == Q.ge(-e, -1)
    assert le*-2 == Q.ge(e*-2, -2)
    assert le*3 == Q.le(e*3, 3)

from sympy import (cse, expand_mul, I, Integral, Rational, solve, sqrt,
                   sqrtdenest, symbols, radsimp, root)
from sympy.simplify.sqrtdenest import denester
from sympy.utilities.pytest import XFAIL

def test_sqrtdenest():
    d = {sqrt(5 + 2 * sqrt(6)): sqrt(2) + sqrt(3),
        sqrt(sqrt(2)): sqrt(sqrt(2)),
        sqrt(5+sqrt(7)): sqrt(5+sqrt(7)),
        sqrt(3+sqrt(5+2*sqrt(7))):
            sqrt(6+3*sqrt(7))/(sqrt(2)*(5+2*sqrt(7))**Rational(1,4)) +
            3*(5+2*sqrt(7))**Rational(1,4)/(sqrt(2)*sqrt(6+3*sqrt(7))),
        sqrt(3+2*sqrt(3)): 3**Rational(1,4)/sqrt(2)+3/(sqrt(2)*3**Rational(1,4))}
    for i in d:
        assert sqrtdenest(i) == d[i] or denester([i])[0] == d[i]

    # this test caused a pattern recognition failure in sqrtdenest
    # nest = sqrt(2) + sqrt(5) - sqrt(7)
    nest = symbols('nest')
    x0, x1, x2, x3, x4, x5, x6 = symbols('x:7')
    l = sqrt(2) + sqrt(5)
    r = sqrt(7) + nest
    s = (l**2 - r**2).expand() + nest**2 # == nest**2
    ok = solve(nest**4 - s**2, nest)[1] # this will change if results order changes
    assert abs((l - r).subs(nest, ok).n()) < 1e-12
    x0 = sqrt(3)
    x2 = root(45*I*x0 - 28, 3)
    x3 = 19/x2
    x4 = x2 + x3
    x5 = -x4 - 14
    x6 = sqrt(-x5)
    ans = -x0*x6/3 + x0*sqrt(-x4 + 28 - 6*sqrt(210)*x6/x5)/3
    assert expand_mul(radsimp(ok) - ans) == 0
    # issue 2554
    eq = sqrt(sqrt(sqrt(2) + 2) + 2)
    assert sqrtdenest(eq) == eq

# more complex example:
@XFAIL # this fails on amd64
def test_sqrtdenest2():
    assert sqrtdenest(sqrt(16-2*sqrt(29)+2*sqrt(55-10*sqrt(29)))) == \
            sqrt(5) + sqrt(11-2*sqrt(29))

def test_issue_2758():
    from sympy.abc import x, y
    z = sqrt(1/(4*sqrt(3) + 7) + 1)
    ans = (sqrt(2) + sqrt(6))/(sqrt(3) + 2)
    assert sqrtdenest(z) == ans
    assert sqrtdenest(1 + z) == 1 + ans
    assert sqrtdenest(Integral(z + 1, (x, 1, 2))) == Integral(1 + ans, (x, 1, 2))
    assert sqrtdenest(x + sqrt(y)) == x + sqrt(y)

from itertools import combinations

from sympy import (ask, cos, exp, FiniteSet, Float, Function, I, log, oo,
    pi, Q, Rational, S, sin, sqrt, symbols)
from sympy.simplify import trigsimp
from sympy.printing import sstr, pretty, latex

x,y,z = symbols('x y z')

def test_printing():
    eq = Q.eq(x,y)
    assert sstr(eq) == "x = y"
    assert pretty(eq) == "x = y"
    assert latex(eq) == "x = y"

def test_add():
    eq1 = Q.eq(x,1)
    eq2 = Q.eq(y,2)
    assert eq1 + 2 == Q.eq(x+2, 3)
    assert 3 + eq1 == Q.eq(x+3, 4)
    assert eq1 + eq2 == Q.eq(x+y, 3)

def test_sub():
    eq1 = Q.eq(x,1)
    eq2 = Q.eq(y,2)
    assert eq1 - 2 == Q.eq(x-2, -1)
    assert 3 - eq1 == Q.eq(-x+3, 2)
    assert eq1 - eq2 == Q.eq(x-y, -1)

def test_mul():
    eq1 = Q.eq(x,1)
    eq2 = Q.eq(y,2)
    assert eq1*2 == Q.eq(x*2, 2)
    assert 3*eq1 == Q.eq(3*x, 3)
    assert eq1*eq2 == Q.eq(x*y, 2)

def test_mul_noncomm():
    a,b,c,d = symbols('a b c d', commutative=False)
    f = Function("f", commutative=False)
    eq = Q.eq(f(x), a/b)
    assert eq*a == Q.eq(f(x)*a, a*b**-1*a)
    assert a*eq == Q.eq(a*f(x), a**2*b**-1)

def test_pow():
    eq1 = Q.eq(x,1)
    eq2 = Q.eq(y,2)
    assert eq1**2 == Q.eq(x**2, 1)
    assert 3**eq1 == Q.eq(3**x, 3)
    assert eq1**eq2 == Q.eq(x**y, 1)

def test_truediv():
    eq1 = Q.eq(x,1)
    eq2 = Q.eq(y,2)
    assert eq1/2 == Q.eq(x/2, S.One/2)
    assert 3/eq1 == Q.eq(3/x, 3)
    assert eq1/eq2 == Q.eq(x/y, S.One/2)

def test_diff():
    eq = Q.eq(sin(x)*cos(y), (sin(x+y) + sin(x-y))/2)
    assert eq.diff(x) == Q.eq(
        (sin(x)*cos(y)).diff(x),
        ((sin(x+y) + sin(x-y))/2).diff(x)
    )
    assert eq.diff(y, y) == Q.eq(
        (sin(x)*cos(y)).diff(y, 2),
        ((sin(x+y) + sin(x-y))/2).diff(y, 2)
    )

def test_integrate():
    eq = Q.eq(sin(x)*cos(y), (sin(x+y) + sin(x-y))/2)
    assert eq.integrate(x) == Q.eq(
        (sin(x)*cos(y)).integrate(x),
        ((sin(x+y) + sin(x-y))/2).integrate(x)
    )
    assert eq.integrate((y, 0, 1)) == Q.eq(
        (sin(x)*cos(y)).integrate((y, 0, 1)),
        ((sin(x+y) + sin(x-y))/2).integrate((y, 0, 1))
    )

def test_func():
    eq = Q.eq(x,y)
    assert sin(eq) == Q.eq(sin(x), sin(y))
    assert cos(eq) == Q.eq(cos(x), cos(y))

def test_sqrt():
    eq = Q.eq(x,y)
    assert sqrt(eq) == Q.eq(sqrt(x), sqrt(y))

def test_rewrite():
    eq = Q.eq(y, exp(x*I))
    assert eq.rewrite(cos) == Q.eq(y, I*sin(x)+cos(x))

def test_simplify():
    eq = Q.eq(y, (x+x**2)/(x*sin(y)**2+x*cos(y)**2))
    assert eq.simplify() == Q.eq(x - y, -1)

    assert Q.eq(y, x).simplify() == Q.eq(x, y)
    assert Q.eq(x - 1, 0).simplify() == Q.eq(x, 1)
    assert Q.eq(x - 1, x).simplify() == Q.eq(0, 1)
    assert Q.eq(2*x - 1, x).simplify() == Q.eq(x, 1)
    assert Q.eq(2*x, 4).simplify() == Q.eq(x, 2)
    z = cos(1)**2 + sin(1)**2 - 1  # z.is_zero is None
    assert Q.eq(z*x, 0).simplify() == Q.eq(0, 0)
    assert Q.eq(x*(y + 1) - x*y - x + 1, x).simplify() == Q.eq(x, 1)

def test_doit():
    assert Q.eq(x, 0).doit() == Q.eq(x, 0)
    assert Q.eq(1, 1).doit() == Q.eq(1, 1)

def test_ask():
    assert ask(Q.eq(0, 0)) is True
    assert ask(Q.eq(1, 0)) is False
    assert ask(Q.eq(I, 2)) is False
    a = Float('.000000000000000000001', '')
    b = Float('.0000000000000000000001', '')
    assert ask(Q.eq(pi+a, pi+b)) is False
    assert ask(Q.eq(x, x)) is True

    assert ask(Q.eq(log(cos(2)**2 + sin(2)**2), 0)) is True

    assert ask(Q.eq(x, 1), Q.negative(x)) is False

def test_reversed():
    assert Q.eq(x, y).reversed == Q.eq(y, x)

def test_issue_10304():
    d = cos(1)**2 + sin(1)**2 - 1
    assert d.is_comparable is False  # if this fails, find a new d
    e = 1 + d*I
    assert Q.eq(e, 0).simplify() == Q.eq(0, 1)
    assert ask(Q.eq(e, 0)) is False

def test_issue_18412():
    d = (Rational(1, 6) + z / 4 / y)
    assert Q.eq(x, pi * y**3 * d).replace(y**3, z) == Q.eq(x, pi * z * d)

# XXX : Migrate test_issue_10401 from test_relational after
# inequality handlers are implemented

def test_issue_10633():
    assert ask(Q.eq(True, False)) is False
    assert ask(Q.eq(False, True)) is False
    assert ask(Q.eq(True, True)) is True
    assert ask(Q.eq(False, False)) is True

def test_issue_10927():
    assert ask(Q.eq(x, oo)) is None
    assert ask(Q.eq(x, -oo)) is None

# XXX : Migrate test_binary_symbols from test_relational after
# Equal.binary_symbols is implemented

def test_reversedsign():
    assert Q.eq(x,y).reversedsign == Q.eq(-x, -y)
    assert Q.eq(x,y).reversedsign.reversedsign == Q.eq(x,y)
    assert Q.eq(-x,y).reversedsign.reversedsign == Q.eq(-x,y)
    assert Q.eq(x,-y).reversedsign.reversedsign == Q.eq(x,-y)
    assert Q.eq(-x,-y).reversedsign.reversedsign == Q.eq(-x,-y)

def test_reversed_reversedsign():
    assert Q.eq(x, y).reversed.reversedsign == Q.eq(x, y).reversedsign.reversed
    assert Q.eq(-x, y).reversed.reversedsign == Q.eq(-x, y).reversedsign.reversed
    assert Q.eq(x, -y).reversed.reversedsign == Q.eq(x, -y).reversedsign.reversed
    assert Q.eq(-x, -y).reversed.reversedsign == Q.eq(-x, -y).reversedsign.reversed

def test_improved_canonical():
    def test_different_forms(listofforms):
        for form1, form2 in combinations(listofforms, 2):
            assert form1.canonical == form2.canonical

    def generate_forms(expr):
        return [expr, expr.reversed, expr.reversedsign,
                expr.reversed.reversedsign]

    test_different_forms(generate_forms(Q.eq(x, -y)))

def test_set_equality_canonical():
    a, b, c = symbols('a b c')
    A = Q.eq(FiniteSet(a, b, c), FiniteSet(1, 2, 3))
    assert A.canonical == A.reversed

def test_trigsimp():
    # issue 16736
    s, c = sin(2*x), cos(2*x)
    eq = Q.eq(s, c)
    assert trigsimp(eq) == eq  # no rearrangement of sides
    # simplification of sides might result in
    # an unevaluated Eq
    changed = trigsimp(Q.eq(s + c, sqrt(2)))
    assert ask(changed.subs(x, pi/8)) is True
    # or an evaluated one
    assert ask(trigsimp(Q.eq(cos(x)**2 + sin(x)**2, 1))) is True

def test_multivariate_linear_function_simplification():
    assert Q.eq(2*x + y, 2*x + y - 3).simplify() == Q.eq(0, 3)

def test_nonpolymonial_relations():
    assert Q.eq(cos(x), 0).simplify() == Q.eq(cos(x), 0)

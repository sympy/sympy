from sympy import (diff, expand, Eq, Integral, integrate, Interval, lambdify,
                   log, oo, Piecewise, piecewise_fold, pi, S, symbols, solve,
                   Rational)
from sympy.utilities.pytest import XFAIL, raises

x,y = symbols('x,y')

def test_piecewise():

    # Test canonization
    assert Piecewise((x, x > 1)) == Piecewise((x, x > 1), S.NaN)
    assert Piecewise((x, x > 1), (0, True), S.NaN) == 0
    assert Piecewise((x, x > 1), (0, False)) == S.NaN
    assert Piecewise((x, x > 1), (0, False), 1) == 1
    assert Piecewise((x, True), S.NaN) == x
    raises(TypeError,"Piecewise((x,x**2))")
    raises(TypeError,"Piecewise((1,Interval(0,1,False,True)),(0,x<1))")
    raises(TypeError,"Piecewise((x,1))")
    raises(TypeError,"Piecewise((None, x>0))")
    raises(TypeError,"Piecewise((x, x>0),None)")
    raises(ValueError,"Piecewise((x,))")
    raises(ValueError,"Piecewise((x,x>1,x>0))")
    raises(ValueError,"Piecewise(x,(x,x<1))")
    raises(ValueError,"Piecewise(x,x)")

    # Test properties
    p1 = Piecewise((x, x > 1), x**2)
    p2 = Piecewise((x, x > 1), (x**2, x <=1))
    assert p1.exprcondpairs == ((x, x > 1),)
    assert p1.otherwise == x**2
    assert p2.exprcondpairs == ((x, x > 1),(x**2, x <= 1))
    assert p2.otherwise == S.NaN

    # Test subs
    p = Piecewise((-1, x < -1), (x**2, x < 0), (log(x), x >=0))
    p_x2 = Piecewise((-1, x**2 < -1), (x**4, x**2 < 0), (log(x**2), x**2 >=0))
    assert p.subs(x,x**2) == p_x2
    assert p.subs(x,-5) == -1
    assert p.subs(x,-1) == 1
    assert p.subs(x,1) == log(1)

    # More subs tests
    p2 = Piecewise((1, x < pi), (-1, x < 2*pi), (0, x > 2*pi))
    assert p2.subs(x, 2) == 1
    assert p2.subs(x,4) == -1
    assert p2.subs(x,10) == 0

    # More subs tests
    p3 = Piecewise( (x**2, Interval(-3, 3).contains(x)) )
    assert p3.subs(x,-2) == 4
    assert p3.subs(x,4) == S.NaN

    # Test evalf
    assert p.evalf() == p
    assert p.evalf(subs={x:-2}) == -1
    assert p.evalf(subs={x:-1}) == 1
    assert p.evalf(subs={x:1}) == log(1)

    # Test doit
    f_int = Piecewise((Integral(x,(x,0,1)), x < 1))
    assert f_int.doit() == Piecewise( (1.0/2.0, x < 1) )

    # Test differentiation
    f = x
    fp = x*p
    dp = Piecewise((0, x < -1), (2*x, x < 0), (1/x, x >= 0))
    fp_dx = x*dp + p
    assert diff(p,x) == dp
    assert diff(f*p,x) == fp_dx

    # Test simple arithmetic
    assert x*p == fp
    assert x*p + p == p + x*p
    assert p + f == f + p
    assert p + dp == dp + p
    assert p - dp == -(dp - p)

    # Test _eval_interval
    f1 = x*y + 2
    f2 = x*y**2 + 3
    peval = Piecewise( (f1, x<0), (f2, x>0))
    peval_interval = f1.subs(x,0) - f1.subs(x,-1) + f2.subs(x,1) - f2.subs(x,0)
    assert peval._eval_interval(x, -1, 1) == peval_interval

    # Test integration
    p_int =  Piecewise((-x,x < -1), (x**3/3.0, x < 0), (-x + x*log(x), x >= 0))
    assert integrate(p,x) == p_int
    p = Piecewise((x, x < 1),(x**2, -1 <= x),(x,3<x))
    assert integrate(p,(x,-2,2)) == 5.0/6.0
    assert integrate(p,(x,2,-2)) == -5.0/6.0
    p = Piecewise((0, x < 0), (1,x < 1), (0, x < 2), (1, x < 3), 0)
    assert integrate(p, (x,-oo,oo)) == 2
    p = Piecewise((x, x < -10),(x**2, x <= -1),(x, 1 < x))
    assert integrate(p,(x,-2,2)) == S.NaN

    # Test commutativity
    assert p.is_commutative is True

    # Test leading term
    # This may not be the best behavior, but is necessary to pass tests
    # See XFAIL test_piecewise_leading_term for better behavior
    p = Piecewise((x**3 + x, x < -1), (x**2 - 1, x < 0), log(x))
    assert p.as_leading_term(x) == log(x)

@XFAIL
def test_piecewise_leading_term():
    # Test as_leading_term
    p = Piecewise((x**3 + x, x < -1), (x**2 - 1, x < 0), log(x))
    assert p.as_leading_term(x) == Piecewise((x, x < -1), (-1, x < 0), log(x))

def test_piecewise_free_symbols():
    a = symbols('a')
    f = Piecewise((x , a<0), y)
    assert f.free_symbols == set([x, y, a])

def test_piecewise_integrate():
    f = Piecewise(((x - 2)**2, x >= 0), 1)
    assert integrate(f, (x, -2, 2)) == Rational(14, 3)

    g = Piecewise(((x - 5)**5, x >= 4), f)
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == Rational(43, 6)

    g = Piecewise(((x - 5)**5, x >= 4), (f, x < 4))
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == Rational(43, 6)

    g = Piecewise(((x - 5)**5, x >= 2), (f, x < 2))
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == -Rational(701, 6)

    g = Piecewise(((x - 5)**5, x >= 2), f)
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == -Rational(701, 6)

    g = Piecewise(((x - 5)**5, x >= 2), 2 * f)
    assert integrate(g, (x, -2, 2)) == 2 * Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == -Rational(673, 6)

    g = Piecewise((1, x > 0), (0, Eq(x, 0)), (-1, x < 0))
    assert integrate(g, (x, -1, 1)) == 0

    g = Piecewise((1, x - y < 0), 0)
    assert integrate(g, (y, -oo, oo)) == oo

def test_piecewise_solve():
    abs2 = Piecewise((-x, x <= 0), (x, x > 0))
    f = abs2.subs(x, x - 2)
    assert solve(f, x) == [2]
    assert solve(f - 1,x) == [1, 3]

    f = Piecewise(((x - 2)**2, x >= 0), 1)
    assert solve(f, x) == [2]

    g = Piecewise(((x - 5)**5, x >= 4), f)
    assert solve(g, x) == [2, 5]

    g = Piecewise(((x - 5)**5, x >= 4), (f, x < 4))
    assert solve(g, x) == [2, 5]

    g = Piecewise(((x - 5)**5, x >= 2), (f, x < 2))
    assert solve(g, x) == [5]

    g = Piecewise(((x - 5)**5, x >= 2), f)
    assert solve(g, x) == [5]

    g = Piecewise(((x - 5)**5, x >= 2), (f, True), (10, False), evaluate=False)
    assert solve(g, x) == [5]

    assert solve(Piecewise((x, x < 0), x), x) == [0]

# See issue 1253 (enhance the solver to handle inequalities).
@XFAIL
def test_piecewise_solve2():
    f = Piecewise(((x - 2)**2, x >= 0), 0)
    assert solve(f, x) == [2, Interval(0, oo, True, True)]

def test_piecewise_fold():

    p = Piecewise((x, x < 1), (1, 1 <= x))

    assert piecewise_fold(x*p) == Piecewise((x**2, x < 1), (x, 1 <= x))
    assert piecewise_fold(p+p) == Piecewise((2*x, x < 1), (2, 1 <= x))
    assert piecewise_fold(Piecewise((1, x < 0), 2) \
                          + Piecewise((10, x < 0), -10)) == \
           Piecewise((11, x < 0), -8)

    p1 = Piecewise((0, x < 0), (x, x <= 1), 0)
    p2 = Piecewise((0, x < 0), (1 - x, x <=1), 0)

    p = 4*p1 + 2*p2
    assert integrate(piecewise_fold(p),(x,-oo,oo)) == integrate(2*x + 2, (x, 0, 1))

def test_piecewise_fold_expand():
    p1 = Piecewise((1,Interval(0,1,False,True).contains(x)),0)

    #p2 = piecewise_fold(expand((1-x)*p1))
    #assert p2 == Piecewise((1 - x, Interval(0,1,False,True).contains(x)), \
    #    Piecewise((-x, Interval(0,1,False,True).contains(x)), 0))

    p2 = expand(piecewise_fold((1-x)*p1))
    assert p2 == Piecewise((1 - x, Interval(0,1,False,True).contains(x)), 0)

def test_piecewise_duplicate():
    p = Piecewise((x, x < -10),(x**2, x <= -1),(x, 1 < x))
    assert p == Piecewise(*p.args)

def test_doit():
    p1 = Piecewise((x, x < 1), (x**2, -1 <= x), (x, 3 < x))
    p2 = Piecewise((x, x < 1), (Integral(2 * x), -1 <= x), (x, 3 < x))
    assert p2.doit() == p1
    assert p2.doit(deep = False) == p2

def test_piecewise_interval():
    p1 = Piecewise((x, Interval(0,1).contains(x)), 0)
    assert p1.subs(x, -0.5) == 0
    assert p1.subs(x, 0.5) == 0.5
    assert p1.diff(x) == Piecewise((1, Interval(0, 1).contains(x)), 0)
    assert integrate(p1, x) == Piecewise((x**2/2, Interval(0, 1).contains(x)), 0)

def test_piecewise_collapse():
    p1 = Piecewise((x, x<0),(x**2,x>1))
    p2 = Piecewise((p1,x<0),(p1,x>1))
    assert p2.doit() == Piecewise((x, x < 0), (x**2, 1 < x))

    p1 = Piecewise(Piecewise((x,x<0),1))
    assert p1 == Piecewise(Piecewise((x,x<0),1))

def test_piecewise_lambdify():
    p = Piecewise(
        (x**2, x < 0),
        (x, Interval(0, 1, False, True).contains(x)),
        (2 - x, x >= 1),
        0
    )
    f = lambdify(x, p)
    assert f(-2.0) == 4.0
    assert f(0.0) == 0.0
    assert f(0.5) == 0.5
    assert f(2.0) == 0.0

def test_piecewise_series():
    from sympy import sin, cos, O
    p1 = Piecewise((sin(x), x<0),(cos(x),x>0))
    p2 = Piecewise((x+O(x**2), x<0),(1+O(x**2),x>0))
    assert p1.nseries(x,n=2) == p2

def test_piecewise_evaluate():
    assert Piecewise(x) != x
    assert Piecewise(x).is_Piecewise
    assert Piecewise(x, evaluate=True) == x

    assert Piecewise((x, True), S.NaN) == x
    assert Piecewise((x, True), S.NaN, evaluate=False) == Piecewise(x)

    assert Piecewise((x, x < 1), (0, True), S.NaN) == 0
    assert Piecewise((x, x < 1), (0, True), S.NaN, evaluate=False) == Piecewise((x, x < 1), 0)

    assert Piecewise((x, 1 > 2), (-x, False)) == S.NaN
    p1 = Piecewise((x, x > 0), (2*x, x <= 0), 3*x)
    p2 = Piecewise((p1, x > 0), (2*x**2, x <= 0))
    p3 = Piecewise((x**2, x > 0), (2*x**2, x <= 0), p1)
    assert p2 == Piecewise((x, x > 0), (2*x**2, x <=0))
    assert p3 == Piecewise((x**2, x > 0), (2*x**2, x<=0), 3*x)
    p2 = Piecewise((p1, x > 0), (2*x**2, x <= 0), evaluate=False)
    p3 = Piecewise((x**2, x > 0), (2*x**2, x <= 0), p1, evaluate=False)
    assert p2.exprcondpairs[0][0].is_Piecewise
    assert p3.otherwise.is_Piecewise

@XFAIL
def test_piecewise_evaluate_false():
    # Requires Issue 3025
    assert Piecewise((x, True), S.NaN, evaluate=False).is_Piecewise
    assert Piecewise((x, 1 < 2), S.NaN, evaluate=False).is_Piecewise
    assert Piecewise((x, 1 > 2), S.NaN, evaluate=False).is_Piecewise
    p = Piecewise((x, x < 1), (0, True), S.NaN, evaluate=False)
    assert p.is_Piecewise
    assert len(p.args) == 3
    assert Piecewise((x, 1 > 2), (-x, False), evaluate=False).is_Piecewise

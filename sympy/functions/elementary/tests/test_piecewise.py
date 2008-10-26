from sympy import diff, Integral, integrate, log, oo, Piecewise, \
    piecewise_fold, raises, symbols

x,y = symbols('xy')

def test_piecewise():

    # Test canonization
    assert Piecewise((x, x < 1), (0, True)) == Piecewise((x, x < 1), (0, True))
    assert Piecewise((x, x < 1), (0, False), (-1, 1>2)) == Piecewise((x, x < 1))
    assert Piecewise((x, True)) == x
    raises(TypeError,"Piecewise(x)")
    raises(TypeError,"Piecewise((x,x**2))")

    # Test subs
    p = Piecewise((-1, x < -1), (x**2, x < 0), (log(x), x >=0))
    p_x2 = Piecewise((-1, x**2 < -1), (x**4, x**2 < 0), (log(x**2), x**2 >=0))
    assert p.subs(x,x**2) == p_x2
    assert p.subs(x,-5) == -1
    assert p.subs(x,-1) == 1
    assert p.subs(x,1) == log(1)

    # More subs test
    p2 = Piecewise((1, x < pi), (-1, x < 2*pi), (0, x > 2*pi))
    assert p2.subs(x,2) == 1
    assert p2.subs(x,4) == -1
    assert p2.subs(x,10) == 0

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
    p = Piecewise((0, x < 0), (1,x < 1), (0, x < 2), (1, x < 3), (0, True))
    assert integrate(p, (x,-oo,oo)) == 2
    p = Piecewise((x, x < -10),(x**2, x <= -1),(x, 1 < x))
    raises(ValueError, "integrate(p,(x,-2,2))")

def test_piecewise_fold():

    p = Piecewise((x, x < 1), (1, 1 <= x))

    assert piecewise_fold(x*p) == Piecewise((x**2, x < 1), (x, 1 <= x))
    assert piecewise_fold(p+p) == Piecewise((2*x, x < 1), (2, 1 <= x))

    p1 = Piecewise((0, x < 0), (x, x <= 1), (0, True))
    p2 = Piecewise((0, x < 0), (1 - x, x <=1), (0, True))

    p = 4*p1 + 2*p2
    assert integrate(piecewise_fold(p),(x,-oo,oo)) == integrate(2*x + 2, (x, 0, 1))

from sympy import S, cos, sin, pi, sqrt, Function
from sympy.physics.continuum_mechanics.column import Column
from sympy.utilities.pytest import XFAIL, raises
from sympy.core.symbol import Symbol, symbols


E, I  = symbols('E, I', positive=True)
l = Symbol('l')
P = Symbol('P', positive=True)

def test_column():
    y = Function('y')
    x = Symbol('x')
    M = Symbol('M')
    d = Symbol('d')
    F = Symbol('F')
    C1, C2 = symbols('C1, C2')

    # test for pinned-pinned end-condition
    c1 = Column(l, E, I, P, top="pinned", bottom="pinned")
    c1.solve_slope_deflection()
    assert c1.moment() == P*y(x)
    assert c1.deflection() == C1*sin(sqrt(P)*x/(sqrt(E)*sqrt(I)))
    assert c1.slope() == C1*sqrt(P)*cos(sqrt(P)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I))
    assert c1.critical_load() == pi**2*E*I/l**2

    # test for fixed-fixed end-condition
    c2 = Column(l, E, I, P, top="fixed", bottom="fixed")
    c2.solve_slope_deflection()
    assert c2.moment() == -M + P*y(x)
    assert c2.deflection() == -M*cos(sqrt(P)*x/(sqrt(E)*sqrt(I)))/P + M/P
    assert c2.slope() == M*sin(sqrt(P)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I)*sqrt(P))
    assert c2.critical_load() == 4*pi**2*E*I/l**2

    # test for free-fixed end-condition
    c3 = Column(l, E, I, P, top="free", bottom="fixed")
    c3.solve_slope_deflection()
    assert c3.moment() == -P*d + P*y(x)
    assert c3.deflection() == -d*cos(sqrt(P)*x/(sqrt(E)*sqrt(I))) + d
    assert c3.slope() == sqrt(P)*d*sin(sqrt(P)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I))
    assert c3.critical_load() == pi**2*E*I/(4*l**2)

    # test for pinned-fixed end-condition
    c4 = Column(l, E, I, P, top="pinned", bottom="fixed")
    c4.solve_slope_deflection()
    assert c4.moment() == -F*(l - x) + P*y(x)
    assert c4.deflection() == sqrt(E)*F*sqrt(I)*sin(sqrt(P)*x/(sqrt(E)*sqrt(I)))/P**(3/S(2)) - F*l*cos(sqrt(P)*x/(sqrt(E)*sqrt(I)))/P + F*l/P - F*x/P
    assert c4.slope() == F*cos(sqrt(P)*x/(sqrt(E)*sqrt(I)))/P - F/P + F*l*sin(sqrt(P)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I)*sqrt(P))

    raises(ValueError, lambda: Column(l, E, I, P, top="free", bottom="free"))

@XFAIL
def test_critical_load_pinned_fixed():
    # the deflction equation of pinned-fixed end-condition
    # comes out to be of the form `a*cos(x) - b*sin(x)` for which
    # solve should return in the form `x = atan(a/b)`. solve() currently
    # gives the answer in other form. Either we can get a way to convert it
    # into desired form or make solve() return it in the required form.
    c = Column(l, E, I, P, top="pinned", bottom="fixed")
    c.solve_slope_deflection()
    c.critical_load() == 2*pi**2*E*I/l**2

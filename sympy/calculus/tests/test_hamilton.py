from sympy import Symbol, Function, Derivative as D, Eq, cos, sin
from sympy.utilities.pytest import raises
from sympy.calculus.hamilton import hamilton_equations as hamilton


def test_hamilton_errors():
    # Test by giving the Hamiltonian Interface bad-functions

	# ERROR CASE: 1
	# =============
    x = Function('x')
    p = Function('p')
    t = Symbol('t')
    h1 = x(t)**2 / 2 + p(t)**2 / 2
    h2 = x(p)**2 / 2 + p(t)**2 / 2
    raises(TypeError, lambda: hamilton(h1, x(t), p, t))
    raises(TypeError, lambda: hamilton(h2, x(t), p(t), t))
	
	# ERROR CASE: 2
	# =============
    x = Function('x')
    p = Function('p')
    t = Symbol('t')
    h = x(t) ** 2 / 2 + p(t) ** 2 / 2
    # Giving function instead of symbol as time
    raises(TypeError, lambda: hamilton(h, [x(t)], [p(t)], x(t)))

	# ERROR CASE: 3
	# =============
    x = Function('x')
    p = Function('p')
    k = Symbol('k')
    t = Symbol('t')
    h = x(k) ** 2 / 2 + p(k) ** 2 / 2
    # Giving function instead of symbol as time
    raises(ValueError, lambda: hamilton(h, [x(k)], [p(k)], t))
	
	# ERROR CASE: 4
	# =============
    x = Symbol('x')
    p = Symbol('p')
    t = Symbol('t')
    h1 = x ** 2 / 2 + p ** 2 / 2
    # Giving symbol instead of function as the Hamiltonian
    raises(TypeError, lambda: hamilton(h1, [x], [p], t))


def test_hamilton_shm1():
    # Test by solving for SHM in 1-Dimension
    x = Function('x')
    p = Function('p')
    t = Symbol('t')
    h1 = x(t)**2 / 2 + p(t)**2 / 2
    assert hamilton(h1, x(t), p(t), t) == [
        Eq(x(t) + D(p(t), t), 0), Eq(p(t) - D(x(t), t), 0)
    ]


def test_hamilton_shm2():
    # Test by solving for SHM in 2-Dimensions
    x = Function('x')
    p = Function('p')
    t = Symbol('t')
    y = Function('y')
    q = Function('q')
    h2 = x(t)**2 / 2 + y(t)**2 / 2 + p(t)**2 / 2 + q(t)**2 / 2
    assert hamilton(h2, [x(t), y(t)], [p(t), q(t)], t) == [
        Eq(x(t) + D(p(t), t), 0), Eq(p(t) - D(x(t), t), 0),
        Eq(y(t) + D(q(t), t), 0), Eq(q(t) - D(y(t), t), 0)
    ]


def test_hamilton_shm3():
    # Test by solving for SHM over a cylinder
    t = Symbol('t')
    qr = Function('R')
    qt = Function('Theta')
    qz = Function('Z')
    pr = Function('pR')
    pt = Function('pTheta')
    pz = Function('pZ')
    h3 = pt(t)**2/2 + pz(t)**2/2 + qz(t)**2/2
    assert hamilton(h3, [qr(t), qt(t), qz(t)], [pr(t), pt(t), pz(t)], t) == [
        Eq(D(pr(t), t), 0), Eq(-D(qr(t), t), 0), Eq(D(pt(t), t), 0),
        Eq(pt(t) - D(qt(t), t), 0), Eq(qz(t) + D(pz(t), t), 0),
        Eq(pz(t) - D(qz(t), t), 0)
    ]

import sympy.physics.mechanics.models as models
from sympy import (cos, sin, Matrix, symbols)
from sympy.physics.mechanics import (dynamicsymbols)


def test_multi_mass_spring_damper_inputs():

    c0, k0, m0 = symbols("c0 k0 m0")
    g = symbols("g")
    v0, x0, f0 = dynamicsymbols("v0 x0 f0")

    kane1 = models.multi_mass_spring_damper(1)
    assert kane1.rhs() == Matrix([[v0], [(-c0*v0 - k0*x0)/m0]])

    kane2 = models.multi_mass_spring_damper(1, True)
    assert kane2.rhs() == Matrix([[v0], [(-c0*v0 + g*m0 - k0*x0)/m0]])

    kane3 = models.multi_mass_spring_damper(1, True, True)
    assert kane3.rhs() == Matrix([[v0], [(-c0*v0 + g*m0 - k0*x0 + f0)/m0]])

    kane4 = models.multi_mass_spring_damper(1, False, True)
    assert kane4.rhs() == Matrix([[v0], [(-c0*v0 - k0*x0 + f0)/m0]])


def test_multi_mass_spring_damper_higher_order():
    c0, k0, m0 = symbols("c0 k0 m0")
    c1, k1, m1 = symbols("c1 k1 m1")
    c2, k2, m2 = symbols("c2 k2 m2")
    v0, x0 = dynamicsymbols("v0 x0")
    v1, x1 = dynamicsymbols("v1 x1")
    v2, x2 = dynamicsymbols("v2 x2")

    kane1 = models.multi_mass_spring_damper(3)
    assert kane1.rhs() == Matrix([[v0], [v1], [v2], [(-c0*v0 - k0*x0 - m2*(-c2*v2 - k2*x2 - m2*(-c0*v0 - k0*x0)/(m0 + m1 + m2) - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)*(-c1*v1 - k1*x1 - (m1 + m2)*(-c0*v0 - k0*x0)/(m0 + m1 + m2))/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2)))/(-m2**2/(m0 + m1 + m2) + m2 - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)**2/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2))) - (m1 + m2)*(-c1*v1 - k1*x1 - (m1 + m2)*(-c0*v0 - k0*x0)/(m0 + m1 + m2) - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)*(-c2*v2 - k2*x2 - m2*(-c0*v0 - k0*x0)/(m0 + m1 + m2) - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)*(-c1*v1 - k1*x1 - (m1 + m2)*(-c0*v0 - k0*x0)/(m0 + m1 + m2))/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2)))/(-m2**2/(m0 + m1 + m2) + m2 - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)**2/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2))))/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2)))/(m0 + m1 + m2)], [(-c1*v1 - k1*x1 - (m1 + m2)*(-c0*v0 - k0*x0)/(m0 + m1 + m2) - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)*(-c2*v2 - k2*x2 - m2*(-c0*v0 - k0*x0)/(m0 + m1 + m2) - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)*(-c1*v1 - k1*x1 - (m1 + m2)*(-c0*v0 - k0*x0)/(m0 + m1 + m2))/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2)))/(-m2**2/(m0 + m1 + m2) + m2 - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)**2/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2))))/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2))], [(-c2*v2 - k2*x2 - m2*(-c0*v0 - k0*x0)/(m0 + m1 + m2) - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)*(-c1*v1 - k1*x1 - (m1 + m2)*(-c0*v0 - k0*x0)/(m0 + m1 + m2))/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2)))/(-m2**2/(m0 + m1 + m2) + m2 - (-m2*(m1 + m2)/(m0 + m1 + m2) + m2)**2/(m1 + m2 - (m1 + m2)**2/(m0 + m1 + m2)))]])


def test_n_link_pendulum_on_cart_inputs():
    l0, m0 = symbols("l0 m0")
    m1 = symbols("m1")
    g = symbols("g")
    q0, q1, F, T1 = dynamicsymbols("q0 q1 F T1")
    u0, u1 = dynamicsymbols("u0 u1")

    kane1 = models.n_link_pendulum_on_cart(1)
    assert kane1.rhs() == Matrix([[u0], [u1], [(l0*m1*(-g*l0*m1*cos(q1) + l0*m1*(l0*m1*u1**2*cos(q1) + F)*sin(q1)/(m0 + m1))*sin(q1)/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1) + l0*m1*u1**2*cos(q1) + F)/(m0 + m1)], [(-g*l0*m1*cos(q1) + l0*m1*(l0*m1*u1**2*cos(q1) + F)*sin(q1)/(m0 + m1))/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1)]])

    kane2 = models.n_link_pendulum_on_cart(1, False)
    assert kane2.rhs() == Matrix([[u0], [u1], [(l0*m1*(-g*l0*m1*cos(q1) + l0**2*m1**2*u1**2*sin(q1)*cos(q1)/(m0 + m1))*sin(q1)/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1) + l0*m1*u1**2*cos(q1))/(m0 + m1)], [(-g*l0*m1*cos(q1) + l0**2*m1**2*u1**2*sin(q1)*cos(q1)/(m0 + m1))/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1)]])

    kane3 = models.n_link_pendulum_on_cart(1, False, True)
    assert kane3.rhs() == Matrix([[u0], [u1], [(l0*m1*u1**2*cos(q1) + l0*m1*(-g*l0*m1*cos(q1) + l0**2*m1**2*u1**2*sin(q1)*cos(q1)/(m0 + m1) + T1)*sin(q1)/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1))/(m0 + m1)], [(-g*l0*m1*cos(q1) + l0**2*m1**2*u1**2*sin(q1)*cos(q1)/(m0 + m1) + T1)/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1)]])

    kane4 = models.n_link_pendulum_on_cart(1, True, False)
    assert kane4.rhs() == Matrix([[u0], [u1], [(l0*m1*(-g*l0*m1*cos(q1) + l0*m1*(l0*m1*u1**2*cos(q1) + F)*sin(q1)/(m0 + m1))*sin(q1)/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1) + l0*m1*u1**2*cos(q1) + F)/(m0 + m1)], [(-g*l0*m1*cos(q1) + l0*m1*(l0*m1*u1**2*cos(q1) + F)*sin(q1)/(m0 + m1))/(-l0**2*m1**2*sin(q1)**2/(m0 + m1) + l0**2*m1)]])


def test_n_link_pendulum_on_cart_higher_order():
    l0, m0 = symbols("l0 m0")
    l1, m1 = symbols("l1 m1")
    m2 = symbols("m2")
    g = symbols("g")
    q0, q1, q2 = dynamicsymbols("q0 q1 q2")
    u0, u1, u2 = dynamicsymbols("u0 u1 u2")
    F, T1 = dynamicsymbols("F T1")

    kane1 = models.n_link_pendulum_on_cart(2)
    assert kane1.rhs() == Matrix([[u0], [u1], [u2], [(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + l1*m2*(-g*l1*m2*cos(q2) + l0*l1*m2*(sin(q1)*cos(q2) - sin(q2)*cos(q1))*u1**2 + l1*m2*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)*sin(q2)/(m0 + m1 + m2) - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))*(-g*l0*m1*cos(q1) - g*l0*m2*cos(q1) + l0*l1*m2*(-sin(q1)*cos(q2) + sin(q2)*cos(q1))*u2**2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)/(m0 + m1 + m2))/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2)))*sin(q2)/(-l1**2*m2**2*sin(q2)**2/(m0 + m1 + m2) + l1**2*m2 - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))**2/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2))) - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(-g*l0*m1*cos(q1) - g*l0*m2*cos(q1) + l0*l1*m2*(-sin(q1)*cos(q2) + sin(q2)*cos(q1))*u2**2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)/(m0 + m1 + m2) - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))*(-g*l1*m2*cos(q2) + l0*l1*m2*(sin(q1)*cos(q2) - sin(q2)*cos(q1))*u1**2 + l1*m2*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)*sin(q2)/(m0 + m1 + m2) - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))*(-g*l0*m1*cos(q1) - g*l0*m2*cos(q1) + l0*l1*m2*(-sin(q1)*cos(q2) + sin(q2)*cos(q1))*u2**2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)/(m0 + m1 + m2))/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2)))/(-l1**2*m2**2*sin(q2)**2/(m0 + m1 + m2) + l1**2*m2 - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))**2/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2))))/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2)) + F)/(m0 + m1 + m2)], [(-g*l0*m1*cos(q1) - g*l0*m2*cos(q1) + l0*l1*m2*(-sin(q1)*cos(q2) + sin(q2)*cos(q1))*u2**2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)/(m0 + m1 + m2) - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))*(-g*l1*m2*cos(q2) + l0*l1*m2*(sin(q1)*cos(q2) - sin(q2)*cos(q1))*u1**2 + l1*m2*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)*sin(q2)/(m0 + m1 + m2) - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))*(-g*l0*m1*cos(q1) - g*l0*m2*cos(q1) + l0*l1*m2*(-sin(q1)*cos(q2) + sin(q2)*cos(q1))*u2**2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)/(m0 + m1 + m2))/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2)))/(-l1**2*m2**2*sin(q2)**2/(m0 + m1 + m2) + l1**2*m2 - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))**2/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2))))/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2))], [(-g*l1*m2*cos(q2) + l0*l1*m2*(sin(q1)*cos(q2) - sin(q2)*cos(q1))*u1**2 + l1*m2*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)*sin(q2)/(m0 + m1 + m2) - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))*(-g*l0*m1*cos(q1) - g*l0*m2*cos(q1) + l0*l1*m2*(-sin(q1)*cos(q2) + sin(q2)*cos(q1))*u2**2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))*(l0*m1*u1**2*cos(q1) + l0*m2*u1**2*cos(q1) + l1*m2*u2**2*cos(q2) + F)/(m0 + m1 + m2))/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2)))/(-l1**2*m2**2*sin(q2)**2/(m0 + m1 + m2) + l1**2*m2 - (l0*l1*m2*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) + l1*m2*(-l0*m1*sin(q1) - l0*m2*sin(q1))*sin(q2)/(m0 + m1 + m2))**2/(l0**2*m1 + l0**2*m2 - (-l0*m1*sin(q1) - l0*m2*sin(q1))**2/(m0 + m1 + m2)))]])

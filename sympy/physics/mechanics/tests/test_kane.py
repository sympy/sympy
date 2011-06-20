from sympy import symbols, Matrix, sin, cos, tan, expand
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                     RigidBody, Kane, inertia, Particle)

def test_one_dof():
    q, qd, u, ud = dynamicsymbols('q qd u ud')
    m, c, k = symbols('m c k')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, u * N.x)

    kd = {qd: u}
    FL = [(P, (-k * q - c * u) * N.x)]
    pa = Particle()
    pa.mass = m
    pa.point = P
    BL = [pa]

    KM = Kane(N)
    KM.gen_speeds([u])
    KM.kindiffeq(kd)
    KM.form_frstar(BL)
    KM.form_fr(FL)
    MM = KM.mass_matrix()
    forcing = KM.forcing()
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand(-(q * k + u * c) / m)

def test_two_dof():
    q1, q2, q1d, q2d = dynamicsymbols('q1 q2 q1d q2d')
    u1, u2, u1d, u2d = dynamicsymbols('u1 u2 u1d u2d')
    m, c1, c2, k1, k2 = symbols('m c1 c2 k1 k2')
    N = ReferenceFrame('N')
    P1 = Point('P1')
    P2 = Point('P2')
    P1.set_vel(N, u1 * N.x)
    P2.set_vel(N, (u1 + u2) * N.x)
    kd = {q1d: u1, q2d: u2}

    FL = [(P1, (-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2) * N.x), (P2, (-k2 *
        q2 - c2 * u2) * N.x)]
    pa1 = Particle()
    pa2 = Particle()
    pa1.mass = m
    pa2.mass = m
    pa1.point = P1
    pa2.point = P2
    BL = [pa1, pa2]

    KM = Kane(N)
    KM.gen_speeds([u1, u2])
    KM.kindiffeq(kd)
    KM.form_frstar(BL)
    KM.form_fr(FL)
    MM = KM.mass_matrix()
    forcing = KM.forcing()
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand((-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2)/m)
    assert expand(rhs[1]) == expand((k1 * q1 + c1 * u1 - 2 * k2 * q2 - 2 *
                                    c2 * u2) / m)

def test_pend():
    q, qd, u, ud = dynamicsymbols('q qd u ud')
    m, l, g = symbols('m l g')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, -l * u * sin(q) * N.x + l * u * cos(q) * N.y)
    kd = {qd: u}

    FL = [(P, m * g * N.x)]
    pa = Particle()
    pa.mass = m
    pa.point = P
    BL = [pa]

    KM = Kane(N)
    KM.gen_speeds([u])
    KM.kindiffeq(kd)
    KM.form_fr(FL)
    KM.form_frstar(BL)
    MM = KM.mass_matrix()
    forcing = KM.forcing()
    rhs = MM.inv() * forcing
    rhs.simplify()
    assert expand(rhs[0]) == expand(-g / l * sin(q))


def test_rolling_disc():
    # Rolling Disc Example
    q1, q2, q3, q1d, q2d, q3d = dynamicsymbols('q1 q2 q3 q1d q2d q3d')
    u1, u2, u3, u1d, u2d, u3d = dynamicsymbols('u1 u2 u3 u1d u2d u3d')
    r, m, g = symbols('r m g')

    N = ReferenceFrame('N')
    Y = N.orientnew('Y', 'Simple', q1, 3)
    L = Y.orientnew('L', 'Simple', q2, 1)
    R = L.orientnew('R', 'Simple', q3, 2)
    R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
    R.set_ang_acc(N, R.ang_vel_in(N).dt(R) + (R.ang_vel_in(N) ^ R.ang_vel_in(N)))

    C = Point('C')
    C.set_vel(N, 0)
    Dmc = C.newpoint('Dmc', r * L.z)
    Dmc.v2pt(C, N, R)
    Dmc.a2pt(C, N, R)

    I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)

    kd = dict({q1d: (u3/cos(q3)),
        q2d: (u1),
        q3d: (u2 - u3*tan(q2))})

    ForceList = [(Dmc, - m * g * Y.z)]
    BodyD = RigidBody()
    BodyD.mc = Dmc
    BodyD.inertia = (I, Dmc)
    BodyD.frame = R
    BodyD.mass = m
    BodyList = [BodyD]

    KM = Kane(N)
    KM.gen_speeds([u1, u2, u3])
    KM.kindiffeq(kd)
    fr = KM.form_fr(ForceList)
    frstar = KM.form_frstar(BodyList)
    MM = KM.mass_matrix()
    forcing = KM.forcing()
    rhs = MM.inv() * forcing
    assert rhs.expand() == Matrix([(10*u2*u3*r - 5*u3**2*r*tan(q2) +
        4*g*sin(q2))/(5*r), -2*u1*u3/3, u1*(-2*u2 + u3*tan(q2))]).expand()


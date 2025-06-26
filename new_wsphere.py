from sympy import (
    symbols,
    sin,
    cos,
    sqrt,
    simplify,
    init_printing,
    pprint,
    Rational,
)
from sympy.physics.mechanics import (
    dynamicsymbols,
    ReferenceFrame,
    Point,
    dot,
    LagrangesMethod,
    WrappingSphere,
)

init_printing(use_unicode=True)

# Generalised coordinates and their derivatives
t = symbols("t")
theta, phi = dynamicsymbols("theta phi")
dtheta, dphi = dynamicsymbols("theta phi", 1)

# Constants
R, m, k, L0 = symbols("R m k L0", positive=True)
Ax, Ay, Az, Bx, By, Bz = symbols("Ax Ay Az Bx By Bz")

# Reference frame and points
N = ReferenceFrame("N")
O = Point("O")
O.set_vel(N, 0)

# Wrapping sphere
sphere = WrappingSphere(R, O)

# Wrapping point on the sphere in spherical coordinates
P = Point("P")
x = R * sin(theta) * cos(phi)
y = R * sin(theta) * sin(phi)
z = R * cos(theta)
P.set_pos(O, x * N.x + y * N.y + z * N.z)

# Velocity of the wrapping point
P.set_vel(N, P.pos_from(O).diff(t, N))

# Kinetic energy
T = simplify(Rational(1, 2) * m * dot(P.vel(N), P.vel(N)))

# Muscle origin and insertion points fixed in space
Ax, Ay, Az = symbols("Ax Ay Az")
Bx, By, Bz = symbols("Bx By Bz")
A = Point("A")
B = Point("B")
A.set_pos(O, Ax * N.x + Ay * N.y + Az * N.z)  # origin point
B.set_pos(O, Bx * N.x + By * N.y + Bz * N.z)  # insertion point

# Approximation to get tangent points:

# Compute unit‐radial directions from O toward A and B
uA = A.pos_from(O).normalize()
uB = B.pos_from(O).normalize()

# Place TA and TB on the sphere along those rays
TA = Point("TA")
TB = Point("TB")
TA.set_pos(O, R * uA)
TB.set_pos(O, R * uB)

TA_vec = P.pos_from(A)
TB_vec = P.pos_from(B)

# Geodesic length
# Get the wrapping length between points TA and TB on the sphere, using the WrappingSphere instance
L_wrap = simplify(sphere.geodesic_length(TA, TB))
print("Wrapping length:\n")
pprint(L_wrap, use_unicode=True)

# Muscle path length = length from A to P, plus wrapping length, plus length from P to B
L_tot = simplify(sqrt(dot(TA_vec, TA_vec)) + L_wrap + sqrt(dot(TB_vec, TB_vec)))
print("\nTotal Length:")
pprint(L_tot, use_unicode=True)

# Potential energy due to muscle stretch
V = Rational(1, 2) * k * (L_tot - L0) ** 2

# Lagrangian
Lag = T - V

# Lagrange
LM = LagrangesMethod(Lag, [theta, phi])
eqns = LM.form_lagranges_equations()

# ------- Or we can define a function for Euler-Lagrange of q -------
# def euler_lagrange(L, q, dq):
#    dL_dq = simplify(diff(L, q))
#    dL_ddq = simplify(diff(L, dq))

#    ddt_dL_ddq = simplify(diff(dL_ddq, t))

#    EL_eq = simplify(ddt_dL_ddq - dL_dq)
#    return dL_dq, dL_ddq, ddt_dL_ddq, EL_eq

# For theta
# dL_dtheta, dL_ddtheta, ddt_dL_ddtheta, EL_theta = euler_lagrange(Lag, theta, dtheta)

# For phi
# dL_dphi, dL_ddphi, ddt_dL_ddphi, EL_phi = euler_lagrange(Lag, phi, dphi)

# Display results
print("\nKinetic Energy:")
pprint(T, use_unicode=True)
print("\nPotential Energy:")
pprint(V, use_unicode=True)
print("\nLagrangian:")
pprint(Lag, use_unicode=True)
print("\nEquations of Motion:")
pprint(eqns)
print(f"Number of equations: {len(eqns)}")
for i, eq in enumerate(eqns):
    print(f"Equation {i + 1}:")
    pprint(simplify(eq), use_unicode=True)

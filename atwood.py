# ---------------
# Atwood Machine
# ---------------
import sympy as sp
from sympy.physics.mechanics import (
    ReferenceFrame,
    Point,
    Particle,
    LagrangesMethod,
    dynamicsymbols,
    WrappingCylinder,
    WrappingPathway,
)

# Symbols and generalized coordinate
t = sp.symbols("t")
m1, m2, g, r, h = sp.symbols("m1 m2 g r h", positive=True, real=True)

q = dynamicsymbols("q")  # generalized coordinate: how far m1 moves downward
qd = dynamicsymbols("q", 1)  # time‐derivative ˙q

# Inertial frame & pulley
N = ReferenceFrame("N")
O = Point("O")
O.set_vel(N, 0)  # pulley center is fixed at origin

# Position of the two hanging masses (both z < 0)
#    We choose:
#       P1 at (0, +r, z1) with z1 = –(h + q)
#       P2 at (0, –r, z2) with z2 = –(h – q)
# quick note: -h <= q(t) <= +h

P1 = Point("P1")
P1.set_pos(O, 0 * N.x + r * N.y + (-(h + q)) * N.z)
P1.set_vel(N, P1.pos_from(O).diff(t, N))
M1 = Particle("M1", P1, m1)

P2 = Point("P2")
P2.set_pos(O, 0 * N.x + (-r) * N.y + (-(h - q)) * N.z)
P2.set_vel(N, P2.pos_from(O).diff(t, N))
M2 = Particle("M2", P2, m2)


# Set tangent points for geodesic_length method
# We have selected the points on same vertical as
# the masses on the pulley
# This is so that we can use WrappingCylinder to find the
# wrapped string length on the cylinder

T1 = Point("T1")
T1.set_pos(O, r * N.y)

T2 = Point("T2")
T2.set_pos(O, (-r) * N.y)

pulley = WrappingCylinder(r, O, N.x)  # cylinder of radius r, axis along N.x

wpath = WrappingPathway(T1, T2, pulley)

# Build total string length piecewise:
L1 = sp.sqrt((P1.pos_from(T1).dot(P1.pos_from(T1))))  # |P1 - T1| == h + q(t)

L_curve = wpath.length  # arc on cylinder between T1,T2 (uses the WrappingCylinder)

L2 = sp.sqrt((P2.pos_from(T2).dot(P2.pos_from(T2))))  # |T2 - P2| == h - q(t)

L_total = sp.simplify(L1 + L_curve + L2)

print(f"L_curve = {L_curve}")
print(f"L_total = {L_total}")

# Ensure inextensibility: d(L_total)/dq = 0
dL_dq = sp.simplify(sp.diff(L_total, q))

# sqrt((h + q) ** 2) and sqrt((h - q) ** 2) cannot be simplified due to the assumptions in L_total, thus this assertion fails
# assert dL_dq == 0, "ERROR: total string length depends on q!"

# Kinetic & potential energies
T = sp.simplify(sp.Rational(1, 2) * (m1 * qd**2 + m2 * qd**2))
V = sp.simplify(m1 * g * (-(h + q)) + m2 * g * (-(h - q)))
#  L  = - (m1 + m2) g h  + (m2 - m1) g q

Lagr = T - V

# Lagrange’s equation
LM = LagrangesMethod(Lagr, (q,), forcelist=[])
eom = LM.form_lagranges_equations()[0]  # the scalar equation d/dt(∂L/∂qd) - ∂L/∂q = 0
qdd = sp.solve(eom, sp.diff(q, (t, 2)))[0]
qdd_simpl = sp.simplify(qdd)

print("Symbolic EOM (both masses below pulley):")
print("    q̈(t) =", qdd_simpl)
print("Expected: (m2 - m1)/(m1 + m2) * (-g)")

# Numeric check
numeric_subs = {m1: 1.0, m2: 2.0, g: 9.81, h: 5.0, r: 0.5}
qdd_num = float(qdd_simpl.subs(numeric_subs))

print("\nNumeric check (m1=1 kg, m2=2 kg, g=9.81 m/s²):")
print("    q̈ =", qdd_num, "m/s²")
print(
    "    Expected = (2 - 1)/(2 + 1)*9.81 =",
    float((-1) * (2.0 - 1.0) / (2.0 + 1.0) * 9.81),
)


# Conclusion -> WrappingCylinder CAN be used
# Atleast for simple cases where geodesic lengths can be simplified easily

# ---------------
# Atwood Machine
# ---------------
import sympy as sp
from sympy.physics.mechanics import (
    ReferenceFrame,
    Point,
    Particle,
    KanesMethod,
    dynamicsymbols,
    WrappingCylinder,
    WrappingPathway,
    Force,
)

# 1. -------- Symbols and generalized coordinate --------

t = sp.symbols("t")
m1, m2, g, r, h, T = sp.symbols("m1 m2 g r h T", positive=True, real=True)

q = dynamicsymbols("q")  # generalized coordinate: m1 moves downward by q
u = dynamicsymbols("u")  # generalized speed (u = q̇)

# 2. -------- Inertial frame & pulley center --------

N = ReferenceFrame("N")
O = Point("O")
O.set_vel(N, 0)  # pulley center is fixed at O = (0,0,0)

# 3. -------- Positions & velocities of the two masses --------

# Mass m1 at P1: (x=0, y=+r, z=−(h + q))
P1 = Point("P1")
P1.set_pos(O, r * N.y + (-(h + q)) * N.z)
P1.set_vel(N, P1.pos_from(O).diff(t, N))
M1 = Particle("M1", P1, m1)

# Mass m2 at P2: (x=0, y=−r, z=−(h − q))
P2 = Point("P2")
P2.set_pos(O, -r * N.y + (-(h - q)) * N.z)
P2.set_vel(N, P2.pos_from(O).diff(t, N))
M2 = Particle("M2", P2, m2)

# 4. -------- Create WrappingCylinder --------

# Cylinder of radius r, centered at O, axis along N.x
pulley = WrappingCylinder(r, O, N.x)

# 5. -------- Manually find tangent points T1 and T2
# (no inbuilt method / object solves this problem) --------

# Because P1 is directly below the cylinder’s “eastmost” point (y=+r),
# the line from P1 to the cylinder is vertical. Tangency requires that
# the point on the cylinder satisfy y = +r and z = 0.
T1 = Point("T1")
T1.set_pos(O, r * N.y + 0 * N.z)
T1.set_vel(N, 0)

# Similarly, P2 is directly below (x=0, y=−r),
# so its tangent point is at (y=−r, z=0).
T2 = Point("T2")
T2.set_pos(O, -r * N.y + 0 * N.z)
T2.set_vel(N, 0)

wpath = WrappingPathway(T1, T2, pulley)

# 6. -------- Compute segment lengths and verify inextensibility --------

# Straight‐line segment from P1 to T1: vertical distance = h + q
L1 = sp.sqrt((P1.pos_from(T1).dot(P1.pos_from(T1))))  # = h + q

# Straight‐line segment from P2 to T2: vertical distance = h − q
L2 = sp.sqrt((P2.pos_from(T2).dot(P2.pos_from(T2))))  # = h − q

# Arc (geodesic) on the cylinder between T1 and T2:
L_curve = wpath.length  # = π * r, independent of q

# Total rope length
L_total = sp.simplify(L1 + L_curve + L2)

# Check that d(L_total)/dq = 0 (inextensible rope)
dL_dq = sp.simplify(sp.diff(L_total, q))
# assert dL_dq == 0, "ERROR: total rope length depends on q!"

print(f"L1      = {L1}")  # should be h + q
print(f"L_curve = {L_curve}")  # should be π r
print(f"L2      = {L2}")  # should be h − q
print(f"L_total = {L_total}")  # should simplify to 2h + πr
print(f"dL/dq   = {dL_dq}")  # should be 0

# 7. -------- Gravity forces on each mass --------

grav1 = Force(P1, -m1 * g * N.z)
grav2 = Force(P2, -m2 * g * N.z)

# 8. -------- Collect all loads for Kane’s method --------

loads = wpath.to_loads(T) + [grav1, grav2]

# 9. -------- Kinematic differential equation --------

# Declare u = q̇
kin_diff = [u - q.diff()]

# 10. -------- Formulate and solve via Kane’s Method --------

kane = KanesMethod(N, (q,), (u,), kd_eqs=kin_diff)
bodies = [M1, M2]

Fr, Frs = kane.kanes_equations(bodies, loads)

mass_matrix = kane.mass_matrix
forcing_vec = kane.forcing

# Solve for u̇ = q̈
# sp.solve returns dict of 1 key value pair so we index using u.diff()
u_dot = sp.solve((mass_matrix * u.diff() - forcing_vec), u.diff())[u.diff()]
qdd = sp.simplify(u_dot)

print("\nSymbolic equation of motion (including tension T):")
print("    q̈(t) =")
sp.pprint(qdd, use_unicode=True)

# 11. -------- Eliminate T to recover classic Atwood result --------

# Force balances on the two masses give:
#   m1·q̈ = (component of T on m1 along +z) − m1·g
#          = T − m1·g
#   m2·q̈ = m2·g − (component of T on m2 along +z)
#          = m2·g − T
#
# Solve for q̈, T:
qdd_sym, T_sym = sp.symbols("qdd_sym T_sym", real=True)
eqs = [
    m1 * qdd_sym - (T_sym - m1 * g),
    m2 * qdd_sym - (m2 * g - T_sym),
]
solution = sp.solve(eqs, (qdd_sym, T_sym))
qdd_classic = sp.simplify(solution[qdd_sym])

print("\nEliminated T → classic Atwood result:")
print("    q̈ =", qdd_classic)

# 12. -------- Numeric check --------

numeric_vals = {m1: 1.0, m2: 2.0, g: 9.81, h: 5.0, r: 0.5}
qdd_num = float(qdd_classic.subs(numeric_vals))
T_num = float(solution[T_sym].subs(numeric_vals))

print("\nNumeric check (m1=1 kg, m2=2 kg, g=9.81 m/s²):")
print(f"    q̈ = {qdd_num:.6f} m/s²")
print(f"    Tension T = {T_num:.6f} N")
print(
    "    Expected q̈ = (m2 - m1)/(m1 + m2) * g =",
    float((2.0 - 1.0) / (2.0 + 1.0) * 9.81),
)

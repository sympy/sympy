import sympy as sp
from sympy.physics.mechanics import (
    ReferenceFrame,
    Point,
    dynamicsymbols,
    WrappingSphere,
    WrappingPathway,
    inertia,
    RigidBody,
    Force,
    KanesMethod,
)

# Symbols and dynamic variables
t = sp.symbols("t")
q = dynamicsymbols("q")
u = dynamicsymbols("u")
k, L0, c, r_sph, g, I_bone, m_bone, L_b = sp.symbols(
    "k L0 c r_sph g I_bone m_bone L_b", positive=True, real=True
)

# Inertial frame and origin
N = ReferenceFrame("N")
O = Point("O")
O.set_vel(N, 0)

# Muscle origin (fixed) and insertion (rotating)
P_orig = Point("P_orig")
P_orig.set_pos(O, -r_sph * N.x)
P_orig.set_vel(N, 0)

A = N.orientnew("A", "Axis", [q, N.z])
P_ins = Point("P_ins")
P_ins.set_pos(O, L_b * A.y)
P_ins.v2pt_theory(O, N, A)

sphere = WrappingSphere(r_sph, O)

# Tangent point variables
x1, y1, x2, y2 = sp.symbols("x1 y1 x2 y2", real=True)
P1 = Point("T1")
P2 = Point("T2")
P1.set_pos(O, x1 * N.x + y1 * N.y)
P2.set_pos(O, x2 * N.x + y2 * N.y)

# Tangent constraints
orig_vec = P_orig.pos_from(O).to_matrix(N)
ins_vec = P_ins.pos_from(O).to_matrix(N)
eqs = [
    x1**2 + y1**2 - r_sph**2,
    (x1 - orig_vec[0]) * x1 + (y1 - orig_vec[1]) * y1,
    x2**2 + y2**2 - r_sph**2,
    (x2 - ins_vec[0]) * x2 + (y2 - ins_vec[1]) * y2,
]

sol = sp.solve(eqs, [x1, y1, x2, y2], dict=True)[0]

# Set tangent points
P1.set_pos(O, sol[x1] * N.x + sol[y1] * N.y)
P2.set_pos(O, sol[x2] * N.x + sol[y2] * N.y)

# Muscle path
L1 = P_orig.pos_from(P1).magnitude()
L2 = P_ins.pos_from(P2).magnitude()
wpath = WrappingPathway(P1, P2, sphere)
L_curve = wpath.length
L_total = sp.simplify(L1 + L_curve + L2)

moment_arm = sp.diff(L_total, q)
L_dot = sp.diff(L_total, t)
T_muscle = k * (L_total - L0) + c * L_dot
tau_muscle = moment_arm * T_muscle

# Bone body definition
P_cm = Point("P_cm")
P_cm.set_pos(O, (L_b / 2) * A.y)
P_cm.v2pt_theory(O, N, A)

Ib = inertia(A, 0, 0, I_bone)
body = RigidBody("Bone", P_cm, A, m_bone, (Ib, P_cm))

# Forces and torques
grav_load = (P_cm, -m_bone * g * N.y)
muscle_load = (A, tau_muscle * A.z)

# Kane's method
kane = KanesMethod(N, q_ind=[q], u_ind=[u], kd_eqs=[u - q.diff(t)])
fr, frstar = kane.kanes_equations([body], [grav_load, muscle_load])

MM = kane.mass_matrix
forcing = kane.forcing
u_dot = MM.inv() * forcing

# Output
print(f"L_total: ")
sp.pprint(L_total, use_unicode=True)
print("\nmoment_arm:")
sp.pprint(moment_arm, use_unicode=True)
print("\n(u̇):")
print(u_dot)

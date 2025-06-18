import sympy as sp
from sympy.physics.mechanics import (
    ReferenceFrame,
    Point,
    dynamicsymbols,
    KanesMethod,
    inertia,
    RigidBody,
    WrappingSphere,
    WrappingPathway,
)
import time

# Symbols and dynamic/constant variables
t = sp.symbols("t")
q1, q2 = dynamicsymbols("q1 q2")  # 2 DOF coordinates
u1, u2 = dynamicsymbols("u1 u2")  # generalized speeds
r_sph, d = sp.symbols("r_sph d", positive=True, real=True)  # sphere radius + offset
L1, L2 = sp.symbols("L1 L2", positive=True, real=True)
m1, m2 = sp.symbols("m1 m2", positive=True, real=True)
I1xx, I1yy, I1zz = sp.symbols("I1xx I1yy I1zz", positive=True, real=True)
I2xx, I2yy, I2zz = sp.symbols("I2xx I2yy I2zz", positive=True, real=True)
g = sp.symbols("g", positive=True, real=True)
T = sp.symbols("T", real=True)  # muscle tension

# Kinematic differential equations
kin_diff = [u1 - q1.diff(t), u2 - q2.diff(t)]

# Reference frames and origin
N = ReferenceFrame("N")
A1 = N.orientnew("A1", "Axis", [q1, N.z])
A2 = A1.orientnew("A2", "Axis", [q2, A1.x])
print("A1 ReferenceFrame: ")
print(A1.dcm(N))
print("A2 ReferenceFrame: ")
print(A2.dcm(N))

O = Point("O")
O.set_vel(N, 0)

# Link1 COM and dynamics
P1_cm = O.locatenew("P1_cm", (L1 / 2) * A1.y)
P1_cm.v2pt_theory(O, N, A1)
I1 = inertia(A1, I1xx, I1yy, I1zz)
Body1 = RigidBody("Link1", P1_cm, A1, m1, (I1, P1_cm))

# Link2 COM and dynamics
P2_cm = O.locatenew("P2_cm", L1 * A1.y + (L2 / 2) * A2.y)
P2_cm.v2pt_theory(O, N, A2)
I2 = inertia(A2, I2xx, I2yy, I2zz)
Body2 = RigidBody("Link2", P2_cm, A2, m2, (I2, P2_cm))

# Muscle origin & insertion points
o = O.locatenew("P_orig", (r_sph + d) * N.x)  # o -> P_original
o.set_vel(N, 0)
i = O.locatenew("P_ins", L1 * A1.y + L2 * A2.y)  # i -> P_insertion
i.v2pt_theory(O, N, A2)

# Wrapping sphere at origin
C = Point("C")
C.set_pos(O, 0 * N.x + 0 * N.y + 0 * N.z)
C.set_vel(N, 0)
sphere = WrappingSphere(r_sph, C)  # Point C is just Point O (origin)

# Symbols for tangent points
x1, y1, z1, x2, y2, z2 = sp.symbols("x1 y1 z1 x2 y2 z2", real=True)

# Compute plane normal through C, origin, insertion
v_orig = o.pos_from(C).to_matrix(N)
v_ins = i.pos_from(C).to_matrix(N)

# 1. plane normal
n_vec = v_orig.cross(v_ins)
n_unit = n_vec / sp.sqrt(n_vec.dot(n_vec))

# 2. build in-plane orthonormal basis {e1,e2}
# project v_orig onto plane
v_proj = v_orig - (n_unit.dot(v_orig)) * n_unit
# e1 aligned with projected v_orig
e1 = v_proj / sp.sqrt(v_proj.dot(v_proj))
# e2 completes right-handed basis
e2 = n_unit.cross(e1)

# 3. get plane-coordinates of external points
x_o = e1.dot(v_orig)
y_o = e2.dot(v_orig)

x_i = e1.dot(v_ins)
y_i = e2.dot(v_ins)


# function to compute tangent points in plane for coords (xE,yE)
def plane_tangents(xE, yE):
    u = sp.symbols("u", real=True)
    sol_pts = []
    # handle case where yE == 0 (avoid division by zero)
    if yE.equals(0):
        # orthogonality: u*xE = r^2 -> u = r^2/xE
        u_val = r_sph**2 / xE
        # solve for v: u^2 + v^2 = r^2
        v = sp.symbols("v", real=True)
        eq_v = sp.simplify(u_val**2 + v**2 - r_sph**2)
        v_sols = sp.solve(eq_v, v)
        for idx, v_val in enumerate(v_sols, start=1):
            vec = (
                (u_val * e1[0] + v_val * e2[0]) * N.x
                + (u_val * e1[1] + v_val * e2[1]) * N.y
                + (u_val * e1[2] + v_val * e2[2]) * N.z
            )
            P = C.locatenew(f"P{idx}", vec)
            sol_pts.append(P)
    else:
        # general case
        # orthogonality: u*xE + v*yE = r^2  =>  v = (r^2 - u*xE)/yE
        v_expr = (r_sph**2 - u * xE) / yE
        # circle: u^2 + v^2 = r^2
        eq = sp.simplify(u**2 + v_expr**2 - r_sph**2)
        # solve quadratic in u
        u_sols = sp.solve(eq, u)
        for idx, u_val in enumerate(u_sols, start=1):
            v_val = sp.simplify(v_expr.subs(u, u_val))
            vec = (
                (u_val * e1[0] + v_val * e2[0]) * N.x
                + (u_val * e1[1] + v_val * e2[1]) * N.y
                + (u_val * e1[2] + v_val * e2[2]) * N.z
            )
            P = O.locatenew(f"P{idx}", vec)
            sol_pts.append(P)
    return sol_pts


# compute tangent points
start = time.time()
P_tangents_orig = plane_tangents(x_o, y_o)
P_tangents_ins = plane_tangents(x_i, y_i)
print(f"Computed in {time.time() - start:.4f}s")

# pick the correct pair (e.g., the two real ones)
print(P_tangents_orig)
print("-" * 90)
print(P_tangents_ins)
P1, P2 = P_tangents_orig[0], P_tangents_ins[1]

# build wrapping pathway and length
# Muscle path: straight + wrap + straight
wpath = WrappingPathway(P1, P2, sphere)
WPath = sp.Function("WPath")
wpath_length = WPath(
    L1, L2, q1, q2, r_sph, d
)  # wpath.length # sp.symbols("L_curve", positive=True, real=True)
print(f"Geodesic length on sphere = {wpath_length}")
L_total = o.pos_from(P1).magnitude() + wpath_length + i.pos_from(P2).magnitude()

# Moment arms about each DOF
# is the sign correct?
m_arm_q1 = sp.diff(L_total, q1)
m_arm_q2 = sp.diff(L_total, q2)

# Applied loads: gravity + muscle torques
gravity_loads = [(P1_cm, -m1 * g * N.y), (P2_cm, -m2 * g * N.y)]
muscle_loads = [(A1, m_arm_q1 * T * A1.z), (A2, m_arm_q2 * T * A2.x)]

loads = gravity_loads + muscle_loads

# Kane's method
kane = KanesMethod(N, q_ind=[q1, q2], u_ind=[u1, u2], kd_eqs=kin_diff)
fr, frstar = kane.kanes_equations(bodies=[Body1, Body2], loads=loads)

# Generalized speeds derivatives
[u1, u2, u1_dot, u2_dot] = kane.rhs()
print(f"u1_dot = {u1_dot}")
print("-" * 80)
print(f"u2_dot = {u2_dot}")

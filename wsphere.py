import sympy as sp
from sympy import Matrix
from sympy.physics.mechanics import (
    ReferenceFrame, Point, dynamicsymbols,
    KanesMethod, inertia, RigidBody,
    WrappingSphere, WrappingPathway
)
import time

# Symbols and dynamic/constant variables
t = sp.symbols('t')
q1, q2 = dynamicsymbols('q1 q2')  # 2 DOF coordinates
u1, u2 = dynamicsymbols('u1 u2')  # generalized speeds
r_sph, d = sp.symbols('r_sph d', positive=True, real=True)  # sphere radius + offset
L1, L2 = sp.symbols('L1 L2', positive=True, real=True)
m1, m2 = sp.symbols('m1 m2', positive=True, real=True)
I1xx, I1yy, I1zz = sp.symbols('I1xx I1yy I1zz', positive=True, real=True)
I2xx, I2yy, I2zz = sp.symbols('I2xx I2yy I2zz', positive=True, real=True)
g = sp.symbols('g', positive=True, real=True)
T = sp.symbols('T', real=True)  # muscle tension

# Kinematic differential equations
kin_diff = [u1 - q1.diff(t), u2 - q2.diff(t)]

# Reference frames and origin
N = ReferenceFrame('N')
A1 = N.orientnew('A1', 'Axis', [q1, N.z])
A2 = A1.orientnew('A2', 'Axis', [q2, A1.x])
print(f"A1 ReferenceFrame: ")
print(A1.dcm(N))
print(f"A2 ReferenceFrame: ")
print(A2.dcm(N))

O = Point('O')
O.set_vel(N, 0)

# Link1 COM and dynamics
P1_cm = O.locatenew('P1_cm', (L1/2)*A1.y)
P1_cm.v2pt_theory(O, N, A1)
I1 = inertia(A1, I1xx, I1yy, I1zz)
Body1 = RigidBody('Link1', P1_cm, A1, m1, (I1, P1_cm))

# Link2 COM and dynamics
P2_cm = O.locatenew('P2_cm', L1*A1.y + (L2/2)*A2.y)
P2_cm.v2pt_theory(O, N, A2)
I2 = inertia(A2, I2xx, I2yy, I2zz)
Body2 = RigidBody('Link2', P2_cm, A2, m2, (I2, P2_cm))

# Muscle origin & insertion points
o = O.locatenew('P_orig', (r_sph + d)*N.x) # o -> P_original
o.set_vel(N, 0)
i = O.locatenew('P_ins', L1*A1.y + L2*A2.y) # i -> P_insertion
i.v2pt_theory(O, N, A2)

# Wrapping sphere at origin
C = Point('C'); C.set_pos(O, 0*N.x + 0*N.y + 0*N.z)
C.set_vel(N, 0)
sphere = WrappingSphere(r_sph, C) # Point C is just Point O (origin)

# Symbols for tangent points
x1, y1, z1, x2, y2, z2 = sp.symbols('x1 y1 z1 x2 y2 z2', real=True)

# Compute plane normal through C, origin, insertion
v_orig = o.pos_from(C).to_matrix(N)
v_ins  = i.pos_from(C).to_matrix(N)
n = v_orig.cross(v_ins)

# Build 6 equations for two tangent points
eqs = []
for xi, yi, zi, vec in [(x1, y1, z1, v_orig), (x2, y2, z2, v_ins)]:
    # point on sphere surface
    eqs.append(xi**2 + yi**2 + zi**2 - r_sph**2)
    # radius perpendicular to rope
    eqs.append((Matrix([xi, yi, zi]) - vec).dot(Matrix([xi, yi, zi])))
    # coplanarity: dot with plane normal = 0
    eqs.append(Matrix([xi, yi, zi]).dot(n))

# Solve for tangent coordinates
start = time.time()
sol = sp.solve(eqs, [x1, y1, z1, x2, y2, z2], dict=True)[0] # very CPU intensive and time consuming
print(f"sol = ")
print(sol)
print(f"Time consumed by call to sympy.solve() = {time.time() - start}")

# Define tangent points in N
P1 = O.locatenew('P1', sol[x1]*N.x + sol[y1]*N.y + sol[z1]*N.z)
P2 = O.locatenew('P2', sol[x2]*N.x + sol[y2]*N.y + sol[z2]*N.z)

# Muscle path: straight + wrap + straight
wpath = WrappingPathway(P1, P2, sphere)
L_total = o.pos_from(P1).magnitude() + wpath.length + i.pos_from(P2).magnitude()

# Moment arms about each DOF
# is the sign correct?
m_arm_q1 = sp.diff(L_total, q1)
m_arm_q2 = sp.diff(L_total, q2)

# Applied loads: gravity + muscle torques
gravity_loads = [(P1_cm, -m1*g*N.y), (P2_cm, -m2*g*N.y)]
muscle_loads = [(A1, m_arm_q1*T*A1.z), (A2, m_arm_q2*T*A2.x)]

loads = gravity_loads + muscle_loads

# Kane's method
kane = KanesMethod(N, q_ind=[q1, q2], u_ind=[u1, u2], kd_eqs=kin_diff)
fr, frstar = kane.kanes_equations(
    bodies=[Body1, Body2], loads=loads
)

# Generalized speeds derivatives
u_dot = kane.rhs()
print('Generalized speed derivatives:', u_dot)

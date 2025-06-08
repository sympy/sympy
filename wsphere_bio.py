import sympy as sm
import sympy.physics.mechanics as me
import sympy.physics.biomechanics as bm

# Define generalized coordinates and speeds
q1, q2, q3, q4 = me.dynamicsymbols("q1, q2, q3, q4", real=True)
u1, u2, u3, u4 = me.dynamicsymbols("u1, u2, u3, u4", real=True)

# Define constants and parameters
dx, dy, dz = sm.symbols("dx, dy, dz", real=True, nonnegative=True)
lA, lC, lD = sm.symbols("lA, lC, lD", real=True, positive=True)
mA, mC, mD = sm.symbols("mA, mC, mD", real=True, positive=True)
g, k, c, r = sm.symbols("g, k, c, r", real=True, positive=True)

# Reference frames and points
N, A, B, C, D = sm.symbols("N, A, B, C, D", cls=me.ReferenceFrame)
O, P1, P2, P3, P4 = sm.symbols("O, P1, P2, P3, P4", cls=me.Point)
Ao, Co, Cm, Dm, Do = sm.symbols("Ao, Co, Cm, Dm, Do", cls=me.Point)

# Orient frames
A.orient_axis(N, q1, N.z)
B.orient_axis(N, q2, N.y)
C.orient_axis(B, q3, B.z)
D.orient_axis(C, q4, C.y)

# Angular velocities
A.set_ang_vel(N, u1 * N.z)
B.set_ang_vel(N, u2 * N.y)
C.set_ang_vel(B, u3 * B.z)
D.set_ang_vel(C, u4 * C.y)

# Positions of points
Ao.set_pos(O, dx * N.x)
P1.set_pos(Ao, lA * A.y)
P2.set_pos(O, dy * N.y + dz * N.z)
Co.set_pos(P2, (lC / 2) * C.z)
Cm.set_pos(P2, (lC / 3) * C.z)
P3.set_pos(P2, lC * C.z)
Dm.set_pos(P3, (lD / 3) * D.z)
Do.set_pos(P3, (lD / 2) * D.z)
P4.set_pos(P3, lD * D.z)

# Velocities of points
O.set_vel(N, 0)
Ao.set_vel(N, 0)
P1.v2pt_theory(Ao, N, A)
P2.set_vel(N, 0)
Co.v2pt_theory(P2, N, C)
Cm.v2pt_theory(P2, N, C)
P3.v2pt_theory(P2, N, C)
Dm.v2pt_theory(P3, N, D)
Do.v2pt_theory(P3, N, D)
P4.v2pt_theory(P3, N, D)

# Holonomic constraint: P4 must coincide with P1
holonomic = (P4.pos_from(O) - P1.pos_from(O)).to_matrix(N)

# Inertia dyadics
IA = me.Inertia(me.inertia(A, mA / 12 * lA**2, mA / 2 * lA**2, mA / 12 * lA**2), Ao)
IC = me.Inertia(me.inertia(C, mC / 12 * lC**2, mC / 12 * lC**2, mC / 2 * lC**2), Co)
ID = me.Inertia(me.inertia(D, mD / 12 * lD**2, mD / 12 * lD**2, mD / 2 * lD**2), Do)

# Rigid bodies
lever = me.RigidBody("lever", masscenter=Ao, frame=A, mass=mA, inertia=IA)
u_arm = me.RigidBody("upper_arm", masscenter=Co, frame=C, mass=mC, inertia=IC)
l_arm = me.RigidBody("lower_arm", masscenter=Do, frame=D, mass=mD, inertia=ID)

# Gravity forces
gravC = me.Force(u_arm, mC * g * N.z)
gravD = me.Force(l_arm, mD * g * N.z)

# Resistance torque at base (lever)
lever_resistance = me.Torque(A, (-k * q1 - c * u1) * N.z)

# Biceps: linear (unwrapped) pathway between Cm and Dm
biceps_pathway = me.LinearPathway(Cm, Dm)
biceps_activation = bm.FirstOrderActivationDeGroote2016.with_defaults("biceps")
biceps = bm.MusculotendonDeGroote2016.with_defaults(
    "biceps", biceps_pathway, biceps_activation
)

# Create WrappingSphere at the elbow pin joint
sphere = me.WrappingSphere(r, P3)

tf = me.dynamicsymbols._t

# A custom pathway that handles all segments
class SphericalExtensorPathway(me.PathwayBase):
    def __init__(self, origin, insertion, sphere, joint_axis):
        super().__init__(origin, insertion)
        self.origin = origin
        self.insertion = insertion
        self.sphere = sphere
        self.center = sphere.point
        self.joint_axis = joint_axis

        # Compute fixed distances and tangent angles
        self.d1 = self.center.pos_from(origin).magnitude()
        self.d2 = self.center.pos_from(insertion).magnitude()
        self.phi1 = sm.asin(sphere.radius / self.d1)
        self.phi2 = sm.asin(sphere.radius / self.d2)

        # Define tangent points on sphere surface
        self.T1 = me.Point("T1")
        self.T2 = me.Point("T2")
        ea = (origin.pos_from(self.center)).normalize()
        ei = (insertion.pos_from(self.center)).normalize()
        # Finding tangent points and vectors manually again
        t1_vec = (
            -sphere.radius * sm.cos(self.phi1) * (ea.cross(joint_axis).normalize())
            + sphere.radius * sm.sin(self.phi1) * ea
        )
        t2_vec = (
            +sphere.radius * sm.cos(self.phi2) * (ei.cross(joint_axis).normalize())
            + sphere.radius * sm.sin(self.phi2) * ei
        )
        self.T1.set_pos(self.center, t1_vec)
        self.T2.set_pos(self.center, t2_vec)

        # Underlying wrapped (geodesic) segment
        self.wrap = me.WrappingPathway(self.T1, self.T2, sphere)

    @property
    def length(self):
        # straight-line + geodesic + straight-line
        L1 = self.d1 * sm.cos(self.phi1)
        L2 = self.wrap.length
        L3 = self.d2 * sm.cos(self.phi2)
        return L1 + L2 + L3

    @property
    def extension_velocity(self):
        # time derivative of sphere arc only
        return sphere.radius * (self.phi1 + self.phi2 + self.wrap.length.diff(tf))

    def to_loads(self, Fm):
        # Forces at T1, center, T2
        f1, f_center, f2 = self.wrap.to_loads(Fm)
        return [
            me.Force(self.origin, f1.vector),
            f_center,
            me.Force(self.insertion, f2.vector),
        ]


# Instantiate and plug into muscle
triceps_activation = bm.FirstOrderActivationDeGroote2016.with_defaults("triceps")
triceps_path = SphericalExtensorPathway(Cm, Dm, sphere, N.y)
triceps = bm.MusculotendonDeGroote2016.with_defaults(
    "triceps", triceps_path, triceps_activation
)

# Build combined loads
loads = biceps.to_loads() + triceps.to_loads() + [lever_resistance, gravC, gravD]

# Set up Kane's method
kane = me.KanesMethod(
    N,
    (q1,),
    (u1,),
    kd_eqs=(
        u1 - q1.diff(),
        u2 - q2.diff(),
        u3 - q3.diff(),
        u4 - q4.diff(),
    ),
    q_dependent=(q2, q3, q4),
    configuration_constraints=holonomic,
    velocity_constraints=holonomic.diff(me.dynamicsymbols._t),
    u_dependent=(u2, u3, u4),
)

# Form equations of motion
Fr, Frs = kane.kanes_equations((lever, u_arm, l_arm), loads)

# Extract symbols for lambdification
q, u = kane.q, kane.u
a = biceps.x.col_join(triceps.x)
x = q.col_join(u).col_join(a)
e = biceps.r.col_join(triceps.r)

# Parameter vector p (numeric values later)
p = sm.Matrix(
    [
        dx,
        dy,
        dz,
        lA,
        lC,
        lD,
        mA,
        mC,
        mD,
        g,
        k,
        c,
        r,
        biceps.F_M_max,
        biceps.l_M_opt,
        biceps.l_T_slack,
        triceps.F_M_max,
        triceps.l_M_opt,
        triceps.l_T_slack,
    ]
)

# Create lambdified functions for mass matrix, forcing, and activation ODEs
eval_diffeq = sm.lambdify(
    (q, u, a, e, p),
    (kane.mass_matrix, kane.forcing, biceps.rhs().col_join(triceps.rhs())),
    cse=True,
)
eval_holonomic = sm.lambdify((q, p), holonomic, cse=True)

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

# Numeric values for parameters
p_vals = np.array(
    [
        0.31,  # dx [m]
        0.15,  # dy [m]
        -0.31,  # dz [m]
        0.2,  # lA [m]
        0.3,  # lC [m]
        0.3,  # lD [m]
        1.0,  # mA [kg]
        2.3,  # mC [kg]
        1.7,  # mD [kg]
        9.81,  # g [m/s^2]
        5.0,  # k [Nm/rad]
        0.5,  # c [N·s/m]
        0.03,  # r [m]
        500.0,  # biceps F_M_max
        0.6 * 0.3,  # biceps l_M_opt
        0.55 * 0.3,  # biceps l_T_slack
        500.0,  # triceps F_M_max
        0.6 * 0.3,  # triceps l_M_opt
        0.65 * 0.3,  # triceps l_T_slack
    ]
)

# Initial guess for dependent q2, q3, q4 given q1 = 5°
q_vals = np.array(
    [
        np.deg2rad(5.0),  # q1 [rad]
        np.deg2rad(-10.0),  # q2 [rad]
        np.deg2rad(0.0),  # q3 [rad]
        np.deg2rad(75.0),  # q4 [rad]
    ]
)


def eval_holo_fsolve(x):
    q1_val = q_vals[0]
    q2_val, q3_val, q4_val = x
    return eval_holonomic((q1_val, q2_val, q3_val, q4_val), p_vals).squeeze()


# Solve holonomic constraint for q2, q3, q4
q_vals[1:] = fsolve(eval_holo_fsolve, q_vals[1:])
print("Dependent angles (deg):", np.rad2deg(q_vals))

# Initial speeds and activations
u_vals = np.zeros(4)  # [u1, u2, u3, u4]
a_vals = np.zeros(2)  # [a_biceps, a_triceps]
e_vals = np.zeros(2)  # [e_biceps, e_triceps]

print(eval_diffeq(q_vals, u_vals, a_vals, e_vals, p_vals))

# Define the time derivative function
def eval_r_constant(t):
    # No excitation for either muscle
    return np.array([0.0, 0.0])


def eval_rhs(t, x, r, p):
    q_val = x[0:4]
    u_val = x[4:8]
    a_val = x[8:10]
    e_val = r(t)
    qd = u_val
    m_mat, f_vec, ad = eval_diffeq(q_val, u_val, a_val, e_val, p)
    ud = np.linalg.solve(m_mat, f_vec).squeeze()
    return np.hstack((qd, ud, ad.squeeze()))


# Time span for integration
t0, tf = 0.0, 3.0
ts = np.linspace(t0, tf, 301)
x0 = np.hstack((q_vals, u_vals, a_vals))

# Integrate with zero excitation
sol_zero = solve_ivp(
    lambda t, x: eval_rhs(t, x, eval_r_constant, p_vals), (t0, tf), x0, t_eval=ts
)

# Plotting routine for state trajectories
import matplotlib.pyplot as plt

def plot_traj(t, x, syms):
    fig, axes = plt.subplots(5, 2, sharex=True, figsize=(8, 10))
    for ax, traj, sym in zip(axes.T.flatten(), x.T, syms):
        if not sym.name.startswith("a"):
            traj = np.rad2deg(traj)
        ax.plot(t, traj)
        ax.set_ylabel(sm.latex(sym, mode="inline"))
    for ax in axes[-1, :]:
        ax.set_xlabel("Time [s]")
    fig.tight_layout()
    return axes


# Symbols for plotting
state_syms = list(q) + list(u) + list(biceps.x) + list(triceps.x)

# Plot with zero excitation
plot_traj(ts, sol_zero.y.T, state_syms)


# Now apply a time-varying excitation: biceps active from 0.5s to 1.5s
def eval_r_pulse(t):
    if 0.5 <= t <= 1.5:
        return np.array([0.8, 0.0])
    else:
        return np.array([0.0, 0.0])


sol_pulse = solve_ivp(
    lambda t, x: eval_rhs(t, x, eval_r_pulse, p_vals), (t0, tf), x0, t_eval=ts
)

# Plot with pulsed excitation
plot_traj(ts, sol_pulse.y.T, state_syms)
plt.show()

# Example from Kane 1985 Sec 4.10 COULOMB FRICTION FORCES

import sympy as sp
from sympy import symbols, integrate, pi
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import KanesMethod, Point, ReferenceFrame
from sympy.physics.mechanics import Particle, RigidBody, inertia

# Define generalized coordinates, speeds, and constants:
g, m, M, R, L = symbols('g, m, M, R, L', real=True, positive=True)
beta, theta = symbols('beta, theta', real=True)
nstar, r = symbols('n^* r', real=True, positive=True)
b1, b2 = symbols('b_1:3', real=True, positive=True)
sigma1, sigma2, sigma3 = symbols('sigma_1:4', real=True)
tau2, tau3 = symbols('tau_2:4', real=True)
rho2, rho3 = symbols('rho_2:4', real=True)
mu1, mu2 = symbols(r"\mu'_1:3", real=True, positive=True)
T = dynamicsymbols('T') 
q = q1, q2, q3 = dynamicsymbols('q_1:4')
qd1, qd2, qd3 = dynamicsymbols('q_1:4', level=1)
u1, u2, u3, u4, u5, u6, u7, u8, u9 = dynamicsymbols('u_1:10')

## --- Define ReferenceFrames ---
N = ReferenceFrame('N')
S = N.orientnew('S', 'axis', [q2, N.x])
R = S.orientnew('R', 'axis', [beta, S.z])
E = S.orientnew('E', 'axis', [theta, S.x])

## --- define points ---
v_O   = u3 * S.x #+ u6 * S.y + u7 * S.z
omega = u2 * S.x #+ u4 * S.y + u8 * S.z
R_v_P = u1 * R.x + u9 * R.y + u5 * R.z

R.set_ang_vel(N, omega)
S.set_ang_vel(N, omega)

## --- Define Points and their velocities ---
pO = Point('O')

pO.set_vel(N, v_O)
pO.set_vel(S, v_O)
pO.set_vel(R, v_O)

pRs = pO.locatenew('R^*', L*R.x)
pRs.set_vel(R, 0)
for frame in [S, N]:
    pRs.v1pt_theory(pO, frame, R)

pP = pO.locatenew('P', q1*R.x)
pP.set_vel(R, R_v_P)
for frame in [S, N]:
    pP.v1pt_theory(pO, frame, R)

pP_ = pO.locatenew(r'\bar{P}', q1*R.x)
pP_.set_vel(R, 0)
for frame in [S, N]:
    pP_.v2pt_theory(pO, frame, R)

pQ = pO.locatenew('Q', r*E.y)
pQ.v2pt_theory(pO, N, S)

## --- Kinematic differential equations ---
kde = [qd1 - u1, qd2 - u2]

## --- Velocity constraints ---
vc = None

## --- Define bodies and forces in the system ---
N_P = rho2 * R.y + rho3 * R.z
T_P = - mu1 * N_P.magnitude() * pP.vel(R).normalize()

t23 = - mu2 * nstar * pQ.vel(N).normalize()
dsigma = -nstar * S.x + t23
integral = lambda i: integrate(integrate(i * r, (theta, 0, 2*pi)), (r, b1, b2))

forces = [
    (pP, m*g * S.x),
    (pRs, M*g * S.x),
    (S, T * S.x + tau2 * S.y + tau3 * S.z),
    (pP,   (T_P + N_P)),
    (pP_, -(T_P + N_P)),
    (pQ, dsigma, integral) # the force need integrate
]

bodies = [
    Particle('P', pP, m),
    RigidBody('R', pRs, R, M, (inertia(R, 0, M*L**2/3, M*L**2/3), pRs))
]

## --- Use Kane's method ---

qind = [q1, q2]
uind = [u1, u2]
uaux = [u3, u5, u9]

km = KanesMethod(
    frame=N, q_ind=qind, u_ind=uind, kd_eqs=kde,
    q_dependent=None, configuration_constraints=None,
    u_dependent=None, velocity_constraints=vc,
    u_auxiliary=uaux, bodies=bodies, forcelist=forces)
fr, fr_star = km.kanes_equations()
fr = fr.subs(dict(zip(uaux, len(uaux) * [0])))

# Kane 1985 Sec 4.10 Eq. (39) - (43)
F1_gt = m*g*sp.cos(beta) - mu1*sp.sqrt(rho2**2 + rho3**2)*u1/sp.sqrt(u1**2)
F2_gt = T - 2*pi*nstar/3*(b2**3 - b1**3)*mu2*u2/sp.sqrt(u2**2)
F3_gt = (m+M)*g - pi*nstar*(b2**2-b1**2)
F5_gt = rho3
F9_gt = -m*g*sp.sin(beta) + rho2
assert((fr[0] - F1_gt).simplify() == 0)
assert((fr[1] - F2_gt).simplify() == 0)
assert((fr[2] - F3_gt).simplify() == 0)
assert((fr[3] - F5_gt).simplify() == 0)
assert((fr[4] - F9_gt).simplify() == 0)

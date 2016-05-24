from sympy import symbols
import sympy.physics.mechanics as me

# System state variables
theta = me.dynamicsymbols('theta')
thetad = me.dynamicsymbols('theta', 1)

# Other system variables
m, l, g = symbols('m l g')

# Set up the reference frames
# Reference frame A set up in the plane perpendicular to the page containing
# segment OP
N = me.ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [theta, N.z])

# Set up the points and particles
O = me.Point('O')
P = O.locatenew('P', l * A.x)

Pa = me.Particle('Pa', P, m)

# Set up velocities
A.set_ang_vel(N, thetad * N.z)
O.set_vel(N, 0)
P.v2pt_theory(O, N, A)

# Set up the lagrangian
L = me.Lagrangian(N, Pa)

# Create the list of forces acting on the system
fl = [(P, g * m * N.x)]


def lmethod():
    # Create the equations of motion using lagranges method
    l = me.LagrangesMethod(L, [theta], forcelist=fl, frame=N)
    l.form_lagranges_equations()

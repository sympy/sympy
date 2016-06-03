from sympy import symbols, Matrix, Derivative, simplify, zeros
import sympy.physics.mechanics as me
import sympy.physics.mechanics.eombase as eombase
from sympy.utilities.pytest import raises


# Should contain access to all the property attributes contained in existing EOM
# generators


# System state variables
q = me.dynamicsymbols('q')
qd = me.dynamicsymbols('q', 1)

# Other system variables
m, k, b, t = symbols('m k b t')

# Set up the reference frame
N = me.ReferenceFrame('N')

# Set up the points and particles
P = me.Point('P')
P.set_vel(N, qd * N.x)

Pa = me.Particle('Pa', P, m)

# Define the potential energy and create the Lagrangian
Pa.potential_energy = k * q**2 / 2.0
L = me.Lagrangian(N, Pa)

# Create the list of forces acting on the system
fl = [(P, -b * qd * N.x)]

# Create an instance of Lagranges Method
l = me.LagrangesMethod(L, [q], forcelist=fl, frame=N)

# Create the equations of motion using lagranges method
l.form_lagranges_equations()

# Set up a full mass matrix and forcing vector for a mass spring damper system
mm = Matrix([[m]])
f = Matrix([[-b*Derivative(q, t) - 1.0*k*q]])

# Set up the kinematical diff eqs
k_kqdot = Matrix([1])
k_ku = Matrix([1])
f_k = Matrix([0])

# Determine what inputs should be used
EOM = eombase.EOM([q], [qd], f, mm, k_kqdot, k_ku, f_k,  loads=fl,
                  bodies=[Pa])

def test_eom_properties():
    # assert EOM.auxiliary_eqs == Change
    assert EOM.bodies == [Pa]
    assert simplify(EOM.mass_matrix - Matrix([[m]])) == Matrix([[0]])
    assert simplify(EOM.mass_matrix_full - Matrix([[1, 0],
                                                   [0, m]])) == zeros(2, 2)
    assert EOM.loads == [(P, -b * qd * N.x)]
    assert simplify(EOM.forcing -
                    Matrix([[-b*Derivative(q, t) - 1.0*k*q]])) == zeros(1)
    assert simplify(EOM.forcing_full -
                    Matrix([[Derivative(q, t)],
                            [-b*Derivative(q, t) - 1.0*k*q]])) == zeros(2, 1)
    assert EOM.coordinates == [q]
    assert EOM.speeds == [qd]

    # The properties should not be able to be altered
    Change = "Attempt to set values for the properties"
    with raises(AttributeError):
        EOM.bodies = Change
    with raises(AttributeError):
        EOM.mass_matrix = Change
    with raises(AttributeError):
        EOM.mass_matrix_full = Change
    with raises(AttributeError):
        EOM.loads = Change
    with raises(AttributeError):
        EOM.forcing = Change
    with raises(AttributeError):
        EOM.forcing_full = Change
    with raises(AttributeError):
        EOM.coordinates = Change
    with raises(AttributeError):
        EOM.speeds = Change

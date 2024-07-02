=====================================
Symbolic Systems in Physics/Mechanics
=====================================

The `SymbolicSystem` class in physics/mechanics is a location for the pertinent
information of a multibody dynamic system. In its most basic form it contains
the equations of motion for the dynamic system, however, it can also contain
information regarding the loads that the system is subject to, the bodies that
the system is comprised of and any additional equations the user feels is
important for the system. The goal of this class is to provide a unified output
format for the equations of motion that numerical analysis code can be designed
around.

SymbolicSystem Example Usage
============================

This code will go over the manual input of the equations of motion for the
simple pendulum that uses the Cartesian location of the mass as the generalized
coordinates into `SymbolicSystem`.

The equations of motion are formed in the physics/mechanics/examples_. In that
spot the variables q1 and q2 are used in place of x and y and the reference
frame is rotated 90 degrees.

.. _examples: https://docs.sympy.org/latest/modules/physics/mechanics/examples/lin_pend_nonmin_example.html

::

    >>> from sympy import atan, symbols, Matrix
    >>> from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame,
    ...                                      Particle, Point)
    >>> import sympy.physics.mechanics.system as system
    >>> from sympy.physics.vector import init_vprinting
    >>> init_vprinting(pretty_print=False)

The first step will be to initialize all of the dynamic and constant symbols. ::

    >>> x, y, u, v, lam = dynamicsymbols('x y u v lambda')
    >>> m, l, g = symbols('m l g')

Next step is to define the equations of motion in multiple forms:

    [1] Explicit form where the kinematics and dynamics are combined
        x' = F_1(x, t, r, p)

    [2] Implicit form where the kinematics and dynamics are combined
        M_2(x, p) x' = F_2(x, t, r, p)

    [3] Implicit form where the kinematics and dynamics are separate
        M_3(q, p) u' = F_3(q, u, t, r, p)
        q' = G(q, u, t, r, p)

where

    x : states, e.g. [q, u]
    t : time
    r : specified (exogenous) inputs
    p : constants
    q : generalized coordinates
    u : generalized speeds
    F_1 : right hand side of the combined equations in explicit form
    F_2 : right hand side of the combined equations in implicit form
    F_3 : right hand side of the dynamical equations in implicit form
    M_2 : mass matrix of the combined equations in implicit form
    M_3 : mass matrix of the dynamical equations in implicit form
    G : right hand side of the kinematical differential equations

::

    >>> dyn_implicit_mat = Matrix([[1, 0, -x/m],
    ...                            [0, 1, -y/m],
    ...                            [0, 0, l**2/m]])
    >>> dyn_implicit_rhs = Matrix([0, 0, u**2 + v**2 - g*y])
    >>> comb_implicit_mat = Matrix([[1, 0, 0, 0, 0],
    ...                             [0, 1, 0, 0, 0],
    ...                             [0, 0, 1, 0, -x/m],
    ...                             [0, 0, 0, 1, -y/m],
    ...                             [0, 0, 0, 0, l**2/m]])
    >>> comb_implicit_rhs = Matrix([u, v, 0, 0, u**2 + v**2 - g*y])
    >>> kin_explicit_rhs = Matrix([u, v])
    >>> comb_explicit_rhs = comb_implicit_mat.LUsolve(comb_implicit_rhs)

Now the reference frames, points and particles will be set up so this
information can be passed into `system.SymbolicSystem` in the form of a bodies
and loads iterable. ::

    >>> theta = atan(x/y)
    >>> omega = dynamicsymbols('omega')
    >>> N = ReferenceFrame('N')
    >>> A = N.orientnew('A', 'Axis', [theta, N.z])
    >>> A.set_ang_vel(N, omega * N.z)
    >>> O = Point('O')
    >>> O.set_vel(N, 0)
    >>> P = O.locatenew('P', l * A.x)
    >>> P.v2pt_theory(O, N, A)
    l*omega*A.y
    >>> Pa = Particle('Pa', P, m)

Now the bodies and loads iterables need to be initialized. ::

    >>> bodies = [Pa]
    >>> loads = [(P, g * m * N.x)]

The equations of motion are in the form of a differential algebraic equation
(DAE) and DAE solvers need to know which of the equations are the algebraic
expressions. This information is passed into `SymbolicSystem` as a list
specifying which rows are the algebraic equations. In this example it is a
different row based on the chosen equations of motion format. The row index
should always correspond to the mass matrix that is being input to the
`SymbolicSystem` class but will always correspond to the row index of the
combined dynamics and kinematics when being accessed from the `SymbolicSystem`
class. ::

    >>> alg_con = [2]
    >>> alg_con_full = [4]

An iterable containing the states now needs to be created for the system. The
`SymbolicSystem` class can determine which of the states are considered
coordinates or speeds by passing in the indexes of the coordinates and speeds.
If these indexes are not passed in the object will not be able to differentiate
between coordinates and speeds. ::

    >>> states = (x, y, u, v, lam)
    >>> coord_idxs = (0, 1)
    >>> speed_idxs = (2, 3)

Now the equations of motion instances can be created using the above mentioned
equations of motion formats. ::

    >>> symsystem1 = system.SymbolicSystem(states, comb_explicit_rhs,
    ...                                    alg_con=alg_con_full, bodies=bodies,
    ...                                    loads=loads)
    >>> symsystem2 = system.SymbolicSystem(states, comb_implicit_rhs,
    ...                                    mass_matrix=comb_implicit_mat,
    ...                                    alg_con=alg_con_full,
    ...                                    coord_idxs=coord_idxs)
    >>> symsystem3 = system.SymbolicSystem(states, dyn_implicit_rhs,
    ...                                    mass_matrix=dyn_implicit_mat,
    ...                                    coordinate_derivatives=kin_explicit_rhs,
    ...                                    alg_con=alg_con,
    ...                                    coord_idxs=coord_idxs,
    ...                                    speed_idxs=speed_idxs)

Like coordinates and speeds, the bodies and loads attributes can only be
accessed if they are specified during initialization of the `SymbolicSystem`
class. Lastly here are some attributes accessible from the `SymbolicSystem`
class. ::

    >>> symsystem1.states
    Matrix([
    [     x],
    [     y],
    [     u],
    [     v],
    [lambda]])
    >>> symsystem2.coordinates
    Matrix([
    [x],
    [y]])
    >>> symsystem3.speeds
    Matrix([
    [u],
    [v]])
    >>> symsystem1.comb_explicit_rhs
    Matrix([
    [                          u],
    [                          v],
    [(-g*y + u**2 + v**2)*x/l**2],
    [(-g*y + u**2 + v**2)*y/l**2],
    [m*(-g*y + u**2 + v**2)/l**2]])
    >>> symsystem2.comb_implicit_rhs
    Matrix([
    [                 u],
    [                 v],
    [                 0],
    [                 0],
    [-g*y + u**2 + v**2]])
    >>> symsystem2.comb_implicit_mat
    Matrix([
    [1, 0, 0, 0,      0],
    [0, 1, 0, 0,      0],
    [0, 0, 1, 0,   -x/m],
    [0, 0, 0, 1,   -y/m],
    [0, 0, 0, 0, l**2/m]])
    >>> symsystem3.dyn_implicit_rhs
    Matrix([
    [                 0],
    [                 0],
    [-g*y + u**2 + v**2]])
    >>> symsystem3.dyn_implicit_mat
    Matrix([
    [1, 0,   -x/m],
    [0, 1,   -y/m],
    [0, 0, l**2/m]])
    >>> symsystem3.kin_explicit_rhs
    Matrix([
    [u],
    [v]])
    >>> symsystem1.alg_con
    [4]
    >>> symsystem1.bodies
    (Pa,)
    >>> symsystem1.loads
    ((P, g*m*N.x),)

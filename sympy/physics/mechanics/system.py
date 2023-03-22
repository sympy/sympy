from sympy.core.backend import (eye, zeros, Basic, ImmutableMatrix as Matrix)
from sympy.core.containers import OrderedSet
from sympy.utilities.iterables import iterable
from sympy.physics.vector import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics.functions import (
    find_dynamicsymbols, _validate_coordinates, Lagrangian)
from sympy.physics.mechanics.body_base import BodyBase
from sympy.physics.mechanics.particle import Particle
from sympy.physics.mechanics.loads import Force, Torque, _parse_load
from sympy.physics.mechanics.joint import Joint
from sympy.physics.mechanics.method import _Methods
from sympy.physics.mechanics.kane import KanesMethod
from sympy.physics.mechanics.lagrange import LagrangesMethod
from functools import wraps

__all__ = ['System', 'SymbolicSystem']


def _reset_eom_method(method):
    """Decorator to reset the eom_method if a property is changed."""

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        self._eom_method = None
        return method(self, *args, **kwargs)

    return wrapper


class System(_Methods):
    """Class to define a system and form the equations of motion.

    Explanation
    ===========

    A ``System`` instance stores the different objects associated with a model,
    including bodies, joints, constraints, and other relevant information. With
    all the relationships between components defined, the ``System`` can be used
    to form the equations of motion using a backend, such as ``KanesMethod``.
    The ``System`` has been designed to be compatible with third-party
    libraries for greater flexibility and integration with other tools.

    Attributes
    ==========

    origin : Point
        Global origin of the system.
    frame : ReferenceFrame
        Inertial reference frame of the system.
    x : Vector
        Unit vector in the x direction of the inertial reference frame.
    y : Vector
        Unit vector in the y direction of the inertial reference frame.
    z : Vector
        Unit vector in the z direction of the inertial reference frame.
    q : Matrix
        Matrix of all the generalized coordinates.
    u : Matrix
        Matrix of all the generalized speeds.
    q_ind : Matrix
        Matrix of the independent generalized coordinates.
    q_dep : Matrix
        Matrix of the dependent generalized coordinates.
    u_ind : Matrix
        Matrix of the independent generalized speeds.
    u_dep : Matrix
        Matrix of the dependent generalized speeds.
    kdes : Matrix
        Matrix of the kinematic differential equations.
    bodies : tuple of BodyBase subclasses
        Tuple of all bodies that make up the system.
    joints : tuple of Joint
        Tuple of all joints that connect bodies in the system.
    loads : tuple of LoadBase subclasses
        Tuple of all loads that have been applied to the system.
    holonomic_constraints : Matrix
        Matrix with the holonomic constraints as rows.
    nonholonomic_constraints : Matrix
        Matrix with the nonholonomic constraints as rows.
    eom_method : class
        Backend for forming the equations of motion.

    """

    def __init__(self, origin=None, frame=None):
        if origin is None:
            origin = Point('inertial_origin')
        elif not isinstance(origin, Point):
            raise TypeError('Origin must be an instance of Point.')
        self._origin = origin
        if frame is None:
            frame = ReferenceFrame('inertial_frame')
        elif not isinstance(frame, ReferenceFrame):
            raise TypeError('Frame must be an instance of ReferenceFrame.')
        self._frame = frame
        self._origin.set_vel(self._frame, 0)
        self._q_ind = Matrix(1, 0, []).T
        self._q_dep = Matrix(1, 0, []).T
        self._u_ind = Matrix(1, 0, []).T
        self._u_dep = Matrix(1, 0, []).T
        self._kdes = Matrix(1, 0, []).T
        self._hol_coneqs = Matrix(1, 0, []).T
        self._nonhol_coneqs = Matrix(1, 0, []).T
        self._bodies = []
        self._joints = []
        self._loads = []
        self._eom_method = None

    @classmethod
    def from_newtonian(cls, newtonian):
        """Constructs the system with respect to a Newtonian body."""
        if isinstance(newtonian, Particle):
            raise ValueError('A Particle has no frame so cannot act as '
                             'the Newtonian.')
        system = cls(origin=newtonian.masscenter, frame=newtonian.frame)
        system.add_bodies(newtonian)
        return system

    @property
    def origin(self):
        """Global origin of the system."""
        return self._origin

    @property
    def frame(self):
        """Inertial reference frame of the system."""
        return self._frame

    @property
    def x(self):
        """Unit vector in the x direction of the inertial reference frame."""
        return self._frame.x

    @property
    def y(self):
        """Unit vector in the y direction of the inertial reference frame."""
        return self._frame.y

    @property
    def z(self):
        """Unit vector in the z direction of the inertial reference frame."""
        return self._frame.z

    @property
    def bodies(self):
        """Tuple of all bodies that have been added to the system."""
        return tuple(self._bodies)

    @bodies.setter
    @_reset_eom_method
    def bodies(self, bodies):
        bodies = System._objects_to_list(bodies)
        self._check_objects(bodies, [], BodyBase, 'Bodies', 'bodies')
        self._bodies = bodies

    @property
    def joints(self):
        """Tuple of all joints that have been added to the system."""
        return tuple(self._joints)

    @property
    def loads(self):
        """Tuple of loads that have been applied on the system."""
        return tuple(self._loads)

    @loads.setter
    @_reset_eom_method
    def loads(self, loads):
        loads = System._objects_to_list(loads)
        self._loads = [_parse_load(load) for load in loads]

    @property
    def q(self):
        """Matrix of all the generalized coordinates."""
        return self._q_ind.col_join(self._q_dep)

    @property
    def u(self):
        """Matrix of all the generalized speeds."""
        return self._u_ind.col_join(self._u_dep)

    @property
    def q_ind(self):
        """Matrix of the independent generalized coordinates."""
        return self._q_ind

    @q_ind.setter
    @_reset_eom_method
    def q_ind(self, q_ind):
        self._q_ind, self._q_dep = self._parse_coordinates(
            self._objects_to_list(q_ind), True, [], self.q_dep, True)

    @property
    def q_dep(self):
        """Matrix of the dependent generalized coordinates."""
        return self._q_dep

    @q_dep.setter
    @_reset_eom_method
    def q_dep(self, q_dep):
        self._q_ind, self._q_dep = self._parse_coordinates(
            self._objects_to_list(q_dep), False, self.q_ind, [], True)

    @property
    def u_ind(self):
        """Matrix of the independent generalized speeds."""
        return self._u_ind

    @u_ind.setter
    @_reset_eom_method
    def u_ind(self, u_ind):
        self._u_ind, self._u_dep = self._parse_coordinates(
            self._objects_to_list(u_ind), True, [], self.u_dep, False)

    @property
    def u_dep(self):
        """Matrix of the dependent generalized speeds."""
        return self._u_dep

    @u_dep.setter
    @_reset_eom_method
    def u_dep(self, u_dep):
        self._u_ind, self._u_dep = self._parse_coordinates(
            self._objects_to_list(u_dep), False, self.u_ind, [], False)

    @property
    def kdes(self):
        """Kinematic differential equations, which describe the coupling
        between the generalized coordinates and the generalized speeds."""
        return self._kdes

    @kdes.setter
    @_reset_eom_method
    def kdes(self, kdes):
        kdes = System._objects_to_list(kdes)
        self._kdes = self._parse_expressions(
            kdes, [], 'kinematic differential equations')

    @property
    def holonomic_constraints(self):
        """Matrix with the holonomic constraints as rows."""
        return self._hol_coneqs

    @holonomic_constraints.setter
    @_reset_eom_method
    def holonomic_constraints(self, constraints):
        constraints = System._objects_to_list(constraints)
        self._hol_coneqs = self._parse_expressions(
            constraints, [], 'holonomic constraints')

    @property
    def nonholonomic_constraints(self):
        """Matrix with the nonholonomic constraints as rows."""
        return self._nonhol_coneqs

    @nonholonomic_constraints.setter
    @_reset_eom_method
    def nonholonomic_constraints(self, constraints):
        constraints = System._objects_to_list(constraints)
        self._nonhol_coneqs = self._parse_expressions(
            constraints, [], 'nonholonomic constraints')

    @property
    def eom_method(self):
        """Backend for forming the equations of motion."""
        return self._eom_method

    @staticmethod
    def _objects_to_list(lst):
        """Helper to convert passed objects to a list."""
        if not iterable(lst):  # Only one object
            return [lst]
        return list(lst[:])  # converts Matrix and tuple to flattened list

    @staticmethod
    def _check_objects(objects, obj_lst, expected_type, obj_name, type_name):
        """Helper to check the objects that are being added to the system."""
        seen = set(obj_lst)
        duplicates = set()
        wrong_types = set()
        for obj in objects:
            if not isinstance(obj, expected_type):
                wrong_types.add(obj)
            if obj in seen:
                duplicates.add(obj)
            else:
                seen.add(obj)
        if wrong_types:
            raise TypeError(f'{obj_name} {wrong_types} are not {type_name}.')
        if duplicates:
            raise ValueError(f'{obj_name} {duplicates} have already been added '
                             f'to the system.')

    def _parse_coordinates(self, new_coords, independent, old_coords_ind,
                           old_coords_dep, is_coordinates=True):
        """Helper to parse coordinates and speeds."""
        # Construct lists of the independent and dependent coordinates
        if is_coordinates:
            q_ind, q_dep = old_coords_ind[:], old_coords_dep[:]
            u_ind, u_dep = self.u_ind[:], self.u_dep[:]
            coords_ind, coords_dep = q_ind, q_dep
        else:
            q_ind, q_dep = self.q_ind[:], self.q_dep[:]
            u_ind, u_dep = old_coords_ind[:], old_coords_dep[:]
            coords_ind, coords_dep = u_ind, u_dep
        if not iterable(independent):
            independent = [independent] * len(new_coords)
        for coord, indep in zip(new_coords, independent):
            if indep:
                coords_ind.append(coord)
            else:
                coords_dep.append(coord)
        # Check types and duplicates
        _validate_coordinates(q_ind + q_dep, u_ind + u_dep)
        return (Matrix(1, len(coords_ind), coords_ind).T,
                Matrix(1, len(coords_dep), coords_dep).T)

    @staticmethod
    def _parse_expressions(new_expressions, old_expressions, name,
                           check_negatives=False):
        """Helper to parse expressions like constraints."""
        old_expressions = old_expressions[:]
        new_expressions = list(new_expressions)  # Converts a possible tuple
        if check_negatives:
            check_exprs = old_expressions + [-expr for expr in old_expressions]
        else:
            check_exprs = old_expressions
        System._check_objects(new_expressions, check_exprs, Basic, name,
                              'expressions')
        for expr in new_expressions:
            if expr == 0:
                raise ValueError(f'Parsed {name} are zero.')
        return Matrix(1, len(old_expressions) + len(new_expressions),
                      old_expressions + new_expressions).T

    @_reset_eom_method
    def add_coordinates(self, *coordinates, independent=True):
        """Add generalized coordinate(s) to the system.

        Parameters
        ==========

        *coordinates : dynamicsymbols
            One or more generalized coordinates to be added to the system.
        independent : bool or list of bool, optional
            Boolean whether a coordinate is dependent or independent. The
            default is True, so the coordinates are added as independent by
            default.

        """
        self._q_ind, self._q_dep = self._parse_coordinates(
            coordinates, independent, self.q_ind, self.q_dep, True)

    @_reset_eom_method
    def add_speeds(self, *speeds, independent=True):
        """Add generalized speed(s) to the system.

        Parameters
        ==========

        *speeds : dynamicsymbols
            One or more generalized speeds to be added to the system.
        independent : bool or list of bool, optional
            Boolean whether a speed is dependent or independent. The default is
            True, so the speeds are added as independent by default.

        """
        self._u_ind, self._u_dep = self._parse_coordinates(
            speeds, independent, self.u_ind, self.u_dep, False)

    @_reset_eom_method
    def add_kdes(self, *kdes):
        """Add kinematic differential equation(s) to the system.

        Parameters
        ==========

        *kdes : Expr
            One or more kinematic differential equations.

        """
        self._kdes = self._parse_expressions(
            kdes, self.kdes, 'kinematic differential equations',
            check_negatives=True)

    @_reset_eom_method
    def add_holonomic_constraints(self, *constraints):
        """Add holonomic constraint(s) to the system.

        Parameters
        ==========

        *constraints : Expr
            One or more holonomic constraints.

        """
        self._hol_coneqs = self._parse_expressions(
            constraints, self._hol_coneqs, 'holonomic constraints',
            check_negatives=True)

    @_reset_eom_method
    def add_nonholonomic_constraints(self, *constraints):
        """Add nonholonomic constraint(s) to the system.

        Parameters
        ==========

        *constraints : Expr
            One or more nonholonomic constraints.

        """
        self._nonhol_coneqs = self._parse_expressions(
            constraints, self._nonhol_coneqs, 'nonholonomic constraints',
            check_negatives=True)

    @_reset_eom_method
    def add_bodies(self, *bodies):
        """Add body(ies) to the system.

        Parameters
        ==========

        bodies : Particle or RigidBody
            One or more bodies.

        """
        self._check_objects(bodies, self.bodies, BodyBase, 'Bodies', 'bodies')
        self._bodies.extend(bodies)

    @_reset_eom_method
    def add_loads(self, *loads):
        """Add load(s) to the system.

        Parameters
        ==========

        *loads : Force or Torque
            One or more loads.

        """
        loads = [_parse_load(load) for load in loads]  # Checks the loads
        self._loads.extend(loads)

    @_reset_eom_method
    def apply_force(self, point, force, reaction_point=None):
        """Apply a force to the system.

        Explanation
        ===========

        Applies the force on point. If a reaction point is supplied, an equal
        and oppposite force is applied to the reaction point, i.e. -force.

        Parameters
        ==========

        point: Point or Particle or RigidBody
            The point on which the force is applied. If a body is supplied, then
            the center of mass of that body is used.
        force: Vector
            The force to be applied.
        reaction_point : Point or Particle or RigidBody, optional
            The point on which an equal and opposite force is applied. If a body
            is supplied, then the center of mass of that body is used. The
            default is None.

        Example
        =======

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import System, RigidBody
        >>> system = System()
        >>> g = symbols('g')
        >>> rb = RigidBody('rb')
        >>> system.apply_force(rb.masscenter, - rb.mass * g * system.z)
        >>> system.loads
        ((rb_masscenter, - g*rb_mass*inertial_frame.z),)

        To further demonstrate the use of ``apply_force`` attribute, consider
        two bodies connected by a spring. For ease the bodies centers of mass
        are connected.

        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> x = dynamicsymbols('x')
        >>> rb1 = RigidBody('rb1')
        >>> rb2 = RigidBody('rb2')
        >>> rb2.masscenter.set_pos(rb1.masscenter, x * system.x)
        >>> spring_force = x * rb1.masscenter.pos_from(rb2.masscenter)
        >>> system.apply_force(rb1, spring_force, rb2)
        >>> system.loads
        ((rb_masscenter, - g*rb_mass*inertial_frame.z), (rb1_masscenter, - x(t)**2*inertial_frame.x), (rb2_masscenter, x(t)**2*inertial_frame.x))

        Notes
        =====

        Method as it may change a bit with the introduction of Actuator.

        """
        self._loads.append(Force(point, force))
        if reaction_point is not None:
            self._loads.append(Force(reaction_point, -force))

    @_reset_eom_method
    def apply_torque(self, frame, torque, reaction_frame=None):
        """Apply a torque to the system.

        Explanation
        ===========

        Applies the torque on frame. If a reaction frame is supplied, an equal
        and oppposite torque is applied to the reaction point, i.e. -force.

        Parameters
        ==========

        frame: ReferenceFrame or Particle or RigidBody
            The frame to which the torque is applied. If a body is supplied,
            then the frame associated with that body is used.
        torque: Vector
            The torque to be applied.
        reaction_frame : ReferenceFrame or Particle or RigidBody, optional
            The frame on which an equal and opposite torque is applied. If a
            body is supplied, then the frame associated with that body is used.
            The default is None.

        Example
        =======

        >>> from sympy.physics.mechanics import RigidBody, System
        >>> system = System()
        >>> rb = RigidBody('rb')
        >>> system.apply_torque(rb.frame, -system.z)
        >>> system.loads
        ((rb_frame, - inertial_frame.z),)

        To further demonstrate the use, let us consider two bodies such that
        a torque ``T`` is acting on one body, and ``-T`` on the other.

        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> rb1 = RigidBody('rb1')
        >>> rb2 = RigidBody('rb2')
        >>> v = dynamicsymbols('v')
        >>> T = v * system.y  # Torque
        >>> system.apply_torque(rb1, T, rb2)
        >>> system.loads
         ((rb_frame, - inertial_frame.z), (rb1_frame, v(t)*inertial_frame.y), (rb2_frame, - v(t)*inertial_frame.y))

        Notes
        =====

        Method as it may change a bit with the introduction of Actuator.

        """
        self._loads.append(Torque(frame, torque))
        if reaction_frame is not None:
            self._loads.append(Torque(reaction_frame, -torque))

    @_reset_eom_method
    def remove_load(self, location):
        """Remove the loads applied about a point or frame.

        Parameters
        ==========

        location : Point or ReferenceFrame or Particle or RigidBody
            Location about which the loads should be removed. If a body is
            provided then the loads acting upon its center of mass and
            associated frame are removed.

        Returns
        =======

        tuple of Force and Torque
            The removed loads.

        Examples
        ========

        >>> from sympy.physics.mechanics import RigidBody, System
        >>> system = System()
        >>> rb = RigidBody('rb')
        >>> system.apply_force(rb, system.x)
        >>> system.apply_torque(rb, system.z)
        >>> system.loads
        ((rb_masscenter, inertial_frame.x), (rb_frame, inertial_frame.z))
        >>> system.remove_load(rb.masscenter)
        ((rb_masscenter, inertial_frame.x),)
        >>> system.loads
        ((rb_frame, inertial_frame.z),)

        """
        if isinstance(location, BodyBase):
            removed_loads = []
            removed_loads.extend(self.remove_load(location.masscenter))
            if hasattr(location, 'frame'):  # Particle has no frame
                removed_loads.extend(self.remove_load(location.frame))
            return tuple(removed_loads)
        removed_loads = []
        updated_loads_list = []
        for ld in self._loads:
            if ld.location == location:
                removed_loads.append(ld)
            else:
                updated_loads_list.append(ld)
        self._loads = updated_loads_list
        return tuple(removed_loads)

    @_reset_eom_method
    def clear_loads(self):
        """Remove all loads from the system"""
        self._loads = []

    @_reset_eom_method
    def add_joints(self, *joints):
        """Add joint(s) to the system.

        Explanation
        ===========

        This methods adds one or more joints to the system including its
        associated objects, i.e. generalized coordinates, generalized speeds,
        kinematic differential equations and the bodies.

        Parameters
        ==========

        *joints : subclass of Joint
            One or more joints.

        Notes
        =====

        For the generalized coordinates, generalized speeds and bodies it is
        checked whether they are already known by the system instance. If they
        are, then they are not added. The kinematic differential equations are
        however always added to the system, so do not also add those manually
        beforehand.

        """
        self._check_objects(joints, self.joints, Joint, 'Joints', 'joints')
        self._joints.extend(joints)
        coordinates, speeds, kdes, bodies = (OrderedSet() for _ in range(4))
        for joint in joints:
            coordinates.update(joint.coordinates)
            speeds.update(joint.speeds)
            kdes.update(joint.kdes)
            bodies.update((joint.parent, joint.child))
        coordinates = coordinates.difference(self.q)
        speeds = speeds.difference(self.u)
        kdes = kdes.difference(self.kdes[:] + (-self.kdes)[:])
        bodies = bodies.difference(self.bodies)
        self.add_coordinates(*tuple(coordinates))
        self.add_speeds(*tuple(speeds))
        self.add_kdes(*(kde for kde in tuple(kdes) if not kde == 0))
        self.add_bodies(*tuple(bodies))

    def get_body(self, name):
        """Retrieve a body from the system by name.

        Parameters
        ==========

        name : str
            The name of the body to retrieve.

        Returns
        =======

        RigidBody or Particle
            The body with the given name, or None if no such body exists.

        """
        for body in self._bodies:
            if body.name == name:
                return body

    def _form_eoms(self):
        return self.form_eoms()

    def form_eoms(self, eom_method=KanesMethod):
        """Form the equations of motion of the system.

        Parameters
        ==========

        eom_method : class
            Backend class to be used for forming the equations of motion. The
            default is ``KanesMethod``.

        Returns
        ========

        Matrix
            Vector of equations of motions.

        Examples
        ========

        This is a simple example for a one degree of freedom translational
        spring-mass-damper.

        >>> from sympy import S, symbols
        >>> from sympy.physics.mechanics import (
        ...     LagrangesMethod, dynamicsymbols, PrismaticJoint, System,
        ...     Particle, RigidBody)
        >>> q = dynamicsymbols('q')
        >>> qd = dynamicsymbols('q', 1)
        >>> m, k, b = symbols('m k b')
        >>> wall = RigidBody('W')
        >>> system = System.from_newtonian(wall)
        >>> bob = Particle('P', mass=m)
        >>> bob.potential_energy = S.Half * k * q**2
        >>> system.add_joints(PrismaticJoint('J', wall, bob, q, qd))
        >>> system.apply_force(bob, b * qd * system.x, reaction_point=wall)
        >>> system.form_eoms(LagrangesMethod)
        Matrix([[-b*Derivative(q(t), t) + k*q(t) + m*Derivative(q(t), (t, 2))]])

        We can also solve for the states using the 'rhs' method.

        >>> system.rhs()
        Matrix([
        [               Derivative(q(t), t)],
        [(b*Derivative(q(t), t) - k*q(t))/m]])

        """
        # KanesMethod does not accept empty iterables
        loads = self.loads if self.loads else None
        if issubclass(eom_method, KanesMethod):
            self._eom_method = KanesMethod(
                self.frame, self.q_ind, self.u_ind, kd_eqs=self.kdes,
                q_dependent=self.q_dep, u_dependent=self.u_dep,
                configuration_constraints=self.holonomic_constraints,
                velocity_constraints=self.holonomic_constraints.diff(
                    dynamicsymbols._t).col_join(self.nonholonomic_constraints),
                forcelist=loads, bodies=self.bodies,
                explicit_kinematics=False)
        elif issubclass(eom_method, LagrangesMethod):
            self._eom_method = eom_method(
                Lagrangian(self.frame, *self.bodies), self.q, loads,
                self.bodies, self.frame, self.holonomic_constraints,
                self.nonholonomic_constraints)
        else:
            raise NotImplementedError(f'{eom_method} has not been implemented.')
        return self.eom_method._form_eoms()

    def rhs(self, inv_method=None):
        """Returns equations that can be solved numerically.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrices.MatrixBase.inv`

        Returns
        ========

        Matrix
            Numerically solvable equations.

        See Also
        ========

        sympy.physics.mechanics.kane.KanesMethod.rhs:
            KanesMethod's rhs function.
        sympy.physics.mechanics.lagrange.LagrangesMethod.rhs:
            LagrangesMethod's rhs function.

        """
        return self.eom_method.rhs(inv_method=inv_method)

    @property
    def mass_matrix(self):
        r"""The mass matrix of the system.

        Explanation
        ===========

        The mass matrix $M_d$ and the forcing vector $f_d$ of a system describe
        the system's dynamics according to the following equations:
        $$M_d \dot{u} = f_d$$
        where $\dot{u}$ is the time derivative of the generalized speeds.

        """
        return self.eom_method.mass_matrix

    @property
    def mass_matrix_full(self):
        r"""The mass matrix of the system, augmented by the kinematic
        differential equations in explicit or implicit form.

        Explanation
        ===========

        The full mass matrix $M_m$ and the full forcing vector $f_m$ of a system
        describe the dynamics and kinematics according to the following
        equation:
        $$M_m \dot{x} = f_m$$
        where $x$ is the state vector stacking $q$ and $u$.

        """
        return self.eom_method.mass_matrix_full

    @property
    def forcing(self):
        """The forcing vector of the system."""
        return self.eom_method.forcing

    @property
    def forcing_full(self):
        """The forcing vector of the system, augmented by the kinematic
        differential equations in explicit or implicit form."""
        return self.eom_method.forcing_full


class SymbolicSystem:
    """SymbolicSystem is a class that contains all the information about a
    system in a symbolic format such as the equations of motions and the bodies
    and loads in the system.

    There are three ways that the equations of motion can be described for
    Symbolic System:


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

        Parameters
        ==========

        coord_states : ordered iterable of functions of time
            This input will either be a collection of the coordinates or states
            of the system depending on whether or not the speeds are also
            given. If speeds are specified this input will be assumed to
            be the coordinates otherwise this input will be assumed to
            be the states.

        right_hand_side : Matrix
            This variable is the right hand side of the equations of motion in
            any of the forms. The specific form will be assumed depending on
            whether a mass matrix or coordinate derivatives are given.

        speeds : ordered iterable of functions of time, optional
            This is a collection of the generalized speeds of the system. If
            given it will be assumed that the first argument (coord_states)
            will represent the generalized coordinates of the system.

        mass_matrix : Matrix, optional
            The matrix of the implicit forms of the equations of motion (forms
            [2] and [3]). The distinction between the forms is determined by
            whether or not the coordinate derivatives are passed in. If
            they are given form [3] will be assumed otherwise form [2] is
            assumed.

        coordinate_derivatives : Matrix, optional
            The right hand side of the kinematical equations in explicit form.
            If given it will be assumed that the equations of motion are being
            entered in form [3].

        alg_con : Iterable, optional
            The indexes of the rows in the equations of motion that contain
            algebraic constraints instead of differential equations. If the
            equations are input in form [3], it will be assumed the indexes are
            referencing the mass_matrix/right_hand_side combination and not the
            coordinate_derivatives.

        output_eqns : Dictionary, optional
            Any output equations that are desired to be tracked are stored in a
            dictionary where the key corresponds to the name given for the
            specific equation and the value is the equation itself in symbolic
            form

        coord_idxs : Iterable, optional
            If coord_states corresponds to the states rather than the
            coordinates this variable will tell SymbolicSystem which indexes of
            the states correspond to generalized coordinates.

        speed_idxs : Iterable, optional
            If coord_states corresponds to the states rather than the
            coordinates this variable will tell SymbolicSystem which indexes of
            the states correspond to generalized speeds.

        bodies : iterable of Body/Rigidbody objects, optional
            Iterable containing the bodies of the system

        loads : iterable of load instances (described below), optional
            Iterable containing the loads of the system where forces are given
            by (point of application, force vector) and torques are given by
            (reference frame acting upon, torque vector). Ex [(point, force),
            (ref_frame, torque)]

    Attributes
    ==========

    coordinates : Matrix, shape(n, 1)
        This is a matrix containing the generalized coordinates of the system

    speeds : Matrix, shape(m, 1)
        This is a matrix containing the generalized speeds of the system

    states : Matrix, shape(o, 1)
        This is a matrix containing the state variables of the system

    alg_con : List
        This list contains the indices of the algebraic constraints in the
        combined equations of motion. The presence of these constraints
        requires that a DAE solver be used instead of an ODE solver.
        If the system is given in form [3] the alg_con variable will be
        adjusted such that it is a representation of the combined kinematics
        and dynamics thus make sure it always matches the mass matrix
        entered.

    dyn_implicit_mat : Matrix, shape(m, m)
        This is the M matrix in form [3] of the equations of motion (the mass
        matrix or generalized inertia matrix of the dynamical equations of
        motion in implicit form).

    dyn_implicit_rhs : Matrix, shape(m, 1)
        This is the F vector in form [3] of the equations of motion (the right
        hand side of the dynamical equations of motion in implicit form).

    comb_implicit_mat : Matrix, shape(o, o)
        This is the M matrix in form [2] of the equations of motion.
        This matrix contains a block diagonal structure where the top
        left block (the first rows) represent the matrix in the
        implicit form of the kinematical equations and the bottom right
        block (the last rows) represent the matrix in the implicit form
        of the dynamical equations.

    comb_implicit_rhs : Matrix, shape(o, 1)
        This is the F vector in form [2] of the equations of motion. The top
        part of the vector represents the right hand side of the implicit form
        of the kinemaical equations and the bottom of the vector represents the
        right hand side of the implicit form of the dynamical equations of
        motion.

    comb_explicit_rhs : Matrix, shape(o, 1)
        This vector represents the right hand side of the combined equations of
        motion in explicit form (form [1] from above).

    kin_explicit_rhs : Matrix, shape(m, 1)
        This is the right hand side of the explicit form of the kinematical
        equations of motion as can be seen in form [3] (the G matrix).

    output_eqns : Dictionary
        If output equations were given they are stored in a dictionary where
        the key corresponds to the name given for the specific equation and
        the value is the equation itself in symbolic form

    bodies : Tuple
        If the bodies in the system were given they are stored in a tuple for
        future access

    loads : Tuple
        If the loads in the system were given they are stored in a tuple for
        future access. This includes forces and torques where forces are given
        by (point of application, force vector) and torques are given by
        (reference frame acted upon, torque vector).

    Example
    =======

    As a simple example, the dynamics of a simple pendulum will be input into a
    SymbolicSystem object manually. First some imports will be needed and then
    symbols will be set up for the length of the pendulum (l), mass at the end
    of the pendulum (m), and a constant for gravity (g). ::

        >>> from sympy import Matrix, sin, symbols
        >>> from sympy.physics.mechanics import dynamicsymbols, SymbolicSystem
        >>> l, m, g = symbols('l m g')

    The system will be defined by an angle of theta from the vertical and a
    generalized speed of omega will be used where omega = theta_dot. ::

        >>> theta, omega = dynamicsymbols('theta omega')

    Now the equations of motion are ready to be formed and passed to the
    SymbolicSystem object. ::

        >>> kin_explicit_rhs = Matrix([omega])
        >>> dyn_implicit_mat = Matrix([l**2 * m])
        >>> dyn_implicit_rhs = Matrix([-g * l * m * sin(theta)])
        >>> symsystem = SymbolicSystem([theta], dyn_implicit_rhs, [omega],
        ...                            dyn_implicit_mat)

    Notes
    =====

    m : number of generalized speeds
    n : number of generalized coordinates
    o : number of states

    """

    def __init__(self, coord_states, right_hand_side, speeds=None,
                 mass_matrix=None, coordinate_derivatives=None, alg_con=None,
                 output_eqns={}, coord_idxs=None, speed_idxs=None, bodies=None,
                 loads=None):
        """Initializes a SymbolicSystem object"""

        # Extract information on speeds, coordinates and states
        if speeds is None:
            self._states = Matrix(coord_states)

            if coord_idxs is None:
                self._coordinates = None
            else:
                coords = [coord_states[i] for i in coord_idxs]
                self._coordinates = Matrix(coords)

            if speed_idxs is None:
                self._speeds = None
            else:
                speeds_inter = [coord_states[i] for i in speed_idxs]
                self._speeds = Matrix(speeds_inter)
        else:
            self._coordinates = Matrix(coord_states)
            self._speeds = Matrix(speeds)
            self._states = self._coordinates.col_join(self._speeds)

        # Extract equations of motion form
        if coordinate_derivatives is not None:
            self._kin_explicit_rhs = coordinate_derivatives
            self._dyn_implicit_rhs = right_hand_side
            self._dyn_implicit_mat = mass_matrix
            self._comb_implicit_rhs = None
            self._comb_implicit_mat = None
            self._comb_explicit_rhs = None
        elif mass_matrix is not None:
            self._kin_explicit_rhs = None
            self._dyn_implicit_rhs = None
            self._dyn_implicit_mat = None
            self._comb_implicit_rhs = right_hand_side
            self._comb_implicit_mat = mass_matrix
            self._comb_explicit_rhs = None
        else:
            self._kin_explicit_rhs = None
            self._dyn_implicit_rhs = None
            self._dyn_implicit_mat = None
            self._comb_implicit_rhs = None
            self._comb_implicit_mat = None
            self._comb_explicit_rhs = right_hand_side

        # Set the remainder of the inputs as instance attributes
        if alg_con is not None and coordinate_derivatives is not None:
            alg_con = [i + len(coordinate_derivatives) for i in alg_con]
        self._alg_con = alg_con
        self.output_eqns = output_eqns

        # Change the body and loads iterables to tuples if they are not tuples
        # already
        if not isinstance(bodies, tuple) and bodies is not None:
            bodies = tuple(bodies)
        if not isinstance(loads, tuple) and loads is not None:
            loads = tuple(loads)
        self._bodies = bodies
        self._loads = loads

    @property
    def coordinates(self):
        """Returns the column matrix of the generalized coordinates"""
        if self._coordinates is None:
            raise AttributeError("The coordinates were not specified.")
        else:
            return self._coordinates

    @property
    def speeds(self):
        """Returns the column matrix of generalized speeds"""
        if self._speeds is None:
            raise AttributeError("The speeds were not specified.")
        else:
            return self._speeds

    @property
    def states(self):
        """Returns the column matrix of the state variables"""
        return self._states

    @property
    def alg_con(self):
        """Returns a list with the indices of the rows containing algebraic
        constraints in the combined form of the equations of motion"""
        return self._alg_con

    @property
    def dyn_implicit_mat(self):
        """Returns the matrix, M, corresponding to the dynamic equations in
        implicit form, M x' = F, where the kinematical equations are not
        included"""
        if self._dyn_implicit_mat is None:
            raise AttributeError("dyn_implicit_mat is not specified for "
                                 "equations of motion form [1] or [2].")
        else:
            return self._dyn_implicit_mat

    @property
    def dyn_implicit_rhs(self):
        """Returns the column matrix, F, corresponding to the dynamic equations
        in implicit form, M x' = F, where the kinematical equations are not
        included"""
        if self._dyn_implicit_rhs is None:
            raise AttributeError("dyn_implicit_rhs is not specified for "
                                 "equations of motion form [1] or [2].")
        else:
            return self._dyn_implicit_rhs

    @property
    def comb_implicit_mat(self):
        """Returns the matrix, M, corresponding to the equations of motion in
        implicit form (form [2]), M x' = F, where the kinematical equations are
        included"""
        if self._comb_implicit_mat is None:
            if self._dyn_implicit_mat is not None:
                num_kin_eqns = len(self._kin_explicit_rhs)
                num_dyn_eqns = len(self._dyn_implicit_rhs)
                zeros1 = zeros(num_kin_eqns, num_dyn_eqns)
                zeros2 = zeros(num_dyn_eqns, num_kin_eqns)
                inter1 = eye(num_kin_eqns).row_join(zeros1)
                inter2 = zeros2.row_join(self._dyn_implicit_mat)
                self._comb_implicit_mat = inter1.col_join(inter2)
                return self._comb_implicit_mat
            else:
                raise AttributeError("comb_implicit_mat is not specified for "
                                     "equations of motion form [1].")
        else:
            return self._comb_implicit_mat

    @property
    def comb_implicit_rhs(self):
        """Returns the column matrix, F, corresponding to the equations of
        motion in implicit form (form [2]), M x' = F, where the kinematical
        equations are included"""
        if self._comb_implicit_rhs is None:
            if self._dyn_implicit_rhs is not None:
                kin_inter = self._kin_explicit_rhs
                dyn_inter = self._dyn_implicit_rhs
                self._comb_implicit_rhs = kin_inter.col_join(dyn_inter)
                return self._comb_implicit_rhs
            else:
                raise AttributeError("comb_implicit_mat is not specified for "
                                     "equations of motion in form [1].")
        else:
            return self._comb_implicit_rhs

    def compute_explicit_form(self):
        """If the explicit right hand side of the combined equations of motion
        is to provided upon initialization, this method will calculate it. This
        calculation can potentially take awhile to compute."""
        if self._comb_explicit_rhs is not None:
            raise AttributeError("comb_explicit_rhs is already formed.")

        inter1 = getattr(self, 'kin_explicit_rhs', None)
        if inter1 is not None:
            inter2 = self._dyn_implicit_mat.LUsolve(self._dyn_implicit_rhs)
            out = inter1.col_join(inter2)
        else:
            out = self._comb_implicit_mat.LUsolve(self._comb_implicit_rhs)

        self._comb_explicit_rhs = out

    @property
    def comb_explicit_rhs(self):
        """Returns the right hand side of the equations of motion in explicit
        form, x' = F, where the kinematical equations are included"""
        if self._comb_explicit_rhs is None:
            raise AttributeError("Please run .combute_explicit_form before "
                                 "attempting to access comb_explicit_rhs.")
        else:
            return self._comb_explicit_rhs

    @property
    def kin_explicit_rhs(self):
        """Returns the right hand side of the kinematical equations in explicit
        form, q' = G"""
        if self._kin_explicit_rhs is None:
            raise AttributeError("kin_explicit_rhs is not specified for "
                                 "equations of motion form [1] or [2].")
        else:
            return self._kin_explicit_rhs

    def dynamic_symbols(self):
        """Returns a column matrix containing all of the symbols in the system
        that depend on time"""
        # Create a list of all of the expressions in the equations of motion
        if self._comb_explicit_rhs is None:
            eom_expressions = (self.comb_implicit_mat[:] +
                               self.comb_implicit_rhs[:])
        else:
            eom_expressions = (self._comb_explicit_rhs[:])

        functions_of_time = set()
        for expr in eom_expressions:
            functions_of_time = functions_of_time.union(
                find_dynamicsymbols(expr))
        functions_of_time = functions_of_time.union(self._states)

        return tuple(functions_of_time)

    def constant_symbols(self):
        """Returns a column matrix containing all of the symbols in the system
        that do not depend on time"""
        # Create a list of all of the expressions in the equations of motion
        if self._comb_explicit_rhs is None:
            eom_expressions = (self.comb_implicit_mat[:] +
                               self.comb_implicit_rhs[:])
        else:
            eom_expressions = (self._comb_explicit_rhs[:])

        constants = set()
        for expr in eom_expressions:
            constants = constants.union(expr.free_symbols)
        constants.remove(dynamicsymbols._t)

        return tuple(constants)

    @property
    def bodies(self):
        """Returns the bodies in the system"""
        if self._bodies is None:
            raise AttributeError("bodies were not specified for the system.")
        else:
            return self._bodies

    @property
    def loads(self):
        """Returns the loads in the system"""
        if self._loads is None:
            raise AttributeError("loads were not specified for the system.")
        else:
            return self._loads

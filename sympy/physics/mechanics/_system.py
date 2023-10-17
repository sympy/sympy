"""Implementation of system for book-keeping all objects defining a model.

Notes
=====

This module is experimental and so is named with a leading underscore to
indicate that the API is not yet stabilized and could be subject to breaking
changes.

"""

from functools import wraps

from sympy import Basic
from sympy import ImmutableMatrix as Matrix
from sympy.core.containers import OrderedSet
from sympy.physics.mechanics.actuator import ActuatorBase
from sympy.physics.mechanics.body_base import BodyBase
from sympy.physics.mechanics.functions import Lagrangian, _validate_coordinates
from sympy.physics.mechanics.joint import Joint
from sympy.physics.mechanics.kane import KanesMethod
from sympy.physics.mechanics.lagrange import LagrangesMethod
from sympy.physics.mechanics.loads import _parse_load, gravity
from sympy.physics.mechanics.method import _Methods
from sympy.physics.mechanics.particle import Particle
from sympy.physics.vector import Point, ReferenceFrame, dynamicsymbols
from sympy.utilities.iterables import iterable

__all__ = ['System']


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
        Matrix of all the generalized coordinates, i.e. the independent
        generalized coordinates stacked upon the dependent.
    u : Matrix
        Matrix of all the generalized speeds, i.e. the independent generealized
        speeds stacked upon the dependent.
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
    actuators : tuple of ActuatorBase subclasses
        Tuple of all actuators present in the system.
    holonomic_constraints : Matrix
        Matrix with the holonomic constraints as rows.
    nonholonomic_constraints : Matrix
        Matrix with the nonholonomic constraints as rows.
    eom_method : subclass of KanesMethod or LagrangesMethod
        Backend for forming the equations of motion.

    Examples
    ========

    In the example below a cart with a pendulum is created. The cart moves along
    the x axis of the rail and the pendulum rotates about the z axis. The length
    of the pendulum is ``l`` with the pendulum represented as a particle. To
    move the cart a time dependent force ``F`` is applied to the cart.

    As the ``_system.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``System`` must therefore be
    imported directly from the ``sympy.physics.mechanics._system`` module.

    >>> from sympy.physics.mechanics._system import System

    We first need to import some functions and create some of our variables.

    >>> from sympy import symbols, simplify
    >>> from sympy.physics.mechanics import (
    ...     mechanics_printing, dynamicsymbols, RigidBody, Particle,
    ...     ReferenceFrame, PrismaticJoint, PinJoint)
    >>> mechanics_printing(pretty_print=False)
    >>> g, l = symbols('g l')
    >>> F = dynamicsymbols('F')

    The next step is to create bodies. It is also useful to create a frame for
    locating the particle with respect to the pin joint later on, as a particle
    does not have a body-fixed frame.

    >>> rail = RigidBody('rail')
    >>> cart = RigidBody('cart')
    >>> bob = Particle('bob')
    >>> bob_frame = ReferenceFrame('bob_frame')

    Initialize the system, with the rail as the Newtonian reference. The body is
    also automatically added to the system.

    >>> system = System.from_newtonian(rail)
    >>> print(system.bodies[0])
    rail

    Create the joints, while immediately also adding them to the system.

    >>> system.add_joints(
    ...     PrismaticJoint('slider', rail, cart, joint_axis=rail.x),
    ...     PinJoint('pin', cart, bob, joint_axis=cart.z,
    ...              child_interframe=bob_frame,
    ...              child_point=l * bob_frame.y)
    ... )
    >>> system.joints
    (PrismaticJoint: slider  parent: rail  child: cart,
    PinJoint: pin  parent: cart  child: bob)

    While adding the joints, the associated generalized coordinates, generalized
    speeds, kinematic differential equations and bodies are also added to the
    system.

    >>> system.q
    Matrix([
    [q_slider],
    [   q_pin]])
    >>> system.u
    Matrix([
    [u_slider],
    [   u_pin]])
    >>> system.kdes
    Matrix([
    [u_slider - q_slider'],
    [      u_pin - q_pin']])
    >>> [body.name for body in system.bodies]
    ['rail', 'cart', 'bob']

    With the kinematics established, we can now apply gravity and the cart force
    ``F``.

    >>> system.apply_gravity(-g * system.y)
    >>> system.add_loads((cart.masscenter, F * rail.x))
    >>> system.loads
    ((rail_masscenter, - g*rail_mass*rail_frame.y),
     (cart_masscenter, - cart_mass*g*rail_frame.y),
     (bob_masscenter, - bob_mass*g*rail_frame.y),
     (cart_masscenter, F*rail_frame.x))

    With the entire system defined, we can now form the equations of motion.
    Before forming the equations of motion, one can also run some checks that
    will try to identify some common errors.

    >>> system.validate_system()
    >>> system.form_eoms()
    Matrix([
    [bob_mass*l*u_pin**2*sin(q_pin) - bob_mass*l*cos(q_pin)*u_pin'
     - (bob_mass + cart_mass)*u_slider' + F],
    [                   -bob_mass*g*l*sin(q_pin) - bob_mass*l**2*u_pin'
     - bob_mass*l*cos(q_pin)*u_slider']])
    >>> simplify(system.mass_matrix)
    Matrix([
    [ bob_mass + cart_mass, bob_mass*l*cos(q_pin)],
    [bob_mass*l*cos(q_pin),         bob_mass*l**2]])
    >>> system.forcing
    Matrix([
    [bob_mass*l*u_pin**2*sin(q_pin) + F],
    [          -bob_mass*g*l*sin(q_pin)]])

    The complexity of the above example can be increased if we add a constraint
    to prevent the particle from moving in the horizontal (x) direction. This
    can be done by adding a holonomic constraint. After which we should also
    redefine what our (in)dependent generalized coordinates and speeds are. Note
    that the backend for forming the equations of motion is reset automatically
    because the system properties are changed in this process.

    >>> type(system.eom_method)
    <class 'sympy.physics.mechanics.kane.KanesMethod'>
    >>> system.add_holonomic_constraints(
    ...     bob.masscenter.pos_from(rail.masscenter).dot(system.x)
    ... )
    >>> system.q_ind = system.get_joint('pin').coordinates
    >>> system.q_dep = system.get_joint('slider').coordinates
    >>> system.u_ind = system.get_joint('pin').speeds
    >>> system.u_dep = system.get_joint('slider').speeds
    >>> type(system.eom_method)
    <class 'NoneType'>

    With the updated system the equations of motion can be formed again.

    >>> system.validate_system()
    >>> system.form_eoms()
    Matrix([[-bob_mass*g*l*sin(q_pin)
             - bob_mass*l**2*u_pin'
             - bob_mass*l*cos(q_pin)*u_slider'
             - l*(bob_mass*l*u_pin**2*sin(q_pin)
             - bob_mass*l*cos(q_pin)*u_pin'
             - (bob_mass + cart_mass)*u_slider')*cos(q_pin)
             - l*F*cos(q_pin)]])
    >>> simplify(system.mass_matrix)
    Matrix([
    [bob_mass*l**2*sin(q_pin)**2, -cart_mass*l*cos(q_pin)],
    [               l*cos(q_pin),                       1]])
    >>> simplify(system.forcing)
    Matrix([
    [-l*(bob_mass*g*sin(q_pin) + bob_mass*l*u_pin**2*sin(2*q_pin)/2
     + F*cos(q_pin))],
    [
    l*u_pin**2*sin(q_pin)]])

    """

    def __init__(self, frame=None, origin=None):
        """Initialize the system.

        Parameters
        ==========

        frame : ReferenceFrame, optional
            The inertial frame of the system. If none is supplied, a new frame
            will be created.
        origin : Point, optional
            The origin of the system. If none is supplied, a new origin will be
            created.

        """
        if frame is None:
            frame = ReferenceFrame('inertial_frame')
        elif not isinstance(frame, ReferenceFrame):
            raise TypeError('Frame must be an instance of ReferenceFrame.')
        self._frame = frame
        if origin is None:
            origin = Point('inertial_origin')
        elif not isinstance(origin, Point):
            raise TypeError('Origin must be an instance of Point.')
        self._origin = origin
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
        self._actuators = []
        self._eom_method = None

    @classmethod
    def from_newtonian(cls, newtonian):
        """Constructs the system with respect to a Newtonian body."""
        if isinstance(newtonian, Particle):
            raise TypeError('A Particle has no frame so cannot act as '
                            'the Newtonian.')
        system = cls(frame=newtonian.frame, origin=newtonian.masscenter)
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
        bodies = self._objects_to_list(bodies)
        self._check_objects(bodies, [], BodyBase, 'Bodies', 'bodies')
        self._bodies = bodies

    @property
    def joints(self):
        """Tuple of all joints that have been added to the system."""
        return tuple(self._joints)

    @joints.setter
    @_reset_eom_method
    def joints(self, joints):
        joints = self._objects_to_list(joints)
        self._check_objects(joints, [], Joint, 'Joints', 'joints')
        self._joints = []
        self.add_joints(*joints)

    @property
    def loads(self):
        """Tuple of loads that have been applied on the system."""
        return tuple(self._loads)

    @loads.setter
    @_reset_eom_method
    def loads(self, loads):
        loads = self._objects_to_list(loads)
        self._loads = [_parse_load(load) for load in loads]

    @property
    def actuators(self):
        """Tuple of actuators present in the system."""
        return tuple(self._actuators)

    @actuators.setter
    @_reset_eom_method
    def actuators(self, actuators):
        actuators = self._objects_to_list(actuators)
        self._check_objects(actuators, [], ActuatorBase, 'Actuators',
                            'actuators')
        self._actuators = actuators

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
        kdes = self._objects_to_list(kdes)
        self._kdes = self._parse_expressions(
            kdes, [], 'kinematic differential equations')

    @property
    def holonomic_constraints(self):
        """Matrix with the holonomic constraints as rows."""
        return self._hol_coneqs

    @holonomic_constraints.setter
    @_reset_eom_method
    def holonomic_constraints(self, constraints):
        constraints = self._objects_to_list(constraints)
        self._hol_coneqs = self._parse_expressions(
            constraints, [], 'holonomic constraints')

    @property
    def nonholonomic_constraints(self):
        """Matrix with the nonholonomic constraints as rows."""
        return self._nonhol_coneqs

    @nonholonomic_constraints.setter
    @_reset_eom_method
    def nonholonomic_constraints(self, constraints):
        constraints = self._objects_to_list(constraints)
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
        """Helper to check the objects that are being added to the system.

        Explanation
        ===========
        This method checks that the objects that are being added to the system
        are of the correct type and have not already been added. If any of the
        objects are not of the correct type or have already been added, then
        an error is raised.

        Parameters
        ==========
        objects : iterable
            The objects that would be added to the system.
        obj_lst : list
            The list of objects that are already in the system.
        expected_type : type
            The type that the objects should be.
        obj_name : str
            The name of the category of objects. This string is used to
            formulate the error message for the user.
        type_name : str
            The name of the type that the objects should be. This string is used
            to formulate the error message for the user.

        """
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
        coords_ind, coords_dep = old_coords_ind[:], old_coords_dep[:]
        if not iterable(independent):
            independent = [independent] * len(new_coords)
        for coord, indep in zip(new_coords, independent):
            if indep:
                coords_ind.append(coord)
            else:
                coords_dep.append(coord)
        # Check types and duplicates
        if is_coordinates:
            _validate_coordinates(coords_ind + coords_dep,
                                  self.u_ind[:] + self.u_dep[:])
        else:
            _validate_coordinates(self.q_ind[:] + self.q_dep[:],
                                  coords_ind + coords_dep)
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
            One or more holonomic constraints, which are expressions that should
            be zero.

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
            One or more nonholonomic constraints, which are expressions that
            should be zero.

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
    def apply_gravity(self, acceleration):
        """Apply gravity to all bodies in the system.

        Parameters
        ==========

        acceleration : Vector
            The acceleration due to gravity.

        """
        self.add_loads(*gravity(acceleration, *self.bodies))

    @_reset_eom_method
    def add_actuators(self, *actuators):
        """Add actuator(s) to the system.

        Parameters
        ==========

        *actuators : subclass of ActuatorBase
            One or more actuators.

        """
        self._check_objects(actuators, self.actuators, ActuatorBase,
                            'Actuators', 'actuators')
        self._actuators.extend(actuators)

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
        however always added to the system, so you should not also manually add
        those on beforehand.

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

    def get_joint(self, name):
        """Retrieve a joint from the system by name.

        Parameters
        ==========

        name : str
            The name of the joint to retrieve.

        Returns
        =======

        subclass of Joint
            The joint with the given name, or None if no such joint exists.

        """
        for joint in self._joints:
            if joint.name == name:
                return joint

    def _form_eoms(self):
        return self.form_eoms()

    def form_eoms(self, eom_method=KanesMethod, **kwargs):
        """Form the equations of motion of the system.

        Parameters
        ==========

        eom_method : subclass of KanesMethod or LagrangesMethod
            Backend class to be used for forming the equations of motion. The
            default is ``KanesMethod``.

        Returns
        ========

        Matrix
            Vector of equations of motions.

        Examples
        ========

        As the ``_system.py`` module is experimental, it is not yet part of the
        ``sympy.physics.mechanics`` namespace. ``System`` must therefore be
        imported directly from the ``sympy.physics.mechanics._system`` module.

        >>> from sympy.physics.mechanics._system import System

        This is a simple example for a one degree of freedom translational
        spring-mass-damper.

        >>> from sympy import S, symbols
        >>> from sympy.physics.mechanics import (
        ...     LagrangesMethod, dynamicsymbols, PrismaticJoint, Particle,
        ...     RigidBody)
        >>> q = dynamicsymbols('q')
        >>> qd = dynamicsymbols('q', 1)
        >>> m, k, b = symbols('m k b')
        >>> wall = RigidBody('W')
        >>> system = System.from_newtonian(wall)
        >>> bob = Particle('P', mass=m)
        >>> bob.potential_energy = S.Half * k * q**2
        >>> system.add_joints(PrismaticJoint('J', wall, bob, q, qd))
        >>> system.add_loads((bob.masscenter, b * qd * system.x))
        >>> system.form_eoms(LagrangesMethod)
        Matrix([[-b*Derivative(q(t), t) + k*q(t) + m*Derivative(q(t), (t, 2))]])

        We can also solve for the states using the 'rhs' method.

        >>> system.rhs()
        Matrix([
        [               Derivative(q(t), t)],
        [(b*Derivative(q(t), t) - k*q(t))/m]])

        """
        # KanesMethod does not accept empty iterables
        loads = self.loads + tuple(
            load for act in self.actuators for load in act.to_loads())
        loads = loads if loads else None
        if issubclass(eom_method, KanesMethod):
            disallowed_kwargs = {
                "frame", "q_ind", "u_ind", "kd_eqs", "q_dependent",
                "u_dependent", "configuration_constraints",
                "velocity_constraints", "forcelist", "bodies"}
            wrong_kwargs = disallowed_kwargs.intersection(kwargs)
            if wrong_kwargs:
                raise ValueError(
                    f"The following keyword arguments are not allowed to be "
                    f"overwritten in {eom_method.__name__}: {wrong_kwargs}.")
            velocity_constraints = self.holonomic_constraints.diff(
                dynamicsymbols._t).col_join(self.nonholonomic_constraints)
            kwargs = {"frame": self.frame, "q_ind": self.q_ind,
                      "u_ind": self.u_ind, "kd_eqs": self.kdes,
                      "q_dependent": self.q_dep, "u_dependent": self.u_dep,
                      "configuration_constraints": self.holonomic_constraints,
                      "velocity_constraints": velocity_constraints,
                      "forcelist": loads, "bodies": self.bodies,
                      "explicit_kinematics": False, **kwargs}
            self._eom_method = eom_method(**kwargs)
        elif issubclass(eom_method, LagrangesMethod):
            disallowed_kwargs = {
                "frame", "qs", "forcelist", "bodies", "hol_coneqs",
                "nonhol_coneqs", "Lagrangian"}
            wrong_kwargs = disallowed_kwargs.intersection(kwargs)
            if wrong_kwargs:
                raise ValueError(
                    f"The following keyword arguments are not allowed to be "
                    f"overwritten in {eom_method.__name__}: {wrong_kwargs}.")
            kwargs = {"frame": self.frame, "qs": self.q, "forcelist": loads,
                      "bodies": self.bodies,
                      "hol_coneqs": self.holonomic_constraints,
                      "nonhol_coneqs": self.nonholonomic_constraints, **kwargs}
            if "Lagrangian" not in kwargs:
                kwargs["Lagrangian"] = Lagrangian(kwargs["frame"],
                                                  *kwargs["bodies"])
            self._eom_method = eom_method(**kwargs)
        else:
            raise NotImplementedError(f'{eom_method} has not been implemented.')
        return self.eom_method._form_eoms()

    def rhs(self, inv_method=None):
        """Compute the equations of motion in the explicit form.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrices.MatrixBase.inv`

        Returns
        ========

        Matrix
            Equations of motion in the explicit form.

        See Also
        ========

        sympy.physics.mechanics.kane.KanesMethod.rhs:
            KanesMethod's ``rhs`` function.
        sympy.physics.mechanics.lagrange.LagrangesMethod.rhs:
            LagrangesMethod's ``rhs`` function.

        """
        return self.eom_method.rhs(inv_method=inv_method)

    @property
    def mass_matrix(self):
        r"""The mass matrix of the system.

        Explanation
        ===========

        The mass matrix $M_d$ and the forcing vector $f_d$ of a system describe
        the system's dynamics according to the following equations:

        .. math::
            M_d \dot{u} = f_d

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

        .. math::
            M_m \dot{x} = f_m

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

    def validate_system(self, eom_method=KanesMethod, check_duplicates=False):
        """Validates the system using some basic checks.

        Explanation
        ===========

        This method validates the system based on the following checks:

        - The number of dependent generalized coordinates should equal the
          number of holonomic constraints.
        - All generalized coordinates defined by the joints should also be known
          to the system.
        - If ``KanesMethod`` is used as a ``eom_method``:
            - All generalized speeds and kinematic differential equations
              defined by the joints should also be known to the system.
            - The number of dependent generalized speeds should equal the number
              of velocity constraints.
            - The number of generalized coordinates should be less than or equal
              to the number of generalized speeds.
            - The number of generalized coordinates should equal the number of
              kinematic differential equations.
        - If ``LagrangesMethod`` is used as ``eom_method``:
            - There should not be any generalized speeds that are not
              derivatives of the generalized coordinates (this includes the
              generalized speeds defined by the joints).

        Parameters
        ==========

        eom_method : subclass of KanesMethod or LagrangesMethod
            Backend class that will be used for forming the equations of motion.
            There are different checks for the different backends. The default
            is ``KanesMethod``.
        check_duplicates : bool
            Boolean whether the system should be checked for duplicate
            definitions. The default is False, because duplicates are already
            checked when adding objects to the system.

        Notes
        =====

        This method is not guaranteed to be backwards compatible as it may
        improve over time. The method can become both more and less strict in
        certain areas. However a well-defined system should always pass all
        these tests.

        """
        msgs = []
        # Save some data in variables
        n_hc = self.holonomic_constraints.shape[0]
        n_nhc = self.nonholonomic_constraints.shape[0]
        n_q_dep, n_u_dep = self.q_dep.shape[0], self.u_dep.shape[0]
        q_set, u_set = set(self.q), set(self.u)
        n_q, n_u = len(q_set), len(u_set)
        # Check number of holonomic constraints
        if n_q_dep != n_hc:
            msgs.append(f'The number of dependent generalized coordinates '
                        f'{n_q_dep} should be equal to the number of holonomic '
                        f'constraints {n_hc}.')
        # Check if all joint coordinates and speeds are present
        missing_q = set()
        for joint in self.joints:
            missing_q.update(set(joint.coordinates).difference(q_set))
        if missing_q:
            msgs.append(f'The generalized coordinates {missing_q} used in '
                        f'joints are not added to the system.')
        # Method dependent checks
        if issubclass(eom_method, KanesMethod):
            n_kdes = len(self.kdes)
            missing_kdes, missing_u = set(), set()
            for joint in self.joints:
                missing_u.update(set(joint.speeds).difference(u_set))
                missing_kdes.update(set(joint.kdes).difference(
                    self.kdes[:] + (-self.kdes)[:]))
            if missing_u:
                msgs.append(f'The generalized speeds {missing_u} used in '
                            f'joints are not added to the system.')
            if missing_kdes:
                msgs.append(f'The kinematic differential equations '
                            f'{missing_kdes} used in joints are not added to '
                            f'the system.')
            if n_u_dep != n_hc + n_nhc:
                msgs.append(f'The number of dependent generalized speeds '
                            f'{n_u_dep} should be equal to the number of '
                            f'velocity constraints {n_hc + n_nhc}.')
            if n_q > n_u:
                msgs.append(f'The number of generalized coordinates {n_q} '
                            f'should be less than or equal to the number of '
                            f'generalized speeds {n_u}.')
            if n_u != n_kdes:
                msgs.append(f'The number of generalized speeds {n_u} should be '
                            f'equal to the number of kinematic differential '
                            f'equations {n_kdes}.')
        elif issubclass(eom_method, LagrangesMethod):
            not_qdots = set(self.u).difference(self.q.diff(dynamicsymbols._t))
            for joint in self.joints:
                not_qdots.update(set(
                    joint.speeds).difference(self.q.diff(dynamicsymbols._t)))
            if not_qdots:
                msgs.append(f'The generalized speeds {not_qdots} are not '
                            f'supported by this method. Only derivatives of the'
                            f' generalized coordinates are supported. If these '
                            f'symbols are used in your expressions, then this '
                            f'will result in wrong equations of motion.')
        else:
            raise NotImplementedError(f'{eom_method} has not been implemented.')
        if check_duplicates:  # Should be redundant
            duplicates_to_check = [('generalized coordinates', self.q),
                                   ('generalized speeds', self.u),
                                   ('bodies', self.bodies),
                                   ('joints', self.joints)]
            for name, lst in duplicates_to_check:
                seen = set()
                duplicates = {x for x in lst if x in seen or seen.add(x)}
                if duplicates:
                    msgs.append(f'The {name} {duplicates} exist multiple times '
                                f'within the system.')
        if msgs:
            raise ValueError('\n'.join(msgs))

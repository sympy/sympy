from sympy import Symbol
from sympy.physics.mechanics import (RigidBody, Particle, ReferenceFrame,
                                     inertia)
from sympy.physics.vector import Point, Vector

__all__ = ['Body']


class Body(RigidBody, Particle):
    """
    Body is a common representation of RigidBody or a Particle.

    A Body represents either a rigid body or particle in classical mechanics.
    Bodies have a body-fixed reference frame, a mass, a mass center and
    possibly a body-fixed inertia.

    Parameters
    ----------
    name: String
        Defines the name of the body. It is used as the base for defining body
        specific properties.
    masscenter : Point, optional
        A point that represents the center of mass of the body or particle. If no
        point is given, a point is generated.
    frame : ReferenceFrame (optional)
        The ReferenceFrame that represents the reference frame of the body. If
        no frame is given, a frame is generated.
    mass : Sympifyable, optional
        A Sympifyable object which represents the mass of the body. if no mass
        is passed, one is generated.
    body_inertia : Dyadic
        Central inertia dyadic of the body. If none is passed while creating
        RigidBody, a default inertia is generated.

    Examples
    --------
    1. Default behaviour. It creates a RigidBody after defining mass,
     mass center, frame and inertia.

    >>> from sympy.physics.mechanics import Body
    >>> body = Body('name_of_body')

    2. Passing attributes of Rigidbody. All the arguments needed to create a
     RigidBody can be passed while creating a Body too.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, Point, inertia
    >>> from sympy.physics.mechanics import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> frame = ReferenceFrame('frame')
    >>> ixx = Symbol('ixx')
    >>> body_inertia = inertia(frame, ixx, 0, 0)
    >>> body = Body('name_of_body',masscenter,mass,frame,body_inertia)

    3. Creating a Particle. If masscenter and mass are passed, and inertia is
     not then a Particle is created.

    >>> from sympy import Symbol
    >>> from sympy.physics.vector import Point
    >>> from sympy.physics.mechanics import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> body = Body('name_of_body',masscenter,mass)

    Similarly, A frame can also be passed while creating a Particle.

    """
    def __init__(self, name, masscenter=None, mass=None, frame=None,
                 central_inertia=None):

        self.name = name
        self.loads = []

        if frame is None:
            frame = ReferenceFrame(name + '_frame')

        if masscenter is None:
            masscenter = Point(name + '_masscenter')

        if central_inertia is None and mass is None:
            ixx = Symbol(name + '_ixx')
            iyy = Symbol(name + '_iyy')
            izz = Symbol(name + '_izz')
            izx = Symbol(name + '_izx')
            ixy = Symbol(name + '_ixy')
            iyz = Symbol(name + '_iyz')
            _inertia = (inertia(frame, ixx, iyy, izz, ixy, iyz, izx),
                        masscenter)
        else:
            _inertia = (central_inertia, masscenter)

        if mass is None:
            _mass = Symbol(name + '_mass')
        else:
            _mass = mass

        masscenter.set_vel(frame, 0)

        # If user passes masscenter and mass then a particle is created
        # otherwise a rigidbody. As a result a body may or may not have inertia.
        if central_inertia is None and mass is not None:
            self.frame = frame
            self.masscenter = masscenter
            Particle.__init__(self, name, masscenter, _mass)
        else:
            RigidBody.__init__(self, name, masscenter, frame, _mass, _inertia)

    def apply_force(self, vec, point=None):
        """
        Adds the force to the point (center of mass by default) on the body.

        Parameters
        ----------
        vec: Vector
            Defines the force vector. Can be any vector w.r.t any frame or
            combinations of frame.
        point: Point, optional
            Defines the point on which the force must be applied. Default is
            Body's center of mass.

        Example
        -------
        To apply a unit force in x direction of body's frame to body's
        center of mass.

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Body
        >>> body = Body('body')
        >>> g = Symbol('g')
        >>> body.apply_force(body.mass * g * body.frame.x)

        To apply force to any other point than center of mass, pass that point
        as well.

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Body
        >>> parent = Body('parent')
        >>> child = Body('child')
        >>> g = Symbol('g')
        >>> frame = parent.frame
        >>> l = Symbol('l')
        >>> point = child.masscenter.locatenew('force_point', l * body.frame.y)
        >>> gravity = child.mass * g
        >>> body.apply_force(gravity * body.frame.x, point)

        """
        if not isinstance(point, Point):
            if point is None:
                point = self.masscenter  # masscenter
            else:
                raise TypeError("A Point must be supplied to apply force to.")
        if not isinstance(vec, Vector):
            raise TypeError("A Vector must be supplied to apply force.")

        self.loads.append((point, vec))

    def apply_torque(self, vec):
        """
        Adds torque to the body.

        Parameters
        ----------
        vec: Vector
            Defines the force vector. Can be any vector w.r.t any frame or
            combinations of frame.
        """
        if not isinstance(vec, Vector):
            raise TypeError("A Vector must be supplied to add torque.")
        self.loads.append((self.frame, vec))

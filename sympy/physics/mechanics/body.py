from sympy import Symbol
from sympy.physics.mechanics import RigidBody, Particle, ReferenceFrame, \
    inertia
from sympy.physics.vector import Point, Vector

__all__ = ['Body']


class Body(RigidBody, Particle):
    """
    A Body which can be connected by joints.

    Instances of Body can be connected using specific joints and can further
    be used to generate equations of motion. It can be seen as an extension
    of already present RigidBody and Particle classes and can inherit from
    either of them based on the arguments passed.

    Parameters
    ---------
    name: string
        defines the name of the body. It is used as the base for defining body
        specific properties.
    masscenter : Point (optional)
        The point which represents the center of mass of the rigid body.
    frame : ReferenceFrame (optional)
        The ReferenceFrame which the rigid body is fixed in.
    mass : Sympifyable (optional)
        The body's mass.
    body_inertia : Dyadic (instance of inertia)
        The body's inertia about center of mass.

    Example:
    --------
    1. Default behaviour. It creates a RigidBody after defining mass,
     masscenter, frame and inertia.

    >>> from pydy.bodies import Body
    >>> body = Body('name_of_body')

    2. Passing attributes of Rigidbody. All the arguments needed to create a
     RigidBody can be passed while creating a Body too.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, Point, inertia
    >>> from sympy.physics.mechanics import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> frame = ReferenceFrame('frame')
    >>> body_inertia = inertia(frame, 1, 0, 0)
    >>> body = Body('name_of_body', masscenter, mass, frame, body_inertia)

    3. Creating a Particle. If masscenter and mass are passed, and inertia is
     not then a Particle is created.

    >>> from sympy import Symbol
    >>> from sympy import Point
    >>> from sympy.physics.mechanics import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> body = Body('name_of_body', masscenter, mass)

    Similarly, A frame can also be passed while creating a Particle.

    """
    def __init__(self, name, masscenter=None, mass=None, frame=None,
                 body_inertia=None):

        self.parent = None
        self.child = None
        self.force_list = []
        self._counter = 0

        if masscenter is None:
            self._masscenter = Point(name + '_masscenter')
        else:
            self._masscenter = masscenter

        if mass is None:
            _mass = Symbol(name + '_mass')
        else:
            _mass = mass

        if frame is None:
            self._frame = ReferenceFrame(name + '_frame')
        else:
            self._frame = frame

        if body_inertia is None and mass is None:
            _inertia = (inertia(self._frame, 1, 1, 1), self._masscenter)
        else:
            _inertia = (body_inertia, self._masscenter)

        self._masscenter.set_vel(self._frame, 0)

        # If user passes masscenter and mass then a particle is created
        # otherwise a rigidbody. As a result a body may or may not have inertia.
        if body_inertia is None and mass is not None:
            Particle.__init__(self, name, self._masscenter, _mass)
        else:
            RigidBody.__init__(self, name, self._masscenter, self._frame, _mass, _inertia)

    def add_force(self, force_vector, point=None, frame=None):
        """
        Adds force to the body by adding a force vector in form of a tuple
        of length 3 and a point on which the force will be applied.

        Parameters
        ---------
        force_vector: A 3-Tuple
            Defines the force vector w.r.t frame passed (or Body's frame). The
            tuple defines the values of vector in x, y and z directions of the
            frame.
        point: Point
            Defines the point on which the force must be applied.
        frame: ReferenceFrame
            Defines the frame w.r.t which force vector must be defined.

        Example
        --------
        To add a force of magnitude 1 in y direction on a point at distance of
        1 in x direction. All the directions are w.r.t body's frame.

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Body
        >>> body = Body('body')
        >>> g = Symbol('g')
        >>> body.add_force((body.get_mass() * g, 0, 0), body.get_masscenter())

        To define the force_vector w.r.t some other frame, pass the frame as well.

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Body
        >>> parent = Body('parent')
        >>> child = Body('child')
        >>> g = Symbol('g')
        >>> frame = parent.get_frame()
        >>> masscenter = child.get_masscenter()
        >>> gravity = child.get_mass() * g
        >>> body.add_force((gravity, 0, 0), masscenter, frame)

        """
        if point is None:
            point = self._masscenter  # masscenter

        if frame is None:
            frame = self._frame

        if not isinstance(force_vector, tuple):
            raise TypeError("Force vector must be a tuple of length 3")
        else:
            force_vector = self._convert_tuple_to_vector(force_vector, frame)
        self.force_list.append((point, force_vector))
        self._counter += 1

    def _convert_tuple_to_vector(self, pos_tuple, frame):
        """
        converts 3-Tuple into a vector in given frame. Values of the Tuple
        corresponds with x, y and z axis of the frame.
        """
        if len(pos_tuple) != 3:
            raise TypeError('position tuple must be of length 3')
        else:
            unit_vectors = [frame.x, frame.y, frame.z]
            vector = Vector(0)
            for i in range(3):
                vector += pos_tuple[i] * unit_vectors[i]
            return vector

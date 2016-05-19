from sympy.core.basic import Basic
from sympy.vector.scalar import BaseScalar
from sympy import eye, trigsimp, ImmutableMatrix as Matrix, Symbol
from sympy.core.compatibility import string_types, range
from sympy.core.cache import cacheit
from sympy.vector.orienters import (Orienter, AxisOrienter, BodyOrienter,
                                    SpaceOrienter, QuaternionOrienter)
from sympy.vector.coordsysrect import CoordSysCartesian
import sympy.vector


class CoordSysSpherical(Basic):
    """
    Represents a coordinate system in 3-D space.
    """

    def __new__(cls, name, location=None, rotation_matrix=None,
                parent=None, vector_names=None, variable_names=None):
        """
        The orientation/location parameters are necessary if this system
        is being defined at a certain orientation or location wrt another.

        Parameters
        ==========

        name : str
            The name of the new CoordSysSpherical instance.

        location : Vector
            The position vector of the new system's origin wrt the parent
            instance.

        rotation_matrix : SymPy ImmutableMatrix
            The rotation matrix of the new coordinate system with respect
            to the parent. In other words, the output of
            new_system.rotation_matrix(parent).

        parent : CoordSysSpherical
            The coordinate system wrt which the orientation/location
            (or both) is being defined.

        vector_names, variable_names : iterable(optional)
            Iterables of 3 strings each, with custom names for base
            vectors and base scalars of the new system respectively.
            Used for simple str printing.

        """

        name = str(name)
        Vector = sympy.vector.Vector
        BaseVector = sympy.vector.BaseVector
        Point = sympy.vector.Point
        if not isinstance(name, string_types):
            raise TypeError("name should be a string")

        # If orientation information has been provided, store
        # the rotation matrix accordingly
        if rotation_matrix is None:
            parent_orient = Matrix(eye(3))
        else:
            if not isinstance(rotation_matrix, Matrix):
                raise TypeError("rotation_matrix should be an Immutable" +
                                "Matrix instance")
            parent_orient = rotation_matrix

        # If location information is not given, adjust the default
        # location as Vector.zero
        if parent is not None:
            if not isinstance(parent, CoordSysSpherical):
                raise TypeError("parent should be a " +
                                "CoordSysSpherical/None")
            if location is None:
                location = Vector.zero
            else:
                if not isinstance(location, Vector):
                    raise TypeError("location should be a Vector")
                # Check that location does not contain base
                # scalars
                for x in location.free_symbols:
                    if isinstance(x, BaseScalar):
                        raise ValueError("location should not contain" +
                                         " BaseScalars")
            origin = parent.origin.locate_new(name + '.origin',
                                              location)
        else:
            location = Vector.zero
            origin = Point(name + '.origin')

        # All systems that are defined as 'roots' are unequal, unless
        # they have the same name.
        # Systems defined at same orientation/position wrt the same
        # 'parent' are equal, irrespective of the name.
        # This is true even if the same orientation is provided via
        # different methods like Axis/Body/Space/Quaternion.
        # However, coincident systems may be seen as unequal if
        # positioned/oriented wrt different parents, even though
        # they may actually be 'coincident' wrt the root system.
        if parent is not None:
            obj = super(CoordSysSpherical, cls).__new__(
                cls, Symbol(name), location, parent_orient, parent)
        else:
            obj = super(CoordSysSpherical, cls).__new__(
                cls, Symbol(name), location, parent_orient)
        obj._name = name

        # Initialize the base vectors
        if vector_names is None:
            vector_names = (name + '.e_r', name + '.e_theta', name + '.e_phi')
            latex_vects = [(r'\mathbf{\hat{e_r}_{%s}}' % name),
                           (r'\mathbf{\hat{e_theta}_{%s}}' % name),
                           (r'\mathbf{\hat{e_phi}_{%s}}' % name)]
            pretty_vects = (name + '_e_r', name + '_e_theta', name + '_e_phi')
        else:
            _check_strings('vector_names', vector_names)
            vector_names = list(vector_names)
            latex_vects = [(r'\mathbf{\hat{%s}_{%s}}' % (x, name)) for
                           x in vector_names]
            pretty_vects = [(name + '_' + x) for x in vector_names]

        obj._e_r = BaseVector(vector_names[0], 0, obj,
                            pretty_vects[0], latex_vects[0])
        obj._e_theta = BaseVector(vector_names[1], 1, obj,
                            pretty_vects[1], latex_vects[1])
        obj._e_phi = BaseVector(vector_names[2], 2, obj,
                            pretty_vects[2], latex_vects[2])

        # Initialize the base scalars
        if variable_names is None:
            variable_names = (name + '.r', name + '.theta', name + '.phi')
            latex_scalars = [(r"\mathbf{{r}_{%s}}" % name),
                             (r"\mathbf{{theta}_{%s}}" % name),
                             (r"\mathbf{{phi}_{%s}}" % name)]
            pretty_scalars = (name + '_r', name + '_theta', name + '_phi')
        else:
            _check_strings('variable_names', vector_names)
            variable_names = list(variable_names)
            latex_scalars = [(r"\mathbf{{%s}_{%s}}" % (x, name)) for
                             x in variable_names]
            pretty_scalars = [(name + '_' + x) for x in variable_names]

        obj._r = BaseScalar(variable_names[0], 0, obj,
                            pretty_scalars[0], latex_scalars[0])
        obj._theta = BaseScalar(variable_names[1], 1, obj,
                            pretty_scalars[1], latex_scalars[1])
        obj._phi = BaseScalar(variable_names[2], 2, obj,
                            pretty_scalars[2], latex_scalars[2])


        # Assign a Del operator instance
        from sympy.vector.deloperator import Del
        obj._delop = Del(obj)

        # Assign params
        obj._parent = parent
        if obj._parent is not None:
            obj._root = obj._parent._root
        else:
            obj._root = obj

        obj._parent_rotation_matrix = parent_orient
        obj._origin = origin

        # Return the instance
        return obj

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__

    def __iter__(self):
        return iter([self.i, self.j, self.k])

    @property
    def origin(self):
        return self._origin

    @property
    def delop(self):
        return self._delop

    @property
    def e_r(self):
        return self._e_r

    @property
    def e_theta(self):
        return self._e_theta

    @property
    def e_phi(self):
        return self._e_phi

    @property
    def r(self):
        return self._r

    @property
    def theta(self):
        return self._theta

    @property
    def phi(self):
        return self._phi

    def base_vectors(self):
        return self._e_r, self._e_theta, self._e_phi

    def base_scalars(self):
        return self._r, self._theta, self._phi

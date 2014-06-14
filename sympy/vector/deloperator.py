from sympy.core import Basic
from sympy import diff
from sympy.vector.vector import Vector, i, j, k
from sympy.vector.scalar import x, y, z


class Del(Basic):
    """
    Represents the vector differential operator, usually represented in
    mathematical expressions as the 'nabla' symbol.
    """

    def __init__(self):
        self._x, self._y, self._z = x, y, z
        self._i, self._j, self._k = i, j, k

    def __call__(self, scalar_field):
        """
        Represents the gradient of the given scalar field.

        Parameters
        ==========

        scalar_field : SymPy expression
            The scalar field to calculate the gradient of.

        Examples
        ========

        >>> from sympy.vector import x, y, z, delop
        >>> delop(x*y*z)
        y*z*i + x*z*j + x*y*k

        """

        vx = diff(scalar_field, self._x)
        vy = diff(scalar_field, self._y)
        vz = diff(scalar_field, self._z)

        return vx*self._i + vy*self._j + vz*self._k

    def dot(self, vect):
        """
        Represents the dot product between this operator and a given
        vector - equal to the divergence of the vector field.

        Parameters
        ==========

        vect : Vector
            The vector whose divergence is to be calculated.

        Examples
        ========

        >>> from sympy.vector import i, j, k, x, y, z, delop
        >>> v = x*y*z * (i + j + k)
        >>> delop & v
        x*y + x*z + y*z
        >>> delop.dot(i)
        0

        """

        vx = diff(vect.dot(self._i), self._x)
        vy = diff(vect.dot(self._j), self._y)
        vz = diff(vect.dot(self._k), self._z)

        return vx + vy + vz

    __and__ = dot

    def cross(self, vect):
        """
        Represents the cross product between this operator and a given
        vector - equal to the curl of the vector field.

        Parameters
        ==========

        vect : Vector
            The vector whose curl is to be calculated.

        Examples
        ========

        >>> from sympy.vector import i, j, k, x, y, z, delop
        >>> v = x*y*z * (i + j + k)
        >>> delop ^ v
        (-x*y + x*z)*i + (x*y - y*z)*j + (-x*z + y*z)*k
        >>> delop.cross(i)
        0

        """

        vectx = vect.dot(self._i)
        vecty = vect.dot(self._j)
        vectz = vect.dot(self._k)
        outvec = Vector.Zero
        outvec += (diff(vectz, self._y) - diff(vecty, self._z)) * self._i
        outvec += (diff(vectx, self._z) - diff(vectz, self._x)) * self._j
        outvec += (diff(vecty, self._x) - diff(vectx, self._y)) * self._k

        return outvec

    __xor__ = cross

delop = Del()

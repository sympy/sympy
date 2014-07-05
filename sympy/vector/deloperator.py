from sympy.core import Basic
from sympy import diff
from sympy.vector.vector import Vector
from sympy.vector.functions import express
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.core import S


class Del(Basic):
    """
    Represents the vector differential operator, usually represented in
    mathematical expressions as the 'nabla' symbol.
    """

    def __new__(cls, system):
        if not isinstance(system, CoordSysCartesian):
            raise TypeError("system should be a CoordSysCartesian")
        obj = super(Del, cls).__new__(cls, system)
        obj._x, obj._y, obj._z = system.x, system.y, system.z
        obj._i, obj._j, obj._k = system.i, system.j, system.k
        obj._system = system
        obj._name = system.__str__() + ".del"
        return obj

    @property
    def system(self):
        return self._system

    def __call__(self, scalar_field):
        """
        Represents the gradient of the given scalar field.

        Parameters
        ==========

        scalar_field : SymPy expression
            The scalar field to calculate the gradient of.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> C.delop(C.x*C.y*C.z)
        C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k

        """

        scalar_field = express(scalar_field, self.system,
                               variables = True)
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

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> v = C.x*C.y*C.z * (C.i + C.j + C.k)
        >>> C.delop & v
        C.x*C.y + C.x*C.z + C.y*C.z
        >>> C.delop.dot(C.i)
        0

        """

        vx = _diff_conditional(vect.dot(self._i), self._x)
        vy = _diff_conditional(vect.dot(self._j), self._y)
        vz = _diff_conditional(vect.dot(self._k), self._z)

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

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> v = C.x*C.y*C.z * (C.i + C.j + C.k)
        >>> C.delop ^ v
        (-C.x*C.y + C.x*C.z)*C.i + (C.x*C.y - C.y*C.z)*C.j + (-C.x*C.z + C.y*C.z)*C.k
        >>> C.delop.cross(C.i)
        0

        """

        vectx = express(vect.dot(self._i), self.system)
        vecty = express(vect.dot(self._j), self.system)
        vectz = express(vect.dot(self._k), self.system)
        outvec = Vector.zero
        outvec += (diff(vectz, self._y) - diff(vecty, self._z)) * self._i
        outvec += (diff(vectx, self._z) - diff(vectz, self._x)) * self._j
        outvec += (diff(vecty, self._x) - diff(vectx, self._y)) * self._k

        return outvec

    __xor__ = cross

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__



def _diff_conditional(expr, base_scalar):
    """
    First re-expresses expr in the system that base_scalar belongs to.
    If base_scalar appears in the re-expressed form, differentiates
    it wrt base_scalar.
    Else, returns S(0)
    """

    new_expr = express(expr, base_scalar.system, variables = True)
    if base_scalar in new_expr.atoms():
        return diff(new_expr, base_scalar)
    return S(0)

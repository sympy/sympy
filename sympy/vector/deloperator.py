from sympy.core import Basic
from sympy.core.function import Derivative
from sympy.vector.vector import Vector
from sympy.vector.functions import express
from sympy.vector.coordsys import CoordinateSystem, CartesianCoordinateSystem
from sympy.vector.scalar import BaseScalar
from sympy.core import S


class Del(Basic):
    """
    Represents the vector differential operator, usually represented in
    mathematical expressions as the 'nabla' symbol.
    """

    def __new__(cls, system):
        if not isinstance(system, CoordinateSystem):
            raise TypeError("system should be a CoordinateSystem")
        obj = super(Del, cls).__new__(cls, system)
        obj._e1, obj._e2, obj._e3 = system._e1, system._e2, system._e3
        obj._e1cap, obj._e2cap, obj._e3cap = system._e1cap, system._e2cap, system._e3cap
        obj._system = system
        obj._name = system.__str__() + ".delop"
        return obj

    @property
    def system(self):
        return self._system

    def gradient(self, scalar_field, doit=False):
        """
        Returns the gradient of the given scalar field, as a
        Vector instance.

        Parameters
        ==========

        scalar_field : SymPy expression
            The scalar field to calculate the gradient of.

        doit : bool
            If True, the result is returned after calling .doit() on
            each component. Else, the returned expression contains
            Derivative instances

        Examples
        ========

        >>> from sympy.vector import CartesianCoordinateSystem
        >>> C = CartesianCoordinateSystem('C')
        >>> C.delop.gradient(9)
        (Derivative(9, C.x))*C.i + (Derivative(9, C.y))*C.j +
            (Derivative(9, C.z))*C.k
        >>> C.delop(C.x*C.y*C.z).doit()
        C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k

        """

        scalar_field = express(scalar_field, self.system,
                               variables=True)
        lame_params = self.system.lame_parameters()
        vx = 1/lame_params[0]*Derivative(scalar_field, self._e1)
        vy = 1/lame_params[1]*Derivative(scalar_field, self._e2)
        vz = 1/lame_params[2]*Derivative(scalar_field, self._e3)

        if doit:
            return (vx * self._e1cap + vy * self._e2cap + vz * self._e3cap).doit()
        return vx * self._e1cap + vy * self._e2cap + vz * self._e3cap

    __call__ = gradient
    __call__.__doc__ = gradient.__doc__

    def dot(self, vect, doit=False):
        """
        Represents the dot product between this operator and a given
        vector - equal to the divergence of the vector field.

        Parameters
        ==========

        vect : Vector
            The vector whose divergence is to be calculated.

        doit : bool
            If True, the result is returned after calling .doit() on
            each component. Else, the returned expression contains
            Derivative instances

        Examples
        ========

        >>> from sympy.vector import CartesianCoordinateSystem
        >>> C = CartesianCoordinateSystem('C')
        >>> C.delop.dot(C.x*C.i)
        Derivative(C.x, C.x)
        >>> v = C.x*C.y*C.z * (C.i + C.j + C.k)
        >>> (C.delop & v).doit()
        C.x*C.y + C.x*C.z + C.y*C.z

        """
        lame_params = self.system.lame_parameters()
        dprod = lame_params[0] * lame_params[1] * lame_params[2]
        vx = 1/dprod*_diff_conditional(dprod/lame_params[0]*vect.dot(self._e1cap),
                    self._e1)
        vy = 1/dprod*_diff_conditional(dprod/lame_params[1]*vect.dot(self._e2cap),
                    self._e2)
        vz = 1/dprod*_diff_conditional(dprod/lame_params[2]*vect.dot(self._e3cap),
                    self._e3)

        if doit:
            return (vx + vy + vz).doit()
        return vx + vy + vz

    __and__ = dot
    __and__.__doc__ = dot.__doc__

    def cross(self, vect, doit=False):
        """
        Represents the cross product between this operator and a given
        vector - equal to the curl of the vector field.

        Parameters
        ==========

        vect : Vector
            The vector whose curl is to be calculated.

        doit : bool
            If True, the result is returned after calling .doit() on
            each component. Else, the returned expression contains
            Derivative instances

        Examples
        ========

        >>> from sympy.vector import CartesianCoordinateSystem
        >>> C = CartesianCoordinateSystem('C')
        >>> v = C.x*C.y*C.z * (C.i + C.j + C.k)
        >>> C.delop.cross(v, doit = True)
        (-C.x*C.y + C.x*C.z)*C.i + (C.x*C.y - C.y*C.z)*C.j +
            (-C.x*C.z + C.y*C.z)*C.k
        >>> (C.delop ^ C.i).doit()
        0

        """
        lame_params = self.system.lame_parameters()
        dprod = lame_params[0] * lame_params[1] * lame_params[2]
        vectx = express(vect.dot(self._e1cap), self.system, variables=True) * \
                lame_params[0]
        vecty = express(vect.dot(self._e2cap), self.system, variables=True) * \
                lame_params[1]
        vectz = express(vect.dot(self._e3cap), self.system, variables=True) * \
                 lame_params[2]
        outvec = Vector.zero
        outvec += lame_params[0]/dprod*(Derivative(vectz, self._e2) -
                   Derivative(vecty, self._e3)) * self._e1cap
        outvec += lame_params[1]/dprod*(Derivative(vectx, self._e3) -
                   Derivative(vectz, self._e1)) * self._e2cap
        outvec += lame_params[2]/dprod*(Derivative(vecty, self._e1) -
                   Derivative(vectx, self._e2)) * self._e3cap

        if doit:
            return outvec.doit()
        return outvec

    __xor__ = cross
    __xor__.__doc__ = cross.__doc__

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

    new_expr = express(expr, base_scalar.system, variables=True)
    if base_scalar in new_expr.atoms(BaseScalar):
        return Derivative(new_expr, base_scalar)
    return S(0)

from sympy.core import Basic
from sympy.core.function import Derivative
from sympy.vector.vector import Vector
from sympy.vector.functions import express
from sympy.vector.coordsysrect import CoordSystem3D
from sympy.core import S


class Del(Basic):
    """
    Represents the vector differential operator, usually represented in
    mathematical expressions as the 'nabla' symbol.
    """

    def __new__(cls, system):
        if not isinstance(system, CoordSystem3D):
            raise TypeError("system should be a CoordSystem3D")
        obj = super(Del, cls).__new__(cls, system)
        obj._x, obj._y, obj._z = system.x, system.y, system.z
        obj._i, obj._j, obj._k = system.i, system.j, system.k
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

        >>> from sympy.vector import CoordSystem3D
        >>> C = CoordSystem3D('C')
        >>> C.delop.gradient(9)
        (Derivative(9, C.x))*C.i + (Derivative(9, C.y))*C.j + (Derivative(9, C.z))*C.k
        >>> C.delop(C.x*C.y*C.z).doit()
        C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k

        """

        scalar_field = express(scalar_field, self.system,
                               variables=True)
        lame_params = self.system.coordinate_relations()
        vx = 1/lame_params[0]*Derivative(scalar_field, self._x)
        vy = 1/lame_params[1]*Derivative(scalar_field, self._y)
        vz = 1/lame_params[2]*Derivative(scalar_field, self._z)

        if doit:
            return (vx*self._i + vy*self._j + vz*self._k).doit()
        return vx*self._i + vy*self._j + vz*self._k

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

        >>> from sympy.vector import CoordSystem3D
        >>> C = CoordSystem3D('C')
        >>> C.delop.dot(C.x*C.i)
        Derivative(C.x, C.x)
        >>> v = C.x*C.y*C.z * (C.i + C.j + C.k)
        >>> (C.delop & v).doit()
        C.x*C.y + C.x*C.z + C.y*C.z

        """

        lame_params = self.system.coordinate_relations()
        dprod = lame_params[0] * lame_params[1] * lame_params[2]
        vx = 1/dprod*_diff_conditional(dprod/lame_params[0]*vect.dot(self._i),
                                       self._x)
        vy = 1/dprod*_diff_conditional(dprod/lame_params[1]*vect.dot(self._j),
                                       self._y)
        vz = 1/dprod*_diff_conditional(dprod/lame_params[2]*vect.dot(self._k),
                                       self._z)

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

        >>> from sympy.vector import CoordSystem3D
        >>> C = CoordSystem3D('C')
        >>> v = C.x*C.y*C.z * (C.i + C.j + C.k)
        >>> C.delop.cross(v, doit = True)
        (-C.x*C.y + C.x*C.z)*C.i + (C.x*C.y - C.y*C.z)*C.j + (-C.x*C.z + C.y*C.z)*C.k
        >>> (C.delop ^ C.i).doit()
        0

        """

        lame_params = self.system.coordinate_relations()
        dprod = lame_params[0] * lame_params[1] * lame_params[2]
        vectx = express(vect.dot(self._i), self.system, variables=True) * \
                lame_params[0]
        vecty = express(vect.dot(self._j), self.system, variables=True) * \
                lame_params[1]
        vectz = express(vect.dot(self._k), self.system, variables=True) * \
                lame_params[2]
        outvec = Vector.zero
        # assert dprod == 1
        outvec += lame_params[0]/dprod*(Derivative(vectz, self._y) -
                   Derivative(vecty, self._z)) * self._i
        outvec += lame_params[1]/dprod*(Derivative(vectx, self._z) -
                   Derivative(vectz, self._x)) * self._j
        outvec += lame_params[2]/dprod*(Derivative(vecty, self._x) -
                   Derivative(vectx, self._y)) * self._k

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

    new_expr = express(expr, base_scalar.system, variables = True)
    if base_scalar in new_expr.atoms():
        return Derivative(new_expr, base_scalar)
    return S(0)

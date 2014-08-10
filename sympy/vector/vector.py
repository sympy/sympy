from sympy.core.assumptions import StdFactKB
from sympy.core import S, Pow
from sympy.core.expr import AtomicExpr
from sympy import diff as df, sqrt, ImmutableMatrix as Matrix
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.basisdependent import BasisDependent, \
     BasisDependentAdd, BasisDependentMul, BasisDependentZero
from sympy.vector.dyadic import BaseDyadic, Dyadic, DyadicAdd
from sympy.core.compatibility import u


class Vector(BasisDependent):
    """
    Super class for all Vector classes.
    Ideally, neither this class nor any of its subclasses should be
    instantiated by the user.
    """

    is_Vector = True
    _op_priority = 12.0

    @property
    def components(self):
        """
        Returns the components of this vector in the form of a
        Python dictionary mapping BaseVector instances to the
        corresponding measure numbers.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> v = 3*C.i + 4*C.j + 5*C.k
        >>> v.components
        {C.i: 3, C.j: 4, C.k: 5}

        """
        #The '_components' attribute is defined according to the
        #subclass of Vector the instance belongs to.
        return self._components

    def magnitude(self):
        """
        Returns the magnitude of this vector.
        """
        return sqrt(self & self)

    def normalize(self):
        """
        Returns the normalized version of this vector.
        """
        return self / self.magnitude()

    def dot(self, other):
        """
        Returns the dot product of this Vector, either with another
        Vector, or a Dyadic, or a Del operator.
        If 'other' is a Vector, returns the dot product scalar (Sympy
        expression).
        If 'other' is a Dyadic, the dot product is returned as a Vector.
        If 'other' is an instance of Del, returns the directional
        derivate operator as a Python function. If this function is
        applied to a scalar expression, it returns the directional
        derivative of the scalar field wrt this Vector.

        Parameters
        ==========

        other: Vector/Dyadic/Del
            The Vector or Dyadic we are dotting with, or a Del operator .

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> C.i.dot(C.j)
        0
        >>> C.i & C.i
        1
        >>> v = 3*C.i + 4*C.j + 5*C.k
        >>> v.dot(C.k)
        5
        >>> (C.i & C.delop)(C.x*C.y*C.z)
        C.y*C.z
        >>> d = C.i.outer(C.i)
        >>> C.i.dot(d)
        C.i

        """

        from sympy.vector.functions import express
        #Check special cases
        if isinstance(other, Dyadic):
            if isinstance(self, VectorZero):
                return Vector.zero
            outvec = Vector.zero
            for k, v in other.components.items():
                vect_dot = k.args[0].dot(self)
                outvec += vect_dot * v * k.args[1]
            return outvec
        from sympy.vector.deloperator import Del
        if not isinstance(other, Vector) and not isinstance(other, Del):
            raise TypeError(str(other) + " is not a vector, dyadic or " +
                            "del operator")

        #Check if the other is a del operator
        if isinstance(other, Del):
            def directional_derivative(field):
                field = express(field, other.system, variables = True)
                out = self.dot(other._i) * df(field, other._x)
                out += self.dot(other._j) * df(field, other._y)
                out += self.dot(other._k) * df(field, other._z)
                if out == 0 and isinstance(field, Vector):
                    out = Vector.zero
                return out
            return directional_derivative

        if isinstance(self, VectorZero) or isinstance(other, VectorZero):
            return S(0)

        v1 = express(self, other._sys)
        v2 = express(other, other._sys)
        dotproduct = S(0)
        for x in other._sys.base_vectors():
            dotproduct += (v1.components.get(x, 0) *
                           v2.components.get(x, 0))

        return dotproduct

    def __and__(self, other):
        return self.dot(other)
    __and__.__doc__ = dot.__doc__

    def cross(self, other):
        """
        Returns the cross product of this Vector with another Vector or
        Dyadic instance.
        The cross product is a Vector, if 'other' is a Vector. If 'other'
        is a Dyadic, this returns a Dyadic instance.

        Parameters
        ==========

        other: Vector/Dyadic
            The Vector or Dyadic we are crossing with.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> C.i.cross(C.j)
        C.k
        >>> C.i ^ C.i
        0
        >>> v = 3*C.i + 4*C.j + 5*C.k
        >>> v ^ C.i
        5*C.j + (-4)*C.k
        >>> d = C.i.outer(C.i)
        >>> C.j.cross(d)
        (-1)*(C.k|C.i)

        """

        #Check special cases
        if isinstance(other, Dyadic):
            if isinstance(self, VectorZero):
                return Dyadic.zero
            outdyad = Dyadic.zero
            for k, v in other.components.items():
                cross_product = self.cross(k.args[0])
                outer = cross_product.outer(k.args[1])
                outdyad += v * outer
            return outdyad
        elif not isinstance(other, Vector):
            raise TypeError(str(other) + " is not a vector")
        elif (isinstance(self, VectorZero) or
              isinstance(other, VectorZero)):
            return Vector.zero

        #Compute cross product
        def _det(mat):
            """This is needed as a little method for to find the determinant
            of a list in python.
            SymPy's Matrix won't take in Vector, so need a custom function.
            The user shouldn't be calling this.

            """

            return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] *
                                 mat[2][1])
                    + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] *
                    mat[2][2]) + mat[0][2] * (mat[1][0] * mat[2][1] -
                    mat[1][1] * mat[2][0]))

        outvec = Vector.zero
        for system, vect in other.separate().items():
            tempi = system.i
            tempj = system.j
            tempk = system.k
            tempm = [[tempi, tempj, tempk],
                     [self & tempi, self & tempj, self & tempk],
                     [vect & tempi, vect & tempj, vect & tempk]]
            outvec += _det(tempm)

        return outvec

    def __xor__(self, other):
        return self.cross(other)
    __xor__.__doc__ = cross.__doc__

    def outer(self, other):
        """
        Returns the outer product of this vector with another, in the
        form of a Dyadic instance.

        Parameters
        ==========

        other : Vector
            The Vector with respect to which the outer product is to
            be computed.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> N = CoordSysCartesian('N')
        >>> N.i.outer(N.j)
        (N.i|N.j)

        """

        #Handle the special cases
        if not isinstance(other, Vector):
            raise TypeError("Invalid operand for outer product")
        elif (isinstance(self, VectorZero) or
              isinstance(other, VectorZero)):
            return Dyadic.zero

        #Iterate over components of both the vectors to generate
        #the required Dyadic instance
        args = []
        for k1, v1 in self.components.items():
            for k2, v2 in other.components.items():
                args.append((v1*v2) * BaseDyadic(k1, k2))

        return DyadicAdd(*args)

    def __or__(self, other):
        return self.outer(other)
    __or__.__doc__ = outer.__doc__

    def to_matrix(self, system):
        """
        Returns the matrix form of this vector with respect to the
        specified coordinate system.

        Parameters
        ==========

        system : CoordSysCartesian
            The system wrt which the matrix form is to be computed

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> from sympy.abc import a, b, c
        >>> v = a*C.i + b*C.j + c*C.k
        >>> v.to_matrix(C)
        Matrix([
        [a],
        [b],
        [c]])

        """

        return Matrix([self.dot(unit_vec) for unit_vec in
                       system.base_vectors()])

    def separate(self):
        """
        The constituents of this vector in different coordinate systems,
        as per its definition.

        Returns a dict mapping each CoordSysCartesian to the corresponding
        constituent Vector.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> R1 = CoordSysCartesian('R1')
        >>> R2 = CoordSysCartesian('R2')
        >>> v = R1.i + R2.i
        >>> v.separate() == {R1: R1.i, R2: R2.i}
        True

        """

        parts = {}
        for vect, measure in self.components.items():
            parts[vect.system] = (parts.get(vect.system, Vector.zero) +
                                  vect*measure)
        return parts


class BaseVector(Vector, AtomicExpr):
    """
    Class to denote a base vector.
    """

    def __new__(cls, name, index, system, pretty_str, latex_str):
        #Verify arguments
        if not index in range(0, 3):
            raise ValueError("index must be 0, 1 or 2")
        if not isinstance(name, str):
            raise TypeError("name must be a valid string")
        if not isinstance(system, CoordSysCartesian):
            raise TypeError("system should be a CoordSysCartesian")
        #Initialize an object
        obj = super(BaseVector, cls).__new__(cls, S(index),
                                             system)
        #Assign important attributes
        obj._base_instance = obj
        obj._components = {obj: S(1)}
        obj._measure_number = S(1)
        obj._name = name
        obj._pretty_form = u(pretty_str)
        obj._latex_form = latex_str
        obj._system = system

        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)

        #This attr is used for re-expression to one of the systems
        #involved in the definition of the Vector. Applies to
        #VectorMul and VectorAdd too.
        obj._sys = system

        return obj

    @property
    def system(self):
        return self._system

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__


class VectorAdd(BasisDependentAdd, Vector):
    """
    Class to denote sum of Vector instances.
    """

    def __new__(cls, *args, **options):
        obj = BasisDependentAdd.__new__(cls, *args, **options)
        return obj

    def __str__(self, printer=None):
        ret_str = ''
        items = list(self.separate().items())
        items.sort(key = lambda x: x[0].__str__())
        for system, vect in items:
            base_vects = system.base_vectors()
            for x in base_vects:
                if x in vect.components:
                    temp_vect = self.components[x]*x
                    ret_str += temp_vect.__str__(printer) + " + "
        return ret_str[:-3]

    __repr__ = __str__
    _sympystr = __str__


class VectorMul(BasisDependentMul, Vector):
    """
    Class to denote products of scalars and BaseVectors.
    """

    def __new__(cls, *args, **options):
        obj = BasisDependentMul.__new__(cls, *args, **options)
        return obj

    @property
    def base_vector(self):
        """ The BaseVector involved in the product. """
        return self._base_instance

    @property
    def measure_number(self):
        """ The scalar expression involved in the defition of
        this VectorMul.
        """
        return self._measure_number


class VectorZero(BasisDependentZero, Vector):
    """
    Class to denote a zero vector
    """

    _op_priority = 12.1
    _pretty_form = u('0')
    _latex_form = '\mathbf{\hat{0}}'

    def __new__(cls):
        obj = BasisDependentZero.__new__(cls)
        return obj


def _vect_div(one, other):
    """ Helper for division involving vectors. """
    if isinstance(one, Vector) and isinstance(other, Vector):
        raise TypeError("Cannot divide two vectors")
    elif isinstance(one, Vector):
        if other == S.Zero:
            raise ValueError("Cannot divide a vector by zero")
        return VectorMul(one, Pow(other, S.NegativeOne))
    else:
        raise TypeError("Invalid division involving a vector")


Vector._expr_type = Vector
Vector._mul_func = VectorMul
Vector._add_func = VectorAdd
Vector._zero_func = VectorZero
Vector._base_func = BaseVector
Vector._div_helper = _vect_div
Vector.zero = VectorZero()
